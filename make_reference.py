import os
import pandas as pd
import numpy as np
import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pysam




def make_fa_and_gtf_reference(probe_input, output_dir, output_file_prefix, 
                              fwd_primer = "NNNNN", rev_primer = "NNNNN"):
    if isinstance(probe_input, pd.DataFrame):
        probe_df = probe_input
    elif os.path.isfile(probe_input):
        probe_df = pd.read_csv(probe_input, delimiter=',')
    else:
        raise ValueError("The first argument need to be a file or pandas.DataFrame")
    
    probe_df = probe_df.sort_values(by="position")
    fasta_string = ""
    gtf_lines = []
    max_position = 0
    for i in range(len(probe_df)):   
        probe = probe_df["probe"].iloc[i]
        sequence = probe_df["sequence"].iloc[i]
        position = probe_df["position"].iloc[i]
        strand = probe_df["strand"].iloc[i]

        if position > max_position:
            fasta_string += _add_primer(sequence, fwd_primer, rev_primer)
            max_position = position + len(sequence) - 1
        else:
            if max_position-position+1 < len(sequence):
                fasta_string = fasta_string[:-len(rev_primer)] + \
                               _add_primer(sequence[max_position-position+1:],"",rev_primer)
                max_position = max_position - position + len(sequence)
                
        s = len(fasta_string) - len(rev_primer) - len(sequence) + 1
        gtf_lines.append("\t".join(["chr1", probe, "exon", str(s), \
                        str(s+len(sequence)-1), ".", strand, "0", \
                        f'gene_id "{probe}"; transcript_id "{probe}";']) + "\n")
    os.makedirs(output_dir, exist_ok=True)
    full_prefix = os.path.join(output_dir,output_file_prefix)
    with open("%s.fasta"%full_prefix, "wt") as fid_fasta:
        fid_fasta.write(">chr1\n")
        fid_fasta.write(fasta_string)
        
    with open("%s.gtf"%full_prefix, "wt") as fid_gtf:
        for line in gtf_lines:
            fid_gtf.write(line)

    gtf = open(f"{full_prefix}.gtf").read().splitlines()[:-1]
    gtf = [line.split("\t") for line in gtf]
    genome = open(f"{full_prefix}.fasta").read().splitlines()[1]

    error = []
    for line in gtf:
        s = int(line[3])
        e = int(line[4])
        if list(probe_df[probe_df["probe"]==line[1]]["sequence"])[0] != genome[s-1:e]:
            error.append((line[1],list(probe_df[probe_df["probe"]==line[1]]["sequence"])[0], genome[s-1:e]))
    
    with open(f"{full_prefix}_error.txt", "w") as file:
        for t in error:
            file.write(str(t) + "\n")

    


def map_probes(probe_file, bwa_path, reference_genome, merge_overlapping_probes=True,
               output_probe_file="mapped_probes.csv"):
    output_dir = os.path.dirname(probe_file)
    ref_dir = os.path.dirname(reference_genome)
    probe_fa_file = os.path.join(ref_dir, "probe.fa")
    
    _write_probe_to_fasta(probe_file, probe_fa_file)
    
    sam_file = os.path.join(ref_dir,"probe_aligned.sam")
    _align_probe_with_bwa(bwa_path, reference_genome, probe_fa_file, sam_file)
    

    with pysam.AlignmentFile(sam_file, "r") as f:
        data = [{"probe": read.query_name,
                 "gene": read.query_name.split("_")[0],
                 "strand": "+" if not read.is_reverse else "-",
                 "position": read.reference_start,
                 "sequence": read.query_sequence,
                 "cigar": read.cigarstring}
                for read in f.fetch()]

    probe_df = pd.DataFrame(data)
    
    probe_df.to_csv(os.path.join(output_dir,"%s_unmerged.csv"%os.path.splitext(output_probe_file)[0]))
    if merge_overlapping_probes:        
        probe_df = _merge_overlapping_probes(probe_df)
    
    print(os.path.join(output_dir,output_probe_file))
    probe_df.to_csv(os.path.join(output_dir,output_probe_file))
    return probe_df
    
    
    

    
def build_reference_index_with_cell_ranger(path_to_cell_ranger, genome_name, 
                                           fa_gtf_prefix, output_dir):
    cur_dir = os.getcwd()
    fa_file = os.path.abspath(fa_gtf_prefix + ".fasta")
    gtf_file = os.path.abspath(fa_gtf_prefix + ".gtf")
    cmd = [path_to_cell_ranger, "mkref", "--genome=%s"%genome_name,
                    "--fasta=%s"%fa_file,
                    "--genes=%s"%gtf_file]
    print(" ".join(cmd))
    os.chdir(os.path.dirname(output_dir))
    subprocess.call(cmd)
    os.chdir(cur_dir)
    
    


def _add_primer(probe_sequence,fwd_primer="",rev_primer=""):
    return f'{fwd_primer}{probe_sequence}{rev_primer}'
    
    
    
    
def _write_probe_to_fasta(probe_file, probe_fa_file):
    output_dir = os.path.dirname(probe_file)
    probe_df = pd.read_csv(probe_file, delimiter=',')
    gene = sorted(set(probe_df["gene"]))

    seqs = []
    for i in range(len(gene)):
        probes = probe_df[probe_df["gene"]==gene[i]]

        for j in range(len(probes)):
            probe_seq = list(probes["sequence"])[j]
            probe_seq = SeqRecord(Seq(probe_seq), id="%s_%d"%(gene[i],j+1), description="%s_%d"%(gene[i],j+1))
            seqs.append(probe_seq)

    with open(probe_fa_file, "w") as fa_handle:
        SeqIO.write(seqs, fa_handle, "fasta")
        
        
        

def _align_probe_with_bwa(bwa_path, reference_genome, probe_fa_file, sam_file):
    bwa_idx_cmd = [bwa_path, "index", 
               os.path.abspath(reference_genome)]
    subprocess.call(bwa_idx_cmd)

    bwa_mem_cmd = [bwa_path, "mem", reference_genome, probe_fa_file]
    with open(sam_file, "w") as sam_output_file:
        subprocess.call(bwa_mem_cmd, stdout=sam_output_file)
        

        
        
def _merge_overlapping_probes(probe_df):
    genes = sorted(set(probe_df["gene"]))
    
    result = pd.DataFrame()
    for gene in genes:
        probe = probe_df[probe_df["gene"] == gene]
        probe = probe.sort_values(by="position")

        position = []
        sequence = []
        pos = list(probe["position"])[0]
        seq = list(probe["sequence"])[0]
        for i in range(1,len(probe)):
            p = list(probe["position"])[i]
            if p - pos < len(seq):
                s = list(probe["sequence"])[i]
                seq = seq + s[pos+len(seq)-p:]
            else:
                position.append(pos)
                sequence.append(seq)
                pos = list(probe["position"])[i]
                seq = list(probe["sequence"])[i]

        position.append(pos)
        sequence.append(seq)
        n = len(position)

        result = pd.concat([result, 
                pd.DataFrame(data={
                    "probe": ["%s_%d"%(gene,i+1) for i in range(n)],
                    "gene": [gene for _ in range(n)],
                    "strand": [list(probe["strand"])[0] for _ in range(n)],
                    "position": position,
                    "sequence": sequence
                })])
    return result


