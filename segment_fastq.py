import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from itertools import islice
import subprocess
import probacseq_toolkit
import regex
import gzip




def make_construct_regex_pattern(construct_df, insertion=1, deletion=1,
                                 substitution=1, error=2):
    
    construct_regex_pattern = ""
    for i in range(len(construct_df)):
        c = construct_df.iloc[i]
        if c["sequence"] == ".":
            length = ",".join(c["length"].split("-"))
            p = "(?P<%s>.{%s})" % (c["segment_name"], length)
        elif "|" in c["sequence"]:
            p = "(?P<%s>(%s))" % (c["segment_name"], c["sequence"])
        else:
            p = "(?P<%s>%s{i<=%d,d<=%d,s<=%d,e<=%d})" % (c["segment_name"], c["sequence"], \
                                                     insertion, deletion, substitution, error)
        construct_regex_pattern += p
    return construct_regex_pattern
    
    


def segment_fastq(construct_regex, r1_fastq, r2_fastq, use_reverse_complement,
                  output_dir, file_prefix):
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs("%s/unidentified"%output_dir, exist_ok=True)
    UMI = "umi" in construct_regex
    
    with gzip.open(r1_fastq, "rt") as r1_fastq_gz, \
         gzip.open(r2_fastq, "rt") as r2_fastq_gz, \
         open("%s/%s_R1_001.fastq"%(output_dir, file_prefix), "w") as r1_file, \
         open("%s/%s_R2_001.fastq"%(output_dir, file_prefix), "w") as r2_file, \
         open("%s/unidentified/unidentified_%s_R1_001.fastq"%(output_dir, file_prefix), 
                  "w") as r1_unidentified_file, \
         open("%s/unidentified/unidentified_%s_R2_001.fastq"%(output_dir, file_prefix), 
              "w") as r2_unidentified_file:

        forward_reads = SeqIO.parse(r1_fastq_gz, "fastq")
        reverse_reads = SeqIO.parse(r2_fastq_gz, "fastq")
        for R1, R2 in zip(forward_reads, reverse_reads):
            s = regex.search(construct_regex, R2.seq._data.decode())
            if s:                
                probe = R2[s.start("probe"):s.end("probe")]
                if use_reverse_complement:
                    probe = probe.reverse_complement(id=True, description=True, name=True)
                
                if UMI:
                    R1.seq = R1[:16].seq + R2[s.start("umi"):s.end("umi")].seq
                SeqIO.write(R1, r1_file, "fastq")
                SeqIO.write(probe, r2_file, "fastq")
            else:
                SeqIO.write(R1, r1_unidentified_file, "fastq")
                SeqIO.write(R2, r2_unidentified_file, "fastq")




def segment_fastq_multiplex(construct_regex, r1_fastq, r2_fastq, use_reverse_complement,
                            sample_identifiers, output_dir, file_prefix):
    os.makedirs(output_dir, exist_ok=True)
    sample_identifiers.append("unidentified")
    file_handles = {}
    UMI = "umi" in construct_regex
    
    for sample_id in sample_identifiers:
        sample_dir = "%s/%s" %(output_dir,sample_id)
        os.makedirs(sample_dir, exist_ok=True)
        for read in ["R1", "R2"]:
            file_handle = open(os.path.join(sample_dir,\
                    "%s_%s_S1_L001_%s_001.fastq"%(file_prefix,sample_id,read)), "w")
            file_handles["%s_%s"%(sample_id,read)] = file_handle
    
    with gzip.open(r1_fastq, "rt") as r1_fastq_gz, \
         gzip.open(r2_fastq, "rt") as r2_fastq_gz:
        forward_reads = SeqIO.parse(r1_fastq_gz, "fastq")
        reverse_reads = SeqIO.parse(r2_fastq_gz, "fastq")
        for R1, R2 in zip(forward_reads, reverse_reads):
            s = regex.search(construct_regex, R2.seq._data.decode())
            if s:                
                probe = R2[s.start("probe"):s.end("probe")]
                sample_id = s.group("sample")
                if use_reverse_complement:
                    probe = probe.reverse_complement(id=True, description=True, name=True)
                    sample_id = str(Seq(sample_id).reverse_complement())
                if sample_id not in sample_identifiers:
                    sample_id = "unidentified"
                
                if UMI:
                    r1 = R1[:16] + R2[s.start("umi"):s.end("umi")]
                    r1.description = R1.description
                    R1 = r1
                
                SeqIO.write(R1, file_handles["%s_R1"%sample_id], "fastq")
                SeqIO.write(probe, file_handles["%s_R2"%sample_id], "fastq")
            else:
                SeqIO.write(R1, file_handles["unidentified_R1"], "fastq")
                SeqIO.write(R2, file_handles["unidentified_R2"], "fastq")

    for sample_id in file_handles:
        file_handles[sample_id].close()

        
        
        
def run_cellranger_alignment(path_to_cell_ranger, id, transcriptome,
                             fastq_dir, localcores=4, localmem=32):
    cmd = [path_to_cell_ranger, "count", "--id=%s"%id, "--chemistry=threeprime",
           "--transcriptome=%s"%transcriptome, "--fastqs=%s"%fastq_dir,
           "--localcores=%d"%localcores, "--localmem=%d"%localmem]
    print("This is the cell ranger command:")
    print(" ".join(cmd))
    subprocess.call(cmd)

