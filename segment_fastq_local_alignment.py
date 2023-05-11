import os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from itertools import islice
import subprocess
from .utils import parse_construct_sequence




def pairwise_alignment(read, query,
                       target_open_gap_score=-2,
                       target_extend_gap_score=-0.5,
                       query_open_gap_score=-2,
                       query_extend_gap_score=-0.5):
    # Create a PairwiseAligner object
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.target_internal_open_gap_score = target_open_gap_score
    aligner.target_internal_extend_gap_score = target_extend_gap_score
    aligner.query_internal_open_gap_score = query_open_gap_score
    aligner.query_internal_extend_gap_score = query_extend_gap_score

    alignment = aligner.align(read, query)[0]
    unmatched = _count_unmatch(alignment)
    start_coord, end_coord = alignment.coordinates[0,]
    return unmatched, start_coord, end_coord




def _count_unmatch(alignment):
    matches = sum(1 for i in range(len(alignment[0])) if alignment[0][i] == alignment[1][i])
    total_length = len(alignment[0])
    return total_length - matches




def _get_quality_string(record):
    quality_scores = record.letter_annotations["phred_quality"]
    quality_string = ''.join(chr(q + 33) for q in quality_scores)
    return quality_string




def get_coord_with_local_align(read,construct):
    coord = []
    try:
        for segment in construct:
            if segment[1] == "anchor":
                aln = pairwise_alignment(read, segment[0])
                if aln[0] <= 1:
                    coord.append((aln[1],aln[2]))
                else:
                    coord.append(None)
            else:
                coord.append(None)

        for i in range(len(construct)):
            segment =construct[i]

            if i>= 1:
                segment_before = construct[i-1]
            else:
                segment_before = (None,None)
            if i < (len(construct)-1):
                segment_after = construct[i+1]
            else:
                segment_after = (None,None)

            if segment[1] != "anchor":
                if segment_before[1] == "anchor" and segment_after[1] == "anchor":
                    coord[i] = (coord[i-1][1], coord[i+1][0])
                elif segment_before[1] == "anchor":
                    coord[i] = (coord[i-1][1], coord[i-1][1]+len(segment[0]))
                elif segment_after[1] == "anchor":
                    coord[i] = (coord[i+1][0]-len(segment[0]), coord[i+1][0])
                else:
                    coord[i] = (coord[i-1][1], coord[i-1][1]+len(segment[0]))
    except:
        coord = None
    return coord




def segment_fastq_local_alignment(construct_file, r1_fastq, r2_fastq, output_dir, sample_identifier, file_prefix):
    os.makedirs(output_dir, exist_ok=True)
    sample_identifier.append("unidentified")
    file_handles = {}
    construct_sequence = parse_construct_sequence(construct_file)

    for sample_id in sample_identifier:
        sample_dir = "%s/%s" %(output_dir,sample_id)
        os.makedirs(sample_dir, exist_ok=True)
        for read in ["R1", "R2"]:
            file_handle = open(os.path.join(sample_dir,\
                    "%s_%s_S1_L001_%s_001.fastq"%(file_prefix,sample_id,read)), "w")
            file_handles["%s_%s"%(sample_id,read)] = file_handle


    forward_reads = SeqIO.parse(r1_fastq, "fastq")
    reverse_reads = SeqIO.parse(r2_fastq, "fastq")
    for R1, R2 in zip(forward_reads, reverse_reads):
        if "R1" in construct_sequence:
            pass

        R2_construct = construct_sequence["R2"]
        UMI = ""
        sample_id = ""
        probe = R2

        coord = get_coord_with_local_align(R2,R2_construct)
        if coord:
            for i in range(len(R2_construct)):
                if R2_construct[i][1] == "UMI":
                    UMI = R2[coord[i][0]:coord[i][1]].seq._data.decode()
                elif R2_construct[i][1] == "probe":
                    probe = R2[coord[i][0]:coord[i][1]]
                elif R2_construct[i][1] == "sample_identifier":
                    sample_id = R2[coord[i][0]:coord[i][1]].seq._data.decode()
            read_id = "%s-UMI:%s" % (R2.id, UMI)
            R1.id = read_id
            probe.id = read_id
            R1.description = R1.id
            probe.description = probe.id
        if sample_id in sample_identifier:
            SeqIO.write(R1, file_handles["%s_R1"%sample_id], "fastq")
            SeqIO.write(probe, file_handles["%s_R2"%sample_id], "fastq")
        else:
            SeqIO.write(R1, file_handles["unidentified_R1"], "fastq")
            SeqIO.write(probe, file_handles["unidentified_R2"], "fastq")

    for sample_id in file_handles:
        file_handles[sample_id].close()

        
        
        
def run_cellranger_alignment(path_to_cell_ranger, fq_prefix, transcriptome,
                             fastq_dir, localcores=4, localmem=32):
    cmd = [path_to_cell_ranger, "count", "--id=%s"%fq_prefix, "--chemistry=threeprime",
           "--transcriptome=%s"%transcriptome, "--fastqs=%s"%fastq_dir,
           "--localcores=%d"%localcores, "--localmem=%d"%localmem]
    try:
        subprocess.call(cmd)
    except:
        pass        

        
        