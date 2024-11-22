#imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import openpyxl
import gzip
from Bio import SeqIO

# input information
input_organism_name = 'GRCh38.p14'
refseq_Assembly_ID = 'GCF_000001405.40'
input_element_name = ['LINE','SINE','LTR','DNA','RC/Helitron','tRNA','scRNA','Retroposons','RNA','srpRNA'] #study genomic elements
element_name_gff = ['gene']#['gene','pseudogene']
chromosomes_names = ['10']#['1','2','3','4','5','6','7','8','9' ,'10', '11', '12', '13', '14', '15', '16','17', '18', '19', '20', '21', '22','X','Y']

input_plot_len = 8000
annotation_file_path_out_UCSC = r'A:\PhD\genomes\GRCh38.p14\UCSC\hg38.p14.fa.out'
annotation_file_path_out_refseq = r'A:\PhD\genomes\GRCh38.p14\ncbi_dataset\data\GCF_000001405.40\GCF_000001405.40_GRCh38.p14_rm.out'
fasta_path_refseq = r'A:\PhD\genomes\GRCh38.p14\ncbi_dataset\data\GCF_000001405.40\GCF_000001405.40_GRCh38.p14_genomic.fna'
annotation_file_path_GFF = r'A:\PhD\genomes\GRCh38.p14\ncbi_dataset\data\GCF_000001405.40\genomic.gff'
fasta_path_UCSC = r'A:\PhD\genomes\GRCh38.p14\UCSC\hg38.p14.fa.gz'



from my_functions import  check_and_creat_refseq_file
from functions_for_primats import open_fasta_file_multi_record_UCSC_and_refseq
from functions_for_primats import portion_of_TE_elements_in_geneme_with_UCSC_and_refseq
from functions_for_primats import nucleotide_composition_in_seqenence_ATCG
from functions_for_primats import nucleotide_composition_in_seqenence_ATCG_2
from functions_for_primats import test
from functions_for_primats import merge_seq_seq_share_with_1234


refseq_records_lists = check_and_creat_refseq_file(refseq_Assembly_ID,input_organism_name,chromosomes_names)
print(refseq_records_lists)

#refseq

#analysis for genes

for Chr_number,Chr_ID in zip(chromosomes_names,refseq_records_lists):
    Chr = [Chr_number,Chr_ID]

    Chr_sequence_refseq = open_fasta_file_multi_record_UCSC_and_refseq(fasta_path_refseq,Chr)
    plus,comp,both = test(Chr_sequence_refseq,annotation_file_path_GFF,Chr,element_name_gff)
    #nucleotide_composition_in_seqenence_atcg(plus)
    #nucleotide_composition_in_seqenence_atcg(comp)
    #nucleotide_composition_in_seqenence_atcg(both)



    #TEs **********************************************
    Plus_strand,Complement_strand,Both = portion_of_TE_elements_in_geneme_with_UCSC_and_refseq(Chr_sequence_refseq,annotation_file_path_out_refseq,Chr,input_element_name)
    
    #nucleotide_composition_in_seqenence_atcg(Plus_strand)
    #nucleotide_composition_in_seqenence_atcg(Complement_strand)
    #nucleotide_composition_in_seqenence_atcg(Both)

    merged_sequence = merge_seq_seq_share_with_1234(Both,both)

    nucleotide_composition_in_seqenence_ATCG_2(merged_sequence)



















#UCSC

'''for Chr_number,Chr_ID in zip(chromosomes_names,refseq_records_lists):
    Chr = [Chr_number,Chr_ID]
    chr_sequence_UCSC = open_fasta_file_UCSC_multi_record(fasta_path_UCSC,Chr)
    P,C,B = portion_of_TE_elements_in_geneme_with_UCSC(chr_sequence_UCSC,annotation_file_path_out_UCSC,Chr,input_element_name)
    nucleotide_composition_in_seqenence_atcg(P)
    nucleotide_composition_in_seqenence_atcg(C)
    nucleotide_composition_in_seqenence_atcg(B)'''












    







'''for Chr in chromosomes_names:
    open_fasta_file_multi_record_from_refseq(fasta_path_refseq,Chr)'''

















