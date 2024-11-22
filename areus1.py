#imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import openpyxl
import gzip
from Bio import SeqIO

# input information
input_organism_name = 'areus'
refseq_Assembly_ID = 'GCF_000013425.1'
element_name_gff = ['gene']#['gene','pseudogene']

input_plot_len = 8000
tion_file_path_out_refseq = r'A:\PhD\genomes\GRCh38.p14\ncbi_dataset\data\GCF_000001405.40\GCF_000001405.40_GRCh38.p14_rm.out'
fasta_path_refseq = r'A:\PhD\genomes\areus\ncbi_dataset\data\GCF_000013425.1\GCF_000013425.1_ASM1342v1_genomic.fna'
annotation_file_path_GFF = r'A:\PhD\genomes\areus\ncbi_dataset\data\GCF_000013425.1\genomic.gff'
fasta_path_UCSC = r'A:\PhD\genomes\areus\ncbi_dataset\data\GCA_000013425.1\GCA_000013425.1_ASM1342v1_genomic.fna'


from func_for_bac import open_simple_fasta_file_and_prepare_it
from func_for_bac import counting_of_elements
from func_for_bac import nucleotide_composition_in_seqenence_ATCG
from func_for_bac import plus_minis_both_1
#refseq

#analysis for genes



Chr_sequence_refseq = open_simple_fasta_file_and_prepare_it(fasta_path_refseq)
#plus,comp,both = counting_of_elements(Chr_sequence_refseq,annotation_file_path_GFF,element_name_gff,input_plot_len)
plus,comp,both = plus_minis_both_1(Chr_sequence_refseq,annotation_file_path_GFF,element_name_gff)
print('++++++++++++++')
nucleotide_composition_in_seqenence_ATCG(plus)
print('-----------------------')
nucleotide_composition_in_seqenence_ATCG(comp)
print('-+-+-+-+-+-+-+-+-+-+')
nucleotide_composition_in_seqenence_ATCG(both)



