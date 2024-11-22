import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import openpyxl
import gzip
from Bio import SeqIO,Entrez
#from memory_profiler import profile

def check_and_creat_refseq_file(organism_refseq_Assembly_id,input_organism_name,chromosomes_names,):
    import os
    refseq_records = []
    if os.path.exists(f'{input_organism_name}_refseq_records.txt'):
        print('yes refseq_records file exist')

        file = open(f'{input_organism_name}_refseq_records.txt','r')
        for i in file:
            refseq_records.append(i.strip())
        len_file  = len(refseq_records)
        print('len_file',len_file,'and len chromosomes_names',len(chromosomes_names))
        #counter = [i.strip() for i in file]


        if len_file == len(chromosomes_names):
            print('they have equal len Pass!!!!')
            print('refseq_records',refseq_records)

        else:
            print('they have not equal len')
            refseq_records = []
            for Chr in chromosomes_names:
                A = chromosome_number_to_refseq_record(organism_refseq_Assembly_id,Chr)
                refseq_records.append(A)

            with open(f'{input_organism_name}_refseq_records.txt','w') as file:
                for i in refseq_records:
                    file.write(f'{i}\n')
            print('finnally',refseq_records)
    else:
        print('i dont find refseq list file')
        refseq_records = []
        for Chr in chromosomes_names:
            A = chromosome_number_to_refseq_record(organism_refseq_Assembly_id,Chr)
            refseq_records.append(A)


        with open(f'{input_organism_name}_refseq_records.txt','w') as file:
            for i in refseq_records:
                file.write(f'{i}\n')
            print('finnally',refseq_records)
    return refseq_records


def chromosome_number_to_refseq_record(assembly_accession, chromosome_number):
    Entrez.email = 'reza.sotoudeh.van@gmail.com'

    # Construct a query to search for the chromosome in GenBank
    query = f"{assembly_accession}[Assembly] AND Chr{chromosome_number}[All Fields] AND refseq[Filter]"

    try:
        # Perform the search
        handle = Entrez.esearch(db="nucleotide", term=query)
        record = Entrez.read(handle)

        if record['IdList']:
            # Fetch the first accession number from the search results
            accession_id = record['IdList'][0]

            # Fetch the detailed information for the accession
            handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text")
            genbank_record = handle.read()
            #print(genbank_record)

            # Extract accession number from the GenBank record
            accession_number = genbank_record.split("ACCESSION")[1].split()[0].strip()
            VERSION_number = genbank_record.split("VERSION")[1].split()[0].strip()
            print(accession_number)
            print(VERSION_number)
            print(genbank_record.split("DEFINITION")[1].split(',')[0].strip())
            return accession_number
        else:
            print(f"No results found for chromosome {chromosome_number} in assembly {assembly_accession}")
            return None

    except Exception as e:
        print(f"An error occurred: {e}")

        return None

def get_chromosome_definition_from_accession(accession):
    # Provide your email address to NCBI (required for Entrez)
    Entrez.email = "reza.sotoudeh.van@gmail.com"

    # Define the database and the accession number
    database = "nucleotide"
    accession_number = accession


    # Fetch information from NCBI using Entrez
    handle = Entrez.efetch(db=database, id=accession_number, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")

    # Extract chromosome definition from the GenBank record
    chromosome_definition = record.description
    chromosome_definition = chromosome_definition.split()
    chromosome_name =chromosome_definition[3].replace(',','')
    print(chromosome_name)
    return chromosome_name


def open_simple_fasta_file_and_prepare_it(self):
    fasta_path = self.chromosome_sequence
    ch_number = self.number_of_chromosome
    fasta_file = open(fasta_path)
    fasta_file.readline()
    sequence = fasta_file.read()
    sequence = sequence.replace('\n', '')
    fasta_file.close()
    print(ch_number)
    len_chromosome = len(sequence)
    return sequence , len_chromosome


def open_fasta_file_UCSC_multi_record(seq_path,ch_number0):
    sequence = ''
    fasta_path = seq_path

    ch_number = int(ch_number0)
    if ch_number == 23 :
        ch_number = 'X'
    if ch_number == 24 :
        ch_number = 'Y'
    if fasta_path[-2:] == 'gz':
        fasta_path = gzip.open(fasta_path,'rt')
    print('for chromosoeme',ch_number)
    for seq_record in SeqIO.parse(fasta_path, "fasta"):
        if seq_record.id == rf'chr{ch_number}':
            print('sequence id is',seq_record.id)
            sequence = str(seq_record.seq)
            if fasta_path[-2:] == 'gz':
                fasta_path.close()
            l  = len(sequence)
            print('len chromosome is',l,'bp')

            break

    return sequence


def open_fasta_file_multi_record_from_refseq(path,ch_number):
    sequence = ''
    fasta_path = path
    ch_number = ch_number
    print('ch_ref',ch_number)
    for seq_record in SeqIO.parse(fasta_path, "fasta"):
        sequence_id = seq_record.id
        #print(sequence_id)
        if ch_number in sequence_id :
            print('sequence_id',sequence_id)
            ch_integer_number = get_chromosome_definition_from_accession(ch_number)
            sequence = seq_record.seq
            print('*************','Chr',ch_number,sequence_id,ch_integer_number,'**************')
            print('Chr length',len(seq_record.seq)/1000000,'Mbp')
            break



    return sequence


def open_annotation_file_and_prepare_it(anotation_path,ch_number0):
    chr_annotation = anotation_path
    ch_number = ch_number0
    if chr_annotation[-2:] == 'gz':
        file_gff = gzip.open(chr_annotation,'rt')
    else:
        file_gff = open(chr_annotation)
    print('mapping start', ch_number)
    return file_gff


def counting_of_elements( sequence, gff,element, plot_length):

    counting_of_plus_elements = 0
    counting_of_comp_elements = 0
    counting_of_both_elements = 0
    sequences_of_plus_elements = ''
    sequences_of_complement_elements = ''
    sequences_of_both_elements = ''
    #elements_overlap = list(sequence)
    counting_pluses_in_difference_length = [0] *plot_length
    counting_complements_in_difference_length = [0] * plot_length
    counting_both_in_difference_length = [0] *plot_length

    compposintion_in_diffrence_length_plus = [''] * plot_length
    compposintion_in_diffrence_length_comp = [''] * plot_length
    compposintion_in_diffrence_length_both = [''] * plot_length
    # element_location_d = [0] * len(sequence)
    # Alus_location_c = [0] * len(sequence)

    '''A_portion_in_length = [0] * 4000                                                                                                                                                         
    T_portion_in_length = [0] * 4000                                                                                                                                                            
    C_portion_in_length = [0] * 4000                                                                                                                                                            
    G_portion_in_length = [0] * 4000'''

    chromosome_number = self.number_of_chromosome
    element_name = element


    #sequence2_for_overlabing = sequence
    countery  = 0
    for annotation in gff:
        if annotation[0] == '#':
            continue
        annotation.strip()
        annotation = annotation.split()
        if  annotation[3] == '1':
            #print(annotation)
            continue


        repeat_control = []
        if (element_name in annotation[2]) and f'NC_0000{chromosome_number}' in annotation[0] : # or annotation[4] == f'chr{chromosome_number}_'):
            #print(annotation)
            if int(annotation[4]) - int(annotation[3]) == 0:
                print('len was 0')
                continue

            countery +=1
            if annotation[3] in repeat_control:
                print(annotation[3],annotation[4])
            repeat_control.append(annotation[3])
            if annotation[6] == '+':
                sequence_of_element_p = sequence[int(annotation[3]) - 1:int(annotation[4])]

                #print('+',sequence_of_element_p)
                sequences_of_plus_elements += sequence_of_element_p
                #elements_overlap[int(annotation[5]) - 1:int(annotation[6]) - 1] = sequence_of_element_p.lower()
                #element_location_d[int(line[5])] += 10
                len_element_p = len(sequence_of_element_p)
                if len_element_p == 0:
                    print(0)
                    print(sequence_of_element_p)



                if len_element_p > plot_length:
                    print('Not studied Element it is very big',len_element_p)
                if len_element_p < plot_length:
                    counting_of_plus_elements += 1
                    counting_pluses_in_difference_length[len_element_p] += 1
                    compposintion_in_diffrence_length_plus[len_element_p] += sequence_of_element_p
                if int(annotation[4]) - int(annotation[3]) == 0:
                    print('len was 0 not continu')

                '''A_portion_in_length[len_element] += Eplus.count('A')                                                                                                                         
                T_portion_in_length[len_element] += Eplus.count('T')                                                                                                                            
                C_portion_in_length[len_element] += Eplus.count('C')                                                                                                                            
                G_portion_in_length[len_element] += Eplus.count('G')'''
            else:
                #print('-')
                sequence_of_element_c = sequence[int(annotation[3]) - 1:int(annotation[4])]
                sequences_of_complement_elements += sequence_of_element_c
                #print('-',sequence_of_element_c)
                #elements_overlap[int(annotation[5]) - 1:int(annotation[6]) - 1] = sequence_of_element_c.lower()
                #Alus_location_c[int(line[5])] += 10
                len_element_c = len(sequence_of_element_c)
                if len_element_c > plot_length:
                    print('Not studied Element it is very big',len_element_c)
                if len_element_c < plot_length:
                    counting_of_comp_elements += 1
                    counting_complements_in_difference_length[len_element_c] += 1
                    compposintion_in_diffrence_length_comp[len_element_c]+=sequence_of_element_c
                if int(annotation[4]) - int(annotation[3]) == 0:
                    print('len was 0 not continu')
                    #print(counting_pluses_in_difference_length)
    #both making#########################################################
    counting_of_both_elements = counting_of_plus_elements + counting_of_comp_elements
    counting_both_in_difference_length = [sum(i) for i in zip(counting_pluses_in_difference_length,counting_complements_in_difference_length)]
    sequences_of_both_elements = sequences_of_plus_elements + sequences_of_complement_elements
    compposintion_in_diffrence_length_both = [i+j for i,j in zip(compposintion_in_diffrence_length_plus,compposintion_in_diffrence_length_comp)]

    print('element counting finished')
    print('both   plus    comp')
    print(counting_of_both_elements,'  ',counting_of_plus_elements,'  ',counting_of_comp_elements)
    plus_elements_portion_in_sequence = round(len(sequences_of_plus_elements)/len(sequence)*100,2)
    complement_elements_portion_in_sequence = round(len(sequences_of_complement_elements)/len(sequence)*100,2)
    both_elements_portion_in_sequence = round(len(sequences_of_both_elements)/len(sequence)*100,2)
    print(' portion of elements in sequence :\nboth plus comp\n',both_elements_portion_in_sequence,plus_elements_portion_in_sequence,'%',complement_elements_portion_in_sequence,'%')
    len_both_seq = len(sequences_of_both_elements)
    A_b = round(sequences_of_both_elements.count('a') / len_both_seq * 100, 2)
    T_b = round(sequences_of_both_elements.count('t') / len_both_seq * 100, 2)
    C_b = round(sequences_of_both_elements.count('c') / len_both_seq * 100, 2)
    G_b = round(sequences_of_both_elements.count('g') / len_both_seq * 100, 2)
    len_sum_of_pluses = len(sequences_of_plus_elements)
    A_d = round(sequences_of_plus_elements.count('a') / len_sum_of_pluses * 100, 2)
    T_d = round(sequences_of_plus_elements.count('t') / len_sum_of_pluses * 100, 2)
    C_d = round(sequences_of_plus_elements.count('c') / len_sum_of_pluses * 100, 2)
    G_d = round(sequences_of_plus_elements.count('g') / len_sum_of_pluses * 100, 2)
    len_sum_of_complements = len(sequences_of_complement_elements)
    A_c = round(sequences_of_complement_elements.count('a') / len_sum_of_complements * 100, 2)
    T_c = round(sequences_of_complement_elements.count('t') / len_sum_of_complements * 100, 2)
    C_c = round(sequences_of_complement_elements.count('c') / len_sum_of_complements * 100, 2)
    G_c = round(sequences_of_complement_elements.count('g') / len_sum_of_complements * 100, 2)
    print('both composition','+composition', '-composition')
    print('A ',A_b,'    ','A ', A_d, '   ', ' A ', A_c)
    print('T ',T_b,'    ','T ', T_d, '   ', ' T ', T_c)
    print('C ',C_b,'    ','C ', C_d, '   ', ' C ', C_c)
    print('G ',G_b,'    ','G ', G_d, '   ', ' G ', G_c)
    #print(counting_pluses_in_difference_length)
    both_stand_information = [counting_of_both_elements,counting_both_in_difference_length,sequences_of_both_elements,len(sequence),compposintion_in_diffrence_length_both]
    plus_stand_information = [counting_of_plus_elements, counting_pluses_in_difference_length, sequences_of_plus_elements, len(sequence),compposintion_in_diffrence_length_plus]
    comp_strand_information = [counting_of_comp_elements,  counting_complements_in_difference_length, sequences_of_complement_elements, len(sequence),compposintion_in_diffrence_length_comp]
    del sequences_of_plus_elements
    del sequences_of_complement_elements
    del sequences_of_both_elements
    del counting_of_plus_elements
    del counting_of_comp_elements
    del counting_of_both_elements
    del counting_pluses_in_difference_length
    del counting_complements_in_difference_length
    del counting_both_in_difference_length
    del sequence
    del sequence_of_element_p
    del sequence_of_element_c

    gff.close()
    return both_stand_information,plus_stand_information, comp_strand_information, #length_of_overlaping


def segments_to_composition(self,both_element_inf_4,plus_element_inf_4,comp_element_inf_4):
    Ab = [k.count('a') for k in both_element_inf_4]
    Tb = [k.count('t') for k in both_element_inf_4]
    Cb = [k.count('c') for k in both_element_inf_4]
    Gb = [k.count('g') for k in both_element_inf_4]
    B_group = [Ab,Tb,Cb,Gb]

    Ap = [k.count('a') for k in plus_element_inf_4]
    Tp = [k.count('t') for k in plus_element_inf_4]
    Cp = [k.count('c') for k in plus_element_inf_4]
    Gp = [k.count('g') for k in plus_element_inf_4]
    P_group = [Ap,Tp,Cp,Gp]
    Ac = [h.count('a') for h in comp_element_inf_4]
    Tc = [h.count('t') for h in comp_element_inf_4]
    Cc = [h.count('c') for h in comp_element_inf_4]
    Gc = [h.count('g') for h in comp_element_inf_4]
    C_group = [Ac,Tc,Cc,Gc]
    return B_group,P_group,C_group


def scatering_plot_for_frequency_in_differrence_length(self,Both_frequencey,plus_frequencey,comp_frequencey):

    X1 = np.arange(len(plus_frequencey))
    X2 = np.arange(len(comp_frequencey))

    fig,axs = plt.subplots(2)
    axs[0].plot(X1,Both_frequencey, color = 'black')
    axs[1].plot(X1,plus_frequencey,color= 'green')
    axs[1].plot(X2,comp_frequencey,color = 'red')
    plt.show()


def scatering_plot_for_Nucleotid_composition_in_differrence_length(self,B_group,P_group,C_group):

    B_ATCG = [sum(i) for i in zip(B_group[0],B_group[1],B_group[2],B_group[3])]
    Ab = [round(v/(B_ATCG[i]+0.001)*100,2) for i,v in enumerate(B_group[0],start=0)]
    Tb = [round(v/(B_ATCG[i]+0.001)*100,2) for i,v in enumerate(B_group[1],start=0)]
    Cb = [round(v/(B_ATCG[i]+0.001)*100,2) for i,v in enumerate(B_group[2],start=0)]
    Gb = [round(v/(B_ATCG[i]+0.001)*100,2) for i,v in enumerate(B_group[3],start=0)]


    P_ATCG = [sum(i) for i in zip(P_group[0],P_group[1],P_group[2],P_group[3])]
    Ap = [round(v/(P_ATCG[i]+0.001)*100,2) for i,v in enumerate(P_group[0],start=0)]
    Tp = [round(v/(P_ATCG[i]+0.001)*100,2) for i,v in enumerate(P_group[1],start=0)]
    Cp = [round(v/(P_ATCG[i]+0.001)*100,2) for i,v in enumerate(P_group[2],start=0)]
    Gp = [round(v/(P_ATCG[i]+0.001)*100,2) for i,v in enumerate(P_group[3],start=0)]

    C_ATCG = [sum(i) for i in zip(C_group[0],C_group[1],C_group[2],C_group[3])]
    Ac = [round(v/(C_ATCG[i]+0.001)*100,2) for i,v in enumerate(C_group[0],start=0)]
    Tc = [round(v/(C_ATCG[i]+0.001)*100,2) for i,v in enumerate(C_group[1],start=0)]
    Cc = [round(v/(C_ATCG[i]+0.001)*100,2) for i,v in enumerate(C_group[2],start=0)]
    Gc = [round(v/(C_ATCG[i]+0.001)*100,2) for i,v in enumerate(C_group[3],start=0)]
    Ap_Ac = [(Ap[i]+ Ac[i])/2 for i in range(len(Ap))]
    Tp_Tc = [(Tp[i]+ Tc[i])/2 for i in range(len(Ap))]
    Cp_Cc = [(Cp[i]+ Cc[i])/2 for i in range(len(Ap))]
    Gp_Gc = [(Gp[i]+ Gc[i])/2 for i in range(len(Ap))]
    import matplotlib.patches as mpatches

    plt.figure(figsize=(13, 4))
    X1 = np.arange(len(Ap))
    plt.scatter(X1,Ab,color= 'green',s = 10)
    plt.scatter(X1,Tb,color= 'blue',s = 10)
    plt.scatter(X1,Cb,color= 'yellow',s = 10)
    plt.scatter(X1,Gb,color= 'red',s = 10)
    plt.title('number 3 plot')
    plt.show()


    fig, axs = plt.subplots(3,figsize=(13, 4))
    X1 = np.arange(len(Ap))
    axs[0].scatter(X1,Ap,color= 'green',s = 10)
    axs[0].scatter(X1,Tp,color= 'blue',s = 10)
    axs[0].scatter(X1,Cp,color= 'lime',s = 10)
    axs[0].scatter(X1,Gp,color= 'cyan',s = 10)
    #axs[0].title('plot number 4_P')
    axs[1].scatter(X1,Ac,color= 'yellow',s = 10)
    axs[1].scatter(X1,Tc,color= 'red',s = 10)
    axs[1].scatter(X1,Cc,color= 'deeppink',s = 10)
    axs[1].scatter(X1,Gc,color= 'navy',s = 10)

    axs[2].scatter(X1,Ap,color= 'green',s = 10)
    axs[2].scatter(X1,Tp,color= 'blue',s = 10)
    axs[2].scatter(X1,Cp,color= 'lime',s = 10)
    axs[2].scatter(X1,Gp,color= 'cyan',s = 10)
    #axs[0].title('plot number 4_P')
    axs[2].scatter(X1,Ac,color= 'yellow',s = 10)
    axs[2].scatter(X1,Tc,color= 'red',s = 10)
    axs[2].scatter(X1,Cc,color= 'deeppink',s = 10)
    axs[2].scatter(X1,Gc,color= 'navy',s = 10)



    #axs[1].title('plot number 4_C')
    #plt.xlabel('Elements length')
    #plt.ylabel('Nucleotide composition %')

    plt.show()

def calculat_portion_of_elements(sequence, gff,chnumber,element):

    sequence = sequence.upper()
    sequence_plus = sequence
    sequence_plus = list(sequence_plus)
    sequence_comp = sequence
    sequence_comp = list(sequence_comp)
    sequence_both = sequence
    sequence_both = list(sequence_both)
    sequence_lower = sequence.lower()


    chromosome_number = chnumber
    element_name = element
    if chnumber == 23:
        chromosome_number = 'X'
    if chnumber == 24:
        chromosome_number = 'Y'
    gff = open(gff)
    gff.readline()
    gff.readline()
    gff.readline()
    gff.readline()
    gff.readline()
    gff.readline()
    gff.readline()
    gff.readline()
    gff.readline()
    gff.readline()

    for annotation in gff:
        if annotation[0] == '#':
            continue
        annotation.strip()
        annotation = annotation.split()
        if  annotation[3] == '1':
            #print(annotation)
            continue

        if (element_name in annotation[2]) and f'NC_0000{chromosome_number}' in annotation[0] : # or annotation[4] == f'chr{chromosome_number}_'):
            print(annotation)
            start = int(annotation[3])
            end = int(annotation[4])
            if end - start == 0:
                print('len was 0')
                continue


            if annotation[6] == '+':
                sequence_plus[(start - 1):end] = sequence_lower[(start - 1):end]
                sequence_both[(start - 1):end] = sequence_lower[(start - 1):end]

            else:
                sequence_comp[start - 1:end] = sequence_lower[start - 1:end]
                sequence_both[(start - 1):end] = sequence_lower[(start - 1):end]
    sequence_plus = ''.join(sequence_plus)
    sequence_comp = ''.join(sequence_comp)
    sequence_both = ''.join(sequence_both)

    plus_portion = (sequence_plus.count('a') + sequence_plus.count('t') + sequence_plus.count('c') + sequence_plus.count('g'))/len(sequence_plus) * 100
    minis_portion = (sequence_comp.count('a') + sequence_comp.count('t') + sequence_comp.count('c') + sequence_comp.count('g'))/len(sequence_comp)*100
    both_portion = (sequence_both.count('a')+sequence_both.count('t')+sequence_both.count('c')+sequence_both.count('g'))/len(sequence_both) * 100
    print('Plus portion',plus_portion)
    print('Negative portion',minis_portion)
    print('Both portion',both_portion)
    print('p portion',len(sequence_plus)/len(sequence))
    return sequence_plus,sequence_comp,sequence_both

def calculat_portion_of_elements_overlap_sharing(sequence, gff_path,chnumber,element):

    sequence = sequence.upper()#kdnn
    sequence_number_base = sequence.replace('A','1')
    sequence_number_base = sequence_number_base.replace('T','2')
    sequence_number_base = sequence_number_base.replace('C','3')
    sequence_number_base = sequence_number_base.replace('G','4')
    sequence_number_base = list(sequence_number_base)

    sequence_plus = sequence.upper()
    sequence_plus = list(sequence_plus)






    chromosome_number = chnumber
    element_name = element
    if chnumber == 23:
        chromosome_number = 'X'
    if chnumber == 24:
        chromosome_number = 'Y'
    gff = open(gff_path)
    gff.readline()
    gff.readline()
    gff.readline()
    gff.readline()
    gff.readline()
    gff.readline()
    gff.readline()
    gff.readline()
    gff.readline()
    gff.readline()
    #step 1 for plus
    for annotation in gff:
        if annotation[0] == '#':
            continue
        annotation.strip()
        annotation = annotation.split()
        if  annotation[3] == '1':
            #print(annotation)
            continue

        if (element_name in annotation[2]) and f'NC_0000{chromosome_number}' in annotation[0] : # or annotation[4] == f'chr{chromosome_number}_'):
            print(annotation)
            start = int(annotation[3])
            end = int(annotation[4])
            if end - start == 0:
                print('len was 0')
                continue

            if annotation[6] == '-':
                sequence_number_base[(start - 1):end] = sequence_plus[(start - 1):end]



    #step 2 for minis
    gff.close()
    gff = open(gff_path)
    gff.readline()
    gff.readline()
    gff.readline()
    gff.readline()
    gff.readline()
    gff.readline()
    gff.readline()
    gff.readline()
    gff.readline()
    gff.readline() 



    for annotation in gff:
        if annotation[0] == '#':
            continue
        annotation.strip()
        annotation = annotation.split()
        if  annotation[3] == '1':
            #print(annotation)
            continue
        if (element_name in annotation[2]) and f'NC_0000{chromosome_number}' in annotation[0] : # or annotation[4] == f'chr{chromosome_number}_'):
            print(annotation)
            start = int(annotation[3])
            end = int(annotation[4])
            if end - start == 0:
                print('len was 0')
                continue
            if annotation[6] == '+':
                E = sequence_number_base[start - 1:end]
                E = ''.join(E)
                E = E.replace('1','a')
                E = E.replace('2','t')
                E = E.replace('3','c')
                E = E.replace('4','g')
                E = E.replace('A','~')
                E = E.replace('T','!')
                E = E.replace('C','@')
                E = E.replace('G','#')
                E = list(E)
                sequence_number_base[start - 1:end] = E


    sequence_number_base = ''.join(sequence_number_base)

    N = sequence_number_base.count('N')
    print('N',N)
    nonmap_portion = (sequence_number_base.count('1') + sequence_number_base.count('2') + sequence_number_base.count('3') + sequence_number_base.count('4'))/len(sequence_number_base) * 100
    plus_portion = (sequence_number_base.count('A') + sequence_number_base.count('T') + sequence_number_base.count('C') + sequence_number_base.count('G'))/len(sequence_number_base) * 100
    minis_portion = (sequence_number_base.count('a') + sequence_number_base.count('t') + sequence_number_base.count('c') + sequence_number_base.count('g'))/len(sequence_number_base)*100
    share_portion = (sequence_number_base.count('~') + sequence_number_base.count('!') + sequence_number_base.count('@') + sequence_number_base.count('#'))/len(sequence_number_base) * 100
    print('Plus portion',plus_portion)
    print('Negative portion',minis_portion)
    print('share portion',share_portion)
    print('p portion',len(sequence_number_base)/len(sequence))
    print('nonmap',nonmap_portion)
    #print(sequence_number_base)
    return sequence_number_base

def seq_Nucleotide_composition(seq):
    N =  seq.count('N')
    n =  seq.count('n')
    a =  round(seq.count('a')/(len(seq)-N)*100,2)
    t =  round(seq.count('t')/(len(seq)-N)*100,2)
    c =  round(seq.count('c')/(len(seq)-N)*100,2)
    g =  round(seq.count('g')/(len(seq)-N)*100,2)
    A =  round(seq.count('A')/(len(seq)-N)*100,2)
    T =  round(seq.count('T')/(len(seq)-N)*100,2)
    C =  round(seq.count('C')/(len(seq)-N)*100,2)
    G =  round(seq.count('G')/(len(seq)-N)*100,2)


    print('A',A,'T',T,'C',C,'G',G,'N',N)
    print('a',a,'t',t,'c',c,'g',g,'n',n)
    print('all nucleotide',a+t+c+g+A+T+C+G)
    #return A,T,C,G,a,t,c,g


def target_sequences_lowering(seq,start,end):
    start = int(start)
    end = int(end)

    element = list(''.join(seq[(start - 1):end]).lower())

    seq[(start - 1):end] = element
    return seq



def nucleotide_composition_in_seqenence_ATCG_atcg(sequence,plot = True):
    seq_length = len(sequence)
    print('len seq is ',seq_length)
    A = round((sequence.count('A') + sequence.count('a'))/seq_length*100,3)
    T = round((sequence.count('T') + sequence.count('t'))/seq_length*100,3)
    C = round((sequence.count('C') + sequence.count('c'))/seq_length*100,3)
    G = round((sequence.count('G') + sequence.count('g'))/seq_length*100,3)
    N = round((sequence.count('N') + sequence.count('n'))/seq_length*100,3)
    chech_portion = (A + T + C + G + N)
    print('cheked portion',(sequence.count('A') + sequence.count('a')+sequence.count('T') + sequence.count('t')+
                            sequence.count('C') + sequence.count('c')+sequence.count('G') + sequence.count('g')+
                            sequence.count('N') + sequence.count('n'))/seq_length*100)
    print('Aa',A,'%')
    print('Tt',T,'%')
    print('Cc',C,'%')
    print('Gg',G,'%')
    print('Nn',N,'%')
    print('\nnucleotide_composition_in_seqenence Pass!!!!')


    # Sample data
    labels = ['A', 'T', 'C', 'G']
    sizes = [A,T,C,G]  # percentages, should add up to 100%

    # Colors for each category
    colors = ['green', 'blue', 'orange', 'red']

    # Plotting the pie chart
    plt.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)

    # Draw a circle at the center of the pie to make it look like a donut chart
    centre_circle = plt.Circle((0,0),0.70,fc='white')
    fig = plt.gcf()
    fig.gca().add_artist(centre_circle)

    # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.axis('equal')

    # Display the chart
    plt.show()




def counting_of_TE_elements_refseq(sequence, out_mask,chromosome_number,elements, plot_length):
    counting_of_plus_elements = 0
    counting_of_comp_elements = 0
    sequences_of_plus_elements = ''
    sequences_of_complement_elements = ''
    #elements_overlap = list(sequence)
    counting_pluses_in_difference_length = [0] *plot_length
    counting_complements_in_difference_length = [0] * plot_length
    compposintion_in_diffrence_length_plus = [''] * plot_length
    compposintion_in_diffrence_length_comp = [''] * plot_length
    # element_location_d = [0] * len(sequence)
    # Alus_location_c = [0] * len(sequence)
    out_mask.readline()
    out_mask.readline()
    out_mask.readline()

    #sequence2_for_overlabing = sequence
    flaq1 = False
    #'flaq is on for chromosome lins ended'
    for annotation in out_mask:
        annotation = annotation.split()
        if  chromosome_number in annotation[4] :
            flaq1 = True
            for i in elements:
                if i in annotation[10]:
                    if annotation[8] == '+':
                        sequence_of_element_p = sequence[int(annotation[5]) - 1:int(annotation[6]) - 1]
                        print(sequence_of_element_p)
                        sequences_of_plus_elements += sequence_of_element_p
                        #elements_overlap[int(annotation[5]) - 1:int(annotation[6]) - 1] = sequence_of_element_p.lower()
                        #element_location_d[int(line[5])] += 10
                        len_element_p = len(sequence_of_element_p)
                    
                    
                        if len_element_p == 0 :
                            #print(annotation)
                            print('start',int(annotation[5]))
                            print('len',len(sequence))
                    
                    
                        if len_element_p > plot_length:
                            print('Not studied Element it is very big',len_element_p)
                        if len_element_p < plot_length:
                            counting_of_plus_elements += 1
                            counting_pluses_in_difference_length[len_element_p] += 1
                            compposintion_in_diffrence_length_plus[len_element_p] += sequence_of_element_p


                    elif 'C' in annotation:
                        sequence_of_element_c = sequence[int(annotation[5]) - 1:int(annotation[6]) - 1]
                        sequences_of_complement_elements += sequence_of_element_c
                        print(sequence_of_element_c)
                        #elements_overlap[int(annotation[5]) - 1:int(annotation[6]) - 1] = sequence_of_element_c.lower()
                        #Alus_location_c[int(line[5])] += 10
                        len_element_c = len(sequence_of_element_c)
                        if len_element_c > plot_length:
                            print('Not studied Element it is very big',len_element_c)
                        if len_element_c < plot_length:
                            counting_of_comp_elements += 1
                            counting_complements_in_difference_length[len_element_c] += 1
                            compposintion_in_diffrence_length_comp[len_element_c]+=sequence_of_element_c
                    
    print('element counting finished')
    print('  plus    comp')
    print(counting_of_plus_elements,'  ',counting_of_comp_elements)
    plus_elements_portion_in_sequence = round(len(sequences_of_plus_elements)/len(sequence)*100,2)
    complement_elements_portion_in_sequence = round(len(sequences_of_complement_elements)/len(sequence)*100,2)
    print(' portion of elements sequence in + and - strand:\n',plus_elements_portion_in_sequence,'%',complement_elements_portion_in_sequence,'%')
    len_sum_of_pluses = len(sequences_of_plus_elements)
    A_d = round(sequences_of_plus_elements.count('a') / len_sum_of_pluses * 100, 2)
    T_d = round(sequences_of_plus_elements.count('t') / len_sum_of_pluses * 100, 2)
    C_d = round(sequences_of_plus_elements.count('c') / len_sum_of_pluses * 100, 2)
    G_d = round(sequences_of_plus_elements.count('g') / len_sum_of_pluses * 100, 2)
    len_sum_of_complements = len(sequences_of_complement_elements)
    A_c = round(sequences_of_complement_elements.count('a') / len_sum_of_complements * 100, 2)
    T_c = round(sequences_of_complement_elements.count('t') / len_sum_of_complements * 100, 2)
    C_c = round(sequences_of_complement_elements.count('c') / len_sum_of_complements * 100, 2)
    G_c = round(sequences_of_complement_elements.count('g') / len_sum_of_complements * 100, 2)
    print('+composition', '-composition')
    print('A ', A_d, '   ', ' A ', A_c)
    print('T ', T_d, '   ', ' T ', T_c)
    print('C ', C_d, '   ', ' C ', C_c)
    print('G ', G_d, '   ', ' G ', G_c)
    import statistics



    P =   len(sequences_of_plus_elements)/counting_of_plus_elements
    C =   len(sequences_of_complement_elements)/counting_of_comp_elements


    plus_stand_information = [counting_of_plus_elements, counting_pluses_in_difference_length, sequences_of_plus_elements, len(sequence),compposintion_in_diffrence_length_plus,P,len(sequences_of_plus_elements)]
    comp_strand_information = [counting_of_comp_elements,  counting_complements_in_difference_length, sequences_of_complement_elements, len(sequence),compposintion_in_diffrence_length_comp,C,len(sequences_of_complement_elements)]


    print('plus elements length avarage',P,'bp')
    print('compliment elements length avarage',C,'bp')
    print('lensP',len(sequences_of_plus_elements),'lenC',len(sequences_of_complement_elements))

    del sequences_of_plus_elements
    del sequences_of_complement_elements
    del counting_of_plus_elements
    del counting_of_comp_elements
    del counting_pluses_in_difference_length
    del counting_complements_in_difference_length
    del sequence
    del sequence_of_element_p
    del sequence_of_element_c
    del P
    del C

    overlap_information = []
    out_mask.close()
    return plus_stand_information, comp_strand_information, #length_of_overlaping













