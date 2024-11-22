import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import openpyxl
import gzip
from Bio import SeqIO

def open_simple_fasta_file_and_prepare_it(seq_path):

    fasta_file = open(seq_path)
    fasta_file.readline()
    sequence = fasta_file.read()
    sequence = sequence.replace('\n', '')
    fasta_file.close()

    len_chromosome = len(sequence)
    return sequence



def counting_of_elements( sequence,gff_path,element_names, plot_length):

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





    #sequence2_for_overlabing = sequence
    countery  = 0
    gff = open(gff_path)
    for annotation in gff:

        if annotation[0] == '#':
            continue

        annotation = annotation.strip().split()

        if  annotation[3] == '1':
            #print(annotation)
            continue

        for element in element_names:
            repeat_control = []
            if element in annotation[2]  : # or annotation[4] == f'chr{chromosome_number}_'):

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
    A_b = round(sequences_of_both_elements.count('A') / len_both_seq * 100, 2)
    T_b = round(sequences_of_both_elements.count('T') / len_both_seq * 100, 2)
    C_b = round(sequences_of_both_elements.count('C') / len_both_seq * 100, 2)
    G_b = round(sequences_of_both_elements.count('G') / len_both_seq * 100, 2)
    len_sum_of_pluses = len(sequences_of_plus_elements)
    A_d = round(sequences_of_plus_elements.count('A') / len_sum_of_pluses * 100, 2)
    T_d = round(sequences_of_plus_elements.count('T') / len_sum_of_pluses * 100, 2)
    C_d = round(sequences_of_plus_elements.count('C') / len_sum_of_pluses * 100, 2)
    G_d = round(sequences_of_plus_elements.count('G') / len_sum_of_pluses * 100, 2)
    len_sum_of_complements = len(sequences_of_complement_elements)
    A_c = round(sequences_of_complement_elements.count('A') / len_sum_of_complements * 100, 2)
    T_c = round(sequences_of_complement_elements.count('T') / len_sum_of_complements * 100, 2)
    C_c = round(sequences_of_complement_elements.count('C') / len_sum_of_complements * 100, 2)
    G_c = round(sequences_of_complement_elements.count('G') / len_sum_of_complements * 100, 2)
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




def nucleotide_composition_in_seqenence_ATCG(sequence):

    A = sequence.count('A')
    T = sequence.count('T')
    C = sequence.count('C')
    G = sequence.count('G')
    N = sequence.count('N')
    one_masked_count = A + T + C + G + N

    n = sequence.count('n')
    one =  sequence.count('1')
    two =  sequence.count('2')
    three =sequence.count('3')
    four = sequence.count('4')
    Duble_masked_count = one + two + three + four
    a = sequence.count('a')
    t = sequence.count('t')
    c = sequence.count('c')
    g = sequence.count('g')
    unmasked_count = a + t + c + g
    Len_seq = len(sequence)

    unsequenced_region = round((n)/Len_seq*100,2)
    one_masked_portion = round(((A + T + C + G + N )/Len_seq*100),2)
    Duble_masked_portion = round(((one + two + three + four)/Len_seq*100),2)
    unmasked_region = (Len_seq/Len_seq*100) - unsequenced_region - one_masked_portion - Duble_masked_portion

    Ap = round(A/one_masked_count*100,3)
    Tp = round(T/one_masked_count*100,3)
    Cp = round(C/one_masked_count*100,3)
    Gp = round(G/one_masked_count*100,3)
    Np = round(N/one_masked_count*100,3)

    try :
        A1 = round(one/Duble_masked_count*100,2)
    except ZeroDivisionError:
        A1 = 0
    try :
        T2 = round(two/Duble_masked_count*100,2)
    except ZeroDivisionError:
        T2 = 0
    try:
        C3 = round(three/Duble_masked_count*100,2)
    except ZeroDivisionError:
        C3 = 0
    try:
        G4 = round(three/Duble_masked_count*100,2)
    except ZeroDivisionError:
        G4 = 0




    ap = round(a/unmasked_count*100,3)
    tp = round(t/unmasked_count*100,3)
    cp = round(c/unmasked_count*100,3)
    gp = round(g/unmasked_count*100,3)






    #result
    print('len sequence :',Len_seq)
    print('one_masked',one_masked_portion,'%','Duble_masked',Duble_masked_portion,'%','unmasked',unmasked_region,'%'
          'unsequenced',unsequenced_region,'%')
    print('Sum of all',one_masked_portion + Duble_masked_portion + unmasked_region + unsequenced_region)
    print('\none_masked_portion composition')
    print('A',Ap,'%')
    print('T',Tp,'%')
    print('C',Cp,'%')
    print('G',Gp,'%')
    print('N',Np,'%')
    print('sum : ',Ap + Tp + Cp + Gp + Np)
    print('A/T',round(A/T,3))
    print('C/G',round(C/G,3))

    print('Duble masked composition')
    print('A1',A1)
    print('T2',T2)
    print('C3',C3)
    print('G4',G4)
    print('sum : ',A1+T2+C3+G4)

    print('unmasked composition')
    print('a',ap)
    print('t',tp)
    print('c',cp)
    print('g',gp)
    print('sum :',ap+tp+cp+gp)


    print('\nnucleotide_composition_in_seqenence Pass!!!!')
    # Sample data
    labels = ['one_masked_portion', 'Duble_masked_portion', 'unmasked_region', 'unsequenced_region']
    sizes = [one_masked_portion,Duble_masked_portion,unmasked_region,unsequenced_region]  # percentages, should add up to 100%

    explode = (0.1, 0.1, 0.1, 0.1)

    #Plotting the pie chart
    plt.pie(sizes,labels = labels, explode=explode,  autopct='%1.1f%%', shadow=True, startangle=140)

    # Equal aspect ratio ensures that pie is drawn as a circle
    plt.axis('equal')

    # Adding a title
    plt.title('Example Pie Chart')

    # Display the plot
    plt.show()





def plus_minis_both_1(sequence0, gff,elements):
    sequence = sequence0.lower()  # ATCG all
    sequence_plus = sequence
    sequence_plus = list(sequence_plus) #ATCG all
    sequence_comp = sequence
    sequence_comp = list(sequence_comp) #ATCG all
    sequence_both = sequence
    sequence_both = list(sequence_both) #ATCG all
    sequence_lower = sequence.upper() # this is for masking
    gff = open(gff)
    gff.readline()
    gff.readline()
    gff.readline()
    element_names_in_line = []

    counter2 = 0
    for annotation in gff:

        if annotation[0] == '#':
            continue
        annotation.strip()
        annotation = annotation.split()

        for element in elements:
            if element in annotation[2] :

                counter2 += 1
                start = int(annotation[3])
                end = int(annotation[4])
                if end - start == 0:
                    print('len was 0')
                    continue

                if annotation[6] == '+':
                    sequence_plus[(start - 1):end] = sequence_lower[(start - 1):end]
                    sequence_both[(start - 1):end] = sequence_lower[(start - 1):end]
                    #print('sample',len(sequence0[(start - 1):end]),sequence0[(start - 1):end])

                if annotation[6] == '-':
                    sequence_comp[start - 1:end] = sequence_lower[start - 1:end]
                    sequence_both[(start - 1):end] = sequence_lower[(start - 1):end]

    sequence_plus = ''.join(sequence_plus)
    sequence_comp = ''.join(sequence_comp)
    sequence_both = ''.join(sequence_both)
    un_sequenced_length = sequence.count('N')  + sequence.count('n')
    sequenced_length = len(sequence) - un_sequenced_length
    unmasked_portion = (sequence_both.count('a') + sequence_both.count('t') + sequence_both.count('c') + sequence_both.count('g'))/sequenced_length * 100
    plus_portion = (sequence_plus.count('A') + sequence_plus.count('T') + sequence_plus.count('C') + sequence_plus.count('G'))/sequenced_length * 100
    minis_portion = (sequence_comp.count('A') + sequence_comp.count('T') + sequence_comp.count('C') + sequence_comp.count('G'))/sequenced_length*100
    both_portion = (sequence_both.count('A')+sequence_both.count('T')+sequence_both.count('C')+sequence_both.count('G'))/sequenced_length * 100
    print('counter',counter2)
    print('sequenced portion',round(sequenced_length/len(sequence)*100,2))
    print('Plus portion',plus_portion)
    print('Negative portion',minis_portion)
    print('Both Plus and Negative portion',both_portion)
    print('unmasked portion',unmasked_portion)
    print('calculat_portion_of_elements Pass!!')

    return sequence_plus,sequence_comp,sequence_both
















