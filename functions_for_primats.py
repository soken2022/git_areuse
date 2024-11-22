import gzip
from Bio import SeqIO,Entrez
import matplotlib.pyplot as plt

def open_fasta_file_multi_record_UCSC_and_refseq(seq_path,ch_number0):
    sequence = ''
    fasta_path = seq_path


    if ch_number0[0] == 23 :
        ch_number0[0] = 'X'
    if ch_number0[0] == 24 :
        ch_number0[0] = 'Y'
    if fasta_path[-2:] == 'gz':
        fasta_file = gzip.open(fasta_path,'rt')
    else:
        fasta_file = fasta_path
    print('for chromosoeme',ch_number0)
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        if seq_record.id == rf'chr{ch_number0[0]}' or ch_number0[1] in (seq_record.id) :
            print('sequence id is',seq_record.id)
            sequence = str(seq_record.seq)
            if fasta_path[-2:] == 'gz':
                fasta_file.close()
            l  = len(sequence)
            print('len chromosome is',l,'bp')

            break

    return sequence

def portion_of_TE_elements_in_geneme_with_UCSC_and_refseq(sequence0, gff,chnumber,element):
    sequence = sequence0.lower()  # ATCG all
    sequence_plus = sequence
    sequence_plus = list(sequence_plus) #ATCG all
    sequence_comp = sequence
    sequence_comp = list(sequence_comp) #ATCG all
    sequence_both = sequence
    sequence_both = list(sequence_both) #ATCG all
    sequence_lower = sequence.upper() # this is for masking


    chromosome_number = chnumber

    gff = open(gff)
    gff.readline()
    gff.readline()
    gff.readline()
    element_names_in_line = []
    counter = 0
    counter2 = 0
    for annotation in gff:

        if annotation[0] == '#':
            continue
        annotation.strip()
        annotation = annotation.split()


        '''if not annotation[10] in element_names_in_line:
            element_names_in_line.append(annotation[10])
            print(annotation[10]) '''
        '''flaq = False
        for i in element:
            if f'chr{chromosome_number}' == annotation[4]  and   i  in annotation[10]:
                flaq = True '''

        if (f'chr{chromosome_number[0]}' == annotation[4] or chromosome_number[1] in annotation[4]) and ('LINE' in annotation[10] or 'SINE' in annotation[10]
        or 'LTR' in annotation[10] or 'DNA' in annotation[10] or 'RC' in annotation[10]
        or 'tRNA' in annotation[10] or 'scRNA' in annotation[10] or 'srpRNA' in annotation[10]
        or 'RNA' in annotation[10] or 'Retroposons' in annotation[10]): #'Satellite' not in annotation[10] or 'Low_complexity' not in annotation[10]or 'rRNA' not in annotation[10] or 'Simple_repeat' not in annotation[10]):

        #if flaq == True  :
            #if f'chr{chromosome_number}' == annotation[4] and  annotation[10].split('/')[0] not in element:
            #if f'chr{chromosome_number}' == annotation[4]:
            counter2 += 1
            start = int(annotation[5])
            end = int(annotation[6])
            if end - start == 0:
                print('len was 0')
                continue

            if annotation[8] == '+':
                sequence_plus[(start - 1):end] = sequence_lower[(start - 1):end]
                sequence_both[(start - 1):end] = sequence_lower[(start - 1):end]
                #print('sample',len(sequence0[(start - 1):end]),sequence0[(start - 1):end])
                '''if sequence0[(start - 1):end].count('A') > 0 or sequence0[(start - 1):end].count('T') > 0\
                        or sequence0[(start - 1):end].count('C') > 0 or sequence0[(start - 1):end].count('G') > 0:

                    print('Hi!!!!!!!!!!!')
                else:
                    print('Hello')
                    counter +=1
                    print(counter,len(sequence0[(start - 1):end]),sequence0[(start - 1):end])'''

            if annotation[8] == 'C':
                sequence_comp[start - 1:end] = sequence_lower[start - 1:end]
                sequence_both[(start - 1):end] = sequence_lower[(start - 1):end]
                '''if sequence0[(start - 1):end].count('A') > 0 or sequence0[(start - 1):end].count('T') > 0 \
                        or sequence0[(start - 1):end].count('C') > 0 or sequence0[(start - 1):end].count('G') > 0:

                    print('Hi!!!!!!!!!!!')
                else:
                    print('Hello')
                    counter +=1
                    print(counter,len(sequence0[(start - 1):end]),sequence0[(start - 1):end])'''
    sequence_plus = ''.join(sequence_plus)
    sequence_comp = ''.join(sequence_comp)
    sequence_both = ''.join(sequence_both)
    un_sequenced_length = sequence.count('N')  + sequence.count('n')
    sequenced_length = len(sequence) - un_sequenced_length
    unmasked_portion = (sequence_both.count('a') + sequence_both.count('t') + sequence_both.count('c') + sequence_both.count('g'))/sequenced_length * 100
    plus_portion = (sequence_plus.count('A') + sequence_plus.count('T') + sequence_plus.count('C') + sequence_plus.count('G'))/sequenced_length * 100
    minis_portion = (sequence_comp.count('A') + sequence_comp.count('T') + sequence_comp.count('C') + sequence_comp.count('G'))/sequenced_length*100
    both_portion = (sequence_both.count('A') + sequence_both.count('T') + sequence_both.count('C')+sequence_both.count('G'))/sequenced_length * 100
    print('counter',counter2)
    print('sequenced portion',round(sequenced_length/len(sequence)*100,2))
    print('Plus portion',plus_portion)
    print('Negative portion',minis_portion)
    print('Both Plus and Negative portion',both_portion)
    print('unmasked portion',unmasked_portion)
    print('calculat_portion_of_elements Pass!!')
    return sequence_both,sequence_comp,sequence_both

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

    A1 = round(one/Duble_masked_count*100,2)
    T2 = round(two/Duble_masked_count*100,2)
    C3 = round(three/Duble_masked_count*100,2)
    G4 = round(three/Duble_masked_count*100,2)


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

def nucleotide_composition_in_seqenence_ATCG_2(sequence):

    A1 = sequence.count('~')
    T1 = sequence.count('!')
    C1 = sequence.count('@')
    G1 = sequence.count('#')
    N1 = sequence.count('N')
    one_masked_count_type1 = A1 + T1 + C1 + G1 + N1

    A2 = sequence.count('A')
    T2 = sequence.count('T')
    C2 = sequence.count('C')
    G2 = sequence.count('G')
    N2 = sequence.count('N')
    one_masked_count_type2 = A2 + T2 + C2 + G2 + N2

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
    one_masked_portion_type1 = round(((A1 + T1 + C1 + G1 + N1 )/Len_seq*100),2)
    one_masked_portion_type2 = round(((A2 + T2 + C2 + G2 + N2 )/Len_seq*100),2)
    Duble_masked_portion = round(((one + two + three + four)/Len_seq*100),2)
    unmasked_region = round(((Len_seq/Len_seq*100) - unsequenced_region -
                       one_masked_portion_type1 - one_masked_portion_type2 - Duble_masked_portion) ,2)

    A1p = round(A1/one_masked_count_type1*100,3)
    T1p = round(T1/one_masked_count_type1*100,3)
    C1p = round(C1/one_masked_count_type1*100,3)
    G1p = round(G1/one_masked_count_type1*100,3)
    N1p = round(N1/one_masked_count_type1*100,3)

    A2p = round(A2/one_masked_count_type2*100,3)
    T2p = round(T2/one_masked_count_type2*100,3)
    C2p = round(C2/one_masked_count_type2*100,3)
    G2p = round(G2/one_masked_count_type2*100,3)
    N2p = round(N2/one_masked_count_type2*100,3)

    An1 = round(one/Duble_masked_count*100,2)
    Tn2 = round(two/Duble_masked_count*100,2)
    Cn3 = round(three/Duble_masked_count*100,2)
    Gn4 = round(three/Duble_masked_count*100,2)


    ap = round(a/unmasked_count*100,3)
    tp = round(t/unmasked_count*100,3)
    cp = round(c/unmasked_count*100,3)
    gp = round(g/unmasked_count*100,3)






    #result
    print('len sequence :',Len_seq)
    print('one_masked_type1',one_masked_portion_type1,'one_masked_type2',one_masked_portion_type2,'%','Duble_masked',Duble_masked_portion,'%','unmasked',unmasked_region,'%'
                                                                                                                 'unsequenced',unsequenced_region,'%')
    print('Sum of all',one_masked_portion_type1+ one_masked_portion_type2 + Duble_masked_portion + unmasked_region + unsequenced_region)
    print('\none_masked_portion_type1 composition')
    print('A',A1p,'%')
    print('T',T1p,'%')
    print('C',C1p,'%')
    print('G',G1p,'%')
    print('N',N1p,'%')
    print('sum : ',A1p + T1p + C1p + G1p + N1p)
    print('A/T',round(A1/T1,3))
    print('C/G',round(C1/G1,3))

    print('\none_masked_portion_type2 composition')
    print('A',A2p,'%')
    print('T',T2p,'%')
    print('C',C2p,'%')
    print('G',G2p,'%')
    print('N',N2p,'%')
    print('sum : ',A2p + T2p + C2p + G2p + N2p)
    print('A/T',round(A2/T2,3))
    print('C/G',round(C2/G2,3))



    print('Duble masked composition')
    print('A1',An1)
    print('T2',Tn2)
    print('C3',Cn3)
    print('G4',Gn4)
    print('sum : ',An1+Tn2+Cn3+Gn4)
    print('A/T',round(An1/Tn2,3))
    print('C/G',round(Cn3/Gn4,3))

    print('unmasked composition')
    print('a',ap)
    print('t',tp)
    print('c',cp)
    print('g',gp)
    print('sum :',ap+tp+cp+gp)
    print('A/T',round(ap/tp,3))
    print('C/G',round(cp/gp,3))

    print('\nnucleotide_composition_in_seqenence Pass!!!!')
    # Sample data
    labels = ['Genes', 'TEs','Shere-Genes-TEs', 'unmasked_region', 'un_sequenced']
    sizes = [one_masked_portion_type1,one_masked_portion_type2,Duble_masked_portion,unmasked_region,unsequenced_region]  # percentages, should add up to 100%
    print('sizes',sizes)
    explode = (0.1, 0.1, 0.1, 0.1,0.1)

    #Plotting the pie chart
    plt.pie(sizes,labels = labels, explode=explode,  autopct='%1.1f%%', shadow=True, startangle=140)

    # Equal aspect ratio ensures that pie is drawn as a circle
    plt.axis('equal')

    # Adding a title
    plt.title('Es composition')

    # Display the plot
    plt.show()



def counting_of_elements_with_lowering_in_begin(sequence0, gff, ch_id,element_name0,plot_length):
    #ref sequence
    ch_sequence = sequence0.lower()

    #plus sequence
    sequence_plus_mark_1234 = ch_sequence.replace('a','1')
    sequence_plus_mark_1234 = sequence_plus_mark_1234.replace('t','2')
    sequence_plus_mark_1234 = sequence_plus_mark_1234.replace('c','3')
    sequence_plus_mark_1234 = sequence_plus_mark_1234.replace('g','4')
    sequence_plus_mark_1234 = list(sequence_plus_mark_1234)

    #minis sequence
    sequence_plus_mark_5678 = ch_sequence.replace('a','5')
    sequence_plus_mark_5678 = sequence_plus_mark_5678.replace('t','6')
    sequence_plus_mark_5678 = sequence_plus_mark_5678.replace('c','7')
    sequence_plus_mark_5678 = sequence_plus_mark_5678.replace('g','8')
    sequence_plus_mark_5678 = list(sequence_plus_mark_5678)

    ch_sequence = list(ch_sequence)

    print('start function : counting_of_elements function')
    print('chromosome id is',ch_id)
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

    chromosome_number = ch_id
    element_name = element_name0
    gff = open(gff)
    for annotation in gff:

        if annotation[0] == '#':
            continue
        annotation.strip()
        annotation_columns = annotation.split()
        #print(annotation[0])

        if  annotation_columns[3] == '1':
            #print(annotation)
            continue
        if annotation_columns[3] == 'pseudogene' or not annotation_columns[3].isdigit():
            continue

        repeat_control = []
        if (annotation_columns[2] in element_name) and chromosome_number[0] == annotation_columns[0] or chromosome_number[1] in annotation_columns[0] : # or annotation[4] == f'chr{chromosome_number}_'):
            if int(annotation_columns[4]) - int(annotation_columns[3]) == 0:
                #print('len was 0')
                #print(annotation_columns)
                continue


            if annotation_columns[3] in repeat_control:
                print('founded repeat',annotation_columns[3],annotation_columns[4])
                repeat_control.append(annotation_columns[3])
            '''print('GFF information')
            print(annotation_columns)
            seq = sequence0[int(annotation_columns[3]) - 1:int(annotation_columns[4])]
            print('GFF seq',len(seq),seq)
            #gain gene name

            sequcnce_checker(annotation)
            print('\n\n\n\n\n\n')'''


            #get_gene_sequence()

            if annotation_columns[6] == '+':

                candit_seq_for_p = ch_sequence[int(annotation_columns[3]) - 1:int(annotation_columns[4])]

                if '5' or '6' or '7'or '8' or 'A' or 'T' or 'C' or 'G' in candit_seq_for_p:
                    sequence_of_element_p = ''.join(candit_seq_for_p)
                    sequence_of_element_p = sequence_of_element_p.replace('5','A')
                    sequence_of_element_p = sequence_of_element_p.replace('6','T')
                    sequence_of_element_p = sequence_of_element_p.replace('7','C')
                    sequence_of_element_p = sequence_of_element_p.replace('8','G')

                    sequence_of_element_p = sequence_of_element_p.replace('a','1')
                    sequence_of_element_p = sequence_of_element_p.replace('t','2')
                    sequence_of_element_p = sequence_of_element_p.replace('c','3')
                    sequence_of_element_p = sequence_of_element_p.replace('g','4')
                    sequence_of_element_p = list(sequence_of_element_p)
                else:
                    sequence_of_element_p = sequence_plus_mark_1234[int(annotation_columns[3]) - 1:int(annotation_columns[4])]

                ch_sequence[int(annotation_columns[3]) - 1:int(annotation_columns[4])] = sequence_of_element_p
                counting_of_plus_elements += 1
                counting_of_both_elements += 1
                #elements_overlap[int(annotation[5]) - 1:int(annotation[6]) - 1] = sequence_of_element_p.lower()
                #element_location_d[int(line[5])] += 10
                len_element_p = int(annotation_columns[4]) - int(annotation_columns[3])


                if len_element_p >= plot_length:
                    print('Not studied Element it is very big,only dont show in plot',len_element_p)
                if len_element_p < plot_length:
                    counting_pluses_in_difference_length[len_element_p] += 1
                    counting_both_in_difference_length[len_element_p] += 1
                    compposintion_in_diffrence_length_plus[len_element_p] += ''.join(sequence_of_element_p)
                    compposintion_in_diffrence_length_both[len_element_p] += ''.join(sequence_of_element_p)
                    #sequences_of_plus_elements += ''.join(sequence_of_element_p)
                    #sequences_of_both_elements += ''.join(sequence_of_element_p)



            #-----------------------------------------------------------------------

            if annotation_columns[6] == '-':
                candid_for_c = ch_sequence[int(annotation_columns[3]) - 1:int(annotation_columns[4])]
                if '1' or '2' or '3' or '4' or 'A' or 'T' or 'C' or 'G' in  candid_for_c:
                    sequence_of_element_c = ''.join(candid_for_c)
                    sequence_of_element_c = sequence_of_element_c.replace('1','A')
                    sequence_of_element_c = sequence_of_element_c.replace('2','T')
                    sequence_of_element_c = sequence_of_element_c.replace('3','C')
                    sequence_of_element_c = sequence_of_element_c.replace('4','G')

                    sequence_of_element_c = sequence_of_element_c.replace('a','5')
                    sequence_of_element_c = sequence_of_element_c.replace('t','6')
                    sequence_of_element_c = sequence_of_element_c.replace('c','7')
                    sequence_of_element_c = sequence_of_element_c.replace('g','8')

                    sequence_of_element_c = list(sequence_of_element_c)
                else:
                    sequence_of_element_c = sequence_plus_mark_5678[int(annotation_columns[3]) - 1:int(annotation_columns[4])]


                ch_sequence[int(annotation_columns[3]) - 1:int(annotation_columns[4])] = sequence_of_element_c
                counting_of_comp_elements += 1
                counting_of_both_elements += 1
                len_element_c = int(annotation_columns[4]) - int(annotation_columns[3])
                if len_element_c >= plot_length:
                    print('Not studied Element it is very big,only dont show in plot',len_element_c)
                if len_element_c < plot_length:

                    counting_complements_in_difference_length[len_element_c] += 1
                    counting_both_in_difference_length[len_element_c] += 1
                    compposintion_in_diffrence_length_comp[len_element_c] += ''.join(sequence_of_element_c)
                    compposintion_in_diffrence_length_both[len_element_c] += ''.join(sequence_of_element_c)
                    #sequences_of_complement_elements += ''.join(sequence_of_element_c)
                    #sequences_of_both_elements += ''.join(sequence_of_element_c)


    #number
    print('element counting finished\n')
    print('number of elements in chromosome strands')
    print('in both    in plus    in comp')
    print(counting_of_both_elements,'  ',counting_of_plus_elements,'  ',counting_of_comp_elements)
    print('All covered reagion',counting_of_both_elements+counting_of_plus_elements+counting_of_comp_elements)

    #sequence portion

    plus_elements_portion_in_sequence = round((ch_sequence.count('1')+ch_sequence.count('2')+ch_sequence.count('3')
                                               +ch_sequence.count('4'))/len(ch_sequence)*100,2)
    complement_elements_portion_in_sequence = round((ch_sequence.count('5')+ch_sequence.count('6')
                                                     +ch_sequence.count('7')+ch_sequence.count('8'))/len(ch_sequence)*100,2)
    both_elements_portion_in_sequence = round((ch_sequence.count('A')+ch_sequence.count('T')
                                               +ch_sequence.count('C')+ch_sequence.count('G'))/len(ch_sequence)*100,2)

    print(' portion of elements in  chromosome sequence :\nboth plus comp\n',both_elements_portion_in_sequence,'%',plus_elements_portion_in_sequence,'%',complement_elements_portion_in_sequence,'%')
    len_chromosme = len(ch_sequence)
    A_b = round(ch_sequence.count('A') / len_chromosme * 100, 2)
    T_b = round(ch_sequence.count('T') / len_chromosme * 100, 2)
    C_b = round(ch_sequence.count('C') / len_chromosme * 100, 2)
    G_b = round(ch_sequence.count('G') / len_chromosme * 100, 2)
    ATCG_b = (ch_sequence.count('A')+ch_sequence.count('T')+ch_sequence.count('C')+ch_sequence.count('G'))
    len_sum_of_pluses = len(sequences_of_plus_elements)
    A_d = round(ch_sequence.count('1') / len_chromosme * 100, 2)
    T_d = round(ch_sequence.count('2') / len_chromosme * 100, 2)
    C_d = round(ch_sequence.count('3') / len_chromosme * 100, 2)
    G_d = round(ch_sequence.count('4') / len_chromosme * 100, 2)
    ATCG_P = (ch_sequence.count('1')+ch_sequence.count('2')+ch_sequence.count('3')+ch_sequence.count('4'))
    len_sum_of_complements = len(sequences_of_complement_elements)
    A_c = round(ch_sequence.count('5') / len_chromosme * 100, 2)
    T_c = round(ch_sequence.count('6') / len_chromosme * 100, 2)
    C_c = round(ch_sequence.count('7') / len_chromosme * 100, 2)
    G_c = round(ch_sequence.count('8') / len_chromosme * 100, 2)
    ATCG_C = (ch_sequence.count('5')+ch_sequence.count('6')+ch_sequence.count('7')+ch_sequence.count('8'))

    print('both composition','+composition', '-composition')
    print('A ',A_b,'    ','A ', A_d, '   ', ' A ', A_c)
    print('T ',T_b,'    ','T ', T_d, '   ', ' T ', T_c)
    print('C ',C_b,'    ','C ', C_d, '   ', ' C ', C_c)
    print('G ',G_b,'    ','G ', G_d, '   ', ' G ', G_c)

    print('aditional information')
    print('Both composition')
    print(round(ch_sequence.count('A') / ATCG_b * 100, 2))
    print(round(ch_sequence.count('T') / ATCG_b * 100, 2))
    print(round(ch_sequence.count('C') / ATCG_b * 100, 2))
    print(round(ch_sequence.count('G') / ATCG_b * 100, 2))
    print('Plus composition')
    print(round(ch_sequence.count('1') / ATCG_P * 100, 2))
    print(round(ch_sequence.count('2') / ATCG_P * 100, 2))
    print(round(ch_sequence.count('3') / ATCG_P * 100, 2))
    print(round(ch_sequence.count('4') / ATCG_P * 100, 2))
    print('comp composition')
    print(round(ch_sequence.count('5') / ATCG_C * 100, 2))
    print(round(ch_sequence.count('6') / ATCG_C * 100, 2))
    print(round(ch_sequence.count('7') / ATCG_C * 100, 2))
    print(round(ch_sequence.count('8') / ATCG_C * 100, 2))

    #print(counting_pluses_in_difference_length)
    both_stand_information = [counting_of_both_elements,counting_both_in_difference_length,sequences_of_both_elements,len(ch_sequence),compposintion_in_diffrence_length_both]
    plus_stand_information = [counting_of_plus_elements, counting_pluses_in_difference_length, sequences_of_plus_elements, len(ch_sequence),compposintion_in_diffrence_length_plus]
    comp_strand_information = [counting_of_comp_elements,  counting_complements_in_difference_length, sequences_of_complement_elements, len(ch_sequence),compposintion_in_diffrence_length_comp]




    del sequences_of_plus_elements
    del sequences_of_complement_elements
    del sequences_of_both_elements
    del counting_of_plus_elements
    del counting_of_comp_elements
    del counting_of_both_elements
    del counting_pluses_in_difference_length
    del counting_complements_in_difference_length
    del counting_both_in_difference_length
    del ch_sequence
    #del sequence_of_element_p
    #del sequence_of_element_c

    gff.close()
    print('pass function : counting_of_elements_with_lowering_in_begin\n\n')
    return plus_stand_information, comp_strand_information,both_stand_information #length_of_overlaping



def test(sequence0, gff,chnumber,element):
    sequence = sequence0.lower()  # ATCG all
    sequence_plus = sequence
    sequence_plus = list(sequence_plus) #ATCG all
    sequence_comp = sequence
    sequence_comp = list(sequence_comp) #ATCG all
    sequence_both = sequence
    sequence_both = list(sequence_both) #ATCG all
    sequence_lower = sequence.upper() # this is for masking


    chromosome_number = chnumber

    gff = open(gff)
    gff.readline()
    gff.readline()
    gff.readline()
    element_names_in_line = []
    counter = 0
    counter2 = 0
    for annotation in gff:

        if annotation[0] == '#':
            continue
        annotation.strip()
        annotation = annotation.split()

        if (f'chr{chromosome_number[0]}' == annotation[0] or chromosome_number[1] in annotation[0]) and (element[0] in annotation[2] ):
            #if f'chr{chromosome_number}' == annotation[4] and  annotation[10].split('/')[0] not in element:
            #if f'chr{chromosome_number}' == annotation[4]:
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
                '''if sequence0[(start - 1):end].count('A') > 0 or sequence0[(start - 1):end].count('T') > 0\
                        or sequence0[(start - 1):end].count('C') > 0 or sequence0[(start - 1):end].count('G') > 0:

                    print('Hi!!!!!!!!!!!')
                else:
                    print('Hello')
                    counter +=1
                    print(counter,len(sequence0[(start - 1):end]),sequence0[(start - 1):end])'''

            if annotation[6] == '-':
                sequence_comp[start - 1:end] = sequence_lower[start - 1:end]
                sequence_both[(start - 1):end] = sequence_lower[(start - 1):end]
                '''if sequence0[(start - 1):end].count('A') > 0 or sequence0[(start - 1):end].count('T') > 0 \
                        or sequence0[(start - 1):end].count('C') > 0 or sequence0[(start - 1):end].count('G') > 0:

                    print('Hi!!!!!!!!!!!')
                else:
                    print('Hello')
                    counter +=1
                    print(counter,len(sequence0[(start - 1):end]),sequence0[(start - 1):end])'''
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





def merge_seq_seq_share_with_1234(seq1,seq2):
    C = []
    for i,j in zip(seq1,seq2):
        Ni = False
        Nj = False
        if i.isupper() == True:
            Ni = True
        if j.isupper() == True:
            Nj = True
        if Ni == False and Nj == False:
            C.append(i)
        if Ni == True and Nj == False:
            if i == 'A':
                C.append('~')
            if i == 'T':
                C.append('!')
            if i == 'C':
                C.append('@')
            if i == 'G':
                C.append('#')

        if Ni == False and Nj == True:
            C.append(j)
        if Ni == True and Nj == True :
            if i == 'A':
                C.append('1')
            if i == 'T':
                C.append('2')
            if i == 'C':
                C.append('3')
            if i == 'G':
                C.append('4')

    return C


