# Algorithms in Bioinformatics assignment 2 to get k-mers overlap
# Author: 余思克(Sike Yu) HZAU Bioinformatics, Wuhan, Hubei, China
# date: Jun. 2021

# read and store the sequecnes in the list
def get_seq_list(file_path):
    with open(file_path) as file:
        text_list = file.readlines() # read the sequence file
    seq = list()
    for text in text_list:
        text = text.strip() # remove the '\n' character
        if text[0] != '>':
            seq.append(text.upper()) # store each sequence separately in the list
    return seq

# reverse the subsequence and get the complement one of it
def get_rc_subseq(subseq):
    trans_table = str.maketrans('ACGT', 'TGCA') # make the translation pattern
    rc_subseq = subseq[::-1].translate(trans_table) # get the reverse complement subsequence
    return rc_subseq

# get the k-mer subsequences
def get_k_mer(seq, k):
    k_mer = set()
    for i in range(len(seq)-k+1):
        subseq = seq[i:i+k:1] # get the subsequence
        rc_subseq = get_rc_subseq(subseq) # get the reverse complement subsequence
        k_mer.add(subseq if subseq < rc_subseq else rc_subseq) # get the k-mer subsequence set by comparing the two subsequences
    return k_mer

# calculating the jaccard similarity of two sets of k-mer subsequences
def get_jaccard_similarity(k_mer_S, k_mer_R):
    jaccard_similarity = len(k_mer_S & k_mer_R) / len(k_mer_S | k_mer_R) # calculating the jaccard similarity of R and R' in S
    return jaccard_similarity

if __name__ == '__main__':

    import sys
    import time

    output_jaccard_similarity = True
    output_running_time = False # set the prameters to decide whether to output the running time or jaccard similarity in order to test the efficiency of the program or get the result

    file_path_S = sys.argv[1]
    file_path_R = sys.argv[3]
    k = int(sys.argv[2]) # get the command line parameters including the file path of S sequences and R sequence and the parameter k of k-mer subsequence

    start_time = time.time() # start the timer

    seqs_S = get_seq_list(file_path_S)
    seq_R = get_seq_list(file_path_R)[0] # get the sequences

    k_mer_R = get_k_mer(seq_R, k) # get the k-mer subsequences of R
    for seq_S in seqs_S:
        k_mer_S = get_k_mer(seq_S, k) # get the k-mer subsequences of S
        jaccard_similarity = get_jaccard_similarity(k_mer_S, k_mer_R) # calculating the jaccard similarity
        if output_jaccard_similarity:
            print(jaccard_similarity) # output the jaccard similarity

    end_time = time.time() # stop timing
    
    if output_running_time:
        print('Running time: ', end_time - start_time, 's', sep = '') # output the running time of the program