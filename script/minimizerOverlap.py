# Algorithms in Bioinformatics assignment 2 to get minimizers overlap
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

# get the minimizers subsequences
def get_minimizer(seq, k, w):
    minimizer = set()
    length = len(seq)
    for i in range(length-k+1):
        # move the minimizer frame further left
        if i == 0:
            subseqs = list()
            for j in range(w):
                subseq = seq[i+j:i+j+k:1] # get the subsequence
                rc_subseq = get_rc_subseq(subseq) # get the reverse complement subsequence
                subseqs.append(subseq if subseq < rc_subseq else rc_subseq) # get the k-mer subsequence set by comparing the two subsequences
                minimum_subseq = min(subseqs)
                minimizer.add(minimum_subseq) # get the minimizer subsequence and store it in minimizer set
        # move the minimizer frame further right
        elif i > length-k-w+1:
            del subseqs[0]
            minimizer.add(min(subseqs)) # move the minimizer frame by deleting the first subsequence and adding a new minimizer subsequence to the minimizer set
        # move the minimizer by each base normally
        else:
            pre_subseq = subseqs.pop(0) # move the minimizer frame by delete the first subsequence
            subseq = seq[i+w-1:i+w+k-1:1] # get the subsequence
            rc_subseq = get_rc_subseq(subseq) # get the reverse complement subsequence
            subseqs.append(subseq if subseq < rc_subseq else rc_subseq) # get the k-mer subsequence set by comparing the two subsequences
            # if the deleted subsequence is the minimizer subsequence of previous minimizer frame
            if minimum_subseq == pre_subseq:
                minimum_subseq = min(subseqs)
                minimizer.add(minimum_subseq) # get the current minimizer subsequence and move the minimizer by adding it to the minimizer set
            # if the minimizer subsequence of previous minimizer frame is still in the current minimizer frame
            else:
                minimum_subseq = min(subseqs[w-1],minimum_subseq)
                minimizer.add(minimum_subseq) # get the minimizer subsequence by comparing new subsequence with previous minimizer subsequence and move the minimizer by adding it to the minimizer set
    return minimizer

# calculating the jaccard similarity of two sets of k-mer subsequences
def get_jaccard_similarity(minimizer_S, minimizer_R):
    jaccard_similarity = len(minimizer_S & minimizer_R) / len(minimizer_S | minimizer_R) # calculating the jaccard similarity of R and R' in S
    return jaccard_similarity

if __name__ == '__main__':

    import sys
    import time

    output_jaccard_similarity = False
    output_running_time = True # set the prameters to decide whether to output the running time or jaccard similarity in order to test the efficiency of the program or get the result

    file_path_S = sys.argv[1]
    file_path_R = sys.argv[4]
    k = int(sys.argv[2])
    w = int(sys.argv[3]) # get the command line parameters including the file path of S sequences and R sequence, the parameter k of k-mer subsequence and the parameter w of minimizer frame
    
    start_time = time.time() # start the timer

    seqs_S = get_seq_list(file_path_S)
    seq_R = get_seq_list(file_path_R)[0] # get the sequences

    minimizer_R = get_minimizer(seq_R, k, w) # get the minimizer subsequences of S
    for seq_S in seqs_S:
        minimizer_S = get_minimizer(seq_S, k, w) # get the minimizer subsequences of S
        jaccard_similarity = get_jaccard_similarity(minimizer_S, minimizer_R) # calculating the jaccard similarity
        if output_jaccard_similarity:
            print(jaccard_similarity) # output the jaccard similarity

    end_time = time.time() # stop timing
    
    if output_running_time:
        print('Running time: ', end_time - start_time, 's', sep = '') # output the running time of the program