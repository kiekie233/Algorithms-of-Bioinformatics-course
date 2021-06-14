# Algorithms in Bioinformatics on calculating jaccard similarity

In this Assignment, two program will be designed to calculate the similarity of two sequence. One is the k_mer overlap algorithm, the other is the minimizer overlap algorithm. Besides, we will optimize these two algorithms not only for saving the calculation time but also for save the space. The result, performance and discussion of these two algorithms will show below.

## Program steps

### k-mer overlap

The command format is `python kmerOverlap.py <input.fa> <k> <template.fa>` and  you can change the parameter in the program to choose whether output the jaccard similarity and running time.

```shell
cd script/
python kmerOverlap.py ../data/eg1.fa 7 ../data/template1.fa > ../result/out1.txt
python kmerOverlap.py ../data/eg2.fa 16 ../data/template2.fa > ../result/out2.txt
python kmerOverlap.py ../data/eg3.fa 16 ../data/../data/template3.fa > ../result/out3.txt
python kmerOverlap.py ../data/eg4.fa 16 ../data/template4.fa > ../result/out4.txt # running the program and get the output data
```

### minimizer overlap

The command format is `python minimizerOverlap.py <input.fa> <k> <w> <template.fa>` and  you can change the parameter in the program to choose whether output the jaccard similarity and running time.

```shell
cd script/
python minimizerOverlap.py ../data/eg1.fa 7 30 ../data/template1.fa > ../result/mout1.txt
python minimizerOverlap.py ../data/eg2.fa 16 30 ../data/template2.fa > ../result/mout2.txt
python minimizerOverlap.py ../data/eg3.fa 16 30 ../data/template3.fa > ../result/mout3.txt
python minimizerOverlap.py ../data/eg4.fa 16 30 ../data/template4.fa > ../result/mout4.txt # running the program and get the output data
```

## Common steps in these two program

### Get sequence list

A function is defined to read the sequence data. Considering there is no need to use the name of sequences, the algorithms only read the sequences information by identify the sequence identifier '>' and store them in a list.

### Get reverse complement sequences

A function is defined to get the reverse complement sequence of a sequence. The translation pattern is set as 'ACGT' -> 'TGCA' and the sequence will be translate into complement sequence by this pattern. After that, the reverse complement sequence will be obtained by indexing the complement sequence with the step size of -1.

### Get jaccard similarity

A function is defined to calculate the jaccard similarity of two sets of subsequences. The similarity will be calculated by divide the intersection of these two sets of subsequences by the union of that.

### Other details

The 'sys' module is imported to get the command line parameters including the file path of S and R sequence and the parameter of k or k and w.

The 'time' module is imported to time the running time of the program.

## Algorithm of kmerOverlap.py

The unique algorithm of kmerOverlap.py is as follows:

### Get k-mer subsequences

A function is defined to get the k-mer subsequences from a sequence.

Firstly, a 'set' data structure is created to prepare to store the k-mer subsequences.

Secondly, the sequence is traversed and a k-length subsequence and the reverse complement subsequence is got for every base advance.

Thirdly, the minimum subsequence is got by comparing the above two sequences and is added to the k-mer subsequences set.

Finally, the jaccard similarity will be calculated.

### Time & space complexity

For one read, the time complexity of this algorithm is O(L(k+2k+2)) = O(3Lk+2L). The time complexity of calculating of jaccard similarity of two reads is O(n1+n2+min(n1+n2)), n<<L. So the total time complexity for one read is about **O(3Lk)**. L is the length of the sequence, and k is the size of k-mer.

for one read of S and one read R, the space complexity of this algorithm is O(n1+n2). So the space complexity of one read is **O(n)**. n is the number of k-mer subsequences which is duplicated.

### Running time of the program

For the eg1, the running time is 0.002s. 

For the eg2, the running time is 80.5s.

For the eg3, the running time is 136s.

For the eg4, the running time is 8.25s

![image-20210613230324592](C:\Users\kiekie\AppData\Roaming\Typora\typora-user-images\image-20210613230324592.png)

## Algorithm of minimizerOverlap.py

The unique algorithm of minimizerOverlap.py is as follows: https://github.com/kiekie233/Algorithms-of-Bioinformatics-course

### Get minimizer subsequences

A function is defined to get the minimizer subsequences from a sequence.

Firstly, a 'set' data structure is created to prepare to store the minimizer subsequences and a 'list' data structure is created to prepare to store subsequences in the current minimizer frame.

Secondly, when the minimizer frame is at the beginning of the sequence, the length of frame is extended from k to k+w-1. And when the length extends by 1, the subsequences list is filled by compare the subsequence and the reverse complement one of it in every base in the frame and the minimizer subsequence is got by comparing all the subsequences in the subsequences list.

Thirdly, when the frame is moving on, the first subsequence of the subsequences list will be moved out and the new subsequence is added in it. If the subsequence moved out is the minimizer subsequence of last frame, the next minimizer subsequence of  the current frame will be got by comparing all the subsequences in this frame, otherwise the next minimizer subsequence of the current frame will be got by only comparing the last minimizer subsequence with the new subsequence in the subsequences list for saving the computing time.

Fourthly, when the frame is at the end of the sequence, the length of the frame is reduced from k+w-1 to k. And when the length reduces by 1, the first subsequence of the subsequences list will be moved out and the minimizer subsequence will be got by comparing the remaining subsequences in the frame.

Finally, the the jaccard similarity will be calculated.

### Time & space complexity

For one read S and one read R, the time complexity of this algorithm is O(w(k+2k+2+0.5w+1)+L(w+k+2k+2+((w-1)(1+1)+w+1)/w)+w(0.5w+0.5w)) = O((L+w)(3K+w+3)-0.5w²-2L). The time complexity of calculating of jaccard similarity of two reads is O(n1+n2+min(n1+n2)), n<<L. So the total time complexity for one read is about **O((L+w)(3k+w))**. L is the length of the sequence, k is the size of k-mer, and n is the number of k-mer subsequences which is duplicated.

For one read S and one read R, the space complexity of this algorithm is O(n1+n2+w+w),n>>w. So the space complexity of one read is **O(n)**. n is the number of k-mer subsequences which is duplicated.

### Running time of the program

For the eg1, the running time is 0.001s.

For the eg2, the running time is 117s.

For the eg3, the running time is 184s.

For the eg4, the running time is 12.0s.

![image-20210613231621612](C:\Users\kiekie\AppData\Roaming\Typora\typora-user-images\image-20210613231621612.png)

## Discussion

The J(R,R') and J'(R,R') are both the index to describe the similarity of two sequence. Of course we can get the exact similarity of two sequences by aligning them, but this process will take too much time. So we choose J(R,R') and J'(R,R') to get the approximate conclusion of sequence similarity. For J(R,R'), we get a k-mer for each base which means there are still so much subsequence to be used to calculate the similarity. So we also used a new method minimizer to calculate the similarity. In J'(R,R'), the subsequence is the minimum subsequence of a w-length frame so that each the subsequence represents the local sequence information rather than k-mer which only represent the information of a sequence locus. Besides, the minimizer method will reduce the number of subsequences because several adjacent frame may share the same sequence. All in all, k-mer and minimizer are both great methods to calculate the similarity of two sequences in a fast way, we can choose one of them according to different situation of sequence characteristics.

