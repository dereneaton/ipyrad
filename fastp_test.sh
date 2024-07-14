#!/usr/bin/env bash

fastp=/home/deren/mambaforge-pypy3/envs/ipyrad1/bin/fastp

# 2, global trimming at front (--trim_front)
# 4, quality pruning at 5' (--cut_front)
# 5, quality pruning by sliding window (--cut_right)
# 6, quality pruning at 3' (--cut_tail)

cat <<EOT > /tmp/test.fastq
@NB502016:186:H2HNLBGXN:1:11101:4872:1081 1:N:0:0
AACTGGTGCAGTATAATCTATTTTAAATACCTCCATCCTAAATAGAGGAACCCTATGGTCTAACCCCCATTTAATTGCATTACGAATGTCAA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@NB502016:186:H2HNLBGXN:1:11101:20414:1089 1:N:0:0
AACTGGTGCAGGAAATTGTTTGGCCGCAGGAAGAAGAATAAGAATCGTCCGAATTTGTCTTCCAATAATTCTGACTTGCTTTGTTCTGGTTC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@NB502016:186:H2HNLBGXN:1:11101:21961:1096 1:N:0:0
AACTGGTGCAGCTGAGTAACACCCGCACAAACCAACACCATAGGTATTAGGTAGGACCAACTGCTAGCACATTGGCCGAACTAATCCAAATC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@NB502016:186:H2HNLBGXN:1:11101:4227:1163 1:N:0:0
AACTGGTGCAGCCCGTAAGCCGACACCACCGCGTCTAAAAATCCGTAACCATTAGACAGCCCCGCGGACTTAACGACTTTTCCACAAACGGC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOT

# THIS WORKS, 6 bases are trimmed from front so that all new reads start with 'TGCAG'
$fastp -i /tmp/test.fastq -o /tmp/test.trimmed.fastq -h /tmp/test.html -j /tmp/test.json -f 6
cat /tmp/test.fastq
echo "================================================================================================="
cat /tmp/test.trimmed.fastq

# THIS DOESN'T WORK, for some reason now 9 bases are trimmed!
$fastp -i /tmp/test.fastq -o /tmp/test.trimmed.fastq -h /tmp/test.html -j /tmp/test.json -3 -5 -f 6
cat /tmp/test.fastq
echo "================================================================================================="
cat /tmp/test.trimmed.fastq

# this shows that quality pruning on its own is not the problem.
$fastp -i /tmp/test.fastq -o /tmp/test.trimmed.fastq -h /tmp/test.html -j /tmp/test.json -5 -3
cat /tmp/test.fastq
echo "================================================================================================="
cat /tmp/test.trimmed.fastq