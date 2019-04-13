#! /usr/bin/env bash

# Set parameters
ref=FILES/chr19_100w.fa
table=FILES/table.txt
result=FILES/result.sam

work_path=$(dirname $(readlink -f $0))
cd $work_path

# Set PATH variable
export PATH=$PATH:$work_path/TOOLS/bwa-0.7.5a:$work_path/TOOLS/simulation

norsim  -r 0 -X 0 -D 0 -B 0 -I $table $ref nor.sim
readgen -d 500 -s 50 -c 20 -I nor.sim.idx $ref nor.sim left.fq right.fq

# Set ref.fa index
bwa index -a bwtsw $ref

# Seek SA coordinates
bwa aln $ref left.fq>left.sai
bwa aln $ref right.fq>right.sai

# Transform SA coordinates, output sam
bwa sampe -f $result $ref left.sai right.sai left.fq right.fq

# Remove unnecessary files
rm *.log *.sim *.idx *.txt *.fq *.sai
cd $work_path/FILES
rm *.amb *.ann *.bwt *.pac *.sa