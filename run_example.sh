#!/bin/bash

exe=./learnDRMN

indir=example_input

suf=atac_qmotif_chromatin_9marks
k=7
v=0
p1=10
p2=50
p3=40

order=${indir}/order.txt
ogids=${indir}/ogids.txt
tree=${indir}/tree.txt
conf=${indir}/${suf}_k${k}/config.txt
src=esc

out=example_out
mkdir -p ${out}

$exe ${order} ${ogids} null ${k} ${tree} ${conf} none ${out} learnCV:${v}:12345 ${src} uniform 0.80 selfInit LEASTFUSED ${p1} ${p2} ${p3} 
