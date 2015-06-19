#! /bin/bash

#This script aims to get the rank of data via GFOLD. The rank is based the value of GFOLD.


#ARG
cd ../data/AGR
gfold diff -s1 AGR1,AGR2,AGR3,AGR4 -s2 AGR5,AGR6,AGR7,AGR8 -suf .rc -o AGR.diff
awk  'BEGIN{OFS="\t"}!/^#/{print $1,$3}' AGR.diff > ../../results/AGR/gfold_AGR.rnk
awk  'BEGIN{OFS="\t"}!/^#/{print $1,$5}' AGR.diff > ../../results/AGR/gfoldFC_AGR.rnk


#NWU
cd ../NWU
time gfold diff -s1 NWU1,NWU2,NWU3,NWU4,NWU5 -s2 NWU6,NWU7,NWU8,NWU9,NWU10 -suf .rc -o NWU.diff
awk  'BEGIN{OFS="\t"}!/^#/{print $1,$3}' NWU.diff > ../../results/NWU/gfold_NWU.rnk
awk  'BEGIN{OFS="\t"}!/^#/{print $1,$5}' NWU.diff > ../../results/NWU/gfoldFC_NWU.rnk

#BGI
cd ../BGI
time gfold diff -s1 BGI1,BGI2,BGI3,BGI4,BGI5 -s2 BGI6,BGI7,BGI8,BGI9,BGI10 -suf .rc -o BGI.diff
awk  'BEGIN{OFS="\t"}!/^#/{print $1,$3}' BGI.diff > ../../results/BGI/gfold_BGI.rnk
awk  'BEGIN{OFS="\t"}!/^#/{print $1,$5}' BGI.diff > ../../results/BGI/gfoldFC_BGI.rnk

#PSU
cd ../PSU
time gfold diff -s1 PSU1,PSU2,PSU3,PSU4,PSU5 -s2 PSU6,PSU7,PSU8,PSU9,PSU10 -suf .rc -o PSU.diff
awk  'BEGIN{OFS="\t"}!/^#/{print $1,$3}' PSU.diff > ../../results/PSU/gfold_PSU.rnk
awk  'BEGIN{OFS="\t"}!/^#/{print $1,$5}' PSU.diff > ../../results/PSU/gfoldFC_PSU.rnk

