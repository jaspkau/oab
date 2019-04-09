#!/bin/bash

#########################
#
# Simple SparCC wrapper 
#
# Use:
# chmod 755 sparccWrapper.sh
# ./sparccWrapper.sh &
#
# Author: Karoline Faust
#
#########################

python scripts/yonatanf-sparcc-05f4d3f31d77/SparCC.py results/net_pseudo/net.txt -i 100 --cor_file=results/net_pseudo/sim_cor.txt

python scripts/yonatanf-sparcc-05f4d3f31d77/MakeBootstraps.py results/net_pseudo/net.txt -n 100 &&

# compute sparcc on resampled (with replacement) datasets
for i in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99
do

python scripts/yonatanf-sparcc-05f4d3f31d77/SparCC.py results/net_pseudo/net.txt.permuted_$i.txt -i 100 --cor_file=results/net_pseudo/sim_cor_$i.txt >> results/net_pseudo/sparcc.log

done &&

# compute p-value from bootstraps
for i in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99
do

python scripts/yonatanf-sparcc-05f4d3f31d77/PseudoPvals.py results/net_pseudo/sim_cor.txt results/net_pseudo/sim_cor_$i.txt 100 -o results/net_pseudo/pvals_two_sided.txt -t two_sided  >> results/net_pseudo/sparcc.log

done

# visualization requires parsing and thresholding the p-value OTU matrix

#python PseudoPvals.py example/basis_corr/cor_sparcc.out example/pvals/perm_cor_#.txt 5 -o example/pvals/pvals.one_sided.txt -t two_sided

