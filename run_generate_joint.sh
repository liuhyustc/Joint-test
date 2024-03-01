#!/bin/bash  
# args=commandArgs(trailingOnly=TRUE)
# N = as.numeric(args[1])  # population
# nSNP0 = as.numeric(args[2]) # of casual SNPs
# nSNP = as.numeric(args[3])  # of non-casual SNPs
# alphaY = log(as.numeric(args[4])) # Y -> D
# rho = as.numeric(args[5]) #AR1(rho) between casual SNPs and non-causal SNPs
# eta = log(as.numeric(args[6])) # the degree of correlation between E and causal SNPs
# test = args[7] ## null or alter
# # set = args[6]  ## null--joint/main
# eff = log(as.numeric(args[8]))


for rho in {0,0.2,0.6,0.9}
do
    for expdelta in {1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2}
    do
        for test in {null,alter}
        do
            Rscript Generate_sample_joint.R 2000000 5 30 1.5 $rho $expdelta $test
        done
    done
done





# for i in "${arr[@]}"  
# do  
# echo $i  
# done


# for rho in {0,0.2,0.6,0.9} 
# do  
#     for expdelta in {1,2} 
#     do
#      args="$rho $expdelta"
#      echo $args
#     done
# done