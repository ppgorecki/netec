#!/bin/bash

eval "git clone https://github.com/j-paszek/genomicduplicationilp.git"
data_path="genomicduplicationilp/wgd-simulated"

# Simulate true gene trees and sequences
for i in {0..5}; do
    if [ $i -eq 0 ]; then
        s_path="s_tree.nexus"
    else
        s_path="${data_path}/n20-wgd-${i}.nexus"
    fi

    echo ${s_path}
    command1="simphy -sr ${s_path} -rl F:100 -rg 1 -sp F:10 -su F:0.0000000004 -lb F:2e-10 -ld F:2e-10 -ll 4 -hg LN:1.5,1 -ot 0 -oc 1 -sp f:10 -v 3 -cs 85341 -o ${data_path}/n20-wgd-${i}"
    eval "$command1"

    command2="perl INDELIble_wrapper.pl ${data_path}/n20-wgd-${i} parameters_indelible.txt 177 1"
    eval "$command2"
done

# ML Inference
for subdir in $(find ${data_path} -type f -name "*.phy"); do
    outpath="${subdir}_tree_ML"
    ./FastTree -nt -gtr -nosupport $subdir > $outpath
done

python3 cleanup_and_rooting.py
