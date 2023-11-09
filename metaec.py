#!/usr/bin/python3

import argparse

from metatreeop import count_wgd_nodes_combined
from treeop import str2tree, Tree
import time
import re
import os
import sys
    
    
import random

sys.setrecursionlimit(3000)

def main():
    parser = argparse.ArgumentParser(
        description="Run WGD reconciliation algorithm for input species tree and gene trees with ?")
    parser.add_argument("--gene_trees", help="Path to a file with newline separated gene trees", type=str, default="data_sim/wgd-1-gene-trees")
    parser.add_argument("--species_tree", help="Path to a file with a species tree", type=str, default="data_sim/s_tree")
    parser.add_argument("--initial_gene_tree", help="Path to a file with the initial gene tree (already precomputed)", type=str, default=None)
    parser.add_argument("--out_file", help="Path to an output file with results", type=str, default="")
    parser.add_argument("--randomize_from", help="Start randomizing from a given size of binom(n,k) in the main loop", type=int, default=0)
    parser.add_argument("--noimprovement_stop", help="How many times to run DP with no improvement (0 - do not stop)", type=int, default=0)
    parser.add_argument("--reversed_climb", help="Start from fixedwgd and interatively search in larger sets untils solution is found", action='store_true')   
    parser.add_argument("--distribution_maps", help="Add distributions maps", action='store_true')   
    parser.add_argument("--reference_trees", help="A path to a reference gene trees", type=str, default=None)   
    parser.add_argument("--distribution_maps_epi", help="Compute distributions using episode set from initial gene trees", action='store_true')   
    parser.add_argument("--print_distr_maps", help="Print the output tree with distribution maps", action='store_true')
    parser.add_argument("--verbose", help="0 - none, 1 - basic, 2 - print wgd nodes", type=int, default=1)
    parser.add_argument("--distr_counts", help="Do not normalize distr maps", action='store_true')

    args = parser.parse_args()

    gene_trees = open(args.gene_trees).read().split()
    gene_trees = [Tree(str2tree(g_str)) for g_str in gene_trees]
    species_tree = open(args.species_tree).read()
    species_tree = Tree(str2tree(species_tree))    

    t = time.process_time()
    
    setid=re.sub('[A-Za-z/-]','',args.gene_trees)
    setid=re.sub('^_*','',setid)


    initial_gene_tree = None
    if args.initial_gene_tree:
        with open(args.initial_gene_tree) as f:
            initial_gene_tree = f.read()

    reference_trees = None
    if args.reference_trees:
        with open(args.reference_trees) as f:
            reference_trees = [Tree(str2tree(g_str)) for g_str in f.read().split() ]


    cost, used_nodes, exactsolution, outstats = count_wgd_nodes_combined(
            species_tree, 
            gene_trees, 
            wgddebug = False, 
            outfile = args.out_file,
            noimprovement_stop = args.noimprovement_stop,
            randomize_from = args.randomize_from,
            setid = setid,        
            reversed_climb = args.reversed_climb,
            initial_gene_tree = initial_gene_tree,
            distribution_maps = args.distribution_maps,
            reference_trees = reference_trees,
            distribution_maps_epi = args.distribution_maps_epi,
            print_distr_maps = args.print_distr_maps,
            outgroup="o",
            verbose=args.verbose,
            distr_counts=args.distr_counts
            )

    endtime = time.process_time() - t

    if args.out_file:        
        with open(args.out_file,"w") as f:
            f.write(f"gene_trees_file={args.gene_trees}\n")        
            f.write(f"species_tree_file={args.species_tree}\n")                
            f.write(f"{outstats}")
            f.write(f"randomize_from={args.randomize_from}\n")
            f.write(f"noimprovement_stop={args.noimprovement_stop}\n")
            f.write(f"time={endtime}")
        
    if args.verbose:
        print(f"[{setid}] Cost: {cost} Exact:{exactsolution}")

    if args.verbose==2:
        print("Used nodes: ")
        for node in used_nodes:
            print(node)
    
    assert len(used_nodes) == cost


if __name__ == "__main__":
    main()

