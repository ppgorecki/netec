#!/usr/bin/python3

import argparse

from metatreeop import count_wgd_nodes_combined
from treeop import str2tree, Tree
from netop import Network
import time
import re
import os
import sys
    
    
import random

sys.setrecursionlimit(3000)

def main():
    parser = argparse.ArgumentParser(
        description="Run WGD reconciliation algorithm for input network/tree and gene trees with ?")
    
    parser.add_argument("--gene_trees", help="Path to a file with newline separated gene trees", type=str, default=None)
    
    parser.add_argument("--network", help="Path to a file with a species network", type=str, default=None)
    
    parser.add_argument("--species_tree", help="Path to a file with a species tree", type=str, default="data_sim/s_tree")
    
    parser.add_argument("--initial_gene_tree", help="Path to a file with the initial gene tree (already precomputed)", type=str, default=None)
    
    parser.add_argument("--out_file", help="Path to an output file with results", type=str, default="")
    
    parser.add_argument("--randomize_from", help="Start randomizing from a given size of binom(n,k) in the main loop", type=int, default=0)
    
    parser.add_argument("--noimprovement_stop", help="How many times to run DP with no improvement (0 - do not stop)", type=int, default=0)

    parser.add_argument("--reversed_climb", help="Start from fixedwgd and interatively search in larger sets untils solution is found", action='store_true')   
    parser.add_argument("--distribution_maps", help="Add distributions maps", action='store_true')   
    parser.add_argument("--reference_trees", help="A path to a reference gene trees", type=str, default=None)   
    parser.add_argument("--distribution_maps_epi", help="Compute distributions using episode set from initial gene trees (networks not implemented yet)", action='store_true')   

    parser.add_argument("--print_distr_maps", help="Print the output tree with distribution maps (networks not implemented yet)", action='store_true')

    parser.add_argument("--save_embedding", help="Save inferred embedding to [outfile].embedding", action='store_true')

    parser.add_argument("--episummaryfile", help="Save the species/network with the episize attributes", action="store_true")

    parser.add_argument("--verbose", help="0 - none, 1 - basic, 2 - print wgd nodes", type=int, default=1)

    parser.add_argument("--distr_counts", help="Do not normalize distr maps", action='store_true')

    parser.add_argument("--gsestyle", help="Use gse output for attributes, reticulations ids without #; default is newick", action='store_true')

    parser.add_argument("--user_episodes", help="User defined list episodes; a list of node identifiers, e.g., '2 4 10' or use 'all' for all)", type=str, default='')    

    parser.add_argument("--fixed_episodes_ext", help="User defined list of fixed episodes; use if episodes are known to be used in order to optimize computations", type=str, default='')    

    parser.add_argument("--no_fixed_episodes_search", help="Skip fixed episode search using DP (def. False)",   action='store_true')



    args = parser.parse_args()


    user_episodes = []
    if args.user_episodes:
        if args.user_episodes == 'all':
            user_episodes = 'all'
        else:
            try:
                user_episodes = list(map(int,args.user_episodes.split()))
            except:
                print("Incorrect format of user episodes. Example '2 4 10' or all", file=sys.stderr)
                sys.exit(-1)

    fixed_episodes_ext = []
    if args.fixed_episodes_ext:
        try:
            fixed_episodes_ext = list(map(int,args.fixed_episodes_ext.split()))
        except:
            print("Incorrect format of fixed episodes. Example '2 4 10'", file=sys.stderr)
            args.print_help()
            sys.exit(-1)

    if not args.gene_trees:
        print("Gene trees not specified", file=sys.stderr)
        parser.print_help()
        sys.exit(-1)

    gene_trees = open(args.gene_trees).read().split()
    gene_trees = [ Tree(str2tree(g_str)) for g_str in gene_trees]

    # To be conistent with the orignal species tree version
    #if args.network:    
    network, netclass = args.network, Network
    #else:
    #    network, netclass = args.species_tree, Tree

    network = netclass(str2tree(open(network).read()))
    
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
            reference_trees = [ Tree(str2tree(g_str)) for g_str in f.read().split() ]

    out_file = args.out_file          
    out_basefile="" 
    out_dir = "."+os.path.sep
    if out_file:
        if os.path.isdir(out_file):            
            out_dir = out_file
            out_file = out_dir + os.path.sep + "metaec.log" # default
            out_basefile = out_dir + os.path.sep 
        else:
            if len(out_file)>5 and out_file[-3]=='.':
                out_basefile = out_file[:-3]                

    cost, used_nodes, exactsolution, outstats = count_wgd_nodes_combined(
            network, 
            gene_trees, 
            wgddebug = False,             
            out_file = out_file,
            out_basefile = out_basefile,
            noimprovement_stop = args.noimprovement_stop,
            randomize_from = args.randomize_from,
            setid = setid,        
            reversed_climb = args.reversed_climb,
            initial_gene_tree = initial_gene_tree,
            distribution_maps = args.distribution_maps,
            reference_trees = reference_trees,
            distribution_maps_epi = args.distribution_maps_epi,
            print_distr_maps = args.print_distr_maps,
            save_embedding = args.save_embedding,
            outgroup="o",
            verbose=args.verbose,
            distr_counts=args.distr_counts,
            gsestyle=args.gsestyle,
            user_episodes = user_episodes,
            fixed_episodes_ext = fixed_episodes_ext,
            find_fixed_episodes = not args.no_fixed_episodes_search
            )

    endtime = time.process_time() - t
     
    if out_file:
        with open(out_file,"w") as f:
            f.write(f"gene_trees_file={args.gene_trees}\n")        
            f.write(f"network_file={args.network}\n")                
            f.write(f"{outstats}")
            f.write(f"randomize_from={args.randomize_from}\n")
            f.write(f"noimprovement_stop={args.noimprovement_stop}\n")
            f.write(f"time={endtime}")

    # Write output network (gse)
    if args.episummaryfile:
        pos = outstats.find('outspeciestree_wo')
        if pos>=0:            
            with open(out_basefile+"episummary","w") as f:
                f.write(outstats[pos:].split('\n')[0][19:-1])
        
    if args.verbose:
        print(f"[{setid}] Cost: {cost} Exact:{exactsolution}")

    if args.verbose==2:
        print("Used nodes: ")
        for node in used_nodes:
            print(node)
    
    assert len(used_nodes) == cost


if __name__ == "__main__":
    main()


