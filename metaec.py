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
    
    parser.add_argument("--out_file", help="Path to an output file or directory with results", type=str, default="")
    
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

    parser.add_argument("--locked_epi_support", help="For every gene tree and every net node identify locked episodes",  action='store_true')

    parser.add_argument("--mindup", help="Filter gene trees with duplication count < MINDUP (only for networks with no reticulations, i.e., trees; ignored for networks); default is 0", type=int, default=0)

    parser.add_argument("--print_dup_stats", help="Print duplication statistics: for each k=0,1,2,... show how many gene trees have k duplications", action='store_true')

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

    network, netclass = args.network, Network
    network = netclass(str2tree(open(network).read()))

    # Filter gene trees based on mindup parameter (only for trees, not networks)
    if args.mindup > 0:
        if len(network.reticulations) == 0:
            original_count = len(gene_trees)
            filtered_gene_trees = []
            for gt in gene_trees:
                gt.set_lca_mapping(network)
                if gt.dupcost(network) >= args.mindup:
                    filtered_gene_trees.append(gt)
            gene_trees = filtered_gene_trees
            if args.verbose:
                print(f"Filtered gene trees: {original_count} -> {len(gene_trees)} (mindup={args.mindup})")
        else:
            if args.verbose:
                print(f"Warning: --mindup parameter ignored (network has {len(network.reticulations)} reticulation(s))")

    # Print duplication statistics if requested
    if args.print_dup_stats:
        dup_counts = {}
        for gt in gene_trees:
            if not hasattr(gt.root, 'lcamap'):  # Check if LCA mapping already exists
                gt.set_lca_mapping(network)
            dup_cost = gt.dupcost(network)
            dup_counts[dup_cost] = dup_counts.get(dup_cost, 0) + 1

        print("Duplication statistics:")
        if dup_counts:
            max_dup = max(dup_counts.keys())
            for k in range(max_dup + 1):
                count = dup_counts.get(k, 0)
                print(f"  {k} duplications: {count} gene tree(s)")
        else:
            print("  No gene trees to analyze")

        return

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

    if args.locked_epi_support:        

        locked_episodes = []
        for gt in gene_trees:
            net, locked_epi = count_wgd_nodes_combined(
                network, 
                [gt], 
                wgddebug = False,             
                out_file = None,
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
                find_fixed_episodes = not args.no_fixed_episodes_search,
                locked_epi_support = args.locked_epi_support
            )
            locked_episodes.append(locked_epi)            

        d = {}
        with open(out_basefile+"locked_epi",'w') as f:
            for i in locked_episodes:
                f.write(" ".join(str(e.num) for e in i)+"\n")

                for e in i:
                    if e.num not in d: d[e.num]=0
                    d[e.num]+=1

        for n in net.nodes:
            if n.num in d:
                n.lockedepisupport = d[n.num]

        with open(out_basefile+"locked_epi_net",'w') as f:
            f.write(net.root.attrrepr(['num','lockedepisupport'], gsestyle=args.gsestyle))
                
        return


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


