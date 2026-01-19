#!/usr/bin/python3

import argparse

from algorithms import count_wgd_nodes_combined, epiattr, wgdnums
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
        description="Run WGD reconciliation algorithm for input network/tree and gene trees")
    
    parser.add_argument("--gene_trees", help="Path to a file with newline separated gene trees", type=str, default=None)
    
    parser.add_argument("--network", help="Path to a file with a species network", type=str, default=None)
    
    parser.add_argument("--species_tree", help="Path to a file with a species tree", type=str, default="data_sim/s_tree")
    
    parser.add_argument("--initial_gene_tree", help="Path to a file with the initial gene tree (already precomputed)", type=str, default=None)
    
    parser.add_argument("--out_dir", help="Path to directory with results; default results", type=str, default="")
    
    parser.add_argument("--randomize_from", help="Start randomizing from a given size of binom(n,k) in the main loop", type=int, default=0)
    
    parser.add_argument("--noimprovement_stop", help="How many times to run DP with no improvement (0 - do not stop)", type=int, default=0)

    parser.add_argument("--reversed_climb", help="Start from fixedwgd and interatively search in larger sets untils solution is found", action='store_true')   
    parser.add_argument("--distribution_maps", help="Add distributions maps", action='store_true')   
    parser.add_argument("--reference_trees", help="A path to a reference gene trees", type=str, default=None)   
    parser.add_argument("--distribution_maps_epi", help="Compute distributions using episode set from initial gene trees (networks not implemented yet)", action='store_true')   

    parser.add_argument("--print_distr_maps", help="Print the output tree with distribution maps (networks not implemented yet)", action='store_true')

    parser.add_argument("--save_embedding", help="Save inferred embedding to [outfile].embedding", action='store_true')    

    parser.add_argument("--verbose", help="0 - none, 1 - basic, 2 - print wgd nodes", type=int, default=1)

    parser.add_argument("--distr_counts", help="Do not normalize distr maps", action='store_true')

    parser.add_argument("--gsestyle", help="Use gse output for attributes, reticulations ids without #; default is newick", action='store_true')

    parser.add_argument("--wgddebug", help="Print all tabs from DP programming run (debug)", action='store_true')

    parser.add_argument("--user_episodes", help="User defined list episodes as a list of node identifiers, e.g., '2 4 10', 'all' for all, or a file name); the computations are done only for the given set", type=str, default='')    
    parser.add_argument("--fixed_episodes", help="List of precomputed fixed episodes; use if episodes are known to optimize computations with --no_fixed_episodes_search False ", type=str, default='')   

    parser.add_argument("--no_fixed_episodes_search", help="Skip fixed episodes search (def. False)", action='store_true')

    parser.add_argument("--extended_episodes_search", help="Identify additional episodes with large number of duplications after the best episodes were identified; saved as eeepisizepost attribute", action='store_true')    

    parser.add_argument("--extended_episodes_from_fixedepi", help="Identify non-fixed episodes with large number of duplications; saved as eeepisize attribute", action='store_true')    

    parser.add_argument("--locked_epi_support", help="For every gene tree and every net node identify locked episodes; saved as lockedepisupport attribute",  action='store_true')

    parser.add_argument("--mindup", help="Filter gene trees with duplication count < MINDUP (only for networks with no reticulations, i.e., trees; ignored for networks); default is 0", type=int, default=0)

    parser.add_argument("--print_dup_stats", help="Print duplication statistics: for each k=0,1,2,... show how many gene trees have k duplications", action='store_true')

    parser.add_argument("--fixed_episodes_only", help="Stop after fixed episodes phase", action='store_true')


    args = parser.parse_args()

    outstyleext = ".gse" if args.gsestyle else ".newick"

    def is_existing_file(filepath: str) -> bool:
        return os.path.exists(filepath) and os.path.isfile(filepath)

    def parse_episodes(f, name):
        if not f: 
            return []

        if f == 'all':
            return f
            
        try:
            if is_existing_file(f):
                f = open(f).read()

            return list(map(int, f.split()))

        except:                        
            print(f"Incorrect format for {name}. Example '2 4 10', all or a filename", file=sys.stderr)
            sys.exit(-1)

    
    user_episodes = parse_episodes(args.user_episodes, "user_episodes")
    fixed_episodes = parse_episodes(args.fixed_episodes,"fixed_episodes")
    fixed_episodes_search = not args.no_fixed_episodes_search
        
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
    
    

    initial_gene_tree = None
    if args.initial_gene_tree:
        with open(args.initial_gene_tree) as f:
            initial_gene_tree = f.read()

    reference_trees = None
    if args.reference_trees:
        with open(args.reference_trees) as f:
            reference_trees = [ Tree(str2tree(g_str)) for g_str in f.read().split() ]
    
    if args.out_dir:
        out_dir = args.out_dir        
        setid = f"[{out_dir}] "

    else:
        out_dir = "results" # default
        setid = ""

    if os.path.isfile(out_dir):
        print(f"Error: '{out_dir}' is a file, not a directory", file=sys.stderr)
        return 1

    
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)



    
    out_file = out_dir + os.path.sep + "netec.log" # default
    out_basefile = out_dir + os.path.sep 

    if args.locked_epi_support:        

        if args.verbose>=0:
            print(f"{setid}Computing locked episode support")

        locked_episodes = []
        for gt in gene_trees:
            net, locked_epi = count_wgd_nodes_combined(
                network, 
                [gt], 
                wgddebug = args.wgddebug,             
                out_file = None,
                out_basefile = out_basefile,
                noimprovement_stop = args.noimprovement_stop,
                randomize_from = args.randomize_from,
                setid = setid,        
                reversed_climb = args.reversed_climb,
                initial_gene_tree = initial_gene_tree,
                distribution_maps = args.distribution_maps,
                reference_trees = None,
                distribution_maps_epi = False,
                print_distr_maps = False,
                save_embedding = False,
                outgroup="o",
                verbose=args.verbose-1, # decrease
                distr_counts=False,
                gsestyle=args.gsestyle,                
                fixed_episodes_search = True,
                locked_epi_support = True
            )
            locked_episodes.append(locked_epi)            

        locked_epi_support = {}
        
        with open(out_basefile+"locked_epi",'w') as f:
            for i in locked_episodes:
                f.write(" ".join(str(e.num) for e in i)+"\n")

                for e in i:
                    if e.num not in locked_epi_support: locked_epi_support[e.num]=1
                    else: locked_epi_support[e.num]+=1                

    cost, used_nodes, exactsolution, outstats, stroot, fixed_episodes, best_wgd_nodes, exactsolution, artificial_wgdroot = count_wgd_nodes_combined(
            network, 
            gene_trees, 
            wgddebug = args.wgddebug,             
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
            fixed_episodes = fixed_episodes,
            fixed_episodes_search = fixed_episodes_search,
            extended_episodes_search = args.extended_episodes_search,           
            extended_episodes_from_fixedepi = args.extended_episodes_from_fixedepi,
            fixed_episodes_only= args.fixed_episodes_only            
            )



    endtime = time.process_time() - t
     
    with open(out_file, "w") as f:
        f.write(f"gene_trees_file={args.gene_trees}\n")        
        f.write(f"network_file={args.network}\n")                
        f.write(f"{outstats}")
        f.write(f"randomize_from={args.randomize_from}\n")
        f.write(f"noimprovement_stop={args.noimprovement_stop}\n")
        f.write(f"time={endtime}")

    if args.locked_epi_support:
        for nnum, v in locked_epi_support.items():
            for n in stroot.nodes():
                if n.num == nnum:
                    n.lockedepisupport = v

    if fixed_episodes: 
        with open(out_basefile + "fixed_episodes","w") as f:
            f.write(wgdnums(fixed_episodes,'',''))
        if args.verbose>=2:                
            print(f"{setid}Fixed episodes stored in fixed_episodes file")

    if args.fixed_episodes_only:
        return


    # Write output network (gse) and fixedwgd's    
    with open(out_basefile + "episummary" + outstyleext,"w") as f:
        f.write(stroot.attrrepr(epiattr, ignorezeros=False, gsestyle=args.gsestyle))
        if args.verbose>=2:
            print(f"{setid}Network with attributes saved in {out_basefile}episummary{outstyleext}")


    if best_wgd_nodes:         
        with open(out_basefile + "best_wgd_nodes" +("exact" if exactsolution else "approx"),"w") as f:
            f.write(wgdnums(best_wgd_nodes,'',''))
        
    if args.verbose:
        print(f"{setid}Cost: {cost} Exact:{exactsolution}")
        if artificial_wgdroot is not None:
            print(f"{setid}Artificial WGD root id:", artificial_wgdroot.num)

        print(f"{setid}Best episodes: {wgdnums(best_wgd_nodes,'','')}")

        if args.save_embedding:                
            print(f"{setid}Episode sizes:", " ".join(f"{b.num}:{b.episize}" for b in best_wgd_nodes if hasattr(b, "episize")))

        if args.extended_episodes_from_fixedepi:                        
            print(f"{setid}Extended episode sizes (from fixedepi):", " ".join(f"{b.num}:{b.eeepisize}" for b in stroot.nodes() if hasattr(b, "eeepisize")))

        if args.extended_episodes_search:                        
            print(f"{setid}Extended episode sizes (post):", " ".join(f"{b.num}:{b.eeepisizepost}" for b in stroot.nodes() if hasattr(b, "eeepisizepost")))

        if args.locked_epi_support:                        
            print(f"{setid}Locked episode support:", " ".join(f"{b.num}:{b.lockedepisupport}" for b in stroot.nodes() if hasattr(b, "lockedepisupport")))

            
        # if args.extended_episodes_search:


        #print(f"{setid}Extended episode stats:",s)        




    if args.verbose>2:
        print("Used nodes: ")
        for node in used_nodes:
            print(node)
    
    assert len(used_nodes) == cost


if __name__ == "__main__":
    main()


