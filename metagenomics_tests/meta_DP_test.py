#!/usr/bin/python3 

import argparse
from random import randint

from algorithms import count_wgd_nodes_combined
from treeop import str2tree, randtreestr, Tree


def main():
    parser = argparse.ArgumentParser(
        description="Run WGD reconciliation algorithm for random gene trees with ? and random species trees")
    parser.add_argument("--min_leaf", help="Minimum number of leaves for species tree", type=int, default=3)
    parser.add_argument("--max_leaf", help="Maximum number of leaves for species tree", type=int, default=20)
    parser.add_argument("--min_unknown", help="Maximum number of ? leaves for gene tree", type=int, default=1)
    parser.add_argument("--max_unknown", help="Maximum number of ? leaves for gene tree", type=int, default=3)
    parser.add_argument("--gene_tree_nr", help="Number of gene trees", type=int, default=5)
    args = parser.parse_args()
    if args.min_leaf <= 1:
        print("Min_leaf has to be at least 2")
        exit(1)
    if args.min_leaf > args.max_leaf:
        print("Min_leaf is bigger than max_leaf")
        exit(1)
    while True:
        leaves_s = randint(args.min_leaf, args.max_leaf)
        species_tree = Tree(str2tree(randtreestr(leaves_s, use_numbers=True)))
        print(f"Species tree: {species_tree}")
        gene_trees = []
        for _ in range(args.gene_tree_nr):
            unknown_leaves = randint(args.min_unknown, args.max_unknown)
            leaves_g = randint(3, leaves_s)
            gene_tree = Tree(str2tree(randtreestr(leaves_g, use_numbers=True, unknown_labels=unknown_leaves)))
            gene_trees.append(gene_tree)
        print(f"Gene trees: {gene_trees}")
        cost, used_nodes = count_wgd_nodes_combined(species_tree, gene_trees)
        print(f"Cost: {cost}")
        print("Used nodes: ")
        for node in used_nodes:
            print(node)
        assert len(used_nodes) == cost


if __name__ == "__main__":
    main()
