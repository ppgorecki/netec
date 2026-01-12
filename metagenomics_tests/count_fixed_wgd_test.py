import os

from fixedec import fixedec
from treeop import str2tree, Tree
from algorithms import add_outgroup, combine_gene_trees


def main():
    directory = "preprocessed_data/"
    for file in os.listdir(directory):
        if file.startswith("wgd-1"):
            gene_trees_str = open(directory + file).read().split()
            gene_tree = combine_gene_trees([Tree(str2tree(g_str)) for g_str in gene_trees_str])
            species_tree_str = open(directory + "s_tree").read()
            species_tree = add_outgroup(Tree(str2tree(species_tree_str)))
            wgd_nodes = fixedec(gene_tree, species_tree)
            print(f"File: {file}, number of fixed wgd nodes: {len(wgd_nodes)}/{len(species_tree.nodes)}")


if __name__ == "__main__":
    main()
