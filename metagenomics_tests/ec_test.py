import random
from random import randint

from fixedec import fixedec
from metatreeop import combine_gene_trees, add_outgroup, is_reconciled_using_wgd
from rec import rec
from treeop import str2tree, randtreestr, Tree

min_leaf = 5
max_leaf = 5
min_unknown = 1
max_unknown = 1
gene_tree_nr = 1
wgd_nodes_size = 2*min_leaf + 1

while True:
    leaves_s = randint(min_leaf, max_leaf)
    species_tree = Tree(str2tree(randtreestr(leaves_s, use_numbers=True)))
    gene_trees = []
    for _ in range(gene_tree_nr):
        unknown_leaves = randint(min_unknown, max_unknown)
        leaves_g = randint(3, leaves_s)
        gene_tree = Tree(str2tree(randtreestr(leaves_g, use_numbers=True, unknown_labels=unknown_leaves)))
        gene_trees.append(gene_tree)
    species_tree = add_outgroup(species_tree)
    gene_tree = combine_gene_trees(gene_trees)
    fixed_wgd_nodes = fixedec(gene_tree, species_tree)
    wgd_nodes = random.sample(list(set(species_tree.nodes) - fixed_wgd_nodes), wgd_nodes_size - len(fixed_wgd_nodes))
    is_feasible, used_nodes, inferred_tree = is_reconciled_using_wgd(species_tree, gene_tree, set(wgd_nodes) | fixed_wgd_nodes)
    if not is_feasible:
        continue
    else:
        ec, ec_nodes = rec([Tree(str2tree(inferred_tree))], species_tree)
        if not ec == len(used_nodes):
            print("-----NEW EXAMPLE-----")
            print(f"Species tree: {species_tree}")
            print(f"Gene tree: {gene_tree}")
            print(f"Gene tree inferred: {inferred_tree}")
            print(f"Fixed WGD nodes:")
            print(*fixed_wgd_nodes, sep="\n")
            print(f"Unfixed nodes used by DP:")
            print(*used_nodes, sep="\n")
            print(f"Unfixed nodes used by EC:")
            print(*ec_nodes, sep="\n")
