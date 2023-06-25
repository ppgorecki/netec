from treeop import Tree, str2tree
import random
import os


def clear_taxa(g: Tree, p: float) -> Tree:
    """ Randomly choose whether to replace each leaf of the input tree with ? independently with the given probability"""
    for leaf in g.leaves():
        if random.random() <= p:
            leaf.label = "?"
    return g


thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
repeats = 10
tree_files = ["../data_yeast/gtrees"] + [f"../data_sim/wgd-{i}-gene-trees" for i in [1, 2, 3, 4, 5]]
for threshold in thresholds:
    for i in range(repeats):
        for tree_file in tree_files:
            fn = f"{tree_file}_{threshold}_{i + 1}"
            if os.path.exists(fn):
                print(f"File {fn} exists (skip)")
                continue
            with open(fn, "w") as out:
                for line in open(tree_file).readlines():
                    tree = Tree(str2tree(line.strip()))
                    cleared = str(clear_taxa(tree, threshold))
                    out.write(f"{cleared}\n")
