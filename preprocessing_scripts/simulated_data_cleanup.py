def clean_tree(newick: str) -> str:
    """ Returns newick tree with deleted edge lengths and changed name from X_i_j to X"""
    i = 0
    cleaned = ""
    while i < len(newick):
        if newick[i] == ")":
            cleaned += newick[i]
            if i == len(newick) - 1:
                return cleaned
            i += 1
            while newick[i].isalnum() or newick[i] in ".-+:":
                i += 1
        elif newick[i] == ":":
            while newick[i].isnumeric() or newick[i] in ".-+:":
                i += 1

        elif newick[i] == "_":
            while newick[i].isnumeric() or newick[i] in "_":
                i += 1
        else:
            cleaned += newick[i]
            i += 1
    return cleaned


for wgd in range(1, 6):
    with open(f"preprocessed_data/wgd-{wgd}-gene-trees", "w") as f:
        for gene_tree_nr in range(1, 101):
            missing_digits = 3 - len(str(gene_tree_nr))
            nr_prefix = missing_digits * "0"
            gene_tree = open(
                f"wgd-simulated/n20-wgd-{wgd}/1/g_trees{nr_prefix}{gene_tree_nr}.trees").read().strip().strip(";")
            cleaned_gene_tree = clean_tree(gene_tree)
            end = "\n" if gene_tree_nr != 100 else ""
            f.write(f"{cleaned_gene_tree}{end}")
with open(f"../preprocessed_data/s_tree", "w") as f:
    stree = open("wgd-simulated/s_tree.trees").read().strip().strip(";")
    cleaned_stree = clean_tree(stree)
    f.write(cleaned_stree)
