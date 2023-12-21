import os
import subprocess


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


simdir = "genomicduplicationilp/wgd-simulated"

with open(f"{simdir}/s_tree", "w") as f:
    stree = open(f"{simdir}/s_tree.trees").read().strip().strip(";")
    cleaned_stree = clean_tree(stree)
    f.write(cleaned_stree)

os.mkdir("wgd-summary")
for wgd in range(0, 6):
    with open(f"wgd-summary/wgd-{wgd}-gene-trees", "w") as f, open(f"wgd-summary/wgd-{wgd}-rooted-gene-trees", "w") as rf:
        for gene_tree_nr in range(1, 101):
            missing_digits = 3 - len(str(gene_tree_nr))
            nr_prefix = missing_digits * "0"
            gene_tree_path = f"{simdir}/n20-wgd-{wgd}/1/dataset_{nr_prefix}{gene_tree_nr}_TRUE.phy_tree_ML"
            gene_tree = open(gene_tree_path).read().strip().strip(";")
            cleaned_gene_tree = clean_tree(gene_tree)
            f.write(f"{cleaned_gene_tree}\n")

            # root gene trees
            rooting_query = f'urec -g "{cleaned_gene_tree}" -s "{stree}" -um -rmp1'
            try:
                process = subprocess.Popen(rooting_query, stdout=subprocess.PIPE, shell=True)
                (output, error) = process.communicate()
                rooted_tree = output.decode().split()[0]
                rf.write(f"{rooted_tree}\n")
            except:
                rf.write(f"Rooting error: {error}")
