# NetEC: Duplication Episodes in Phylogenetic Networks

NetEC is a tool for inferring whole genome duplication (WGD) episodes by reconciling gene trees with species trees or networks. It uses dynamic programming to find the minimum number of duplication episodes that explain the observed gene tree duplications.

## Requirements

- Python 3.x
- No external dependencies required

## Input File Formats

### Species Tree / Network

A Newick-formatted file containing the species tree or network topology.

**Species tree example** (`data_sim/s_tree`):
```
((a,b),(c,d))
```

**Network with reticulations** - use `#` labels to mark reticulation nodes:
```
((c)#A,((#A,b),a))
```
This represents a network where node `#A` has two parents. The reticulation is defined by matching internal node `(c)#A` with leaf `#A`.

### Gene Trees

A text file with one gene tree per line in Newick format. Gene labels should match species labels.

**Example** (`gene_trees.txt`):
```
((a,a),(b,c))
((a,b),(c,d))
(((a,a),b),(c,d))
```

## Basic Usage

### Minimal Example

```bash
python3 metaec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree
```

**Output:**
```
[] Cost: 1 Exact:True
```
This shows that 1 duplication episode explains all gene tree duplications.

### With Output File

```bash
mkdir -p results
python3 metaec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --out_file results/output.log
```

### With Output Dir (recommended)

```bash
mkdir -p results
python3 metaec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --out_file results
```

### With Verbose Output

```bash
python3 metaec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --verbose 2
```

Use `--verbose 2` to print the WGD nodes used in the solution.

## Working with Real Data (Yeast Dataset)

```bash
mkdir -p results
python3 metaec.py \
    --gene_trees data_yeast/gtrees \
    --network data_yeast/s_tree \
    --out_file results
```

## Working with Networks (Reticulations)

When using phylogenetic networks with hybridization events:

```bash
# Create network file
echo "((c)#A,((#A,b),a))" > network.txt

# Create gene trees file
echo "((a,a),(b,c))" > gtrees.txt
echo "((a,b),(c,c))" >> gtrees.txt

# Run analysis
python3 metaec.py \
    --gene_trees gtrees.txt \
    --network network.txt
```

## Advanced Options

### Filtering Gene Trees by Duplication Count

Filter out gene trees with fewer than N duplications:

```bash
python3 metaec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --mindup 2
```

Note: `--mindup` only works for trees (networks without reticulations).

### View Duplication Statistics

```bash
python3 metaec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --print_dup_stats
```

**Output:**
```
Duplication statistics:
  0 duplications: 15 gene tree(s)
  1 duplications: 42 gene tree(s)
  2 duplications: 28 gene tree(s)
  ...
```

### User-Defined Episodes

Specify which nodes to consider as potential duplication episodes:

```bash
# Consider specific nodes (by node numbers)
python3 metaec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --user_episodes "2 4 10"

# Consider all nodes
python3 metaec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --user_episodes all
```

### Save Episode Summary

Generate a species tree/network with episode size attributes:

```bash
mkdir -p results
python3 metaec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --out_file results \
    --episummaryfile
cat results/episummary
```

### Save Embedding

Save the inferred gene-species mapping:

```bash
mkdir -p results
python3 metaec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --out_file results \
    --save_embedding
cat results/embedding
```

### Distribution Maps

Add distribution maps to the output:

```bash
python3 metaec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --distribution_maps \
    --print_distr_maps
```

### Performance Optimization

For large datasets, use randomization to speed up computation:

```bash
python3 metaec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --randomize_from 1000 \
    --noimprovement_stop 50
```

- `--randomize_from N`: Start randomizing when binomial coefficient exceeds N
- `--noimprovement_stop N`: Stop after N iterations with no improvement (0 = don't stop)

### Reversed Climb Search

Alternative search strategy starting from fixed WGD and expanding:

```bash
python3 metaec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --reversed_climb
```

### GSE-Style Output

Use GSE (Gene tree Species tree Embedding) format for output attributes:

```bash
python3 metaec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --gsestyle
```

## Complete Example Workflow

```bash
# Create output directory
mkdir -p results/yeast

# 1. Analyze gene trees from simulation data
python3 metaec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --out_file results/sim_analysis.log \
    --print_dup_stats \
    --verbose 1

# 2. Analyze with episode summary output
python3 metaec.py \
    --gene_trees data_yeast/gtrees \
    --network data_yeast/s_tree \
    --out_file results/yeast/ \
    --episummaryfile \
    --save_embedding

# 3. Quick analysis filtering low-duplication trees
python3 metaec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --mindup 1 \
    --verbose 1
```

## Output Files

When using `--out_file`, the following files may be generated:

- `<out_file>` or `<out_dir>/metaec.log` - Main log with parameters and results
- `<out_basefile>.embedding` - Gene-species embedding (with `--save_embedding`)
- `<out_basefile>episummary` - Species tree with episode attributes (with `--episummaryfile`)

## Batch Processing

For processing multiple gene tree files, see the example at:
https://bitbucket.org/pgor17/pandanales

Example shell script pattern:
```bash
mkdir -p results
for file in data_sim/wgd-1-gene-trees_*; do
    python3 metaec.py \
        --gene_trees "$file" \
        --network data_sim/s_tree \
        --out_file "results/$(basename $file).log"
done
```

## Questions

For questions contact: <gorecki@mimuw.edu.pl>

## Funding

This work was supported by National Science Centre grants:
- #2019/33/B/ST6/00737 (tree version)
- #2023/51/B/ST6/02792 (network)
