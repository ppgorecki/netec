# NetEC: Duplication Episodes in Phylogenetic Networks

NetEC is a tool for inferring whole genome duplication (WGD) episodes by reconciling gene trees with species trees or phylogenetic networks. It uses dynamic programming to find the minimum number of duplication episodes that explain the observed gene tree duplications.

Article (submitted): Górecki, P., Rutecka, N., Mykowiecka, A., Paszek, J., Episode Clustering in Phylogenetic Networks

NetEC is partially based on the MetaEC tool https://bitbucket.org/pgor17/metaec for the inference of WGD events, jointly with missing gene-species assignments, but only for species trees. MetaEC references:
- Górecki, P., Rutecka, N., Mykowiecka, A. et al. Unifying duplication episode clustering and gene-species mapping inference. Algorithms Mol Biol 19, 7 (2024). https://doi.org/10.1186/s13015-024-00252-8
- Górecki, P., Rutecka, N., Mykowiecka, A., Paszek, J., Simultaneous reconstruction of duplication episodes and gene-species mappings, WABI 2023.

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
python3 netec.py \
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
python3 netec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --out_file results/output.log
```

### With Output Dir (recommended)

```bash
mkdir -p results
python3 netec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --out_file results
```

### With Verbose Output

```bash
python3 netec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --verbose 2
```

Use `--verbose 2` to print the WGD nodes used in the solution.

## Working with Real Data (Yeast Dataset)

```bash
mkdir -p results
python3 netec.py \
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
python3 netec.py \
    --gene_trees gtrees.txt \
    --network network.txt
```

## Advanced Options

### Filtering Gene Trees by Duplication Count

Filter out gene trees with fewer than N duplications:

```bash
python3 netec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --mindup 2
```

Note: `--mindup` only works for trees (networks without reticulations).

### View Duplication Statistics

```bash
python3 netec.py \
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
python3 netec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --user_episodes "2 4 10"

# Consider all nodes
python3 netec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --user_episodes all
```

### Save Episode Summary

Generate a species tree/network with episode size attributes:

```bash
mkdir -p results
python3 netec.py \
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
python3 netec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --out_file results \
    --save_embedding
cat results/embedding
```

### Extended Episode Analysis

Identify non-fixed episode nodes that are good candidate to be episode.

```bash
python3 netec.py \  
    --gene_trees example_extended_episodes/gene_trees  \
    --network example_extended_episodes/network  \    
    --extended_episodes_search
```

Note that some the candidatte episode combinations have no feasible solution. Such a sitution is rare in practise.

### Locked Episode Support Analysis

Identify which network nodes are required (locked) for reconciling each individual gene tree:

```bash
mkdir -p results
python3 netec.py \
    --gene_trees example_locked_epi/gene_trees  \
    --network example_locked_epi/species_tree  \
    --locked_epi_support \
    --out_file results
```

This generates two output files:

**File: `results/locked_epi`** - Lists locked episode node numbers for each gene tree (one line per gene tree):
```
2
5
3 2

2
0
```

**File: `results/locked_epi_net.newick`** - The species network annotated with support counts:
```
((A[num=2;lockedepisupport=3],B[num=3;lockedepisupport=1])[num=1],(C[num=5;lockedepisupport=1],D[num=6])[num=4])[num=0;lockedepisupport=1]
```

**Interpretation:**
- `lockedepisupport=3` means 3 gene trees require a WGD at that node
- Nodes with high support are strong candidates for WGD events
- Each line in `locked_epi` shows which nodes are mandatory for that specific gene tree

**Example Use Cases:**
- Identify consensus WGD locations across gene trees
- Find outlier gene trees with unusual duplication patterns
- Understand which gene trees impose strict constraints on WGD placement
- Low lockedepisupport may suggest removal gene trees that induce support at that node

### Distribution Maps

Add distribution maps to the output:

```bash
python3 netec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --distribution_maps \
    --print_distr_maps
```

Note: use with species trees only (no reticulations in networks).

### Performance Optimization

For large datasets, use randomization to speed up computation:

```bash
python3 netec.py \
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
python3 netec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --reversed_climb
```

### GSE-Style Output

Use GSE format for output attributes:

```bash
python3 netec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --out_file results
    --gsestyle
```

## Complete Example Workflow

```bash
# Create output directory
mkdir -p results/sim results/yeast

# 1. Analyze gene trees from simulation data
python3 netec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --out_file results/sim \
    --print_dup_stats \
    --verbose 1

# 2. Analyze with episode summary output
python3 netec.py \
    --gene_trees data_yeast/gtrees \
    --network data_yeast/s_tree \
    --out_file results/yeast/ \
    --episummaryfile \
    --save_embedding

# 3. Filtering of low-duplication trees (for non-reticulation networks)
for mindup in 1 2 3 4; do
python3 netec.py \
    --gene_trees data_sim/wgd-1-gene-trees \
    --network data_sim/s_tree \
    --mindup $mindup \
    --verbose 1 \
    --out_file results/sim/mindup"$mindup".log
done
```

## Output Files

When using `--out_file`, the following files may be generated:

- `<out_file>` or `<out_dir>/netec.log` - Main log with parameters and results
- `<out_basefile>.embedding` - Gene-species embedding (with `--save_embedding`)
- `<out_basefile>episummary` - Species tree with episode attributes (with `--episummaryfile`)

## Batch Processing

For processing multiple gene tree files, see the example at:
https://bitbucket.org/pgor17/pandanales

Example shell script pattern:
```bash
for file in data_sim/wgd-1-gene-trees_*; do
	DIR="results/$(basename $file)"
	mkdir -p $DIR
    python3 netec.py \
        --gene_trees "$file" \
        --network data_sim/s_tree \
        --out_file $DIR
done
```

## Fasturec Integration

[Fasturec](https://bitbucket.org/pgor17/fasturec) is a tool for inferring species trees from gene trees using duplication cost minimization. You can use Fasturec to infer a species tree, then analyze duplication episodes with NetEC.

> **Note:** This workflow requires a species tree (no reticulations/networks).

### Installation

```bash
git clone git@bitbucket.org:pgor17/fasturec.git
cd fasturec && make
```

### Workflow

**1. Infer species tree from gene trees:**

```bash
fasturec -q10 -Y -G gt.txt -oft
cut -f1 fu.txt > st.txt
```

- `-q10`: Quick search with 10 random restarts
- `-Y`: Use Yale (YDC) duplication cost
- `-G gt.txt`: Input gene trees file
- `-oft`: Output full tree info to `fu.txt`
- The `cut` command extracts the species tree from Fasturec output

**2. View duplication statistics:**

```bash
python3 netec.py \
    --network st.txt \
    --gene_trees gt.txt \
    --print_dup_stats
```

**3. Run episode analysis with filtered trees:**

```bash
python3 netec.py \
    --network st.txt \
    --gene_trees gt.txt \
    --mindup 4
```

Using `--mindup 4` filters out gene trees with fewer than 4 duplications, focusing the analysis on trees with more duplication signal.

## Questions

For questions contact: <gorecki@mimuw.edu.pl>

## Funding

This work was supported by National Science Centre grant:
- #2023/51/B/ST6/02792 (networks)
