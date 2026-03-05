# NetEC: Duplication Episodes in Phylogenetic Networks

NetEC is a tool for inferring whole genome duplication (WGD) episodes by reconciling gene trees with species trees or phylogenetic networks. It uses dynamic programming to find the minimum number of duplication episodes that explain the observed gene tree duplications.

Article (submitted): Górecki, P., Rutecka, N., Mykowiecka, A., Paszek, J., Episode Clustering in Phylogenetic Networks

For examples and datasets from the above article, refer to:
- https://github.com/ppgorecki/NetEC-Pandanales.git
- https://github.com/ppgorecki/NetEC-Simulations

NetEC is partially based on the MetaEC tool https://bitbucket.org/pgor17/metaec for the inference of WGD events, jointly with missing gene-species assignments, but only for species trees. MetaEC references:
- Górecki, P., Rutecka, N., Mykowiecka, A. et al. Unifying duplication episode clustering and gene-species mapping inference. Algorithms Mol Biol 19, 7 (2024). https://doi.org/10.1186/s13015-024-00252-8
- Górecki, P., Rutecka, N., Mykowiecka, A., Paszek, J., Simultaneous reconstruction of duplication episodes and gene-species mappings, WABI 2023.

Financial support was provided by the NCN grant 2023/51/B/ST6/02792.

## Requirements

- Python 3.x
- No external dependencies required

## Input File Formats

### Species Tree / Network

A Newick-formatted file containing the species tree or network topology.

**Species tree example:**
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
    --gene_trees example1/gene_trees \
    --network example1/network
```

**Output:**
```
Cost: 1 Exact:True
Best episodes: 2
```
This shows that 1 duplication episode explains all gene tree duplications.

The result is saved in results dir, it will be created if not exists.

### With Non-default Output Directory

```bash
python3 netec.py \
    --gene_trees example1/gene_trees \
    --network example1/network \
    --out_dir smp
```

Results are saved to the specified directory (created if needed). Default output directory is `smp/`.

### With Verbose Output

```bash
python3 netec.py \
    --gene_trees example1/gene_trees \
    --network example1/network \
    --verbose 2
```

Use `--verbose 2` to print more detailed info.

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
    --gene_trees example1/gene_trees \
    --network example1/network \
    --mindup 2
```

Note: `--mindup` only works for trees (networks without reticulations).

### View Duplication Statistics

```bash
python3 netec.py \
    --gene_trees example1/gene_trees \
    --network example1/network \
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

Specify which nodes to consider as potential duplication episodes. 
Program adds needed fixed episodes, but full optimization is not executed (only single optimize mode - one DP run).

```bash
# Consider specific nodes (by node numbers)
python3 netec.py \
    --gene_trees example1/gene_trees \
    --network example1/network \
    --user_episodes "2 4 10"

# Consider all nodes
python3 netec.py \
    --gene_trees example1/gene_trees \
    --network example1/network \
    --user_episodes all
```

### Fixed Episodes

Provide precomputed fixed episodes to optimize computations:

```bash
python3 netec.py \
    --gene_trees example1/gene_trees \
    --network example1/network \
    --fixed_episodes "2 4" \
    --no_fixed_episodes_search \
    --verbose 2
```

Use `--no_fixed_episodes_search` to skip automatic fixed episode detection when providing known episodes.

### Save Embedding

Save the inferred gene-species mapping:

```bash
python3 netec.py \
    --gene_trees example1/gene_trees \
    --network example1/network \
    --out_dir results \
    --save_embedding
cat results/embedding
```

Episode summary is automatically saved to `<out_dir>/episummary.newick` (or `.gse` with `--gsestyle`).

Important: always use --save_embedding if the size of episodes is needed (episize attribute).

### Extended Episode Analysis

Identify non-fixed episode nodes that are good candidates to be episodes.

```bash
python3 netec.py \
    --gene_trees example2/gene_trees \
    --network example2/network.newick \
    --extended_episodes_search
```

Extended episode sizes are saved as `eeepisizepost` attribute.

Use `--extended_episodes_from_fixedepi` to identify non-fixed episodes with large number of duplications (saved as `eeepisize` attribute).

Note that some candidate episode combinations have no feasible solution. Such a situation is rare in practice.

### Locked Episode Support Analysis

Identify which network nodes are required (locked) for reconciling each individual gene tree:

```bash
python3 netec.py \
    --gene_trees example1/gene_trees \
    --network example1/network \
    --locked_epi_support \
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
    --gene_trees example1/gene_trees \
    --network example1/network \
    --distribution_maps \
    --print_distr_maps
```

Use `--distr_counts` to output raw counts instead of normalized values.

Note: use with species trees only (no reticulations in networks).

TODO: implement with networks.

### Performance Optimization

For large datasets, use randomization to speed up computation:

```bash
python3 netec.py \
    --gene_trees example1/gene_trees \
    --network example1/network \
    --randomize_from 1000 \
    --noimprovement_stop 50
```

- `--randomize_from N`: Start randomizing when binomial coefficient exceeds N
- `--noimprovement_stop N`: Stop after N iterations with no improvement (0 = don't stop)

### Reversed Climb Search

Alternative search strategy starting from fixed WGD and expanding:

```bash
python3 netec.py \
    --gene_trees example1/gene_trees \
    --network example1/network \
    --reversed_climb
```

### Optimization Modes

Control the optimization algorithm behavior:

```bash
# Full optimization (default)
python3 netec.py \
    --gene_trees example1/gene_trees \
    --network example1/network \
    --optimize full

# Single DP run on fixed + user episodes
python3 netec.py \
    --gene_trees example1/gene_trees \
    --network example1/network \
    --optimize single

# Skip DP algorithms entirely
python3 netec.py \
    --gene_trees example1/gene_trees \
    --network example1/network \
    --optimize no
```

Use `--fixed_episodes_only` to stop after the fixed episodes phase.

### GSE-Style Output

Use GSE format for output attributes:

```bash
python3 netec.py \
    --gene_trees example1/gene_trees \
    --network example1/network \    
    --gsestyle
```

## Complete Example Workflow

This section demonstrates a typical analysis workflow. For a complete real-world example, see the [Pandanales dataset](https://github.com/ppgorecki/NetEC-Pandanales.git).

### Step 1: Identify Fixed Episodes

Fixed episodes are WGD events that must occur at specific nodes to explain the gene tree duplications. Identifying them first optimizes all subsequent runs.

```bash
python3 netec.py \
    --network network.nwk \
    --gene_trees gene_trees.nwk \
    --out_dir fixed \
    --fixed_episodes_only \
    --verbose 2
```

### Step 2: Discovery Mode

#### Infer Episodes

Run after fixed episodes are identified:

```bash
python3 netec.py \
    --network network.nwk \
    --gene_trees gene_trees.nwk \
    --out_dir disco \
    --save_embedding \
    --no_fixed_episodes_search \
    --fixed_episodes "$(cat fixed/fixed_episodes)"
```

#### Search for Extended Episodes

Identify additional candidate WGD nodes beyond the fixed episodes:

```bash
python3 netec.py \
    --network network.nwk \
    --gene_trees gene_trees.nwk \
    --out_dir discoext \
    --extended_episodes_from_fixedepi \
    --save_embedding \
    --no_fixed_episodes_search \
    --fixed_episodes "$(cat disco/best_wgd_nodesexact)" \
    --optimize single
```

#### Extended Episodes Inference

Add specific candidate episodes (e.g., node 21) to the analysis:

```bash
python3 netec.py \
    --network network.nwk \
    --gene_trees gene_trees.nwk \
    --out_dir discoext21 \
    --locked_epi_support \
    --extended_episodes_from_fixedepi \
    --save_embedding \
    --no_fixed_episodes_search \
    --fixed_episodes "$(cat disco/best_wgd_nodesexact) 21" \
    --optimize single \
    --verbose 2
```

To test multiple candidate nodes (e.g., 19, 21, and both), run in parallel:

```bash
parallel python3 netec.py \
    --network network.nwk \
    --gene_trees gene_trees.nwk \
    --out_dir discoext"{1}" \
    --extended_episodes_from_fixedepi \
    --save_embedding \
    --no_fixed_episodes_search \
    --fixed_episodes "$(cat disco/best_wgd_nodesexact) {1}" \
    --optimize single \
    --verbose 2 \
    ::: 19 21 "19 21"
```

### Step 3: Hypothesis-Driven Mode

When you have prior hypotheses about WGD locations (e.g., from literature), test them directly.

#### Define WGD Hypotheses

```bash
# Example: hypothesized WGD nodes from prior studies
WGD="3 4 12 19 28"
```

#### Infer Episodes with Hypotheses

```bash
python3 netec.py \
    --network network.nwk \
    --gene_trees gene_trees.nwk \
    --out_dir hypo \
    --save_embedding \
    --no_fixed_episodes_search \
    --fixed_episodes "$(cat fixed/fixed_episodes) $WGD"
```

#### Search for Extended Episodes

```bash
python3 netec.py \
    --network network.nwk \
    --gene_trees gene_trees.nwk \
    --out_dir hypoext \
    --extended_episodes_from_fixedepi \
    --save_embedding \
    --no_fixed_episodes_search \
    --fixed_episodes "$(cat hypo/best_wgd_nodesexact)" \
    --optimize single
```

#### Extended Episodes Inference

Add additional candidate episodes to the hypothesis:

```bash
python3 netec.py \
    --network network.nwk \
    --gene_trees gene_trees.nwk \
    --out_dir hypoext21 \
    --locked_epi_support \
    --extended_episodes_from_fixedepi \
    --save_embedding \
    --no_fixed_episodes_search \
    --fixed_episodes "$(cat hypo/best_wgd_nodesexact) 21" \
    --optimize single \
    --verbose 2
```

## Output Files

When using `--out_dir`, the following files are generated in the specified directory:

- `netec.log` - Main log with parameters and results
- `episummary.newick` (or `episummary.gse` with `--gsestyle`) - Species tree/network with episode attributes (always generated)
- `embedding` - Gene-species embedding (with `--save_embedding`)
- `fixed_episodes` - List of fixed episode node IDs
- `best_wgd_nodesexact` or `best_wgd_nodesapprox` - Best WGD nodes (exact or approximate solution)
- `locked_epi` - Locked episodes per gene tree (with `--locked_epi_support`)
- `locked_epi_net.newick` - Network with locked episode support attributes (with `--locked_epi_support`)


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
    --mindup 4 \
    --out_dir results/fasturec
```

Using `--mindup 4` filters out gene trees with fewer than 4 duplications, focusing the analysis on trees with more duplication signal.

## Questions

For questions contact: <gorecki@mimuw.edu.pl>

## Funding

This work was supported by National Science Centre grant:
- #2023/51/B/ST6/02792 (networks)
