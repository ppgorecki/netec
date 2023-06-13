# MetaEC: Simultaneous reconstruction of duplication episodes and gene-species mappings

Inferring gene-species assignments and genomic duplication events.

See WABI 2023 article:

Simultaneous reconstruction of duplication episodes and gene-species mappings 
Paweł Górecki, Natalia Rutecka, Agnieszka Mykowiecka and Jarosław Paszek

The package originated from the embretnet repository https://bitbucket.org/pgor17/embretnet/src/meta_DP.

### Single file processing 

Run metaec.py for a single inference.


### Multiple files processing 

Run metaec_multiset.sh to perform MetaEC multiset analysis and output results in the specified output directory as a "wgd.csv" file.

Requirements for metaec_multiset.sh:
- csvmanip is downloaded automatically from https://github.com/ppgorecki/csvmanip using git
- gnu parallel

To process yeast dataset:
```
metaec_multiset.sh -R50 -N50 -j10 -o yeast_out -d data_yeast data_yeast/gtrees_0.[0-6]*
```

To process simulated dataset:
```
metaec_multiset.sh -R50 -N100 -j6 -o sim_out -d data_sim data_sim/wgd-[1-5]-*_0*
```

Here: 
- '-R 50' - randomize from binom(n,k)>50 in the main loop
- '-N 100' - stop if DP will not improve after 100 executions
- '-j 10' - run 10 jobs in parallel
- '-o DIR' - store results in directory DIR
- '-d DIR' - input directory 
- `data_sim/wgd-[1-5]-*_0*` - files to process

If the analysis is interrupted, running the command again will resume processing from where it left off.

To get help execute:
```
metaec_multiset.sh
```

### Questions ###

For questions contact:

<gorecki@mimuw.edu.pl>

### Funding

The support was provided by National Science Centre grant #2019/33/B/ST6/00737.
