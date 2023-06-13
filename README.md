# MetaEC: Simultaneous reconstruction of duplication episodes and gene-species mappings

The software package, which is based on the embretnet repository, is currently available as a branch of embretnet at \url{https://bitbucket.org/pgor17/embretnet/src/meta_DP}. We are planning to detach the repository and develop it as a separate project in the nearest future.


### Usage ###

Run metaec.py for a single inference.

For processing multisets run metaec_multiset.sh

Requirements for metaec_multiset.sh:
* csvmanip (https://github.com/ppgorecki/csvmanip)
* gnu parallel

To process yeast dataset:
```
metaec_multiset.sh -R50 -N50 -j10 -o yeast_out -d data_yeast data_yeast/gtrees_0.[0-6]*
```

To rocess simulated dataset:
```
metaec_multiset.sh -R50 -N100 -j60 -o sim_out data_sim/wgd--_0*
```

### Questions ###

For questions, please contact:

<gorecki@mimuw.edu.pl>