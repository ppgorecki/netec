#!/usr/bin/python3
import os
import sys
from metatreeop import count_wgd_nodes, is_reconciled_using_wgd
from treeop import Tree,str2tree
from rec import ecfeasbible, metaecfeasible



if __name__ == '__main__':

    if len(sys.argv)!=5:
        print(f"""        
Usage: {sys.argv[0]} SpeciesTree PartialGeneTree Repeats EpisodeCandSize

Repeats - how many times to run DP vs Naive
EpisodeCandSize - the size of candidate episode sets; the set is randomly chosen in each run

Naive vs DP tester

""")
        exit(-1)

    st, gt, repeats, episodecandsize = sys.argv[1:5]

    # gt="((?,(b,((?,?),(?,?)))),(e,((d,c),a)))"
    # st="((e,b),(d,(c,a)))"
    # episodes = [st.root, st.root.c[0].c[0] ]

    st = Tree(str2tree(st))
    gt = Tree(str2tree(gt))
    repeats = int(repeats)
    episodecandsize = int(episodecandsize)

    from random import sample

    tcnt = 0
    
    for i in range(repeats):
        episodes = sample(st.nodes, episodecandsize) # random cand. episodes

        naive = metaecfeasible(gt, st, episodes)
        dp, wgdset, feasgt = is_reconciled_using_wgd(st, gt, episodes)
        if naive!=dp:
            # naive and dp disagree, report
            f = open("err.log","a")
            f.write(f"naive={naive} dp={dp} gt={gt} st={st} episodes={episodes}\n")
            f.close()
            print(f"naive={naive} dp={dp} gt={gt} st={st} episodes={episodes}")            
        elif naive:
            tcnt += 1
    
    print(f"Truecnt={tcnt}/{repeats} episodecandsize={episodecandsize}")


   
