#!/usr/bin/python3
import os
import sys
from metatreeop import count_wgd_nodes, is_reconciled_using_wgd_withcounts, is_reconciled_using_wgd
from treeop import Tree,str2tree
from rec import ecfeasbible, metaecfeasible



if __name__ == '__main__':

    if len(sys.argv)<5:
        print(f"""        
Usage: {sys.argv[0]} SpeciesTree PartialGeneTree Repeats EpisodeCandSize [-H]

Repeats - how many times to run DP vs Naive
EpisodeCandSize - the size of candidate episode sets; the set is randomly chosen in each run

Naive vs DP tester

-H - optional, gen heat maps

""")
        exit(-1)
    
    st, gt, repeats, episodecandsize = sys.argv[1:5]
    leafdistr = len(sys.argv)>5
         
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

        print (",".join(str(e.num) for e in episodes))

        if leafdistr:
            naive, distrcnt, leafmap = metaecfeasible(gt, st, episodes, leafdistr = leafdistr)
            print(leafmap)
        else: 
            naive, _, _ = metaecfeasible(gt, st, episodes)

        dp = is_reconciled_using_wgd_withcounts(st, gt, episodes)
        dp2,_,_ = is_reconciled_using_wgd(st, gt, episodes)
        if naive!=dp or dp!=dp2:
            # naive and dp disagree, report
            f = open("err.log","a")
            f.write(f"naive={naive} dp={dp} dp2={dp2} gt={gt} st={st} episodes={episodes}\n")
            f.close()
            print(f"naive={naive} dp={dp} dp2={dp2} gt={gt} st={st} episodes={episodes}")            
        elif naive:
            tcnt += 1
    
    print(f"Truecnt={tcnt}/{repeats} episodecandsize={episodecandsize}")


   
