#!/usr/bin/python3
import os
import sys
from algorithms import count_wgd_nodes, is_reconciled_using_wgd_withcounts, is_reconciled_using_wgd
from treeop import Tree,str2tree
from rec import ecfeasbible, netecfeasible


def compd(d1, d2, slf):
    res = ""
    dif = False
    for s in slf:
        l = s.clusterleaf
        v1 = d1.get(l,0)
        v2 = d2.get(l,0)
        if v1==v2:
            res+=f"{l}:{v1} "
        else:
            res+=f"{l}:{v1},{v2} "
            dif = True
    return res, dif


if __name__ == '__main__':

    if len(sys.argv)<5:
        print(f"""        
Usage: {sys.argv[0]} SpeciesTree PartialGeneTree Repeats EpisodeCandSize EpiNums

Repeats - how many times to run DP vs Naive

EpisodeCandSize - the size of candidate episode sets; the set is randomly chosen in each run

Naive vs DP tester

-H - optional, gen heat maps

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

    epinums = []
    if len(sys.argv)>5:
        epinums = list(map(int, sys.argv[5:]))

    from random import sample

    tcnt = 0
    slf = st.leaves()    
    wgddebug = False
    
    errcnt = 0
    for i in range(repeats):
        if episodecandsize>0:
            episodes = sample(st.nodes, episodecandsize) # random cand. episodes
        else:
            episodes = [ st.nodes[i] for i in epinums ]

        #print (",".join(str(e.num)+str(e) for e in episodes))        
        naive, distrcnt, leafmap = netecfeasible(gt, st, episodes, leafdistr = True)        
        dp, dpdistr = is_reconciled_using_wgd_withcounts(st, gt, episodes,  wgddebug=wgddebug)        

        if naive!=dp: #or dp!=dp2:
            # naive and dp disagree, report
            print ("ERR!")
            f = open("err.log","a")
            f.write(f"naive={naive} dp={dp} gt={gt} st={st} episodes={episodes}\n")
            f.close()
            print(f"naive={naive} dp={dp} gt={gt} st={st} episodes={episodes}")            
            errcnt += 1 
        elif naive:
            tcnt += 1
            if wgddebug: 
                print (f"NAIVE vs DP {naive} {dp}")
        if dp == naive and dp: 
            difx = False
            if wgddebug:
                for i in gt.nodes:
                    print(i,i.num,[c.num for c in i.c])
            for i,(guleaf,dv) in enumerate(dpdistr.items()):
                dstr, dif = compd(leafmap[guleaf], dv, slf)
                if wgddebug:
                    print (f"?={i}.{guleaf.num} {dstr}")
                difx = difx or dif            
            if difx: 
                print("Bad counts! See err.log for details")
                f = open("err.log","a")
                for i in gt.nodes:
                    f.write(f"{i},{i.num},{[c.num for c in i.c]}")
                    for i,(guleaf,dv) in enumerate(dpdistr.items()):
                        dstr, dif = compd(leafmap[guleaf], dv, slf)                    
                        f.write(f"?={i}.{guleaf.num} {dstr}")
                f.close()                
                break        
    
    print(f"Truecnt={tcnt}/{repeats} episodecandsize={episodecandsize}")

    if errcnt:
        print("Check err.log for details")
        sys.exit(-1)


   
