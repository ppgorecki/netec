import itertools
import math
from treeop import Tree, str2tree
from leafdistr import emptyleafcounts, leafcountsupdate

# return left child
def get_left(n):
    return n.c[0]  # n.c[0]  # n.l


# return right child
def get_right(n):
    return n.c[1]  # n.c[1]  # n.r


# return True if 'n' is a duplication
def is_dup(n):
    return n.lcamap == get_left(n).lcamap or n.lcamap == get_right(n).lcamap


def expandintervals(gt, st):
    for n in gt.nodes:
        if n.interval:
            c = n.parent
            while c and c.interval:
                c = c.parent
            if not c:
                n.interval[1] = st.root
                continue
            top = n.interval[0]
            while top.parent != c.lcamap:
                top = top.parent
            n.interval[1] = top


# generates intervals according to PG model
def generate_intervals(gt, st):
    for n in gt.nodes:
        n.interval = None
        if n.leaf():
            continue
        if is_dup(n):
            n.interval = [n.lcamap, None]
    expandintervals(gt, st)


# generates intervals for gene tree nodes
def init(gtrees, st):
    for i, gt in enumerate(gtrees):
        gt.set_lca_mapping(st)
        generate_intervals(gt, st)
        # set ids of gene nodes
        gt.gtid = i
        for g in gt.nodes:
            g.gtid = i
            

# input : set of gene trees (gtrees), species tree (st)
# output : minimal EC score, list of nodes with episodes for that score
def rec(gtrees, st):
    init(gtrees, st)
    dup = []  # duplication nodes == nodes with defined interval
    stnodespostorder = st.root.get_nodes()
    gtnodespostorder = list(itertools.chain.from_iterable(gt.root.get_nodes() for gt in gtrees))

    # limit gene tree nodes only to duplications
    for g in gtnodespostorder:
        if g.interval:
            dup.append(g)

    for s in stnodespostorder:
        s.topinterval = []
        s.botinterval = []
        s.allintervals = []
        s.active = False
        s.min_len_int = math.inf
        s.episize = 0
        s.epigtcount = 0

    for g in dup:
        g.interval[0].botinterval.append(g)
        g.interval[1].topinterval.append(g)

    # set min_len_int value for every node s
    # which is the distance to the nearest top of an interval that pass thru s
    for s in stnodespostorder:
        if s.topinterval:
            s.min_len_int = 0   # node s is an ending of some interval
            s.active = True     # this is a place for an episode, all dup nodes with passing intervals are set to here
            if s.allintervals:
                for i in s.allintervals:
                    i.interval[1].topinterval.remove(i)  # passing interval is removed from interval top list
            if s.botinterval:
                for i in s.botinterval:
                    i.interval[1].topinterval.remove(i)  # starting interval is removed from interval top list
            s.episize = len(s.allintervals) + len(s.botinterval)
            gtids = set(g.gtid for g in s.allintervals).union(g.gtid for g in s.botinterval )
            s.epigtcount = len(gtids) # the number of gene trees participating in the episode
        else:
            if s.botinterval:       # add to parent all intervals that start here
                for i in s.botinterval:
                    s.parent.allintervals.append(i)
            if s.allintervals:      # add to parent all intervals that pass thru here
                for i in s.allintervals:
                    s.parent.allintervals.append(i)

    # Counts selected episode node
    ec_score = 0
    ep_nodes = []
    for s in stnodespostorder:
        if s.active:
            ec_score += 1
            ep_nodes.append(s)

    return ec_score, ep_nodes


def ecfeasbible(gt, st, episodes) -> bool:
    """
    Given a gene tree and a species tree 
    check if duplications from gt can be placed at episodes
    """

    gt.set_lca_mapping(st)
    generate_intervals(gt, st)

    for n in gt.nodes:        
        if n.interval:            
            s,e = n.interval
            e = e.parent
            while s!=e:
                if s in episodes:
                    break # dupl n is in episodes
                s = s.parent
            else:
                break # not feasbile (no breaks)
    else:
        return True # feasible (no not feasible breaks)

    return False

def metaecfeasible(gt, st, episodes, leafdistr=False) -> (bool, int, dict):
    """
    Given a partial gene tree and a species tree 
    check if there is a gene tree gt' that extends gt 
    s.t. duplications from gt' can be placed at episodes 

    Returns triple (feasible, cnt, leafdistributions):
        feasible: bool - true if the scenario is feasible
        cnt: int - if leafdistr is True, the number of feasible leaf-mapping reconstructions 
        leafdistributions: dict - if leafdistr is True, per each unknown leaf the distribution of species leaf reconstruction in feasible mappings 
    """

    gts = str(gt)
    cnt = gts.count("?")

    guls = [ g for g in gt.leaves() if g.clusterleaf == '?' ]

    if leafdistr:
        feascnt = 0        
        leafcounts = emptyleafcounts(guls, st.leaves())

    if not guls:
        res =  ecfeasbible(gt,st,episodes)                
        return res, int(res), {}


    for p in itertools.product(st.leaves(), repeat=cnt):
        cgt = gts
        
        for lb in p:
            cgt = cgt.replace("?", lb.clusterleaf, 1)                    

        res = ecfeasbible(t:=Tree(str2tree(cgt)),st,episodes)        

        if leafdistr:
            if res:
                # aggregate 
                leafcountsupdate(leafcounts, p)
                feascnt += 1

        elif res:
            return res, 1, {} # return the first

    if leafdistr and feascnt:
        return True, feascnt, dict(zip(guls,leafcounts))

    return False, 0, {}


def metaecfeasibleleafmap(gt, st, episodes) -> bool:
    """
    Given a partial gene tree and a species tree 
    check if there is a gene tree gt' that extends gt 
    s.t. duplications from gt' can be placed at episodes 
    """
    gts = str(gt)
    cnt = gts.count("?")


    if not cnt:
        return ecfeasbible(gt,st,episodes)            

    for p in itertools.product(st.leaves(),repeat=cnt):
        cgt = gts
        for lb in p:
            cgt = cgt.replace("?",lb.clusterleaf,1)        
        if ecfeasbible(Tree(str2tree(cgt)),st,episodes):            
            return True

    return False
        






        
