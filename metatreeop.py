import itertools
from typing import Set, Tuple, List
from rec import rec
from fixedec import fixedec
from optTrees_gdscore import opttrees_gdscore
from treeop import Tree, Node, str2tree
import math
from os import getppid
from random import sample, choice
from copy import copy, deepcopy
import statistics 

Unknown = None

def conjuction(a, b):
    if a is False or b is False: return False
    if a is True and b is True: return True
    return Unknown


def disjunction(a, b):
    if a is True or b is True: return True
    if a is False and b is False: return False
    return Unknown


def Lop(a):
    if a is True: return True
    return False


def Mop(a):
    if a is False: return False
    return True


def logic3str(v):
    if type(v) == str: return v
    if v is Unknown: return "U"
    if v is False: return "F"
    return "T"


def ppn(v):
    return "".join(sorted(f"{l.clusterleaf}{l.num}" for l in v.leaves()))+f"#{v.num}"


def pptabs(nm, v, us):
    for s, g in v:
        print(nm, ppn(s), ppn(g), logic3str(v[s, g]))


def ppdistr(d):
    if d is None:
        return "None"
    s = ''
    for gul,distr in d.items():
        s+=f"{gul}{gul.num}:"+" ".join( sl+str(va)  for sl,va in sorted(distr.items()) if va>0)+" "
    return s

def pptabsc(nm, v,  dT, dU):

    def ppc(d, UT):        
        if d is None:
            return f"{UT}=None"
        return f"{UT}=[{ppdistr(d)}]"                

    for s, g in v:
        print(nm, ppn(s), ppn(g), logic3str(v[s, g]), ppc(dT.get((s,g),{}),"t"), " ;;; ",ppc(dU.get((s,g),{}),"u"), )

def pptabs3(st, gt, de, dd, si, deu, ddu, siu):
    def eps(s, g):
        if (s, g) in de and (s, g) in si:
            return logic3str(disjunction(de[s, g], si[s, g]))
        return "-"

    def ppu(u):
        return "{"+" ".join(str(v.num) for v in u)+"}"

    print(f"   {'SN':20}    {'GN':20} De Dd Si Ep | deu ddu siu")
    for g, s in itertools.product(gt.nodes, st.nodes):
        print(
            f"{s.num:2} {ppn(s):15} {g.num:2} {ppn(g):15} {logic3str(de.get((s, g), '-'))}  {logic3str(dd.get((s, g), '-'))}  {logic3str(si.get((s, g), '-'))}  {eps(s, g)}  | {ppu(deu[s,g])} {ppu(ddu[s,g])} {ppu(siu[s,g])}")
        # print(nm,ppn(g),ppn(s),logic3str(v[s,g]))



def is_reconciled_using_wgd(st: Tree, gt: Tree, wgd_nodes: Set[Node], wgddebug=False, excludedoutgroup="") -> Tuple[bool, Set[Node], str]:
    """
    Checks whether S and G can be reconciled using WGD events from a given set nodes of S

    Args:
        st: Tree - a species tree
        gt: Tree - a gene tree with "?" nodes
        wgd_nodes: Set[Node] - allowed duplication epidodes
        wgddebug: Bool - if True print additional debug info
        excludedoutgroup: str - a label excluded from leaf mapping reconstructions

    Returns:
        a tuple contating
        - Bool: if there is a feasible scenario using wgd_nodes
        - Set[Node]: the set of used episodes from wgd_nodes
        - str: a reconstructed gene tree (no ?) if exists
    """

    deltav = {}
    delta_usage = {}
    delta_leafmap = {}
    deltadownv = {}
    deltadown_usage = {}
    deltadown_leafmap = {}
    sigmav = {}
    sigma_usage = {}
    sigma_leafmap = {}

    def delta(s: Node, g: Node) -> Tuple[bool , Set[Node], Tuple[Tuple[Node, str]]]:
        """ g maps to s and g is a duplication
            g not a leaf
            Return F U T
        """
        if (s, g) in deltav:
            return deltav[s, g], delta_usage[s, g], delta_leafmap[s, g]

        is_reconciled, usage, leafmap = False, set(), dict()



        if not g.leaf():
            # if g.num in (2,4): print(ppn(g),ppn(s),"D - sigma",is_reconciled)

            for left, right in (g.c, (g.c[1], g.c[0])):

                is_reconciled_left, usage_left, leafmap_left = deltaexact(s, left)                


                # optimize if is_reconciled_left is False
                if is_reconciled_left is False: continue

                is_reconciled_right, usage_right, leafmap_right = delta_down(s, right)                

                is_reconciled = conjuction(is_reconciled_left, is_reconciled_right)

                if s in wgd_nodes:
                    is_reconciled = Mop(is_reconciled)  # Mop is True or False
                    #print("       X->rr", logic3str(is_reconciled))
                
                    if is_reconciled is True:
                        usage = usage_left | usage_right | {s}
                        leafmap = leafmap_left + leafmap_right
                        break
                else:
                    # s not in wgd_nodes
                    is_reconciled = conjuction(is_reconciled, Unknown)  # with Unknown
                    #print("        ->rr", logic3str(is_reconciled))

                    if is_reconciled is Unknown:
                        usage = usage_left | usage_right
                        leafmap = leafmap_left + leafmap_right
                        break

        deltav[s, g] = is_reconciled
        delta_usage[s, g] = usage
        delta_leafmap[s, g] = leafmap
        return is_reconciled, usage, leafmap

    def deltaexact(s, g):
        """ g maps to s
            Return F U T
        """

        is_reconciled, usage, leafmap = sigma(s, g)  # T or F
        if is_reconciled is True:
            return is_reconciled, usage, leafmap

        return delta(s, g)

    def delta_down(s: Node, g: Node) -> Tuple[bool, Set[Node], Tuple[Tuple[Node, str]]]:
        """ g maps to s or below
            Return F U T
        """

        if (s, g) in deltadownv:
            return deltadownv[s, g], deltadown_usage[s, g], deltadown_leafmap[s, g]

        is_reconciled, usage, leafmap = deltaexact(s, g)        

        usageadd = {}

        if is_reconciled is not True:

            for c in s.c:

                is_reconciled2, usage2, leafmap2 = delta_down(c, g)
                
                if s in wgd_nodes:

                    if is_reconciled2 is Unknown:
                        usage2 = usage2 | {s} # duplication is digested here

                    is_reconciled2 = Mop(is_reconciled2)

                if is_reconciled2 is True or is_reconciled2 is Unknown and is_reconciled is False:
                    is_reconciled, usage, leafmap = is_reconciled2, usage2, leafmap2
                    if is_reconciled is True:
                        break

        deltadownv[s, g] = is_reconciled
        deltadown_usage[s, g] = usage
        deltadown_leafmap[s, g] = leafmap
        return is_reconciled, usage, leafmap

    def sigma(s: Node, g: Node) -> Tuple[bool, Set[Node], Tuple[Tuple[Node, str]]]:
        """
        g maps to s and (g speciation or g a leaf)
        Return F T
        """
        if (s, g) in sigmav:
            return sigmav[s, g], sigma_usage[s, g], sigma_leafmap[s, g]
        elif s.leaf() and g.leaf():
            is_reconciled, usage, leafmap = g.label == s.label or g.label[0] == "?", set(), tuple()
            if g.label[0] == "?":
                if s.label == excludedoutgroup:
                    is_reconciled = False
                else:
                    leafmap = ((g, s.label),)                    
                    

        elif not s.leaf() and not g.leaf():
            is_reconciled, usage, leafmap = False, set(), tuple()
            for left, right in ((s.c[0], s.c[1]), (s.c[1], s.c[0])):
                is_reconciled_left, usage_left, leafmap_left = delta_down(left, g.c[0])

                if is_reconciled_left is False: continue

                is_reconciled_right, usage_right, leafmap_right = delta_down(right, g.c[1])

                is_reconciled = Lop(conjuction(is_reconciled_right, is_reconciled_left))
                if is_reconciled:
                    usage = usage_left | usage_right
                    leafmap = leafmap_left + leafmap_right                    
                    break
        else:
            is_reconciled, usage, leafmap = False, set(), tuple()
        sigmav[s, g] = is_reconciled
        sigma_usage[s, g] = usage
        sigma_leafmap[s, g] = leafmap
        return is_reconciled, usage, leafmap


    
    for (s,g) in itertools.product(st.nodes,gt.nodes):
        delta_down(s,g)    
        delta(s,g)    

    is_valid, node_usage, leafmap = delta_down(st.root, gt.root)

    #wgddebug=True
    if wgddebug:
        pptabs3(st, gt, deltav, deltadownv, sigmav, delta_usage, deltadown_usage, sigma_usage)
        # pptabs("d ",deltav,delta_usage)
        # pptabs("dd",deltadownv,deltadown_usage)
        # pptabs("si",sigmav,sigma_usage)

    if is_valid is True:
        return is_valid, node_usage, gt.nodemaprepr(dict(leafmap))
    return False, set(), ""

def totalcount(d):
    if not d: return 0
    for k in d:
        return sum(d[k].values())

def combinecounts(d1,d2):
    # disjont keys (gene nodes)
    if d1 is None or d2 is None: return None

    if not d1: return d2.copy()
    if not d2: return d1.copy()

    s1 = totalcount(d1)
    s2 = totalcount(d2)    

    d = {}
    for k,dv in d1.items():
        d[k] = { sn:(s2*v) for sn,v in dv.items() }
    for k,dv in d2.items():
        d[k] = { sn:(s1*v) for sn,v in dv.items() }        
    #print(d1,d2,s1,s2,d)
    return d

def adddict(dest,d):
    for k,v in d.items():
        if k in dest:
            dest[k]+=v
        else:
            dest[k]=v

def sumcounts(*dicts):
    dest = {}    
    present = False
    for d in dicts:
        if d is not None:
            present = True
            for k,vd in d.items():
                if k in dest:
                    adddict(dest[k],vd) 
                else:
                    dest[k]=vd.copy()
    
    if present: return dest
    return None


def is_reconciled_using_wgd_withcounts(st: Tree, gt: Tree, wgd_nodes: Set[Node], wgddebug=False, excludedoutgroup="") -> [bool, dict]:
    """
    Checks whether S and G can be reconciled using WGD events from a given set nodes of S with g->s map counts

    Args:   
        st: Tree - a species tree
        gt: Tree - a gene tree with "?" nodes
        wgd_nodes: Set[Node] - allowed duplication epidodes
        wgddebug: Bool - if True print additional debug info
        excludedoutgroup: str - a label excluded from leaf mapping reconstructions

    Returns:
        a tuple contating
        - Bool: if there is a feasible scenario using wgd_nodes
        - Dict: distr. a map witk key=unknownleaf and a value dictionary with key=sleaf value=#occurences
    """

    deltav = {}    
    deltadownv = {}        
    sigmav = {}
    
    delta_countsT = {}
    delta_countsU = {}

    deltadown_countsT = {}
    deltadown_countsU = {}
    deltaexact_countsU = {}
    deltaexact_countsT = {}
    sigma_counts = {} # always T

    gul = [ gl for gl in gt.leaves() if gl.clusterleaf == '?' ] 

    def delta_wc(s: Node, g: Node) -> Tuple[bool , Set[Node], Tuple[Tuple[Node, str]]]:
        """ g maps to s and g is a duplication
            g not a leaf
            Return F U T
        """
        if (s, g) in deltav:
            return deltav[s, g]

        is_reconciled = False

        delta_countsT[s, g] = None
        delta_countsU[s, g] = None

    
        if not g.leaf():
            l,r = g.c

            # avoid double exact
            is_reconciled_exl = deltaexact_wc(s, l)                                
            is_reconciled_exr = deltaexact_wc(s, r)            

            is_reconciled = conjuction(is_reconciled_exl, is_reconciled_exr)

            c0 = combinecounts(deltaexact_countsT[s, l], deltaexact_countsT[s, r]) 
            c1 = combinecounts(deltaexact_countsT[s, l], deltaexact_countsU[s, r])   
            c2 = combinecounts(deltaexact_countsU[s, l], deltaexact_countsT[s, r]) 
            c3 = combinecounts(deltaexact_countsU[s, l], deltaexact_countsU[s, r])            
        
            delta_countsU[s, g] = sumcounts(c0, c1, c2, c3)

            if s in wgd_nodes:                
                delta_countsT[s, g] = c0                                        
                                                        
            else:
                is_reconciled = conjuction(is_reconciled, Unknown)

        
            for is_reconciled_e, left, right in ((is_reconciled_exl, l, r), (is_reconciled_exr,r,l)):
                #if is_reconciled_e is False: continue
                
                for sc in s.c:                     
                    is_reconciled_right = delta_down_wc(sc, right)

                    _is_reconciled = conjuction(is_reconciled_e, is_reconciled_right)

                    c0 = combinecounts(deltaexact_countsT[s, left], deltadown_countsT[sc, right])# 
                    c1 = combinecounts(deltaexact_countsU[s, left], deltadown_countsT[sc, right]) 
                    c2 = combinecounts(deltaexact_countsT[s, left], deltadown_countsU[sc, right])# 
                    c3 = combinecounts(deltaexact_countsU[s, left], deltadown_countsU[sc, right])

                    delta_countsU[s, g] = sumcounts(delta_countsU[s, g], c0, c1, c2, c3)                                       

                    if s in wgd_nodes:

                        _is_reconciled = Mop(_is_reconciled)  # Mop is True or False                        

                        if _is_reconciled is True:                             
                            is_reconciled = True            
                                                                            
                        delta_countsT[s, g] = delta_countsU[s, g]
                    
                    else:
                        # s not in wgd_nodes
                        _is_reconciled = conjuction(_is_reconciled, Unknown)  # with Unknown -> U F
                        
                        if _is_reconciled is Unknown and is_reconciled is not True:
                            is_reconciled = _is_reconciled                            

        deltav[s, g] = is_reconciled        
        return is_reconciled

    def deltaexact_wc(s, g) -> bool:
        """ g maps to s (epsilon from the article)
            Return F U T
        """

        is_reconciled = sigma_wc(s, g)  # T or F
        

        deltaexact_countsU[s, g] = None        
        deltaexact_countsT[s, g] = sigma_counts[s, g]
                                        
        is_reconciled2 = delta_wc(s, g)
        
        if is_reconciled2 is Unknown: # cond. needed
            deltaexact_countsU[s, g] = delta_countsU[s, g]  

        # if is_reconciled2 is True: 
        #     if delta_countsU[s, g] is not None:
        #         raise Exception(delta_countsU[s, g]) 
        # some examples exists
        
        deltaexact_countsT[s, g] = sumcounts(deltaexact_countsT[s, g],  delta_countsT[s, g])

        if is_reconciled is True: # T F
            return True

        return is_reconciled2 # T U F


    def delta_down_wc(s: Node, g: Node) -> bool:
        """ g maps to s or below
            Return F U T
        """

        if (s, g) in deltadownv:
            return deltadownv[s, g]

        # (4)
        is_reconciled = deltaexact_wc(s, g)    
    
        deltadown_countsT[s, g] = deltaexact_countsT[s, g]
        deltadown_countsU[s, g] = deltaexact_countsU[s, g]


        for c in s.c:
            is_reconciled2 = delta_down_wc(c, g)

            
            if s in wgd_nodes:
                # (5)
                #if is_reconciled2 is True:
                deltadown_countsT[s, g] = sumcounts(deltadown_countsT[s, g], deltadown_countsT[c, g], deltadown_countsU[c, g]) # T->T

                is_reconciled2 = Mop(is_reconciled2) # T or F

                if is_reconciled2 is True: 
                    is_reconciled = True

            else:
                # (6)
                # if is_reconciled2 is Unknown and is_reconciled is False:
                # is_reconciled = is_reconciled2
                deltadown_countsT[s, g] = sumcounts(deltadown_countsT[s, g], deltadown_countsT[c, g]) # T->T
                deltadown_countsU[s, g] = sumcounts(deltadown_countsU[s, g], deltadown_countsU[c, g]) # T->T

                if is_reconciled2 is True:
                    is_reconciled = True

                    #deltadown_countsU[s, g] = sumcounts(deltadown_countsU[s, g], deltadown_countsU[c, g]) # U->U
                if is_reconciled2 is Unknown:
                    if is_reconciled is not True:
                        is_reconciled = is_reconciled2                    

        deltadownv[s, g] = is_reconciled
        
        return is_reconciled

    def sigma_wc(s: Node, g: Node) -> Tuple[bool, Set[Node], Tuple[Tuple[Node, str]]]:
        """
        g maps to s and (g speciation or g a leaf)
        Return F T
        """

        if (s, g) in sigmav:
            return sigmav[s, g]

        

        if s.leaf() and g.leaf():
            is_reconciled = g.label == s.label or g.label[0] == "?"
            if is_reconciled == False:
                sigma_counts[s, g] = None
            else:
                sigma_counts[s, g] = {}
            if g.label[0] == "?":
                if s.label == excludedoutgroup:
                    is_reconciled = False
                    sigma_counts[s, g] = None
                else:
                    sigma_counts[s, g] = { g: {s.clusterleaf:1}} 
                                    

        elif not s.leaf() and not g.leaf():
            is_reconciled = False
            sigma_counts[s, g] = None 
            for left, right in ((s.c[0], s.c[1]), (s.c[1], s.c[0])):
                is_reconciled_left = delta_down_wc(left, g.c[0])

                if is_reconciled_left is False: continue # optional

                is_reconciled_right = delta_down_wc(right, g.c[1])

                _is_reconciled = Lop(conjuction(is_reconciled_right, is_reconciled_left)) # True, False

                sigma_counts[s, g] = sumcounts(sigma_counts[s, g], combinecounts(deltadown_countsT[left, g.c[0]], deltadown_countsT[right, g.c[1]] )) 

                if _is_reconciled:                                
                    is_reconciled = _is_reconciled
        else:
            sigma_counts[s, g] = None
            is_reconciled = False

        sigmav[s, g] = is_reconciled
        
        return is_reconciled

    is_valid = delta_down_wc(st.root, gt.root)

    for (s,g) in itertools.product(st.nodes,gt.nodes):
        delta_down_wc(s,g)    
        delta_wc(s,g)    

    if wgddebug:
        pptabs3(st, gt, deltav, deltadownv, sigmav)

        pptabsc("d ",deltav, delta_countsT, delta_countsU)
        pptabsc("dd",deltadownv, deltadown_countsT, deltadown_countsU)
        pptabsc("si",sigmav, sigma_counts, {})
        deltaexactv = { (s,g):deltaexact_wc(s,g) for (s,g) in itertools.product(st.nodes,gt.nodes)}
        pptabsc("dx",deltaexactv, deltaexact_countsT, deltaexact_countsU)

        # ddv = { (s,g):delta_down_wc(s,g) for (s,g) in itertools.product(st.nodes,gt.nodes)}
        # pptabsc("dd2",ddv, deltadown_countsT, deltadown_countsU)

    if is_valid is True:        
        return is_valid, deltadown_countsT[st.root, gt.root]
    return False,None


def split_outgrouped_tree(tree: Tree, outgroup: str = "outgroup") -> List[Tree]:    
    res = []
    for i in tree.nodes:        
        if i.leaf() and i.clusterleaf == outgroup:
            res.append(Tree(str2tree(str(i.parent.c[0]))))
    return res


def combine_gene_trees(gtrees: List[Tree], outgroup: str = "outgroup") -> Tree:
    """ Adds an outgroup species to all input gene trees and merges them into one tree,
    eg. (a, b) and (b, c) -> (((a, b), outgroup), ((b, c), outgroup))
    """
    gtrees = [f"({str(gtree)},{outgroup})" for gtree in gtrees]

    # balanced output to avoid recursion 
    while len(gtrees)>1:
        dest = []
        while len(gtrees)>1:
            dest.append("("+gtrees.pop()+","+gtrees.pop()+")")
        if len(gtrees)==1:
            dest.append(gtrees.pop())
        gtrees = dest

    #combined_gtree = gtrees[0]
    #for i in range(1, len(gtrees)):
    #    combined_gtree = f"({combined_gtree},{gtrees[i]})"
    return Tree(str2tree(gtrees[0]))


def add_outgroup(stree: Tree, outgroup: str) -> Tree:
    return Tree(str2tree(f"({str(stree)},{outgroup})"))


def random_labelling(gtree: Tree, stree: Tree, outgroup: str) -> str:
    streeleaves = [ l for l in stree.leaves() if l.clusterleaf!=outgroup ]
    def _m(g: Node):
        if g.leaf():
            if g.clusterleaf[0] == '?':
                return choice(streeleaves).clusterleaf        
            return g.clusterleaf
        return "("+ ",".join(_m(c) for c in g.c) + ")"

    return _m(gtree.root)

def count_wgd_nodes_combined(
        stree: Tree, 
        gtrees: List[Tree], 
        outgroup: str = "outgroup", 
        wgddebug = False, 
        outfile = None,
        noimprovement_stop=0,
        randomize_from=0,
        setid=None,
        reversed_climb = 0,
        initial_gene_tree = None,
        distribution_maps = False,
        reference_trees = None,
        distribution_maps_epi = False,
        print_distr_maps = False,
        verbose = 1,
        distr_counts = False
        ) -> Tuple[float, Set[Node]]:
    """ 
    Returns a minimal number of nodes in a species tree S that need to contain WGD events
     in order to reconcile S and a set of gene trees with ?
    
    """
    

    for g in gtrees:
        if not g.is_binary():
            raise ValueError(f"Found a non-binary gene tree {g}")

    gtree = combine_gene_trees(gtrees, outgroup)
    stree = add_outgroup(stree, outgroup)

    if reference_trees:
        reference_tree = combine_gene_trees(reference_trees)
    else:
        reference_tree = None

    return count_wgd_nodes(stree, gtree,  outgroup, 
        wgddebug=wgddebug, 
        outfile=outfile,
        noimprovement_stop=noimprovement_stop, 
        randomize_from=randomize_from, 
        setid=setid, 
        reversed_climb=reversed_climb,
        initial_gene_tree = initial_gene_tree,
        distribution_maps = distribution_maps,
        reference_tree=reference_tree,
        distribution_maps_epi = distribution_maps_epi,
        print_distr_maps = print_distr_maps,
        verbose = verbose,
        distr_counts=distr_counts)

def gtwithdistrmaps(st, gt, dpdistr, reference_tree=None, outgroup: str = "outgroup", distr_counts=False ) -> str:

    def traverse_outgroup(n):
        """
        Set _outgroup=True iff the subtree contains outgroup
        """
        if n.leaf():            
            n._outgroup = n.clusterleaf == outgroup
        else: 
            n._outgroup = False
            for c in n.c:
                if traverse_outgroup(c):
                    n._outgroup = True

        return n._outgroup

    def traverse(n):
        
        if n.leaf():
            #if n._outgroup: 
            #   return '==='

            if n in dpdistr:
                s="d="+"{"+",".join(f"'{k}':{v:.4f}" for k, v in dpdistr[n].items())+"}"
                s+=f" s={len(dpdistr[n])}"
                sv = ",".join(f"{v:.4f}" for v in sorted(set(dpdistr[n].values())))
                s+=f" valset=[{sv}]"                                
                s+=f" lenvalset={len(set(dpdistr[n].values()))}"
                snode = st.root.findnode( frozenset(dpdistr[n].keys()))
                mean = statistics.mean(dpdistr[n].values())

                if len(dpdistr[n])>1:
                    pstdev = statistics.pstdev(dpdistr[n].values()) # population stdev
                else: 
                    pstdev = 0

                cv = pstdev/mean # coefficient of variation

                if snode is not None:
                    s += f" stdistr={snode.num}"                    
                    snode.lfmapcnt+=1 # increase the count in s
                    snode.lfmapmean.append(mean)
                    snode.lfmappstdev.append(pstdev)
                    snode.lfmapcv.append(cv)
                else:
                    s+=f" stdistr=-1"
                    st.root.trueroot.lfmapcnt_none += 1 # increase the count in s                    
                    st.root.trueroot.lfmapmean_none.append(mean)
                    st.root.trueroot.lfmappstdev_none.append(pstdev)
                    st.root.trueroot.lfmapcv_none.append(cv)
                if n.refleaf is not None:
                    s+=" r="+n.refleaf.clusterleaf
                    refval = dpdistr[n].get(n.refleaf.clusterleaf,None)
                    if refval is not None:
                        s+=f" rv={refval:.4f}"
                    else:
                        s+=" rv=None"

                return s
            return n.clusterleaf

        #if n._outgroup:
            #return ";".join( sc for c in n.c if (sc:=traverse(c)) )

        return "("+",".join( traverse(c) for c in n.c) +")"

    if reference_tree:
        for l1, l2 in zip(gt.leaves(), reference_tree.leaves()):
            l1.refleaf = l2
    else:
        for l in gt.leaves():
            l.refleaf = None


    distrsum = None
    if not hasattr(gt.root, 'distrsum'): # prevent double normalization        
        for k,d in dpdistr.items():
                s = sum(d.values())                
                if not distr_counts:
                    for slab in d:
                        d[slab] = d[slab]/s

                if distrsum is None:
                    distrsum = s
                else:
                    assert s == distrsum            

    gt.root.distrsum = distrsum # store in the root only

    for sn in st.nodes:
        sn.lfmapcnt = 0 
        sn.lfmapmean = []
        sn.lfmappstdev = []
        sn.lfmapcv = []
    st.root.trueroot.lfmapcnt_none = 0 
    st.root.trueroot.lfmapmean_none = [] 
    st.root.trueroot.lfmappstdev_none = [] 
    st.root.trueroot.lfmapcv_none = []
    
    traverse_outgroup(gt.root)

    res = traverse(gt.root)    

    return res, distrsum 

def randcombinations(X,k):    
    while True:
        yield sample(X, k)

def count_wgd_nodes(
        st: Tree, 
        gt: Tree, 
        outgroup: str = 'outgroup', 
        wgddebug=False, 
        outfile=None,
        noimprovement_stop=0,
        randomize_from=0,
        setid=None,
        reversed_climb=0,
        initial_gene_tree = None,
        distribution_maps = False,
        reference_tree = None,
        distribution_maps_epi = False,   # tree for distribution validation; must have the same topology as gt
        print_distr_maps = False,        
        verbose = 1,
        distr_counts = False
        ) -> Tuple[float, Set[Node]]:
    """ 
    Returns a minimal number of nodes in a species tree S that need to contain WGD events
    in order to reconcile S and a gene tree with ?

    Outgroup is requred in initial_gene_tree, reference_tree, gt and st
        outfile - appends report results 

    """

    #upper_bound = opttrees_gdscore(st, gt) if count_upper_bound else len(st.nodes) - 1
    
    outgrouped = st.root.c[1].clusterleaf == outgroup

    epiattr = ['episize', 'epigtcount', 'epibestwgd', 'num', 'lfmapcnt', 'lfmapcnt_none', 'lfmapcv', 'lfmapcv_none','distrsum' ] 

    # root.c[0] - skip outgroup
    

    if initial_gene_tree:
        # initialize using initial gene tree
        gt_inferred_str = initial_gene_tree
    else:        
        # initialize upper bound using random gene tree
        gt_inferred_str = random_labelling(gt, st, outgroup)        

    best_cost, best_wgd_nodes = rec([Tree(str2tree(gt_inferred_str))], st)

    if outfile:
        with open(outfile+".genetree","w") as f:
            f.write(gt_inferred_str)

    outstats = f"initialgenetree={initial_gene_tree}\n"
    outstats+= f"initialgenetreecost={best_cost}\n"

    fixed_wgd_nodes = fixedec(gt, st)

    maxec = len(st.root.nodes()) 
    potential_wgd_nodes = list(set(st.root.nodes()) - fixed_wgd_nodes)
    samplingsets = False
    
    if wgddebug:                
        print(gt)
        print(st.root.markrepr(fixed_wgd_nodes))

    unklabs = len(gt.unknownlabels())

    if "_" in setid:
        # in case of more complex setids split
        outstats+="".join(f"setid{i+1}={v}\n" for i,v in enumerate(setid.split("_")))

    outstats+=f"setid=\"{setid}\"\n"
    outstats+=f"genetree=\"{gt}\"\n"
    outstats+=f"speciestree=\"{st}\"\n"
    outstats+=f"speciestreefixedwgd=\"{st.root.markrepr(fixed_wgd_nodes)}\"\n"
    outstats+=f"fixedwgd={len(fixed_wgd_nodes)}\n"
    outstats+=f"reversed_climb={reversed_climb}\n"
    outstats+=f"unknownlabels={unklabs}\n"

    if outgrouped:        
        gts_split = split_outgrouped_tree(gt, outgroup)
        stroot = strootnooutgroup = st.root.c[0]
        eccorrection = 1 if len(gts_split)>1 else 0 # additional dupliaction of two trees are present

        outstats+=f"#with no outgroup (wo)\n"
        outstats+=f"genetrees_wo=\"{';'.join(str(g) for g in gts_split)}\"\n"
        outstats+=f"speciestree_wo=\"{stroot}\"\n"
        outstats+=f"speciestreefixedwgd_wo=\"{stroot.markrepr(fixed_wgd_nodes)}\"\n"
        outstats+=f"fixedwgd_wo={len(fixed_wgd_nodes)-eccorrection}\n"
        outstats+=f"#end (wo)\n"

    else:
        stroot = st.root
        gts_split = [gt]
        eccorrection = 0

    st.root.trueroot=stroot # needed in distr counting
        

    climbs=""
    dpcalls=0
    exactsolution = True
    stop = False

    cur_cost_search = len(fixed_wgd_nodes) # only in reversed

    sampling_occured = False
    dpfromlastimprovement = 0

    outstats+="\n\n"

    if not unklabs:
        exactsolution = True



    def reporterror(outfile, info, outstats):            
        f = open("metaec.err.log","a")
        f.write(f"\n========={outfile}============\n")
        f.write(outstats)        
        f.close()
        raise Exception(f"[{outfile}] {info}. See metaec.err.log for details")



    def getdistrmaps(st, gt, best_wgd_nodes, reference_tree, outgroup):

        if verbose:
            print(f"[{setid}] Computing gene-species distribution maps")

        _, dpdistr = is_reconciled_using_wgd_withcounts(st, gt, best_wgd_nodes, excludedoutgroup=outgroup)

        if dpdistr is None:
            reporterror(outfile, "No distribution map (dpdistr is None)", outstats)

        tr, distrsum = gtwithdistrmaps(st, gt, dpdistr, outgroup=outgroup, distr_counts=distr_counts)
        stats="gtdistrmaps="+tr+"\n"

        stats+=f"distrsum={distrsum}\n"

        if reference_tree:
            stats+="leafmapinferrence="+gtwithdistrmaps(gt, dpdistr, reference_tree, outgroup=outgroup, distr_counts=distr_counts)[0]+"\n"

        return stats

    if distribution_maps_epi:
        # calculate only distributions based on the episode set
        # outstats += getdistrmaps(st, gt, best_wgd_nodes, reference_tree, outgroup)
        dpcalls = 1

        exactsolution = best_cost == len(fixed_wgd_nodes)        
       
    else:

        # perform the search
        
        while unklabs:

            all_explored_stop = False

            if reversed_climb:
                if cur_cost_search == best_cost:                     
                    break
                k = cur_cost_search - len(fixed_wgd_nodes)
            else:
                if best_cost == len(fixed_wgd_nodes):   
                    if verbose:
                        print(f"[{setid}] Best cost = fixed wgds. Stop: {best_cost}/{maxec}")                                     
                    break                    
                k = best_cost - len(fixed_wgd_nodes) - 1

            comb = math.comb(len(potential_wgd_nodes),k)

            samplingsets = not (not randomize_from or randomize_from>comb)


      
            if verbose:
                print(f"[{setid}] EC:{best_cost}/{maxec} Test:{k+len(fixed_wgd_nodes)} FxdWgd:{len(fixed_wgd_nodes)} PotentialEpi:{len(potential_wgd_nodes)} K:{k} Comb:{comb} RndSmpl:{samplingsets} StopAfter:{noimprovement_stop} UnknwnLbls:{unklabs}")

            if samplingsets:
                wgd_node_sets = randcombinations(potential_wgd_nodes, k)
            else:
                wgd_node_sets = itertools.combinations(potential_wgd_nodes, k)        
            

            cnt = 0
            for wgd_nodes in wgd_node_sets:                        


                # if (len(wgd_nodes)==2):
                #     wgd_nodes = [ v for v in st.nodes if v.num in (5,1)]

                wgd_node_set = set(wgd_nodes) | fixed_wgd_nodes
                # print ("!!WN", [v.num for v in wgd_nodes])
                # print ("!!WS", [v.num for v in wgd_node_set])
                # print ("!!FI", [v.num for v in fixed_wgd_nodes])
                

                is_feasible, used_wgd_nodes, gt_inferred_str_cur = is_reconciled_using_wgd(st, gt, wgd_node_set, excludedoutgroup=outgroup)                                
                cnt+=1
                dpcalls+=1
                dpfromlastimprovement+=1                
        
                if is_feasible:
                    gt_inferred_str = gt_inferred_str_cur                    
                    gt_inferred = Tree(str2tree(gt_inferred_str))                    
                    
                    ec, used_ec_nodes = rec({gt_inferred}, st)

                    if verbose:
                        print(f"[{setid}] Feasible solution ec={ec} lenec={len(used_ec_nodes)} wgd={len(used_wgd_nodes)} bc={best_cost}")

                    if ec>len(used_wgd_nodes):
                        raise Exception("Incorrect EC>DP cost")

                    # if len(used_wgd_nodes)<ec:
                    #     print (f"[{setid}]","fx",sorted(v.num for v in fixed_wgd_nodes))
                    #     print (f"[{setid}]","wgd",sorted(v.num for v in wgd_nodes))
                    #     print (f"[{setid}]","uec",sorted(v.num for v in used_ec_nodes))
                    #     print (f"[{setid}]","uwg",sorted(v.num for v in used_wgd_nodes))
                    #     print(f"[{setid}]",gt)

                    #     print(f"[{setid}]",st.root.attrrepr(epiattr))
                    #     print(f"[{setid}]",gt_inferred)
                        
                    #     print(f"[{setid}] Warning!!!!!!!!!!!!!!!!")
                    
                    #if ec < len(used_wgd_nodes):
                    best_cost, best_wgd_nodes = ec, used_ec_nodes
                    #else:
                    #    best_cost, best_wgd_nodes = len(used_wgd_nodes), used_wgd_nodes                                
                    dpfromlastimprovement = 0

                    if verbose:
                        print(f"[{setid}] Feasible solution with new best_cost:{best_cost}/{maxec}")


                    if reversed_climb:
                        # solution found
                        exactsolution = not sampling_occured 
                        stop = True    

                    if outfile:
                        with open(outfile+".genetree","w") as f:
                            f.write(gt_inferred_str)                    
            
                    break            
                
                if noimprovement_stop and dpfromlastimprovement>=noimprovement_stop:                    
                    exactsolution = False # unknown
                    stop = True
                    if verbose:
                        print(f"[{setid}] Stopping criterion reached (no exact solution). Stop with {best_cost}/{maxec}")                                     
                    break
            else:   
                # all combinations explored
                
                if reversed_climb: # no solution located; search in larger
                    cur_cost_search+=1
                else:
                    exactsolution = True  # no solution located; accept current (exact)
                    stop = True
                    if verbose:
                        print(f"[{setid}] All combinations explored. Stop with the current best cost: {best_cost}/{maxec}")

            sampling_occured |= samplingsets # important in reverse climb
            
            if climbs: climbs+=";"

            climbs+=f"{k};{cnt};{comb}"
            
            if stop:          
                break
        

        outstats+=f"climbs={climbs}\n"
        outstats+=f"samplingsets={samplingsets}\n"
        outstats+=f"outgenetree=\"{gt_inferred_str}\"\n"
        
    
    
    for n in best_wgd_nodes:
        n.epibestwgd=1
    
    outstats+=f"bestcost={best_cost}\n"
    outstats+=f"dpcalls={dpcalls}\n"
    outstats+=f"outspeciestree=\"{st.root.attrrepr(epiattr)}\"\n"
    outstats+=f"exactsolution={exactsolution}\n"
    outstats+=f"unknownlabels={unklabs}\n"

    if outgrouped:

        # recompute ec to get epicounts for wo

        gt_wo = split_outgrouped_tree(Tree(str2tree(gt_inferred_str)), outgroup)

        st_wo = Tree(str2tree(str(stroot)))
        worec, _ = rec(gt_wo, st_wo)

        # check correctness
        for stn, st_won in zip(stroot.nodes(), st_wo.root.nodes()):
            a='episize'
            if hasattr(stn, a) and hasattr(st_won, a):
                if getattr(stn,a)!=getattr(st_won, a):
                    raise Exception(f"Incorrect EC attribute {a} {getattr(stn,a)} {getattr(st_won,a)}")
                continue

            if not hasattr(stn,a) and not hasattr(st_won, a):
                continue
            else: 
                raise Exception(f"Only one EC attribute present {a} {hasattr(stn,a)} {hasattr(st_won,a)} in {stn.attrrepr(epiattr)} {st_won.attrrepr(epiattr)}")

        gt_inferred_str_wo = ";".join(str(t) for t in gt_wo )

        if distribution_maps:        
            outstats+=getdistrmaps(st, gt, best_wgd_nodes, reference_tree, outgroup)
            

        if print_distr_maps:            
            _, dpdistr = is_reconciled_using_wgd_withcounts(st, gt, best_wgd_nodes, excludedoutgroup=outgroup)
            tr, distrsum = gtwithdistrmaps(st, gt, dpdistr, outgroup=outgroup, distr_counts=distr_counts)
            print(tr + f" distrsum={distrsum}")

        outstats+=f"outgenetrees_wo=\"{gt_inferred_str_wo}\"\n"
        outstats+=f"bestcost_worec={worec}\n"
        outstats+=f"bestcost_wo={best_cost-eccorrection}\n"
        outstats+=f"outspeciestree_wo=\"{stroot.attrrepr(epiattr)}\"\n"
        outstats+=f"outspeciestree_worec=\"{st_wo.root.attrrepr(epiattr)}\"\n"

        if worec != best_cost-eccorrection:
            reporterror(outfile, "Inconsistency in EC cost computation {worec} {best_cost-eccorrection}.", outstats)
    else:
        if distribution_maps:        
            outstats+=getdistrmaps(st, gt, best_wgd_nodes, reference_tree, outgroup)

    
    return best_cost, best_wgd_nodes, exactsolution, outstats
