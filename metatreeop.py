import itertools
from typing import Set, Tuple, List
from rec import rec
from fixedec import fixedec
from optTrees_gdscore import opttrees_gdscore
from treeop import Tree, Node, str2tree
import math
from os import getppid
from random import sample, choice

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
    return "".join(sorted(l.clusterleaf for l in v.leaves()))


def pptabs(nm, v, us):
    for s, g in v:
        print(nm, ppn(g), ppn(s), logic3str(v[s, g]))


def pptabs3(st, gt, de, dd, si):
    def eps(s, g):
        if (s, g) in de and (s, g) in si:
            return logic3str(disjunction(de[s, g], si[s, g]))
        return "-"

    print(f"   {'GN':5}    {'SN':5} De Dd Si Ep")
    for g, s in itertools.product(gt.nodes, st.nodes):
        print(
            f"{g.num:2} {ppn(g):5} {s.num:2} {ppn(s):5} {logic3str(de.get((s, g), '-'))}  {logic3str(dd.get((s, g), '-'))}  {logic3str(si.get((s, g), '-'))}  {eps(s, g)}")
        # print(nm,ppn(g),ppn(s),logic3str(v[s,g]))


def is_reconciled_using_wgd(st: Tree, gt: Tree, wgd_nodes: Set[Node], wgddebug=False, excludedoutgroup="") -> Tuple[bool, Set[Node], str]:
    """
    Checks whether S and G can be reconciled using WGD events from a given set nodes of S

    Args:
        st: Tree - a species tree
        gt: Tree - a gene tree
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
                    if is_reconciled is True:
                        usage = usage_left | usage_right | {s}
                        leafmap = leafmap_left + leafmap_right
                        break
                else:
                    # s not in wgd_nodes
                    is_reconciled = conjuction(is_reconciled, Unknown)  # with Unknown

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

        if is_reconciled is not True:
            for c in s.c:
                is_reconciled2, usage2, leafmap2 = delta_down(c, g)
                if s in wgd_nodes:
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

    is_valid, node_usage, leafmap = delta_down(st.root, gt.root)

    if wgddebug:
        pptabs3(st, gt, deltav, deltadownv, sigmav)
        # pptabs("d ",deltav,delta_usage)
        # pptabs("dd",deltadownv,deltadown_usage)
        # pptabs("si",sigmav,sigma_usage)

    if is_valid is True:
        return is_valid, node_usage, gt.nodemaprepr(dict(leafmap))
    return False, set(), ""


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

def count_wgd_nodes_combined(stree: Tree, gtrees: List[Tree], outgroup: str = "outgroup", 
        wgddebug = False, 
        outfile = None,
        noimprovement_stop=0,
        randomize_from=0,
        setid=None,
        reversed_climb = 0,
        initial_gene_tree = None
        ) -> Tuple[float, Set[Node]]:
    """ Returns a minimal number of nodes in a species tree S that need to contain WGD events
     in order to reconcile S and a set of gene trees with ?"""
    # print(len(gtrees))
    for g in gtrees:
        if not g.is_binary():
            raise ValueError(f"Found a non-binary gene tree {g}")
    gtree = combine_gene_trees(gtrees)
    stree = add_outgroup(stree, outgroup)
    return count_wgd_nodes(stree, gtree,  outgroup, wgddebug=wgddebug, outfile=outfile,
        noimprovement_stop=noimprovement_stop, randomize_from=randomize_from, setid=setid, 
        reversed_climb=reversed_climb,
        initial_gene_tree = initial_gene_tree)

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
        initial_gene_tree = None      
        ) -> Tuple[float, Set[Node]]:
    """ 
    Returns a minimal number of nodes in a species tree S that need to contain WGD events
    in order to reconcile S and a gene tree with ?

        outfile - appends report results 
    """

    #upper_bound = opttrees_gdscore(st, gt) if count_upper_bound else len(st.nodes) - 1

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

    outstats = f"initialgenetree={initial_gene_tree}\ninitialgenetreecost={best_cost}\n"

    fixed_wgd_nodes = fixedec(gt, st)
    maxec = len(st.root.nodes())
    potential_wgd_nodes = list(set(st.root.nodes()) - fixed_wgd_nodes)
    samplingsets = False
    
    if wgddebug:                
        print(gt)
        print(st.root.markrepr(fixed_wgd_nodes))

    unklabs = len(gt.unknownlabels())
    outstats+=f"setid=\"{setid}\"\ngenetree=\"{gt}\"\nspeciestree=\"{st}\"\nspeciestreefixedwgd=\"{st.root.markrepr(fixed_wgd_nodes)}\"\nfixedwgd={len(fixed_wgd_nodes)}\n\nreversed_climb={reversed_climb}\n\nunknownlabels={unklabs}\n\n"
        

    climbs=""
    dpcalls=0
    exactsolution = True
    stop = False

    cur_cost_search = len(fixed_wgd_nodes) # only in reversed

    sampling_occured = False
    dpfromlastimprovement = 0

    if not unklabs:
        exactsolution = True



    while unklabs:

        if reversed_climb:
            if cur_cost_search == best_cost: 
                break
            k = cur_cost_search - len(fixed_wgd_nodes)
        else:
            if best_cost == len(fixed_wgd_nodes):
                break
            k = best_cost - len(fixed_wgd_nodes) - 1

        comb = math.comb(len(potential_wgd_nodes),k)

        samplingsets = not (not randomize_from or randomize_from>comb)
  
        print(f"[{setid}] EC:{best_cost}/{maxec} Test:{k+len(fixed_wgd_nodes)} FixedWgd:{len(fixed_wgd_nodes)} PotentialEpi:{len(potential_wgd_nodes)} K:{k} Combinations:{comb} RndSampling:{samplingsets} StopAfter:{noimprovement_stop} UnknownLabels:{unklabs}")

        if samplingsets:
            wgd_node_sets = randcombinations(potential_wgd_nodes, k)
        else:
            wgd_node_sets = itertools.combinations(potential_wgd_nodes, k)        

        cnt = 0
        for wgd_nodes in wgd_node_sets:                        
            wgd_node_set = set(wgd_nodes) | fixed_wgd_nodes
            is_feasible, used_wgd_nodes, gt_inferred_str = is_reconciled_using_wgd(st, gt, wgd_node_set, excludedoutgroup=outgroup)
            cnt+=1
            dpcalls+=1
            dpfromlastimprovement+=1
            
            

            if is_feasible:
                
                gt_inferred = Tree(str2tree(gt_inferred_str))
                ec, used_ec_nodes = rec({gt_inferred}, st)
                if ec < len(used_wgd_nodes):
                    best_cost, best_wgd_nodes = ec, used_ec_nodes
                else:
                    best_cost, best_wgd_nodes = len(used_wgd_nodes), used_wgd_nodes                                
                dpfromlastimprovement = 0

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
                break
        else:   
            # all combinations explored
            if reversed_climb: # no solution located; search in larger
                cur_cost_search+=1
            else:
                exactsolution = True  # no solution located; accept current (exact)
                break  

        sampling_occured |= samplingsets # important in reverse climb
        
        if climbs: climbs+=";"

        climbs+=f"{k};{cnt};{comb}"
        if stop: 
            break


    outstats+=f"bestcost={best_cost}\ndpcalls={dpcalls}\nclimbs={climbs}\noutgenetree=\"{gt_inferred_str}\"\noutspeciestree=\"{st.root.markrepr(best_wgd_nodes)}\"\noutspeciestreeepicount=\"{st.root.attrrepr('episize')}\"\nsamplingsets={samplingsets}\nexactsolution={exactsolution}\nunknownlabels={unklabs}\n"
        
    return best_cost, best_wgd_nodes, exactsolution, outstats
