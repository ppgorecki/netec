import itertools
import sys
from itertools import product
from collections import Counter
from random import randint, shuffle
from timeit import default_timer as timer
from typing import Tuple, Dict, Set, List
from treeop import Tree, Node, str2tree, getlabs, randtreestr, compcostsmp
import queue

class Network(Tree):
    def __init__(self, tup):
        Tree.__init__(self, tup)

        # recognize reticulations
        dlf = {}
        din = {}
        err = 0
        for i, n in enumerate(self.nodes):
            if n.label and n.label[0] == "#":
                retid = n.label[1:]
                n.reticulation = 1
                if n.leaf():
                    if retid in dlf:
                        print("Reticulation id <%s> already defined" % retid)
                        err = 1
                    dlf[retid] = n
                    n.reticulationleaf = 1
                else:
                    if retid in din:
                        print("Reticulation id <%s> already defined" % retid)
                        err = 1
                    din[retid] = n
                n.retid = retid
            else:
                n.reticulation = 0
                n.retid = 0

        if set(dlf.keys()) != set(din.keys()) or err:
            for k in din:
                if k not in dlf:
                    print("Missing reticulation label <%s> in leaves" % k)
            for k in dlf:
                if k not in din:
                    print("Missing reticulation label <%s> in internal nodes" % k)

            sys.exit(-1)

        self.reticulations = []
        for retnum, retid in enumerate(dlf, 1):
            l = dlf[retid]
            i = din[retid]
            lpar = l.parent
            ipar = i.parent

            i.lparpos = lpar.c.index(l) 
            # lpar.c.remove(l)
            #i.iparpos = ipar.c.index(i) 
            # ipar.c.remove(i)

            lpar.c[i.lparpos] = i  # inserted at i
            #lpar.c.insert(0, i)  # inserted at 0
            #ipar.c.insert(0, i)  # inserted at 0

            i.lftparent = ipar
            i.rghparent = lpar
            i.retnum = retnum

            if i.branchlengthset and l.branchlengthset:  # set only if both are defined
                i.branchlengthset = True
                i.branchlength = (i.branchlength, l.branchlength)
            else:
                i.branchlengthset = False

            self.reticulations.append(i)
            self.nodes.remove(l)

        for n in self.nodes:
            n._setcluster()  # reconstruct clusters

    def isdag(self):
        return sorttop([(n.num, c.num) for n in self.nodes for c in n.c])

    def istimeconsistent(self):
        lftparent_nums = []
        edges = []

        # glue reticulation parents
        for r in self.reticulations:
            lftparent_nums.append(r.lftparent.num)
            r.lftparent.num = r.rghparent.num

        # create edge representation of graph
        for n in self.nodes:
            if n == self.root:
                pass
            elif n.reticulation:
                edges.append((n.lftparent.num, n.num))
            else:
                edges.append((n.parent.num, n.num))

        # unglue reticulation parents
        for i, n in enumerate(self.reticulations):
            n.lftparent.num = lftparent_nums[i]

        return bool(sorttop(edges))

    def ptab(self):
        print("=" * 80)
        for n in self.nodes:
            print(n.num, n.reticulation, n.retid, end='')
            if n.reticulation:
                if hasattr(n, "lftparent"):
                    print("rtp:%d,%d" % (n.lftparent.num, n.rghparent.num), end='')
            elif n.parent:
                print("par:%d" % n.parent.num, end='')
            print('==<<', n.netrepr(), ">>", end='')
            if n.leaf():
                print("LEAF", end='')
            else:
                print("Children=", end='')
                for c in n.c:
                    print(c.num, end=' ')
            print()

    def todotfile(self, f, nodeprefix="", addnodeidcomments=True):
        for n in self.nodes:
            comments = "\n".join(n.comments)
            if comments:
                comments = "\n" + comments
            f.write(nodeprefix)
            if addnodeidcomments:
                nodeidcom = f" {n.num}" + comments
            else:
                nodeidcom = ""
            if n.leaf():
                f.write(f'{n.num} [label="{n}{nodeidcom}"]; #leaf\n')
            else:
                if n.reticulation:
                    f.write(f'{n.num} [shape=box,label="#{n.retid}{nodeidcom}"]; #retic\n')
                else:
                    f.write(f'{n.num} [label="{nodeidcom}"]; #inner\n')
        for n in self.nodes:
            if n.reticulation:
                if hasattr(n, "lftparent"):
                    f.write(
                        f'{nodeprefix}{n.lftparent.num} -> {nodeprefix}{n.num} [color=green,label="{n.retid}l"]; #retedge\n')
                    f.write(
                        f'{nodeprefix}{n.rghparent.num} -> {nodeprefix}{n.num} [color=blue,label="{n.retid}r"]; #retedge\n')
            elif n.parent:
                f.write(f'{nodeprefix}{n.parent.num} -> {nodeprefix}{n.num}; #edge\n')
                # f.write("%s%d -> %s%d; #edge\n" % (nodeprefix, n.parent.num, nodeprefix, n.num))

    def get_leaves(self):
        """Get nodes of 1-indegree and 0-outdegree."""
        return [node for node in self.nodes if node.leaf()]

    def get_labels(self):
        """Get labels of all leaves."""
        return [node.label for node in self.get_leaves()]

    def get_inner_nodes(self):
        """Get all nodes (tree and reticulation) except of leaves."""
        return [node for node in self.nodes if not node.leaf()]

    def get_inner_tree_nodes(self):
        """Get nodes of 1-indegree and 2-outdegree plus root."""
        return [node for node in self.get_inner_nodes() if not node.reticulation]

    def valid_binary(self):
        """
        Sanity check whether network is DAG, bijective, binary nodes and
        reticulations, no parallel edges, unique root. Some of these properties
        are partially assured by the class constructor.

        Returns:
            bool for satisfying conditions
        """

        # check if network directed acyclic
        if not self.isdag():
            return False

        # check if leaves are bijective
        labels = self.get_labels()
        if len(labels) != len(set(labels)):
            return False

        # check if binary
        reticulations = self.reticulations
        for node in reticulations:
            if len(node.c) != 1:
                return False

        inner_tree_nodes = self.get_inner_tree_nodes()
        for node in inner_tree_nodes:
            if len(node.c) != 2:
                return False

        # check for parallel edges
        for node in reticulations:
            if node.lftparent == node.rghparent:
                return False

        for node in inner_tree_nodes:
            if node.c[0] == node.c[1]:
                return False

        return True

    def unfold(self):

        def _unfold(n):
            if n.leaf():
                return ([], [n.label, f"[n{n.num}]", {"netsrc": n}])
            return (tuple(_unfold(c) for c in n.c),
                    [f"[n{n.num} {('#' + str(n.reticulation)) if n.reticulation else ''}]", {"netsrc": n}])

        return _unfold(self.root)

    def treechild(self):
        """
        Check if network is a valid binary network and belongs to the
        Tree-Child class.

        Returns:
            bool for satisfying conditions
        """

        # check if valid network

        if not self.valid_binary():
            return False

        # check if each inner node has a tree or leaf child node

        for node in self.nodes:
            if not node.leaf() and all(child in self.reticulations for child in node.c):
                return False

        return True

    def type1net(self):
        """
        Check if network is a valid binary network such that
        no node has >= two reticulation parents

        Returns:
            bool for satisfying conditions
        """

        # check if valid network

        if not self.valid_binary():
            return False

        # check if each inner node has a tree or leaf child node

        for node in self.nodes:
            if not node.leaf() and len(node.c) > 1 and all(child in self.reticulations for child in node.c):
                return False

        return True

    def __repr__(self):
        return self.root.netrepr()

    def __str__(self):
        return self.root.netrepr()



    def displayedtreebyid(self, displayedtreeid, addnoderef=False):
        """
        Return display tree via id
        """
        ign = 0
        fmt = "{0:0%db}" % len(self.reticulations)  # format with leading zeros

        try:
            return self._bltree(
                *self._displtree(dict(zip(self.reticulations, map(int, fmt.format(displayedtreeid)))), self.root, None, addnoderef))
        except Exception as e:
            print(e)

            ign += 1

        if ign:
            print(f"{ign} tree(s) ignored. Is your network tree-child?", file=sys.stderr)

    # insert bl if present
    def _bltree(self, t, blset, branchlength):
        return t + (f":{branchlength}" if blset else "")

    # generic display tree generator via dictionary of ret. usages
    # given vector of reticulation usages b, a node n and parent node par from which n is reached
    # returns a tuple: (display tree encoded in string, branchlengthset_flag, branchlength ) 
    def _displtree(self, b, n, par, addnoderef=False):

        if n.branchlengthset:
            bl = n.branchlength
        else:
            bl = 0

        noderef = f"[ref={self.nodes.index(n)}]" if addnoderef else ""

        if n.leaf():
            return n.label + noderef, n.branchlengthset, bl

        if n in b:  # reticulation;  aggregate branch lengths
            t = None
            if b[n]:
                if n.lftparent == par:
                    r = self._displtree(b, n.c[0], n, addnoderef=addnoderef)
                    if not r:
                        return None
                    t, blsetc, blc = r
                    if n.branchlengthset:
                        bl = bl[0] + blc
            else:
                if n.rghparent == par:
                    r = self._displtree(b, n.c[0], n, addnoderef=addnoderef)
                    if not r:
                        return None
                    t, blsetc, blc = r
                    if n.branchlengthset:
                        bl = bl[1] + blc

            if t:
                # skip node ref
                if n.branchlengthset and blsetc:
                    return t, True, bl
                return t, False, 0 

            return None
        else:
            l = [self._displtree(b, c, n, addnoderef=addnoderef) for c in n.c]
            l = [s for s in l if s]

            if not l:
                return ''

            if len(l) == 1:
                t, blsetc, blc = l[0]
                if n.branchlengthset and blsetc:
                    return t , True, bl + blc
                return t , False, 0
            s = ",".join(self._bltree(*t) for t in l)

            return "(" + s + ")" + noderef, n.branchlengthset, bl

    def displayedtrees(self):

        ign = 0

        for b in product([0, 1], repeat=len(self.reticulations)):
            try:
                yield self._bltree(*self._displtree(dict(zip(self.reticulations, b)), self.root, None))
            except Exception:
                ign += 1
        if ign:
            print(f"{ign} tree(s) ignored. Is your network tree-child?", file=sys.stderr)

    def _add_reticulations_usages(self, usage1: Dict[Node, int], usage2: Dict[Node, int]) -> Dict[Node, int]:
        usage = {}
        for retnode in self.reticulations:
            usage[retnode] = usage1.get(retnode, 0) | usage2.get(retnode, 0)
        return usage

    def contractret(self, reticulation, preserveleftedge=True):
        """
        Return contracted network str by removing reticulation edge
        """
        # print()
        # print("CONTRACT", reticulation, preserveleftedge, self)
        r = self.root.contractret(reticulation, preserveleftedge, None)
        # print("=======>", reticulation, preserveleftedge, r)

        return r

    def contractbyleaves(self, labels, keep=False):
        """
        Return contracted network str by removing leaves
        """
        return self.root.contractbyleaves(labels, keep, None)

    def findredundantreticulation(self):
        """
        Return one redundant reticulation id or None if such a reticulation is not present
        """
        return self.root.findredundantreticulation()

    def sortedrepr(self):
        d = {(n, None) for n in self.nodes}
        self.root._setsortedrepr(d)


# input: list of directed edges, must be non-empty
# return: None if cycle or empty graph
#         top. sort otherwise
def sorttop(e):
    l = []
    d1 = {}
    d2 = {}
    for x, y in e:
        d1.setdefault(x, []).append(y)
        d2.setdefault(y, []).append(x)
    roots = list(set(d1.keys()).difference(d2.keys()))
    while roots:
        n = roots.pop()
        l.append(n)
        if n not in d1:
            continue
        ms = d1[n][:]
        for m in ms:
            d1[n].remove(m)
            d2[m].remove(n)
            if not d2[m]:
                roots.append(m)
                d2.pop(m)
            if not d1[n]:
                d1.pop(n)
    if d1 or d2:
        return None
    return l


def addretstr(s, reticulations, skip=0, networktype=0, time_consistent=False):
    """
    Exhaustively add given number of reticulations to a tree or tree
    representation of network.

    networktype=0 -> treechild
    networktype=1 -> nontreechild type 1
    networktype=2 -> general
    
    skip is how many reticulations label to skip, e.g. skip=2 omits 'A' and 'B'
    time_consistent=True for networks suitable for HGT model

    Returns None if the network cannot be constructed.
    """

    treechild = networktype == 0

    t = Tree(str2tree(s))
    v = t.nodes.copy()

    reticulations = ['#' + i for i in getlabs(ord('A'), ord('Z'), reticulations + skip)]
    reticulations = reticulations[skip:]
    inserted = []

    # quite ugly code
    while reticulations:

        r = reticulations.pop()
        inserted.append(r)

        # insert a leaf labelled <r>

        while True:
            n = v[randint(0, len(v) - 1)]  # root is allowed
            if not treechild or n.label not in inserted:
                break

        np = n.parent
        ap = Node(([], []), np)
        a = Node(([], [r]), ap)

        if np:
            np.c.remove(n)
            np.c.append(ap)
            v.append(ap)
        else:
            t.root = ap
            v.insert(0, ap)  # new root

        n.parent = ap
        ap.c = [a, n]
        v.append(a)

        vshuffled = v[1:]  # skip the root
        shuffle(vshuffled)

        while vshuffled:

            # insert internal node with the label <r>; avoid cycles
            m = vshuffled.pop()

            mp = m.parent
            b = Node(([], [r]), mp)
            mp.c.remove(m)
            mp.c.append(b)
            b.c = [m]
            b.parent = mp
            m.parent = b
            v.append(b)  # not the last

            # check if dag
            s = str(t)
            net = Network(str2tree(s))  # can be done without using network (todo)
            if net.isdag():
                # dag is found
                # check if treechild if needed
                if networktype == 2 or treechild and net.treechild() or networktype == 1 \
                        and net.type1net():
                    # check if time_consistent if needed
                    if not time_consistent or time_consistent and net.istimeconsistent():
                        break  # OK; next reticulation is OK, end while

            # clean and try again
            v.remove(b)
            m.parent = mp
            mp.c.remove(b)
            mp.c.append(m)

        else:

            # clean; this position of leaf labelled <r> wasn't sucessful
            v.remove(a)
            v.remove(ap)
            if np:
                np.c.append(n)
                np.c.remove(ap)
            else:
                t.root = v[0]
            n.parent = np

            # try again with different position of leaf reticulation
            reticulations.append(r)

    return s  # return tree representation of the dag


def randdagstr(leaves, reticulations, networktype=0):
    """
    Return a network with given number of leaves and reticulations.
    If treechild==0 the network has tree-child property.
    networktype=1 -> nontreechild type 1
    networktype=2 -> general

    Returns None if the network cannot be constructed
    """


    if networktype == 0 and reticulations >= leaves:
        return None

    s = randtreestr(leaves)

    if not reticulations:
        return s
    return addretstr(s, reticulations, networktype=networktype)
