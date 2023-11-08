import sys
from functools import reduce
from random import randint, shuffle


def compcostsmp(g, s, c):
    g.set_lca_mapping(s)
    if c == 'D':
        return g.dupcost(s)
    if c == 'L':
        return g.losscost(s)
    if c == 'U':
        return g.duplosscost(s)
    if c == 'C':
        return g.dccost2(s)
    if c == 'c':
        return g.dccost(s)
    if c == 'R':
        return g.rfcost(s)
    raise Exception("Unknown cost " + c)


def ftread(a, oclass):
    try:
        return [oclass(str2tree(l)) for l in open(a, 'r') if l.strip() and l.strip()[0] != '#']
    except Exception:
        return [oclass(str2tree(x.strip())) for x in a.split(";") if x.strip()]


# string -> list with labels
def tokenizer(s):
    tok = []
    i = 0
    while i < len(s):
        if s[i] in "(),:":
            tok.append(s[i])
            i += 1
        elif s[i].isspace():
            i += 1
        elif s[i].isalnum() or s[i] in "#.+-?_":
            st = i
            while i < len(s) and (s[i].isalnum() or s[i] in "#.-+?_"):
                i += 1
            tok.append(s[st:i])
        elif s[i] == '[':  # comment
            com = ''
            st = i
            while i < len(s) and s[i] != ']':
                i += 1
            i += 1
            tok.append(s[st:i])
        else:
            raise Exception("Non-alphanumeric character: <%s>" % s[i])
    return tok


def str2tree(s):    
    tok = tokenizer(s)

    def _st(tok):
        if not tok:
            raise Exception("String too short..." % s)
        if tok[0] == '(':
            tok.pop(0)
            c = [_st(tok)]

            while tok and tok[0] == ",":
                tok.pop(0)
                c.append(_st(tok))

            if tok[0] != ')':
                raise Exception(") expected found <%s>" % tok[0])
            tok.pop(0)

        else:
            c = []

        # label and :...
        l = []
        while tok and tok[0] not in "(),":
            l.append(tok.pop(0))
        return c, l

    if not tok:
        raise Exception("Syntax error. Tokens left: <%s>" % tok)
    return _st(tok)


def node2label(n):
    return str(n).replace("(", " ").replace(")", " ").replace(",", " ")


class Node:
    def __init__(self, tup, par):

        self.src, self.labels = tup
        self.parent = par

        if self.parent:
            self.height = self.parent.height + 1
        else:
            self.height = 0
        if self.src:
            self.c = [Node(t, self) for t in self.src]
        else:
            self.c = []
        self.interval = None
        self._setcluster()

        self.comments = [c[1:-1] for c in self.labels if c[0] == '[']  # newick comments

        self.branchlengthset = False

        # branch lengths
        if ":" in self.labels:
            p = self.labels.index(":")
            if p + 1 < len(self.labels):
                try:
                    self.branchlength = float(self.labels[p + 1])
                except ValueError:
                    print(f"A number expected. Found: <{self.labels[p + 1]}>", file=sys.stderr)
                    sys.exit(-1)

                self.branchlengthset = True

    def _setcluster(self):
        if self.labels and self.labels[0][0] != '[':
            self.label = self.labels[0]
        else:
            self.label = ''

        if self.leaf():
            self.cluster = frozenset([self.label])
            self.clusterleaf = self.label
        else:
            self.clusterleaf = None
            self.cluster = frozenset().union(*(t.cluster for t in self.c))

    def leaf(self):
        return not self.c

    def setinterval(self, s, distornode):
        if isinstance(distornode, int):
            self.interval = [s, s.ancestor(distornode)]
        else:
            self.interval = [s, distornode]

    def ancestor(self, dist):
        if dist == 0 or not self.parent:
            return self
        return self.parent.ancestor(dist - 1)

    def __str__(self):
        if self.leaf():
            return self.label
        return "(" + ",".join(str(c) for c in self.c) + ")" + self.label

    def netrepr(self):
        if self.leaf():
            return self.clusterleaf
        s = ''
        if self.reticulation:
            s = '#' + self.retid
            # Make sure child of reticulation vertex is processed once
            if hasattr(self, 'visited'):
                del self.visited
                return s
            self.visited = True
        return '(' + ",".join(c.netrepr() for c in self.c) + ')' + s

    def __repr__(self):
        return str(self)

    def nodes(self):  # postorder
        if self.leaf():
            return [self]
        return sum((t.nodes() for t in self.c), [self])

    def get_nodes(self):   # postorder for binary trees
        if self.leaf():
            return [self]
        return self.c[0].get_nodes() + self.c[1].get_nodes() + [self]

    def leaves(self):
        if self.leaf():
            return [self]
        return sum((t.leaves() for t in self.c), [])

    # x = self or x below sel.
    def geq(self, x):
        while x:
            if x == self:
                return True
            x = x.parent
        return False

    def leq(self, x):
        return x.geq(self)

    # added so Python3 can sort nodes
    def __lt__(self, x):
        return x.geq(self) and x != self

    def comparable(self, x):
        return self.geq(x) or self.leq(x)

    def lca(self, y):
        a, b = self, y
        if a.height > b.height:
            a, b = b, a
        # a.h <= b.h
        while b.height != a.height:
            b = b.parent

        # a.h == b.h
        while True:
            if a == b:
                return a
            a = a.parent
            b = b.parent

    def findnode(self, cluster):
        if self.cluster == cluster:
            return self
        if self.leaf():
            return None
        for c in self.c:
            r = c.findnode(cluster)
            if r:
                return r
        return None

    # Helper functions used in DP calculation
    def left_reticulation_used(self):
        return 2 ** (self.retnum * 2 - 1)

    def right_reticulation_used(self):
        return 2 ** (self.retnum * 2)

    def get_node_retusage(self, retusage):
        val = (retusage & self.left_reticulation_used()) + (retusage & self.right_reticulation_used())
        val //= self.left_reticulation_used()
        return val

    def nodemaprepr(self, nodemap):
        """
        Returns str representation where leaf labels are replaced using Node to str map
        """            
        if self in nodemap:
            return nodemap[self]

        if self.leaf():
            return self.label
        return "(" + ",".join(c.nodemaprepr(nodemap) for c in self.c) + ")" + self.label

    def markrepr(self, marknodes):
        """
        Returns str representation where leaf labels are replaced using Node to str map
        """  
        s=""          
        if self in marknodes:
            s=" mark=1"
        
        if self.leaf():
            return self.label+s
        return "(" + ",".join(c.markrepr(marknodes) for c in self.c) + ")" + self.label + s

    def attrrepr(self, attrlabels: list, ignorezeros=True):
        """
        Returns str representation where leaf labels are replaced using Node to str map
        """  
        s=""          
        for attrlabel in attrlabels: 
            if hasattr(self, attrlabel):
                val = getattr(self, attrlabel)
                if not (ignorezeros and not val):            
                    s+=f" {attrlabel}={val}"

        if self.leaf():
            return self.label+s
        return "(" + ",".join(c.attrrepr(attrlabels, ignorezeros) for c in self.c) + ")" + self.label + s





class Tree:
    def __init__(self, tup):
        self.srclist = None

        if isinstance(tup, list):
            self.srclist = tup
            tup = self.srclist[0]
            for t in self.srclist[1:]:
                tup = (tup, t)

        self.root = Node(tup, None)
        self.nodes = self.root.nodes()
        for i, n in enumerate(self.nodes):
            n.num = i
            n.artificial = False
        self.src = tup

        if self.srclist:
            l = self.root
            for i in range(len(self.srclist) - 1):
                l.artificial = True
                l = l.l

    def leaves(self):
        return self.root.leaves()

    def is_binary(self):
        for node in self.nodes:
            if len(node.c) not in (0, 2):
                return False
        return True

    def dccost(self, stree):
        c = 0
        for n in self.nodes:
            if n != self.root:
                c += n.lcamap.height - n.parent.lcamap.height - 1
        return c

    def dccost2(self, stree):
        c = 0
        for n in self.nodes:
            if n != self.root:
                c += n.lcamap.height - n.parent.lcamap.height
        return c

    def rfcost(self, stree):
        return len(self.clusters().symmetric_difference(stree.clusters()))

    def unknownlabels(self):
        return [ l for l in self.leaves() if l.clusterleaf[0]=='?' ]

    def dupcost(self, stree):
        d = 0
        for n in self.nodes:
            if n.leaf():
                continue
            for c in n.c:
                if c.lcamap == n.lcamap:
                    d += 1
                    break
        return d

    def splits(self):
        return [(g.l.cluster, g.r.cluster) for g in self.nodes if not g.leaf()]

    def gencost(self, stree, lambdaf, gammaf):
        sleaves = stree.root.cluster
        gsplits = self.splits()
        ssplits = stree.splits()
        c = sum(lambdaf(a, b, frozenset(s)) for s in sleaves for a, b in gsplits)
        i = sum(gammaf(a, b, x, y) for a, b in gsplits for x, y in ssplits)
        return c + i

    def duplosscost(self, tree):
        return self.dccost(tree) + 3 * self.dupcost(tree)

    def losscost(self, tree):
        return self.dccost(tree) + 2 * self.dupcost(tree)

    def findnodeplus(self, cluster):
        if len(cluster) > 2 and cluster[-2] == "+":
            up = int(cluster[-1])
            cluster = cluster[:-2]
        elif len(cluster) > 1 and cluster[-1][0] == "+":
            up = int(cluster[-1])
            cluster = cluster[:-1]
        else:
            up = 0

        s = self.root.findnode(frozenset(cluster))
        if not s:
            return s

        while up and s.parent:
            s = s.parent
            up = up - 1
        return s

    def __str__(self):
        return str(self.root)

    def set_lca_mapping(self, st):
        def comp_lca_map(n):
            if n.leaf():
                clu = n.clusterleaf
                n.lcamap = None
                while clu and not n.lcamap:
                    n.lcamap = st.root.findnode(frozenset([clu]))
                    clu = clu[:-1]

                if not n.lcamap:
                    raise Exception("Lca mapping not found for ", n)
            else:
                # internal
                n.lcamap = reduce(lambda x, y: x.lca(y), (comp_lca_map(c) for c in n.c))
            return n.lcamap

        comp_lca_map(self.root)

    def __repr__(self):
        return repr(self.root)

    def lcacluster(self, cluster):
        c = set(cluster)
        for n in self.nodes:
            if c.issubset(set(n.cluster)):
                return n
        return None

    def clusters(self):
        return set(n.cluster for n in self.nodes)

    def todotfile(self, f, nodeprefix='g'):
        for n in self.nodes:
            comments = " ".join(n.comments)
            if comments:
                comments = "\n" + comments
            s = ''
            if n.leaf():
                s = str(n)
            f.write("%s%d [label=\"%s %s%s\"];\n" % (nodeprefix, n.num, s, n.num, comments))
        for n in self.nodes:
            if n.parent:
                f.write("%s%d -> %s%d;\n" % (nodeprefix, n.parent.num, nodeprefix, n.num))

    def nodemaprepr(self, nodemap):
        return self.root.nodemaprepr(nodemap)



def getlabs(a, z, labels):
    alfsize = z - a + 1
    letters = [chr(a + i) for i in range(alfsize)]
    return (letters[i % len(letters)] + ("" if i < alfsize else str(i // alfsize)) for i in range(labels))


def randtreestrfromlist(p):
    # gen random tree with all subtrees
    while len(p) > 1:
        a = p.pop(randint(0, len(p) - 1))
        b = p.pop(randint(0, len(p) - 1))
        p.append("(" + a + "," + b + ")")

    return p[0]


def randtreestr(labels, unknown_labels=0, use_numbers=False):
    if use_numbers:
        p = [str(i) for i in range(1, labels + 1)]
    else:
        p = list(getlabs(ord('a'), ord('z'), labels))
    for i in range(unknown_labels):
        p.append("?")
    return randtreestrfromlist(p)


def compatiblesets(x, y):
    return (x.issubset(y) or y.issubset(x)) or not x.intersection(y)


def treestrfromclusters(clusters, binarytree=1):
    """
    Return a tree string from a set of compatible clusters
    Assumption: all trivial clusters are present (root, leaves)
    Multifurcations will be randomly resolved if binarytree=1

    Args:
        clusters: a set of clusters (sets of labels)
        binarytree: flag indicating that the tree must be binary
    """

    clusters = sorted(clusters, key=lambda x: -len(x))
    tree = {}  # cluster -> its subtree
    while clusters:
        n = clusters.pop()

        if len(n) == 1:
            tree[n] = list(n)[0]  # leaf
        else:
            children = [c for c in tree if c.issubset(n)]  # gen children
            subtrees = [tree[c] for c in children]
            if binarytree:
                tree[n] = randtreestrfromlist(subtrees)  # binarize randomly
            else:
                tree[n] = "(" + ",".join(subtrees) + ")"  # allow mutlifurcations

            for c in children:
                del tree[c]  # remove used subtrees

    return tree.pop(n)


def quasiconsensusgen(gtrees, qtreescnt):
    """
    Generate quasi-consesus tree strings from a set of gene trees

    Args:
        gtrees: a collection of gene trees
        qtreescnt: how many quasi-consensus tree to generate
    """

    d = {}
    alllabels = set()
    for g in gtrees:
        for n in g.nodes:
            d[n.cluster] = d.setdefault(n.cluster, 0) + 1
            if len(n.cluster) == 1:
                alllabels = alllabels.union(n.cluster)
    alllabels = frozenset(alllabels)

    counts = [[] for i in range(1 + max(d.values()))]
    for k in d:
        counts[d[k]].append(k)

    for _ in range(qtreescnt):
        curclusters = set([alllabels])
        for q in reversed(counts):
            shuffle(q)  # some randomness
            for c in q:
                if all(compatiblesets(c, cc) for cc in curclusters):
                    curclusters.add(c)

        # compatible clusters are generated
        yield treestrfromclusters(curclusters)

