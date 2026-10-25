"""Memory-bounded reconciliation of independent outgroup-wrapped gene trees.

An evaluator belongs to one search over a fixed network and fixed gene trees.
Only positive witnesses are reused: U remains valid when U <= allowed episodes.
No state is shared between separate searches or public reconciliation calls.
"""
import gc
from treeop import Tree


def measure_episode_sizes(st, embeddings, root_correction=0):
    """Count a stream of embeddings without retaining it or changing node attrs.

    ``embeddings`` is a zero-argument iterator factory; None denotes infeasibility.
    The correction accounts for joins of independent outgroup-wrapped trees.
    """
    missing = object()
    saved = {n: getattr(n, 'episize', missing) for n in st.nodes}
    try:
        for n in st.nodes:
            n.episize = 0
        st.root.episize = root_correction
        for embedding in embeddings():
            if embedding is None:
                return None
            embedding.setepisize()
            del embedding
        return {n: n.episize for n in st.nodes}
    finally:
        for n, value in saved.items():
            if value is missing:
                delattr(n, 'episize')
            else:
                n.episize = value


def _run_dp(function, *args, **kwargs):
    # Avoid repeatedly scanning long-lived input graphs. Collect the short-lived
    # recursive memoization closures after each tree, restoring the caller's GC.
    enabled = gc.isenabled()
    gc.disable()
    try:
        return function(*args, **kwargs)
    finally:
        gc.collect(0)
        if enabled:
            gc.enable()


class BatchEmbedding:
    def __init__(self, shape, items, root):
        self.shape, self.items, self.s = shape, items, root

    def __str__(self):
        def render(shape):
            if isinstance(shape, int):
                return str(self.items[shape])
            return f"({render(shape[0])},{render(shape[1])})+ {self.s.num}"
        return render(self.shape)

    __repr__ = __str__

    def setepisize(self):
        self.s.episize += len(self.items) - 1
        for item in self.items:
            item.setepisize()

    def embeddings(self, stroot):
        if stroot != self.s:
            return [str(self)]
        return [text for item in self.items for text in item.embeddings(stroot)]


class EpisodeBatch:
    def __init__(self, st, gt, outgroup, parts, shape, feasible, embedding):
        self.st, self.gt, self.outgroup = st, gt, outgroup
        self.parts, self.shape = parts, shape
        self._feasible, self._embedding = feasible, embedding
        self._witnesses = [None] * len(parts)

    def _matches(self, st, gt, outgroup, debug):
        return st is self.st and gt is self.gt and outgroup == self.outgroup and not debug

    def feasible(self, st, gt, wgd_nodes, wgddebug=False, excludedoutgroup=""):
        if not self._matches(st, gt, excludedoutgroup, wgddebug):
            return self._feasible(st, gt, wgd_nodes, wgddebug, excludedoutgroup)
        allowed = set(wgd_nodes)
        # Combining >=2 trees creates duplications at the artificial root.
        if st.root not in allowed:
            return False, set(), ""
        used = {st.root}
        for i, part in enumerate(self.parts):
            witness = self._witnesses[i]
            if witness is None or not witness.issubset(allowed):
                ok, witness, _ = _run_dp(
                    self._feasible, st, part, allowed,
                    excludedoutgroup=excludedoutgroup)
                if not ok:
                    return False, set(), ""
                self._witnesses[i] = frozenset(witness)
            used.update(witness)
        return True, used, str(gt)

    def embedding(self, st, gt, wgd_nodes, wgddebug=False, excludedoutgroup=""):
        if not self._matches(st, gt, excludedoutgroup, wgddebug):
            return self._embedding(st, gt, wgd_nodes, wgddebug, excludedoutgroup)
        allowed = set(wgd_nodes)
        if st.root not in allowed:
            return None
        items = []
        for part in self.parts:
            item = _run_dp(self._embedding, st, part, allowed,
                           excludedoutgroup=excludedoutgroup)
            if item is None:
                return None
            items.append(item)
        return BatchEmbedding(self.shape, items, st.root)

    def episode_sizes(self, st, gt, wgd_nodes, wgddebug=False, excludedoutgroup=""):
        """Same counts as embedding().setepisize(), with one tree in memory."""
        if not self._matches(st, gt, excludedoutgroup, wgddebug):
            return measure_episode_sizes(st, lambda: iter([
                _run_dp(self._embedding, st, gt, wgd_nodes, wgddebug,
                        excludedoutgroup)]))
        allowed = set(wgd_nodes)
        if st.root not in allowed:
            return None

        def embeddings():
            for part in self.parts:
                yield _run_dp(self._embedding, st, part, allowed,
                              excludedoutgroup=excludedoutgroup)

        return measure_episode_sizes(st, embeddings, len(self.parts) - 1)


def make_episode_batch(st, gt, outgroup, feasible, embedding):
    """Return None unless both inputs have the standard independent-tree layout.

    Unknown gene labels and nonstandard outgroups use the ordinary DP instead.
    Input topology and labels must remain fixed while the evaluator is in use.
    """
    if (not outgroup or any(label.startswith('?') for label in gt.root.cluster)
            or len(st.root.c) != 2 or not st.root.c[1].leaf()
            or st.root.c[1].label != outgroup
            or any(n.label == outgroup for n in st.root.c[0].leaves())):
        return None
    roots = []

    def split(node):
        if len(node.c) == 2 and all(outgroup in c.cluster for c in node.c):
            return split(node.c[0]), split(node.c[1])
        index = len(roots)
        roots.append(node)
        return index

    shape = split(gt.root)
    if len(roots) < 2 or any(
            len(n.c) != 2 or not n.c[1].leaf() or n.c[1].label != outgroup
            or outgroup in n.c[0].cluster for n in roots):
        return None
    parts = []
    for node in roots:
        # Read-only views preserve node identity and numbering without copying.
        part = Tree.__new__(Tree)
        part.root, part.nodes = node, node.nodes()
        parts.append(part)
    return EpisodeBatch(st, gt, outgroup, parts, shape, feasible, embedding)
