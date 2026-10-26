"""Iterative post-evaluation from Algorithm 2 of btag445 (Section 3.3).

One DAG post-order pass, updating Q immediately after each accepted candidate.
The threshold uses the initial embedding's mean across ALL biological non-leaf
nodes (zero sizes included). Artificial outgroup/root nodes are excluded.
"""
import csv
import json
import math
from pathlib import Path


def validate_iterative_options(tau, optimize='full', fixed_episodes_only=False):
    if not math.isfinite(tau) or tau < 0:
        raise ValueError('extended_episodes_tau must be finite and nonnegative')
    if fixed_episodes_only or optimize not in ('full', 'single'):
        raise ValueError('Iterative extended analysis requires a baseline inference: '
                         'use --optimize full or single without --fixed_episodes_only')


def postorder(root):
    """Deterministic child-before-parent order, visiting reticulations once."""
    done, active, result = set(), set(), []

    def visit(node):
        if node in active:
            raise ValueError('Iterative extended analysis requires an acyclic network')
        if node in done:
            return
        active.add(node)
        for child in sorted(node.c, key=lambda n: n.num):
            visit(child)
        active.remove(node)
        done.add(node)
        result.append(node)

    visit(root)
    return result


def iterative_extended_analysis(network, episodes, measure, tau=2.0, progress=None):
    """Return an auditable report; measure(Q) must return node->size or None.

    Q is the feasible result of the baseline inference (or a supplied scenario).
    Every node outside Q is tested once, including leaves and reticulations.
    Acceptance is strict: size(s) > tau * mean_nonleaf_size(initial Q).
    The mean is computed once and remains fixed throughout the pass.
    The initial episode set is never reduced or optimized again in this phase.
    """
    validate_iterative_options(tau)
    order = postorder(network.root)
    q = set(episodes)
    if not q.issubset(set(order)):
        raise ValueError('Baseline episodes must be nodes of the input network')
    artificial = set()
    if network.outgrouped():
        artificial = {network.root, network.root.c[1]}
    biological = set(order) - artificial
    nonleaves = [n for n in order if n in biological and not n.leaf()]
    initial = set(q)
    current_sizes = measure(q)
    if current_sizes is None:
        raise ValueError('Iterative extended analysis requires a feasible baseline episode set')
    mean = sum(current_sizes.get(n, 0) for n in nonleaves) / len(nonleaves) if nonleaves else 0.0
    threshold = tau * mean
    steps, added = [], []
    candidates = [n for n in order if n in biological and n not in initial]
    for index, node in enumerate(candidates, 1):
        trial = q | {node}
        sizes = measure(trial)
        if sizes is None:
            raise RuntimeError('A superset of a feasible episode set was reported infeasible')
        size = sizes.get(node, 0)
        accepted = size > threshold
        before = len(q)
        if accepted:
            q.add(node)
            added.append(node)
            current_sizes = sizes
        step = dict(step=index, node=node.num, episode_size=size,
                    mean_nonleaf_size=mean, threshold=threshold,
                    accepted=accepted, episodes_before=before, episodes_after=len(q))
        steps.append(step)
        if progress:
            progress(step, len(candidates))
    nums = lambda nodes: sorted(n.num for n in nodes)
    return dict(algorithm='btag445 Algorithm 2', tau=tau, traversal='postorder',
                mean_policy='initial_embedding_all_biological_nonleaf_nodes_including_zeros',
                baseline_mean_nonleaf_size=mean, threshold=threshold,
                nonleaf_count=len(nonleaves), artificial_nodes=nums(artificial),
                baseline_episodes=nums(initial), episodes=nums(q),
                biological_episodes=nums(q & biological),
                added_episodes=[n.num for n in added],
                episode_sizes={str(n.num): current_sizes.get(n, 0) for n in order},
                size_dp_evaluations=1 + len(candidates), steps=steps)


def write_iterative_results(report, network, out_basefile, gsestyle=False):
    """Save the extended result separately, preserving baseline node attributes."""
    base = Path(out_basefile)
    base.mkdir(parents=True, exist_ok=True)
    (base / 'iterative_extended.json').write_text(json.dumps(report, indent=2) + '\n')
    (base / 'iterative_episodes').write_text(' '.join(map(str, report['episodes'])) + '\n')
    (base / 'iterative_added_episodes').write_text(' '.join(map(str, report['added_episodes'])) + '\n')
    fields = ['step', 'node', 'episode_size', 'mean_nonleaf_size', 'threshold',
              'accepted', 'episodes_before', 'episodes_after']
    with (base / 'iterative_steps.tsv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(report['steps'])
    missing = object()
    attributes = ['episize', 'epibestwgd', 'iterative_episode', 'iterative_added_step']
    saved = {(n, a): getattr(n, a, missing) for n in network.nodes for a in attributes}
    accepted_at = {row['node']: row['step'] for row in report['steps'] if row['accepted']}
    try:
        for n in network.nodes:
            n.episize = report['episode_sizes'][str(n.num)]
            n.epibestwgd = int(n.num in report['baseline_episodes'])
            n.iterative_episode = int(n.num in report['episodes'])
            n.iterative_added_step = accepted_at.get(n.num, 0)
        root = network.trueroot if network.outgrouped() else network.root
        suffix = '.gse' if gsestyle else '.newick'
        text = root.attrrepr(['num'] + attributes, ignorezeros=False, gsestyle=gsestyle)
        (base / ('iterative_episummary' + suffix)).write_text(text + '\n')
    finally:
        for (node, attr), value in saved.items():
            if value is missing:
                delattr(node, attr)
            else:
                setattr(node, attr, value)
