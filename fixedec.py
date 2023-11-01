
from typing import Set, Tuple, List
from treeop import Tree, Node, str2tree

def fixedec(genetree: Tree, speciestree: Tree) -> Set[Node]:
	"""
	Returns a set of fixed duplication episode nodes from the species tree
	"""

	if len(speciestree.leaves())==1: return set()

	stleft, stright = speciestree.root.c[0].cluster, speciestree.root.c[1].cluster

	lab2leaf = dict( (s.clusterleaf,s) for s in speciestree.leaves() ) 

	fixedepisodes = set()

	
	def gtraverse(g):
		"""
		Given a gene tree node
		Returns a tuple:
			Has?:bool - True if ? is reachable from g
			cluster - set of leaf labels excluding ?
			stmap - map to the species tree ? is ignored, None if cluster is empty
			dupl - True is node is a duplication (only if ? is not present)
			fixedEC - fixed EC episodes
		"""
		if g.leaf():
			if g.clusterleaf[0] == '?':
				return (True,frozenset(),None, False)
			else:
				return (False, frozenset([g.clusterleaf]), lab2leaf[g.clusterleaf], False)				
		 
		hasunk, labels, stmap, dupl = False, frozenset(), None, False
		chmaps = [] # all maps
		chdupl = [] # maps of dupl. children
		for c in g.c:
			_hasunk, _labels, _stmap, _dupl = gtraverse(c)			
			dupl = dupl or _dupl # child is a duplication
			chmaps.append(_stmap)
			if _dupl: 
				chdupl.append(_stmap) 
			labels = labels.union(_labels)
			hasunk = hasunk or _hasunk			
			if _stmap is not None: 
				if stmap is None: 
					stmap = _stmap
				else:
					stmap = stmap.lca(_stmap)			

		# fixed subtree case
		if not hasunk and dupl:			
			# no ? in the subtree
			# there is a duplicated child 
			# if g speciation directly above dupl. child map
			if stmap not in chmaps: # g is a speciation
				for d in chdupl:
					if d.parent == stmap: 
						fixedepisodes.add(d)
						break

		# root case
		if stmap == speciestree.root and stmap in chmaps:
			fixedepisodes.add(stmap)
		return hasunk, labels, stmap, not hasunk and stmap in chmaps

	gtraverse(genetree.root)	

	return fixedepisodes

if __name__ == "__main__":
	print(fixedec(Tree(str2tree("(((a,a),b),c)")),Tree(str2tree("(c,(a,b))")))) # (a,b)  
	print(fixedec(Tree(str2tree("(((a,a),(b,b)),b),c)")),Tree(str2tree("(c,(a,b))"))) ) # (a,b)
	print(fixedec(Tree(str2tree("(((a,?),(b,b)),b),c)")),Tree(str2tree("(c,(a,b))")))) # empty
	print(fixedec(Tree(str2tree("((a,?),b)")),Tree(str2tree("(a,b)")))) # empty 
	print(fixedec(Tree(str2tree("(((a,?),b),b)")),Tree(str2tree("(a,b)")))) # root episode
	print(fixedec(Tree(str2tree("(((a,?),b),?)")),Tree(str2tree("(a,b)")))) # root episode
	print(fixedec(Tree(str2tree("(((a,?),?),?)")),Tree(str2tree("(a,b)")))) # empty
	