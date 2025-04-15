
from typing import Set, Tuple, List
from treeop import Tree, Node, str2tree
from netop import Network

def findfixedepisodes(g, lab2leaf, stroot, fixedepisodes=None):
	
	if fixedepisodes is None: fixedepisodes = set()
		
	def trav(g):		
		"""
		Given a gene tree node
		Returns a tuple:
			Has?:bool - True if ? is reachable from g
			cluster - set of leaf labels excluding ?
			stmap - map to the species tree ? is ignored, None if cluster is empty
			dupl - True if node is a duplication (only if ? is not present)
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
			_hasunk, _labels, _stmap, _dupl = trav(c)			
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
						# break, bug no break here 					

		# root case
		if stmap == stroot and stmap in chmaps:
			fixedepisodes.add(stmap)

		return hasunk, labels, stmap, not hasunk and stmap in chmaps

	trav(g)
	return fixedepisodes

def fixedec(genetree: Tree, speciestree: Tree) -> Set[Node]:
	"""
	Returns a set of fixed duplication episode nodes from the species tree
	"""

	if len(speciestree.leaves())==1: return set()

	lab2leaf = dict( (s.clusterleaf, s) for s in speciestree.leaves() ) 
	
	return findfixedepisodes(genetree.root, lab2leaf, speciestree.root)		

def fixedecnet(gt, net):

	# Locate subtrees in the network

	fixedepisodes = set()
	# for b in net.root.subnettrees()[1]:
	#  	print("SBTREE", b,"===>", b.cluster)		

	# for sbtree in net.root.subnettrees()[1]:
	# 	#print("Loop",sbtree, gt.root.subtrees(sbtree.cluster))
	# 	lab2leaf = dict( (s.clusterleaf, s) for s in sbtree.leaves() ) 		
	# 	# Locate subtrees matching sbtree in the gene tree
	# 	for gnode in gt.root.subtrees(sbtree.cluster):
	# 		#print("search in ",gnode, sbtree, lab2leaf)
	# 		fixedepisodes = fixedepisodes.union(findfixedepisodes(gnode, lab2leaf, None))

	tripleclusters = {}
	def cluprint(clu):
		if not clu: return '-'
		return "".join(sorted([str(x) for x in clu]))			

	def nprint(clu):
		if not clu: return '-'
		return ".".join(map(str,sorted([x.num for x in clu])))
		
	# Locate all network nodes with no reticulation above
	reticulated = set()
	for s in net.nodes:		
		s.nrc = 0 
		if s.reticulation: 
			reticulated = reticulated.union(s.nodes())			
	nonreticulated = set(net.root.nodes()).difference(reticulated)


	
	for s in nonreticulated:	
		if not s.leaf():
			ln = set(l.clusterleaf for l in s.c[0].leaves())
			rn = set(l.clusterleaf for l in s.c[1].leaves())
			it = ln.intersection(rn)
			ln,rn = ln.difference(rn), rn.difference(ln)			
			tripleclusters[s] = (ln,it,rn)
			#print("TRIPLE",s.num, s, cluprint(ln), cluprint(it), cluprint(rn))
			s.nrtc = 1
	
	
	for g in gt.nodes:
		if g.leaf(): continue
		gl, gr = g.c	
		#print("GNODE", g.num, g, cluprint(g.leaves()))			
		for s in tripleclusters:
			ln, it, rn = tripleclusters[s]
			#print ("Testing", s.num, g.num, cluprint(ln), cluprint(it), cluprint(rn))
			for gc in g.c:
				if gc.leaf(): continue
				if gc.cluster.issubset(ln) and gc.sibling().cluster.issubset(rn):
					# g is s-speciation
					#print(" MATCH g-spec",s.num, g.num, g, cluprint(g.leaves()), cluprint(ln), cluprint(it), cluprint(rn))

					for sc in s.c:
						#print("  SC",sc.num, sc, cluprint(sc.leaves()))
						# check map to s-child (sc)
						if sc.leaf(): # special case
							if gc.c[0].cluster == sc.cluster and gc.c[1].sibling().cluster == sc.cluster:
								fixedepisodes.add(sc)
						elif sc in tripleclusters:
							lnc, itc, rnc = tripleclusters[sc]
							#print("  C-S",sc.num, sc, cluprint(lnc), cluprint(itc), cluprint(rnc))
							for gcc in gc.c:
								#print("  Gtest",gcc, gcc.cluster, gcc.sibling().cluster)
								if gcc.cluster.intersection(lnc) and gcc.cluster.intersection(rnc):
									#gcc maps to sc
									#gc is sc-duplication; gcc maps to sc
									#print("  INSERT!",s.num)
									fixedepisodes.add(sc)
																


			
				
				# fixedepisodes.add(s)
				# break

	for fe in fixedepisodes:
		fe.fixedepi = 1

	#print("&g",gt)
	#print("&s",net.root.attrrepr(['num','nrtc','fixedepi'],gsestyle=1))


	

	# Find cands in gtree for the root episode (true root, not outgrouped)	# if s is the true-root of S, then at least one proper subtree of
	# G contains species (leaf-labels) from both rndiff and lndiff

	
		

	return fixedepisodes


if __name__ == "__main__":

	net1r = "((c)#A,((#A,b),a))"
	net1rb = "((c)#A,((#A,b),(a,d)))"
	gt2 = "(((((a,a),a),d),c),o)"
	gt1 = "(((((a,a),a),d),e),o)"
	gt1 = "(((((a,d),a),d),e),o)"
	gt1 = "((((a,i),a),e),o)"
	net8='(((a,(((d)#B,i), ((h,f),( ((#B,(g,c)),j) )#A) )),((#A,b),e)),o)'
	#net8 = '((a,((d)#B,(((#B,g))#A,(h,f)))),((e,c),(#A,b)))'

	amnet = open('amoutnet.newick').read()
	amgt = open('amoutgt.newick').read()

	def execnetec(g,n):
		fe = fixedecnet(Tree(str2tree(g)),net:=Network(str2tree(n)))
		with open('netec_emb/gtfixedectest.newick','w') as f:
			f.write(g)
		
		with open('netec_emb/netfixedectest.newick','w') as f:
			f.write(net.root.attrrepr(['num','nrtc','fixedepi'],gsestyle=0))
		return fe


	for fixedepi in execnetec(amgt,amnet):
	#for fixedepi in execnetec(gt1, net8):
		print("FIXEDEPI",fixedepi)
	

	#for gtree, net in [(gt1, net1rb)]: # episode at a
	#	print(fixedecnet(Tree(str2tree(gtree)),Network(str2tree(net))))  

	#print(fixedecnet(Tree(str2tree(open('amoutgt.newick').read())),Network(str2tree(open('amoutnet.newick').read()))))  
	
	# print(fixedec(Tree(str2tree("(((a,a),b),c)")),Tree(str2tree("(c,(a,b))")))) # a  
	# print(fixedec(Tree(str2tree("((((a,a),(b,b)),b),c)")),Tree(str2tree("(c,(a,b))"))) ) # a,b, (a,b)
	# print(fixedec(Tree(str2tree("((((a,?),(b,b)),b),c)")),Tree(str2tree("(c,(a,b))")))) # empty
	# print(fixedec(Tree(str2tree("((a,?),b)")),Tree(str2tree("(a,b)")))) # empty 
	# print(fixedec(Tree(str2tree("(((a,?),b),b)")),Tree(str2tree("(a,b)")))) # root episode
	# print(fixedec(Tree(str2tree("(((a,?),b),?)")),Tree(str2tree("(a,b)")))) # root episode
	# print(fixedec(Tree(str2tree("(((a,?),?),?)")),Tree(str2tree("(a,b)")))) # empty




	