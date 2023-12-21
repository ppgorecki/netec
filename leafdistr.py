



def emptyleafcounts(gtunknleaves: list, stleaves: list):
	return [ dict((sl.clusterleaf, 0) for sl in stleaves) for gul in gtunknleaves ]


def leafcountsupdate(lfcounts: list, stmap: tuple):	
	for lfd, sl in zip(lfcounts, stmap):
		lfd[sl.clusterleaf]+=1 
