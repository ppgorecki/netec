
# Tester for network based netec

RETICULATIONS=2
N=5
DESTDIR=netecnettest
GTREES=5

mkdir -p $DESTDIR
src/embnet.py -n rand:$N:$RETICULATIONS -pn > $DESTDIR/net.txt

(( K=N+4 ))
echo $K
for i in `seq $GTREES`
do
	src/embnet.py -g rand:$K:0 -pg 
done | sed "s/[fghi]/?/g" > $DESTDIR/gtrees.txt

./netec.py --gene_trees $DESTDIR/gtrees.txt --network $DESTDIR/net.txt --out_file $DESTDIR/netec.log

