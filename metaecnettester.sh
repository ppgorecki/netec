
# Tester for network based metaec

RETICULATIONS=3
N=5
DESTDIR=metaecnettest
GTREES=3

mkdir -p $DESTDIR
src/embnet.py -n rand:$N:$RETICULATIONS -pn > $DESTDIR/net.txt

(( K=N+4 ))
echo $K
for i in `seq $GTREES`
do
	src/embnet.py -g rand:$K:0 -pg 
done | sed "s/[fghi]/?/g" > $DESTDIR/gtrees.txt

./metaec.py --gene_trees $DESTDIR/gtrees.txt --network $DESTDIR/net.txt --out_file $DESTDIR/metaec.log

