export PYTHONPATH=$PYTHONPATH:..
REP=500
for UNK in 1 2 3
do
	echo Running test 1 rep=$REP unk=$UNK
	python3 naive2dpmetaecmaps.py '(a,((b,e),(c,(d,f))))' '((?,?),(c,(?,a)))' $REP $UNK 
	echo Running test 2 rep=$REP unk=$UNK
	python3 naive2dpmetaecmaps.py '(a,(b,(c,d)))' '((?,?),(?,d))' $REP $UNK
	python3 naive2dpmetaecmaps.py '(a,(b,(c,d)))' '((a,b),(c,d))' 1 $UNK
done

