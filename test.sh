

# To show values of delta's add number of gene tree nodes in -d 

cat << EOF > nets.txt
n1abc
((a, ((b)#Y)#H ), (#H, (c,#Y))) [Unstable network N1]
((a, b), (c, c))
n1ab
((a, ((b)#Y)#H ), (#H, (c,#Y))) [Unstable network N2]
(a, b)
n2
( ((a, (b)#X) , (#X,(c)#Y) ),   (#Y ,d)  ) [Unstable network N3]
((a, b), (c, d))
single1
(((a,#A),(b)#A ),c) [S1]
(a,(b,c))
single2
( (b,(a)#A), (#A,c)) [S2]
(a,(b,c))
EOF

cat nets.txt | while read F && read N && read G
do
	
	[ ${F:0:1} == '#' ] && continue # skip if F[0] = #
	echo "python3 embnet.py -n \"$N\" -g \"$G\"  -pe gG -d \"$F\""


	if python3 embnet.py -n "$N" -g "$G"  -p egG -d "$F:0:1:2:3:4"
	then 
	# #neato -Tpdf 0.dot -o 0.pdf
		dot -Tpdf $F.dot -o $F.pdf && evince $F.pdf & 
		#cat $F.dot
	fi
done

