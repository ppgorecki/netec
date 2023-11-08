export PYTHONPATH=$PYTHONPATH:$PWD
export DT=$PWD
export DATADIR=data_sim
export RANDOMIZE_FROM=100
export NOIMPROVEMENT_STOP=1000
export FORCE=

JOBS=14 

function Usage()
{
cat <<EOF

Perform MetaEC multiset analysis and output results in the specified output directory as a "wgd.csv" file.

Usage: $0 [-r] [-h] [-d DATADIR] [-j NUM] [-o DIR] [-f] [-c] [-R RANDOMIZE_FROM] [-N NOIMPROVEMENT_STOP] GENETREEFILES...

-h Display help information.
-r Perform reverse climb.
-o DIR Output directory (default: resultsTmp_$(basename $DATADIR)).
-d DATADIR Data directory (default: $DATADIR).
-j NUM Number of jobs to use (default: $JOBS).
-f Force recomputation of all data (do not skip already computed).
-c Continue analysis using previously computed gene trees.
-R RANDOMIZE_FROM Randomize starting point (default: $RANDOMIZE_FROM).
-N NOIMPROVEMENT_STOP Number of iterations with no improvement to stop (default: $NOIMPROVEMENT_STOP).

Examples:

Process yeast dataset:
$0 -R50 -N50 -j10 -o yeast_out -d data_yeast data_yeast/gtrees_0.[0-6]*

Process simulated dataset:
$0 -R50 -N100 -j6 -o sim_out -d data_sim data_sim/wgd-[1-5]-*_0*

Output: 5 files sim_out/sim[1-5].txt + wgd.csv

Process simulated dataset (inferred):
$0 -R50 -N100 -j6 -o sim_inferred_out -d data_sim_inferred data_sim_inferred/wgd-[1-5]-*_0.[0-5]*

Output: 5 files sim_inferred_out/sim[1-5].txt + wgd.csv

If the analysis is interrupted, running the command again will resume processing from where it left off.

EOF

} 

[[ $* ]] || { Usage && exit -1 ;  }

set -- $( getopt hro:d:j:fcN:R: $* )

while [ "$1" != -- ]
do
    case $1 in
        -h)   Usage; exit;;
        -r)   export REVERSE_CLIMB="--reversed_climb";;          
		-o)   export RES=$2; shift;;          
		-d)   export DATADIR=$2; shift;;          
		-j)   export JOBS=$2; shift;; 
		-f)   export FORCE=1;;
		-c)   export CONTINUE=1;;
		-R)   export RANDOMIZE_FROM=$2; shift;;          
		-N)   export NOIMPROVEMENT_STOP=$2; shift;;          
    esac
    shift   
done

shift 


[[ $RES ]] || export RES=resultsTmp_$(basename $DATADIR)
export SPECIESTREE=$DATADIR/s_tree


mkdir -p $RES

if [[ $* ]]
then
	FILELIST="$*"
else
	FILELIST=$(  ls $DATADIR/wgd-*-gene-trees_0.*  )
fi
echo $FILELIST



function singledataset()
{
	set $1		
	GENETREES="$2"	
	echo $GENETREES start		
	outfile=$RES/$(basename $GENETREES).dat
	

	if [[ $CONTINUE ]] && [[ -f $outfile.genetree ]]
	then				
		CONT="--initial_gene_tree $outfile.genetree"
		if [[ -f $outfile ]] && grep -q "exactsolution=True" $outfile 
		then
			echo $GENETREES has the exact solution in $outfile with cost $( grep "bestcost" $outfile )
			return  
		fi
	fi


	if [[ ! $CONTINUE ]] && [[ -f $outfile ]] && ! [[ $FORCE ]]
		then
			if grep outspeciestree_worec $outfile
			then
				echo $GENETREES already completed
				return 
			fi
		fi

	#echo "python3 metaec.py --randomize_from $RANDOMIZE_FROM --noimprovement_stop $NOIMPROVEMENT_STOP --gene_trees $GENETREES --species_tree $SPECIESTREE --out_file=$outfile --distribution_maps $CONT $REVERSE_CLIMB"
	python3 metaec.py --randomize_from $RANDOMIZE_FROM --noimprovement_stop $NOIMPROVEMENT_STOP --gene_trees $GENETREES --species_tree $SPECIESTREE --out_file=$outfile --distribution_maps $CONT $REVERSE_CLIMB
	
	echo $GENETREES completed 
	
}

export -f singledataset

for i in $FILELIST; do echo $i; done | nl | parallel --ungroup --progress --no-run-if-empty -j $JOBS singledataset 

#csvmanip/csvmanip.py  -i "genetree,speciestree,outgenetree,outspeciestree,climbs,speciestreefixedwgd" $RES/*.dat > $RES/sim.$GTSET.$$.csv

[[ -d csvmanip ]] || git clone https://github.com/ppgorecki/csvmanip

csvmanip/csvmanip.py  -i "genetree,speciestree,outgenetree,outspeciestree,climbs,speciestreefixedwgd" $RES/*.dat > $RES/wgd.csv

echo $RES/wgd.csv created

if [[ "$DATADIR" =~ data_sim ]] 
then
	for i in 0 1 2 3 4 5
	do
		csvmanip/csvmanip.py  -e "outspeciestree_worec" -i "Id" -q -H $RES/wgd-$i*_0.[0-7]*.dat > $RES/sim_wgd$i.txt
		echo $RES/sim_wgd$i.txt created
	done
else

	#csvmanip/csvmanip.py  -e "outspeciestreeepicount" -i "Id,Source" -q -H $RES/*.dat > $RES/wgd.txt
	csvmanip/csvmanip.py  -e "outspeciestree_worec" -i "Id" -q -H $RES/*.dat > $RES/wgd.txt

	echo $RES/wgd.txt created

fi