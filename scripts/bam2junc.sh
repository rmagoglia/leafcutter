
leafCutterDir=$(dirname $0)/../
echo $leafCutterDir
bamfile=$1
bedfile=$1.bed
juncfile=$2
python $leafCutterDir/scripts/filter_cs.py $bamfile | python $leafCutterDir/scripts/sam2bed.py --use-RNA-strand - $bedfile
$leafCutterDir/scripts/bed2junc.pl $bedfile $juncfile
rm $bedfile
