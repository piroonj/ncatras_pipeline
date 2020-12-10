ref=$1
fq=$2
outBedf=$3
bname=$(basename $fq .fastq)
refbname=$(basename $ref .fa)
alnbname=$bname.$refbname

samtools fqidx $fq
minimap2 -t 35 --secondary=no -cx map-ont $ref $fq > $alnbname.paf
cut -f1,6 $alnbname.paf > $alnbname.read_map_ref.txt
sort -k1,1 -k2,2 $alnbname.read_map_ref.txt | bedtools groupby -g 1 -c 2,2 -o distinct,count_distinct > $alnbname.read_map_ref.group.txt
grep -e "CP055225,transposon" -e "CP055226,transposon" $alnbname.read_map_ref.group.txt | cut -f1 > $alnbname.read_map_ref.group.genomic_tnseq.txt

echo python extCrisprTnseq.py $fq $ref $alnbname.paf $alnbname.read_map_ref.group.genomic_tnseq.txt $outBedf
python extCrisprTnseq.py $fq $ref $alnbname.paf $alnbname.read_map_ref.group.genomic_tnseq.txt $outBedf

sort -k1,1 -k6,6 -k2,2n -k3,3n $outBedf | bedtools groupby -g 1,6,2,3 -c 1 -o count | awk 'BEGIN{OFS="\t"}{print $1,$3,$4,$2$5}' > $(basename $outBedf .bed).bdg

sort -k1,1 -k2,2 $outBedf | bedtools groupby -g 1,2,3 -c 4 -o count | sort -k4,4nr > $(basename $outBedf .bed).uniq.bdg
xs