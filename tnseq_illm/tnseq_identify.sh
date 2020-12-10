bam=$1
bdg=$(basename $bam .bam).bdg
samtools index $bam
bedtools genomecov -ibam $bam -bg > $bdg

samtools view -hb -F 132 $bam | bedtools bamtobed | bedtools flank -g UAMSLAC_transposon.chrom.sizes -i stdin -l 2 -r 0 -s > $bam.Allflank.bed 
bedtools intersect -wa -u -a $bam.Allflank.bed -b UAMSLAC.TA_sties.bed -f 1.0 > $bam.TA.bed

sort -k1,1 -k6,6 -k2,2n -k3,3n $bam.TA.bed | bedtools groupby -g 1,6,2,3 -c 1 -o count | awk 'BEGIN{OFS="\t"}{print $1,$3,$4,$2$5}' > $bam.TA_site.bdg

perl -pe 's/\+|\-//g' $bam.TA_site.bdg | sort -k4,4nr > $bam.TA_site.signal.bed
for f in $bam.TA_site.signal.bed; do sort -k1,1 -k2,2 $f | bedtools groupby -g 1,2,3 -c 4 -o sum | sort -k4,4nr > $bam.TA_site.signal.uniq.bdg; done

sort -k1,1 -k6,6 -k2,2n -k3,3n $bam.Allflank.bed | bedtools groupby -g 1,6,2,3 -c 1 -o count | awk 'BEGIN{OFS="\t"}{print $1,$3,$4,$2$5}' > $bam.Allflank_site.bdg
perl -pe 's/\+|\-//g' $bam.Allflank_site.bdg | sort -k4,4nr > $bam.Allflank_site.signal.bed
for f in $bam.Allflank_site.signal.bed; do sort -k1,1 -k2,2 $f | bedtools groupby -g 1,2,3 -c 4 -o sum | sort -k4,4nr > $bam.Allflank_site.signal.uniq.bdg; done


