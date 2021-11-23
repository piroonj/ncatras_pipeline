# nCATRAs_pipeline

The pipelines for counting TA sites from TnSeq (Transposon sequencing) data for both short and long reads.

---

## [tnseq_illm](https://github.com/piroonj/ncatras_pipeline/tree/master/tnseq_illm)
Count number of transposon insertions from Illumina reads.

#### Data pre-procession
To run the script, the input file for tnseq_identify.sh is BAM format.

#### QC trimming
fastp -i input.fq.gz -o out.fq.gz 

#### Alignment
bwa mem UAMSLAC_transposon.fasta out.fq.gz | samtools sort -o out.bam 

#### Tnseq analysis
tnseq_identify.sh out.bam 

#### Main output file
out.bam.TA_site.signal.uniq.bdg

---

## [tnseq_ont](https://github.com/piroonj/ncatras_pipeline/tree/master/tnseq_ont)

Count number of transposon insertion from Nanopore reads.

#### Data pre-procession
To run the script, the input file for nCATRAs_pipeline.sh is .

#### QC trimming
porechop -i input.fastq -o input.porechop.fastq
cat input.porechop.fastq | NanoFilt -q 7 -l 200 > input.filt.fastq

#### Tnseq analysis
bash nCATRAs_pipeline.sh UAMSLAC_transposon.fasta input.filt.fastq output.target.bed

#### Main output file
output.target.uniq.bdg

---

