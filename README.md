# nCATRAs_pipeline

The pipelines for counting TA sites from TnSeq (Transposon sequencing) data for both short and long reads.

---

## tnseq_illm
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

[tnseq_ont]('tnseq_ont/README.md')

---

