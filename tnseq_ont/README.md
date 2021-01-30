
## tnseq_ont
Count number of transposon insertion from Nanopore reads.

### Data pre-procession
To run the script, the input file for nCATRAs_pipeline.sh is .

#### QC trimming
porechop -i input.fastq -o input.porechop.fastq
cat input.porechop.fastq | NanoFilt -q 7 -l 200 > input.filt.fastq

#### Tnseq analysis
bash nCATRAs_pipeline.sh UAMSLAC_transposon.fasta input.filt.fastq output.target.bed

#### Main output file
output.target.uniq.bdg

---

