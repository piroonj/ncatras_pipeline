ln -s ../0.TnSeq_ONT/0.fastq/TnSeq_ONT.filt.fastq
ln -s ../1.TnSeq_dCas9/0.fastq/TnSeq_dCas9.filt.fastq
ls ../2.TnSeq_mnCATS/0.fastq/*.fastq | xargs -I{} ln -s {}

for f in TnSeq_dCas9 TnSeq_mnCATS_2probes_r1 TnSeq_mnCATS_2probes_r2 TnSeq_mnCATS_2probes_r3 TnSeq_mnCATS_4probes TnSeq_ONT; do qsub -N $f -d $PWD -v ref=UAMSLAC_transposon.fasta,fq=$f.filt.fastq,outBedf=$f.target.bed run_crisprTnSeq.qsub; done

iat TnSeq_mnCATS_2probes_r?.filt.UAMSLAC_transposon.fasta.read_map_ref.txt | grep transposon | cut -f1 | sort -u | wc -l


