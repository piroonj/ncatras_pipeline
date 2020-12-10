import pandas as pd
import sys
# filename = "dCas9_pull_exp1.pass.fastq"
# faname = "old_lac_transposon.fa"
# alnPafF = "dCas9_pull_exp1.pass.old_lac_transposon.paf"
# genomic_tnseqF = "dCas9_pull_exp1.pass.read_map_ref.group.genomic_tnseq.txt"
# outName = "dCas9_pull_exp1.ont.target.bed"

filename = sys.argv[1]
faname = sys.argv[2]
alnPafF = sys.argv[3]
genomic_tnseqF = sys.argv[4]
outName = sys.argv[5]
outNameLabel = "{}.label.txt".format(outName)

readid = set()
with open(genomic_tnseqF) as inf:
    for l in inf:
        readid.add(l.strip())
readid = list(readid)

d = pd.read_csv(alnPafF, sep="\t",usecols=list(range(0,12,1)),header=None)
columns = names = ['qid','qlen','qstart','qend','strand','rid','rlen','rstart','rend',
        'match','mblock','qmap']
d.columns = columns
dFiltSort = d[d['qid'].isin(readid)].sort_values(by=['qid','qstart'])
grouped = dFiltSort.groupby('qid')


i = 1
count = {}
import pysam
import re
def getFqSeq(fq, ref, readid="", pos=[], refid="", hpos=[]):
    if len(pos) == 4:
        pos = sorted(pos)
        qseq = fq.fetch(readid)[pos[1]:pos[2]]
        hstart = int(hpos[2])
        hend = int(hpos[3])
        ovlLen = pos[2]-pos[1]
        if hpos[0]=="3end":
            if hpos[1] == "+":
                rstart = hend-ovlLen
                rend = hend
                rseq = ref.fetch(refid)[rstart:rend]
                ovl_pos = findTA(rseq, 'last')
            elif hpos[1] == "-":
                rstart = hstart
                rend = hstart+ovlLen
                rseq = ref.fetch(refid)[rstart:rend]
                ovl_pos = findTA(rseq, 'first')
        elif hpos[0]=="5end":
            if hpos[1] == "+":
                rstart = hstart
                rend = hstart+ovlLen
                rseq = ref.fetch(refid)[rstart:rend]
                ovl_pos = findTA(rseq, 'first')
            elif hpos[1] == "-":
                rstart = hend-ovlLen
                rend = hend
                rseq = ref.fetch(refid)[rstart:rend]
                ovl_pos = findTA(rseq, 'last')
        else:
            return False
        if ovl_pos:
            TAstart = rstart+ovl_pos[0]
            TAend = rstart+ovl_pos[1]
            TAverified = ref.fetch(refid)[rstart+ovl_pos[0]:rstart+ovl_pos[1]]
            return [refid, TAstart, TAend, readid, 0, hpos[1]]
    return False

def findTA(seq, get='first'):
    allTA = []
    for m in re.finditer("TA",seq):
        allTA.append([m.start(0), m.end(0)])
    if len(allTA) > 0:
        if get == "first":
            out = allTA[0]
        else:
            out = allTA[-1]
        return out
    return False

# filename = "dCas9_pull_exp1.pass.fastq"
fq = pysam.FastaFile(filename)
# faname = "old_lac_transposon.fa"
refFa = pysam.FastaFile(faname)
outf = open(outName,"w")
outLabelf = open(outNameLabel,"w")
for group_name, df_group in grouped:
#     print(group_name, len(df_group))
    if len(df_group) == 2:
        r1qid = df_group.iloc[0,].qid
        r2qid = df_group.iloc[1,].qid
        r1s = int(df_group.iloc[0,].qstart)
        r1e = int(df_group.iloc[0,].qend)
        r2s = int(df_group.iloc[1,].qstart)
        r2e = int(df_group.iloc[1,].qend)
        r1_map =  df_group.iloc[0,].rid
        r2_map =  df_group.iloc[1,].rid
        r1_strand = df_group.iloc[0,].strand 
        r2_strand = df_group.iloc[1,].strand 
        h1s = int(df_group.iloc[0,].rstart)
        h1e = int(df_group.iloc[0,].rend)
        h2s = int(df_group.iloc[1,].rstart)
        h2e = int(df_group.iloc[1,].rend)
        h1_strand = df_group.iloc[0,].strand 
        h2_strand = df_group.iloc[1,].strand 

        groupType = "unk"
        genomeJoinType = "unk"
        if r1s < r2e and r1e > r2s:
            if r2_map == "transposon" and r2_strand == "+" and h2s < 100:
                groupType = "gr1"
                genomeJoinType = "3end"
                out = getFqSeq(fq, refFa, r2qid, [r1s, r1e, r2s, r2e], r1_map,[genomeJoinType,h1_strand,h1s,h1e])
                count[groupType] = count.setdefault(groupType, 0) + 1
                outLabelf.write("{}\t{}\n".format(r2qid,groupType))
            elif r1_map == "transposon" and r1_strand == "-" and h1s < 100:
                groupType = "gr2"
                genomeJoinType = "5end"
                out = getFqSeq(fq, refFa, r2qid, [r1s, r1e, r2s, r2e], r2_map, [genomeJoinType,h2_strand,h2s,h2e])
                count[groupType] = count.setdefault(groupType, 0) + 1
                outLabelf.write("{}\t{}\n".format(r2qid,groupType))
            elif r1_map == "transposon" and r1_strand == "+" and h1e > 3100:
                groupType = "gr3"
                genomeJoinType = "5end"
                out = getFqSeq(fq, refFa, r2qid, [r1s, r1e, r2s, r2e], r2_map, [genomeJoinType,h2_strand,h2s,h2e])
                count[groupType] = count.setdefault(groupType, 0) + 1
                outLabelf.write("{}\t{}\n".format(r2qid,groupType))
            elif r2_map == "transposon" and r2_strand == "-" and h2e > 3100:
                groupType = "gr4"
                genomeJoinType = "3end"
                out = getFqSeq(fq, refFa, r2qid, [r1s, r1e, r2s, r2e], r1_map, [genomeJoinType,h1_strand,h1s,h1e])
                count[groupType] = count.setdefault(groupType, 0) + 1
                outLabelf.write("{}\t{}\n".format(r2qid,groupType))
            else:
                count[groupType] = count.setdefault(groupType, 0) + 1
                outLabelf.write("{}\t{}\n".format(r2qid,groupType))
            if out:
                outf.write("{}\n".format("\t".join(map(str,out))))
        else:
            outLabelf.write("{}\t{}\n".format(group_name,",".join(list(pd.unique(df_group.loc[:,'rid'].sort_values(ascending=True))))))
    else:
        outLabelf.write("{}\t{}\n".format(group_name,",".join(list(pd.unique(df_group.loc[:,'rid'].sort_values(ascending=True))))))

outf.close()
outLabelf.close()
print(count)
