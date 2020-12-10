library(dplyr)
library(GenomicRanges)
library(GenomicAlignments)

shuffleRows <- function(data){
    set.seed(nrow(data))
    data_shuf = data[sample(1:nrow(data)),]
    return(data_shuf)
}

countUniq <- function(data, start, end){
    dataFil <- na.omit(data[seq(start,end),])
    o <- dataFil %>% 
      group_by(paste0(chrom, start)) %>% 
      summarise(uniqCount=n()) %>% 
      summarise(length(uniqCount))
    return(c(nrow(dataFil), o[[1]]))
}

processCounts <- function(data, num, sample, group){
    out <- data.frame(start=c(1),reads=c(0),uniq=c(0),sample=c(sample),group=c(group))
    for (i in seq(1, nrow(data), by=num)){
        start = i
        end = i+num-1
        count = countUniq(data, 1, end)
        out <- rbind(out, data.frame(start=c(1),
                                     reads=c(count[1]),
                                     uniq=c(count[2]),
                                     sample=c(sample),
                                     group=c(group)
                                    )
                    )
    }
    return(out)
}

countGene <- function(genes.gr, data, start, end){
    dataFil <- na.omit(data[seq(start,end),])
    dataFil.gr = GRanges(seqnames=dataFil[,1],
                        ranges=IRanges(start=dataFil[,2],end=dataFil[,3]),
                        strand=dataFil[,6])
    ovl <- summarizeOverlaps(genes.gr, dataFil.gr)
    return(c(nrow(dataFil),sum(assays(ovl)$counts > 0)))
}

processGeneCounts <- function(genes.gr, data, num, sample, group){
    out <- data.frame(start=c(1),reads=c(0),genes=c(0),sample=c(sample),group=c(group))
    for (i in seq(1, nrow(data), by=num)){
        start = i
        end = i+num-1
        count = countGene(genes.gr, data, 1, end)
        print(c(start, end, count))
        out <- rbind(out, data.frame(start=c(1),
                                     reads=c(count[1]),
                                     genes=c(count[2]),
                                     sample=c(sample),
                                     group=c(group)
                                    )
                    )
    }
    return(out)
}

genes.df = read.table("UAMSLAC.ncbi.locus_annotation.bed", header = FALSE,
                     stringsAsFactors=FALSE) 
genes.gr=GRanges(seqnames=genes.df[,1],
                 ranges=IRanges(start=genes.df[,2],end=genes.df[,3]),
                 strand=genes.df[,6])

args <- commandArgs(trailingOnly = TRUE)
file <- args[1] # input bed file
dsample <- args[2] # add sample label
dgroup <- args[3] # add group label
sname <- args[4] # output file prefix

data = read.table(file, header=F, col.names = c("chrom","start","end","name","score","strand"))

d_shuffle <- shuffleRows(data)
print(c(sname, dsample, dgroup))
counts <- processGeneCounts(genes.gr,d_shuffle, 100000, dsample, dgroup)
print(counts)
save(counts, file=paste0(sname,"_genes.RData"))

d_shuffle <- shuffleRows(data)
print(c(sname, dsample, dgroup))
counts <- processCounts(d_shuffle, 50000, dsample, dgroup)
print(counts)
save(counts, file=paste0(sname,"_counts.RData"))

