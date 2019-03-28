# Creates the annotation database that links genes to genomic regions
# This annotation will be saved as a .rda file for the package
# Sources of annotation so far:
# gene exons and introns - from ENSEMBL
# 3' and 5' UTR - from ENSEMBL
# promoters - find overlaps between promoter regions from ENCODE - and regions of 100kb, 50kb, 25kb, 10kb, and 1kb upstream of 5' UTR of each gene

library(GenomicRanges)
library(data.table)
library(pbapply)

# load genes
gene.annotations <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\EWAS\\gene_annotations\\ens_grch37_genes_trans_exons_20-Mar-2017.txt",stringsAsFactors = F)
names(gene.annotations) <- c("Ensembl.Gene.ID","Ensembl.Transcript.ID","Associated.Gene.Name","Chromosome.Name","Strand","Gene.Start..bp.","Gene.End..bp.",
                             "Transcript.Start..bp.","Transcript.End..bp.","Transcription.Start.Site..TSS.","Ensembl.Exon.ID","Exon.Chr.Start..bp.",
                             "Exon.Chr.End..bp.","Exon.Rank.in.Transcript","trans_gene" )
gene.annotations$Strand <- ifelse(gene.annotations$Strand<0,"-","+")
gene.annotations <- gene.annotations[,c("Ensembl.Gene.ID","Chromosome.Name","Gene.Start..bp.","Gene.End..bp.","Strand","Transcription.Start.Site..TSS.")]
gene.annotations <- unique(gene.annotations)

gene.annotation.gr <- makeGRangesFromDataFrame(gene.annotations,
                                               keep.extra.columns = T,strand.field = "Strand",seqnames.field = "Chromosome.Name",
                                               start.field = "Gene.Start..bp.",end.field="Gene.End..bp.")

gene.annotation.gr$Feature <- "Gene_Body"
names(mcols(gene.annotation.gr))[1] <- "genes"

# load promoter
encode.list.p <- list.files("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\EWAS\\Analysis_Peter\\Promoter_annotation_data",full.names = T)

encode.promoters.list <- lapply(encode.list.p,read.delim,header=F)
encode.promoters <- as.data.frame(rbindlist(encode.promoters.list))
encode.promoters$V1 <- sapply(encode.promoters$V1,substring,first=4)
encode.promoters.distal.gr <- makeGRangesFromDataFrame(encode.promoters[grepl("Distal-Prediction",encode.promoters$V4),],
                                                       keep.extra.columns = T,ignore.strand = T,seqnames.field = "V1",
                                                       start.field = "V2",end.field="V3")
encode.promoters.distal.gr <- reduce(encode.promoters.distal.gr)



encode.promoters.prox.gr <- makeGRangesFromDataFrame(encode.promoters[grepl("Proximal-Prediction",encode.promoters$V4),],
                                                     keep.extra.columns = T,ignore.strand = T,seqnames.field = "V1",
                                                     start.field = "V2",end.field="V3")
encode.promoters.prox.gr <- reduce(encode.promoters.prox.gr)

# Add 5' and 3' UTR annotations:

utr <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\EWAS\\Analysis_Peter\\Ensembl_UTR\\ensembl_utrs.txt",stringsAsFactors = F)

utr <- data.frame(Ensembl.Gene.ID=c(utr$Gene.stable.ID[!is.na(utr$X5..UTR.start)],utr$Gene.stable.ID[!is.na(utr$X3..UTR.start)]),
                  CHR=c(utr$Chromosome.scaffold.name[!is.na(utr$X5..UTR.start)],utr$Chromosome.scaffold.name[!is.na(utr$X3..UTR.start)]),
                  Start=c(utr$X5..UTR.start[!is.na(utr$X5..UTR.start)],utr$X3..UTR.start[!is.na(utr$X3..UTR.start)]),
                  End=c(utr$X5..UTR.end[!is.na(utr$X5..UTR.start)],utr$X3..UTR.end[!is.na(utr$X3..UTR.start)]))


# tss.regions.gr <- GRanges(seqnames=seqnames(gene.annotation.gr),
#                           ranges=IRanges(start=ifelse(strand(gene.annotation.gr)=="-",gene.annotation.gr$Transcription.Start.Site..TSS.,gene.annotation.gr$Transcription.Start.Site..TSS.-100000),
#                                          end=ifelse(strand(gene.annotation.gr)=="-",gene.annotation.gr$Transcription.Start.Site..TSS.+100000,gene.annotation.gr$Transcription.Start.Site..TSS.)),
#                           strand=strand(gene.annotation.gr),mcols(complexid.annotation.gr))
#
# promoter.distal.tss.hits <- findOverlaps(encode.promoters.distal.gr,tss.regions.gr)
# split.hits.idx <- split(1:length(promoter.distal.tss.hits),queryHits(promoter.distal.tss.hits),drop=T)
# encode.promoters.distal.gr$genes <- ""
# encode.promoters.distal.gr$genes[as.integer(names(split.hits.idx))] <- pbsapply(split.hits.idx, function(x) {
#   return(tss.regions.gr$Ensembl.Gene.ID[subjectHits(promoter.distal.tss.hits)[x]])
# })


# create sample snp hits nad methy hits
# 50 hits in genes, 50 hits in distal promoters, 50 hits in prox promoters, 10 hits outside those regions
hits <- sample(gene.annotation.gr,50)
mcols(hits) <- NULL
hits <- c(hits,sample(encode.promoters.distal.gr,50),sample(encode.promoters.prox.gr,50))
start(hits) <- end(hits) <- sapply(hits,function(x) sample(start(x):end(x),1))

gr <- GRanges(
  seqnames = sample(1:22,size = 10,replace = T),
  ranges = IRanges(101:110, end = 101:110),
  strand = "*")


hits <- c(hits,gr)

hits$name <- paste("Test",1:length(hits),sep="_")
hits$pheno <- sample(LETTERS[1:8],length(hits),replace=T)
save(hits,file = "D:\\projects\\ComplexID\\exampleHits.Rdata")
hits.pheno <- matrix(c(LETTERS[1:8],seq(from=0.1,to=.8,length.out=8)),ncol=2)

use_data(hits,hits.pheno)
# load PGRS data to test

pgrs.hits <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\data\\results_PGRS_SNPs_pleio_combined_4outcomes_n656_ANNOTATIONS.txt",stringsAsFactors = F)
pgrs.hits <- pgrs.hits[pgrs.hits$index.beta!="0",]
sig.hits <- apply(pgrs.hits,1,function(x) {
  #print(x)
  x <- as.data.frame(t(as.data.frame(x,stringsAsFactors=F)),stringsAsFactors=F)
  #return(x)
  ret <- x[rep(1,length(strsplit(x$index.beta,",")[[1]])),]
  ret$index.beta <- strsplit(x$index.beta,",")[[1]]
  return(ret)
})
library(data.table)
sig.hits <- as.data.frame(rbindlist(sig.hits))
sig.hits$BP <- as.integer(sig.hits$BP)
sig.hits$CHR <- as.integer(sig.hits$CHR)

sig.hits.gr <- makeGRangesFromDataFrame(sig.hits[,c("SNP","index.beta","CHR","BP")],keep.extra.columns = T,ignore.strand = T,seqnames.field = "CHR",
                                        start.field = "BP",end.field="BP")

save(sig.hits.gr,file = "D:\\projects\\ComplexID\\PGRS_Sig_Hits.Rdata")
system.time(test <- runSCICW(sig.hits.gr,sig.hits.gr,phenoSim=1,promoterRange = 10000))

# load STRING PPI
string.ppi <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\string_ppi\\9606.protein.links.v10.5.txt.gz",sep=" ",stringsAsFactors = F)
string.ppi$protein1 <- substring(string.ppi$protein1,first = 6)
string.ppi$protein2 <- substring(string.ppi$protein2,first = 6)

# filter STRING PPI for interactions
string.ppi.filt <- string.ppi[string.ppi$combined_score>700,]

# load ENSEMBL and Entrez Gene Names
gene.names <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\ensembl_gene_to_protein\\Gene_and_protein_names.gz",stringsAsFactors = F)

# annotate the string ppi with uniprot and entrez gene ids in order to combine with CORUM
string.ppi.annot <- merge(string.ppi.filt,gene.names,by.x="protein1",by.y="Protein.stable.ID",all.x=T)
string.ppi.annot <- merge(string.ppi.annot,gene.names,by.x="protein2",by.y="Protein.stable.ID",all.x=T)
string.ppi.annot <- string.ppi.annot[!is.na(string.ppi.annot$EntrezGene.ID.x) & !is.na(string.ppi.annot$EntrezGene.ID.y),]

# load CORUM
corum.annot <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\corum\\coreComplexes.txt",stringsAsFactors = F)
corum.annot <- corum.annot[corum.annot$Organism=="Human",]

# convert CORUM to list connected each complex to all proteins in it
corum.subunits <- strsplit(corum.annot$subunits.Entrez.IDs.,";")
corum.subunits <- lapply(corum.subunits,function(x) x[which(x != "None")])
names(corum.subunits) <- corum.annot$ComplexID
corum.subunits <- corum.subunits[sapply(corum.subunits,length)>1]

library(ggplot2)
ggplot(mapping=aes(sapply(corum.subunits,length)))+geom_histogram()+labs(x="Size of Complex",y="Count of Complexes",title="Distribution of the Sizes of Complexes in CORUM Core Complexes")
corum.subunits[sapply(corum.subunits,length)>70]


# remove duplicate complexes by alphabetizing the subunit lists, and finding the duplicates
corum.subunits <- lapply(corum.subunits, function(x) {
  return(sort(as.integer(x),decreasing = F))
  #paste(sort(as.integer(x),decreasing = F),collapse=";")
})
corum.subunits <- corum.subunits[!duplicated(corum.subunits)]

# remove very large complexes that are outliers (>70 subunits)
corum.subunits <- corum.subunits[sapply(corum.subunits,length)<=70]

# examine fully-overlapping complexes
find.large.complex <- sapply(corum.subunits, function(x) {
  sapply(corum.subunits, function(y) {
    if (length(x) >= length(y))
      if (sum(y %in% x) == length(y))
        return(T)
      return(F)
    if (sum(x %in% y) == length(x))
      return(T)
    return(F)
  })
})
num.smaller.complexes<-apply(find.large.complex,2,sum)
num.larger.complexes<-apply(find.large.complex,1,sum)

ggplot(mapping=aes(num.smaller.complexes-1))+geom_bar()+scale_x_continuous("Number of CORUM sub-complexes in complex", 0:9, 0:9)+labs(title="Distribution of the Number of Sub-Complexes in a Complex")
sum(num.smaller.complexes>1)

ggplot(mapping=aes(num.larger.complexes-1))+geom_bar()+scale_x_continuous("Number of CORUM complexes in another Complex", 0:25, 0:25)+labs(title="Distribution of the Number of Complexes that are part of another Complex")
sum(num.larger.complexes>1)
# remove all complexes that have larger complex that completely incorporate them
corum.subunits <- corum.subunits[num.larger.complexes==1]

# convert STRING's ENSP-ENSP interactions to Entrez-Entrez gene interactions and remove any duplicates
string.ppi.entrez <- string.ppi.annot[,c("EntrezGene.ID.x","EntrezGene.ID.y")]
string.ppi.entrez <- unique(string.ppi.entrez)
string.ppi.entrez <- rbind(string.ppi.entrez,data.frame(EntrezGene.ID.x=string.ppi.entrez$EntrezGene.ID.y,
                                                        EntrezGene.ID.y=string.ppi.entrez$EntrezGene.ID.x))
string.ppi.entrez <- unique(string.ppi.entrez)

# remove interactions between the same gene and itself
string.ppi.entrez <- string.ppi.entrez[string.ppi.entrez$EntrezGene.ID.x!=string.ppi.entrez$EntrezGene.ID.y,]

# create dictionary, linking each gene to all of its neighbors.
geneToGene <- split(string.ppi.entrez$EntrezGene.ID.y,string.ppi.entrez$EntrezGene.ID.x,drop=T)

# check if there are any complexes that do not have a single protein in the PPI
sum(sapply(corum.subunits, function(x) { sum(x %in% as.integer(names(geneToGene))) })==0)

# create network of complexes
C2C.mat <- matrix(0,nrow = length(corum.subunits),ncol=length(corum.subunits))
rownames(C2C.mat) <- colnames(C2C.mat) <- names(corum.subunits)

for (i in 1:nrow(C2C.mat)) {
  row.subunits <- corum.subunits[[i]]
  for (j in 1:i) {
    if (i==j)
      next
    col.subunits <- corum.subunits[[j]]
    intersection <- row.subunits %in% col.subunits
    col.intersection <- col.subunits %in% row.subunits
    denom <- (length(row.subunits) - sum(intersection))*(length(col.subunits) - sum(col.intersection))
    col.subunits <- col.subunits[!col.intersection]
    row.subunits.filt <- row.subunits[!intersection]

    row.interactions <- unlist(geneToGene[as.character(row.subunits.filt)])
    #if(sum(as.integer(col.subunits) %in% unlist(row.interactions)))
    C2C.mat[i,j] <- sum(sapply(col.subunits,function(x) length(which(x==row.interactions))))/denom
    #print(j)
  }
  if (i %% 100 == 0)
    print(i)
}
C2C.mat[upper.tri(C2C.mat)] <- t(C2C.mat)[upper.tri(C2C.mat)]

# create list of individual protein complexes (proteins that are in the STRING DB but not in CORUM)
indivToGene <- geneToGene[!(as.integer(names(geneToGene)) %in% unlist(corum.subunits))]

# create I2C and C2I matrices to connect proteins not in CORUM to complexes in CORUM
I2C.mat <- matrix(0,nrow = length(indivToGene),ncol=length(corum.subunits))
rownames(I2C.mat) <- names(indivToGene)
colnames(I2C.mat) <- names(corum.subunits)

indivGeneDegrees <- sapply(indivToGene,length)
cAproteins <- lapply(corum.subunits,function(x) unlist(geneToGene[as.character(x)]) )

for (i in 1:nrow(I2C.mat)) {
  geneI <- as.integer(names(indivToGene[i]))
  for (j in 1:ncol(I2C.mat)) {
    cA <- cAproteins[[j]]
    I2C.mat[i,j] <- sum(geneI == cA)
  }
  if (i %% 100 == 0)
    print(i)
}

C2I.mat <- t(I2C.mat)

I2C.mat <- I2C.mat/indivGeneDegrees
C2I.mat <- C2I.mat/sapply(corum.subunits,length)

# Create I2I matrix to connect proteins not in CORUM with other proteins not in CORUM
I2I.mat <- matrix(0,nrow = length(indivToGene),ncol=length(indivToGene))
rownames(I2I.mat) <- names(indivToGene)
colnames(I2I.mat) <- names(indivToGene)

for (i in 1:nrow(I2I.mat)) {
  I2I.mat[i,which(as.integer(names(indivToGene)) %in% indivToGene[[i]])] <- 1
  if (i %% 100 == 0)
    print(i)
}

I2I.mat <- I2I.mat/indivGeneDegrees

W <- rbind(cbind(C2C.mat,C2I.mat),cbind(I2C.mat,I2I.mat))

W_norm_t <- t(W)
W_norm_t <- sweep(W_norm_t,2,colSums(W_norm_t),'/')
sum(is.na(W_norm_t))
W_norm_t_spars <- Matrix(W_norm_t, sparse = TRUE)


# Calculate density(Ca) for all complexes
density <- sapply(corum.subunits, function(x) {
  all.interactions <- unlist(geneToGene[as.character(x)])
  edges <- sum(sapply(x,function(y) sum(y == all.interactions)))/2
  return((2*edges)/(length(x)*(length(x)-1)))
})

ggplot(mapping=aes(density))+geom_histogram()+labs(x="Density of Complex",title="Distribution of Density of Multi-Protein Complexes")
table(sapply(corum.subunits,length)[density==1])
table(sapply(corum.subunits,length)[density==0])
min(density[density!=0])

density <- density+min(density[density!=0])
density <- ifelse(density>1,1,density)
ggplot(mapping=aes(density))+geom_histogram()+labs(x="Density of Complex",title="Distribution of Density of Multi-Protein Complexes\nAfter Adjustment")



# Need to change gene annotations to only include genes in the PPI or CORUM AND have an Entrez Gene ID

all.entrez.genes.in.network <- unique(c(string.ppi.entrez$EntrezGene.ID.x,string.ppi.entrez$EntrezGene.ID.y,unlist(corum.subunits)))

gene.annotations.entrez <- merge(gene.annotations,gene.names[,c("Gene.stable.ID","EntrezGene.ID")],by=1)
gene.annotations.entrez <- unique(gene.annotations.entrez)
length(unique(gene.annotations.entrez$Ensembl.Gene.ID[is.na(gene.annotations.entrez$EntrezGene.ID)]))
gene.annotations.entrez <- gene.annotations.entrez[!is.na(gene.annotations.entrez$EntrezGene.ID),]
length(unique(gene.annotations.entrez$EntrezGene.ID[!(gene.annotations.entrez$EntrezGene.ID %in% all.entrez.genes.in.network)]))
gene.annotations.entrez <- gene.annotations.entrez[gene.annotations.entrez$EntrezGene.ID %in% all.entrez.genes.in.network,]
gene.annotations.entrez <- unique(gene.annotations.entrez)
gene.annotations.entrez$Ensembl.Gene.ID <- gene.annotations.entrez$EntrezGene.ID
gene.annotations.entrez$EntrezGene.ID <- NULL
names(gene.annotations.entrez)[1] <- "EntrezGene.ID"

gene.annotation.gr <- makeGRangesFromDataFrame(gene.annotations.entrez,
                                               keep.extra.columns = T,strand.field = "Strand",seqnames.field = "Chromosome.Name",
                                               start.field = "Gene.Start..bp.",end.field="Gene.End..bp.")

gene.annotation.gr$Feature <- "Gene_Body"
names(mcols(gene.annotation.gr))[1] <- "genes"

length(unique(gene.annotation.gr$genes))
sum(width(reduce(gene.annotation.gr)))

# Add 5', 3' UTR regions

utr.entrez <- merge(utr,gene.names[,c("Gene.stable.ID","EntrezGene.ID")],by=1)
utr.entrez <- unique(utr.entrez)
length(unique(utr.entrez$Ensembl.Gene.ID[is.na(utr.entrez$EntrezGene.ID)]))
utr.entrez <- utr.entrez[!is.na(utr.entrez$EntrezGene.ID),]
length(unique(utr.entrez$EntrezGene.ID[!(utr.entrez$EntrezGene.ID %in% all.entrez.genes.in.network)]))
utr.entrez <- utr.entrez[utr.entrez$EntrezGene.ID %in% all.entrez.genes.in.network,]
utr.entrez <- unique(utr.entrez)
utr.entrez$Ensembl.Gene.ID <- utr.entrez$EntrezGene.ID
utr.entrez$EntrezGene.ID <- NULL
names(utr.entrez)[1] <- "EntrezGene.ID"

utr.entrez.gr <- makeGRangesFromDataFrame(utr.entrez,
                                               keep.extra.columns = T,ignore.strand = T,seqnames.field = "CHR",
                                               start.field = "Start",end.field="End")

utr.entrez.gr$Feature <- "UTR"
names(mcols(utr.entrez.gr))[1] <- "genes"

# Create gene-to-complex mapping for the finally scoring (equation 13)

geneToComplex <- lapply(all.entrez.genes.in.network, function(x) {
  corum <- which(sapply(corum.subunits, function(y) x %in% y))
  indiv.names <- rownames(W)[(length(corum.subunits)+1):nrow(W)]
  return(c(corum,which(as.character(x) == indiv.names)+length(corum.subunits)))
})
names(geneToComplex) <- all.entrez.genes.in.network

library(devtools)
.gene.annotation.gr <- gene.annotation.gr
rm(gene.annotation.gr)
.encode.promoters.distal.gr <- encode.promoters.distal.gr
rm(encode.promoters.distal.gr)
.encode.promoters.prox.gr <- encode.promoters.prox.gr
rm(encode.promoters.prox.gr)
.density <- density
rm(density)
.corum.subunits <- corum.subunits
rm(corum.subunits)
.W_norm_t_spars <- W_norm_t_spars
rm(W_norm_t_spars)
.geneToComplex <- geneToComplex
rm(geneToComplex)
.utr.entrez.gr <- utr.entrez.gr
rm(utr.entrez.gr)
use_data(.gene.annotation.gr,.encode.promoters.distal.gr,.encode.promoters.prox.gr,.density,.corum.subunits,.W_norm_t_spars,.geneToComplex,.utr.entrez.gr,internal=T,overwrite=T)

system.time(test <- runSCICW(sig.hits.gr,phenoSim=1,promoterRange = 10000))

# generate phenotype similarity matrix
phenosim.mat <- matrix(c(1,2,3,4,1,1,1,1),ncol=2)

system.time(test <- runComplexID(sig.hits.gr,phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T))

# add hugo gene name and corum complex names
entrezToHugo <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\entrez_to_hugo\\entrez_to_hugo.gz")
entrezToHugo <- entrezToHugo[!is.na(entrezToHugo$EntrezGene.ID),]
hugoNames <- entrezToHugo$HGNC.symbol[match(names(ComplexID:::.geneToComplex),as.character(entrezToHugo$EntrezGene.ID))]

complexNames <- sapply(ComplexID:::.geneToComplex, function(x) {
  if (x[1] <= length(ComplexID:::.corum.subunits))
    return(paste(unique(corum.annot$ComplexName[match(names(x),corum.annot$ComplexID)]),collapse=";"))
  return("None")
})

.hugoNames <- hugoNames
.complexNames <- complexNames

.gene.annotation.gr <- ComplexID:::.gene.annotation.gr
.encode.promoters.distal.gr <- ComplexID:::.encode.promoters.distal.gr
.encode.promoters.prox.gr <- ComplexID:::.encode.promoters.prox.gr
.density <- ComplexID:::.density
.corum.subunits <- ComplexID:::.corum.subunits
.W_norm_t_spars <- ComplexID:::.W_norm_t_spars
.geneToComplex <- ComplexID:::.geneToComplex
.utr.entrez.gr <- ComplexID:::.utr.entrez.gr

use_data(.gene.annotation.gr,.encode.promoters.distal.gr,.encode.promoters.prox.gr,.density,.corum.subunits,.W_norm_t_spars,.geneToComplex,.utr.entrez.gr,.hugoNames,.complexNames,internal=T,overwrite=T)

# Construct igraph object of PPI, add edges to it of genes not in PPI but in CORUM
library(igraph)
complex.size <- sapply(ComplexID:::.corum.subunits,length)
colA <- unlist(mapply(rep,as.integer(names(ComplexID:::.corum.subunits)),complex.size))
corum.edge.list <- matrix(data = c(colA,unlist(ComplexID:::.corum.subunits)),ncol = 2)
corum.edge.list <- lapply(ComplexID:::.corum.subunits,function(x) as.data.frame(t(combn(x,2))))
library(data.table)
corum.edge.list <- as.data.frame(rbindlist(corum.edge.list))

string.ppi.entrez$color <- ifelse(paste(string.ppi.entrez$EntrezGene.ID.x,string.ppi.entrez$EntrezGene.ID.y,sep=";") %in% c(paste(corum.edge.list[,1],corum.edge.list[,2],sep=";"),paste(corum.edge.list[,2],corum.edge.list[,1],sep=";")),"red","black")
corum.edge.list <- cbind(corum.edge.list,
                         matrix(ifelse(paste(corum.edge.list[,1],corum.edge.list[,2],sep=";") %in% c(paste(string.ppi.entrez$EntrezGene.ID.x,string.ppi.entrez$EntrezGene.ID.y,sep=";"),paste(string.ppi.entrez$EntrezGene.ID.y,string.ppi.entrez$EntrezGene.ID.x,sep=";")),"red","green"),ncol = 1))

colnames(corum.edge.list) <- names(string.ppi.entrez)
complete.edge.list <- unique(rbind(string.ppi.entrez,corum.edge.list))
complete.edge.list$EntrezGene.ID.x <- as.character(complete.edge.list$EntrezGene.ID.x)
complete.edge.list$EntrezGene.ID.y <- as.character(complete.edge.list$EntrezGene.ID.y)

ordered <- ifelse(as.integer(complete.edge.list$EntrezGene.ID.x) < as.integer(complete.edge.list$EntrezGene.ID.y),paste(complete.edge.list$EntrezGene.ID.x,complete.edge.list$EntrezGene.ID.y,sep=";"),paste(complete.edge.list$EntrezGene.ID.y,complete.edge.list$EntrezGene.ID.x,sep=";"))
complete.edge.list <- complete.edge.list[!duplicated(ordered),]

#complete.edge.list$color <-

complete.igraph <- graph_from_edgelist(as.matrix(complete.edge.list[,1:2]),directed=F)
E(complete.igraph)$color <- complete.edge.list$color
#geneToGraphVertex <- sapply(all.entrez.genes.in.network, function(x) { which(names(V(complete.igraph)) == as.character(x))  })
#names(geneToGraphVertex) <- all.entrez.genes.in.network

V(complete.igraph)$color <- ifelse(names(V(complete.igraph)) %in% as.character(test$scores$Entrez.Gene.ID)[1:10],"red","orange")
plot.igraph(x = subgraph(complete.igraph,unlist(neighborhood(complete.igraph,order = 0,nodes = as.character(test$scores$Entrez.Gene.ID)[1:10]))))
plot.igraph(x = subgraph(complete.igraph,unlist(neighborhood(complete.igraph,order = 1,nodes = as.character(test$scores$Entrez.Gene.ID)[1:10]))))
plot.igraph(x = subgraph(complete.igraph,unlist(neighborhood(complete.igraph,order = 1,nodes = as.character(test$scores$Entrez.Gene.ID)[1:5]))))
plot.igraph(x = subgraph(complete.igraph,unlist(neighborhood(complete.igraph,order = 1,nodes = as.character(test$scores$Entrez.Gene.ID)[1:3]))))



# try to include complexes underlying genes:
graph.to.plot <- subgraph(complete.igraph,unlist(neighborhood(complete.igraph,order = 2,nodes = as.character(test$scores$Entrez.Gene.ID)[1])))
complexes.to.plot <- ComplexID:::.corum.subunits[sapply(ComplexID:::.corum.subunits,function(x) { sum(as.integer(names(V(graph.to.plot))) %in% x)>0 })]
complexes.to.plot <- lapply(complexes.to.plot,function(x) { as.character(x)[as.character(x) %in% names(V(graph.to.plot))] })

set.seed(1)
plot.igraph(x = subgraph(complete.igraph,unlist(neighborhood(complete.igraph,order = 2,nodes = as.character(test$scores$Entrez.Gene.ID)[1]))))
set.seed(1)
plot.igraph(x=graph.to.plot,mark.groups = complexes.to.plot)
set.seed(1)

# try smaller vertices and remove labels:
plot.igraph(x=graph.to.plot,mark.groups = complexes.to.plot,vertex.size=4,vertex.label=NA)


# try to separate the nodes out:
e <- get.edgelist(graph.to.plot,names = F)
layout <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(graph.to.plot))
plot(graph.to.plot,layout=layout,mark.groups = complexes.to.plot,vertex.size=4,vertex.label=NA)

# plot larger graph:
graph.to.plot <- subgraph(complete.igraph,unlist(neighborhood(complete.igraph,order = 1,nodes = as.character(test$scores$Entrez.Gene.ID)[1:10])))
complexes.to.plot <- ComplexID:::.corum.subunits[sapply(ComplexID:::.corum.subunits,function(x) { sum(as.integer(names(V(graph.to.plot))) %in% x)>0 })]
complexes.to.plot <- lapply(complexes.to.plot,function(x) { as.character(x)[as.character(x) %in% names(V(graph.to.plot))] })
set.seed(1)
plot.igraph(x=graph.to.plot,mark.groups = complexes.to.plot)
set.seed(1)
plot.igraph(x=graph.to.plot,mark.groups = complexes.to.plot,vertex.size=4,vertex.label=NA)
e <- get.edgelist(graph.to.plot,names = F)
layout <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(graph.to.plot))
plot(graph.to.plot,layout=layout,mark.groups = complexes.to.plot,vertex.size=4,vertex.label=NA)

set.seed(1)
plot.igraph(x=graph.to.plot,mark.groups = complexes.to.plot,vertex.size=4,vertex.label.cex=.75)
ent.to.hug <- ComplexID:::.hugoNames
names(ent.to.hug) <- names(ComplexID:::.geneToComplex)
vertex.interest <- c("5802","5789","346007","23213")
graph.to.plot <- subgraph(complete.igraph,unlist(neighborhood(complete.igraph,order = 1,nodes = vertex.interest)))
complexes.to.plot <- ComplexID:::.corum.subunits[sapply(ComplexID:::.corum.subunits,function(x) { sum(as.integer(names(V(graph.to.plot))) %in% x)>0 })]
complexes.to.plot <- lapply(complexes.to.plot,function(x) { as.character(x)[as.character(x) %in% names(V(graph.to.plot))] })
#set.seed(1)
plot.igraph(x=graph.to.plot,mark.groups = complexes.to.plot,vertex.size=4,vertex.label=ent.to.hug[names(V(graph.to.plot))],vertex.label.degree=pi/2,vertex.label.dist=1)



# plot top 50, no neighbors
graph.to.plot <- subgraph(complete.igraph,unlist(neighborhood(complete.igraph,order = 0,nodes = as.character(test$scores$Entrez.Gene.ID)[1:50])))
complexes.to.plot <- ComplexID:::.corum.subunits[sapply(ComplexID:::.corum.subunits,function(x) { sum(as.integer(names(V(graph.to.plot))) %in% x)>0 })]
complexes.to.plot <- lapply(complexes.to.plot,function(x) { as.character(x)[as.character(x) %in% names(V(graph.to.plot))] })
set.seed(1)
plot.igraph(x=graph.to.plot,mark.groups = complexes.to.plot,vertex.size=4,vertex.label=NA)

# plot top 100, no neighbors
graph.to.plot <- subgraph(complete.igraph,unlist(neighborhood(complete.igraph,order = 0,nodes = as.character(test$scores$Entrez.Gene.ID)[1:100])))
complexes.to.plot <- ComplexID:::.corum.subunits[sapply(ComplexID:::.corum.subunits,function(x) { sum(as.integer(names(V(graph.to.plot))) %in% x)>0 })]
complexes.to.plot <- lapply(complexes.to.plot,function(x) { as.character(x)[as.character(x) %in% names(V(graph.to.plot))] })
set.seed(1)
plot.igraph(x=graph.to.plot,mark.groups = complexes.to.plot,vertex.size=4,vertex.label=NA)

#plot top 100, no neighbors, color by presence of hits
graph.to.plot <- subgraph(complete.igraph,unlist(neighborhood(complete.igraph,order = 0,nodes = as.character(test$scores$Entrez.Gene.ID)[1:100])))
V(graph.to.plot)$color <- ifelse(names(V(graph.to.plot)) %in% test$scores$Entrez.Gene.ID[test$scores$Num.Hits>0],"red","orange")
complexes.to.plot <- ComplexID:::.corum.subunits[sapply(ComplexID:::.corum.subunits,function(x) { sum(as.integer(names(V(graph.to.plot))) %in% x)>0 })]
complexes.to.plot <- lapply(complexes.to.plot,function(x) { as.character(x)[as.character(x) %in% names(V(graph.to.plot))] })
set.seed(1)
plot.igraph(x=graph.to.plot,mark.groups = complexes.to.plot,vertex.size=4,vertex.label=NA)



# output use_data with complex network igraph:
.complete.igraph <- complete.igraph
.ent.to.hug <- ent.to.hug


.gene.annotation.gr <- ComplexID:::.gene.annotation.gr
.encode.promoters.distal.gr <- ComplexID:::.encode.promoters.distal.gr
.encode.promoters.prox.gr <- ComplexID:::.encode.promoters.prox.gr
.density <- ComplexID:::.density
.corum.subunits <- ComplexID:::.corum.subunits
.W_norm_t_spars <- ComplexID:::.W_norm_t_spars
.geneToComplex <- ComplexID:::.geneToComplex
.utr.entrez.gr <- ComplexID:::.utr.entrez.gr
.hugoNames <- ComplexID:::.hugoNames
.complexNames <- ComplexID:::.complexNames

use_data(.gene.annotation.gr,.encode.promoters.distal.gr,.encode.promoters.prox.gr,.density,.corum.subunits,.W_norm_t_spars,.geneToComplex,.utr.entrez.gr,.hugoNames,.complexNames,.complete.igraph,.ent.to.hug,internal=T,overwrite=T)


# test plotting functionality
a <- generatePlot(test$scores$Entrez.Gene.ID[1:10])
a <- generatePlot(test$scores$Entrez.Gene.ID[1:10],useHugoNames = F,vertex.size=4,vertex.label=NA)
a <- generatePlot(test$scores$Entrez.Gene.ID[1:10],useHugoNames = F,vertex.size=4,vertex.label=NA,order=1)
generatePlot(test$scores$Entrez.Gene.ID[1:50],useHugoNames = F,vertex.size=4,vertex.label=NA)
generatePlot(test$scores$Entrez.Gene.ID[1:100],useHugoNames = F,vertex.size=4,vertex.label=NA)
generatePlot(test$scores$Entrez.Gene.ID[1:50],useHugoNames = T,vertex.size=4,vertex.label.degree=pi/2,vertex.label.dist=.75,vertex.label.cex=.5)
generatePlot(test$scores$Entrez.Gene.ID[1:100],useHugoNames = T,vertex.size=4,vertex.label.degree=pi/2,vertex.label.dist=.75,vertex.label.cex=.5)

# create database of eqtls:
eqtl.files <- list.files("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\gtex\\GTEx_Analysis_v6p_eQTL\\",full.names=T)
eqtl.files <- eqtl.files[grepl("signif_snpgene_pairs",eqtl.files)]

genenamesmapping <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\ensembl_gene_to_protein\\Gene_names_8_18_17.gz")
genenamesmapping$EntrezGene.ID <- as.character(genenamesmapping$EntrezGene.ID)
genenamesmapping <- genenamesmapping[,c("Gene.stable.ID","EntrezGene.ID")]

eqtl.list <- lapply(eqtl.files, function(x) {
  print(x)
  infile <- read.delim(x,stringsAsFactors = F)
  s <- strsplit(infile[,1],"_")
  infile$snp_chr <- sapply(s,"[[",1)
  infile$snp_pos <- sapply(s,"[[",2)
  return(infile[,c("gene_id","snp_chr","snp_pos")])
})
library(data.table)
eqtl.df <- as.data.frame(rbindlist(eqtl.list))
eqtl.df <- unique(eqtl.df)
eqtl.df$gene_id <- substr(eqtl.df$gene_id,1,15)
eqtl.df <- unique(eqtl.df)

genenamesmapping <- unique(genenamesmapping)

eqtl.df <- merge(eqtl.df,genenamesmapping,by=1)

eqtl.df <- eqtl.df[eqtl.df$EntrezGene.ID %in% as.character(all.entrez.genes.in.network),]
eqtl.df <- eqtl.df[,2:4]
eqtl.df <- unique(eqtl.df)
eqtl.gr <- makeGRangesFromDataFrame(eqtl.df,
                                    keep.extra.columns = T,ignore.strand = T,seqnames.field = "snp_chr",
                                    start.field = "snp_pos",end.field="snp_pos")
eqtl.gr$Feature="eQTL"

# annotate loops
encode.list <- list.files("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\enhancers\\ENCODE_Predictions\\",full.names = T)

encode.enhancers.list <- lapply(encode.list,read.delim,header=F)
library(data.table)
encode.enhancers <- as.data.frame(rbindlist(encode.enhancers.list))
encode.enhancers$V1 <- sapply(encode.enhancers$V1,substring,first=4)
encode.enhancers.gr <- makeGRangesFromDataFrame(encode.enhancers,
                                                       keep.extra.columns = F,ignore.strand = T,seqnames.field = "V1",
                                                       start.field = "V2",end.field="V3")
encode.enhancers.gr <- reduce(encode.enhancers.gr)

vista.fantom.enhancers <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\EWAS\\Analysis_Peter\\Enhancer_annotation_data\\ensembl_hg19_other_regulatory_regions_vista_fantom_enhancers.txt.gz",stringsAsFactors = F)
vista.gr <- makeGRangesFromDataFrame(vista.fantom.enhancers[vista.fantom.enhancers$Feature.Type=="VISTA Enhancers",],
                                     keep.extra.columns = F,ignore.strand = T,seqnames.field = "Chromosome.Name",
                                     start.field = "Start..bp.",end.field="End..bp.")
fantom.robust.gr <- makeGRangesFromDataFrame(vista.fantom.enhancers[vista.fantom.enhancers$Feature.Type=="FANTOM predictions"
                                                                    & vista.fantom.enhancers$Feature.Type.Description=="FANTOM enhancers, robust",],
                                             keep.extra.columns = F,ignore.strand = T,seqnames.field = "Chromosome.Name",
                                             start.field = "Start..bp.",end.field="End..bp.")

all.enhancers.gr <- c(encode.enhancers.gr,vista.gr,fantom.robust.gr)
all.enhancers.gr <- reduce(all.enhancers.gr)

loop.list <- list.files("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\loops\\loop_data\\",full.names = T,pattern="\\.gz$")
loop.list <- lapply(loop.list,read.delim,header=T,stringsAsFactors=F)
loops.df <- as.data.frame(rbindlist(loop.list))
loops.df$chr1 <- ifelse(nchar(loops.df$chr1)>3,sapply(loops.df$chr1,substring,first=4),loops.df$chr1)
loops.df$chr2 <- ifelse(nchar(loops.df$chr2)>3,sapply(loops.df$chr2,substring,first=4),loops.df$chr2)
loops.df <- loops.df[,1:6]
loops.df <- rbind(loops.df,
                  data.frame(chr1=loops.df$chr2,x1=loops.df$y1,x2=loops.df$y2,
                             chr2=loops.df$chr1,y1=loops.df$x1,y2=loops.df$x2))
loops.df <- unique(loops.df)

loops.gr <- makeGRangesFromDataFrame(loops.df,
                                                keep.extra.columns = T,ignore.strand = T,seqnames.field = "chr1",
                                                start.field = "x1",end.field="x2")

enhancer.overlaps <- findOverlaps(loops.gr,all.enhancers.gr,ignore.strand=F)

loops.gr$enhancer.start <- NA
loops.gr$enhancer.end <- NA
loops.gr$enhancer.start[queryHits(enhancer.overlaps)] <- start(all.enhancers.gr)[subjectHits(enhancer.overlaps)]
loops.gr$enhancer.end[queryHits(enhancer.overlaps)] <- end(all.enhancers.gr)[subjectHits(enhancer.overlaps)]
loops.gr <- loops.gr[!is.na(loops.gr$enhancer.start)]

ranges(loops.gr) <- IRanges(start=loops.gr$y1,end=loops.gr$y2)

.encode.promoters.distal.gr <- ComplexID:::.encode.promoters.distal.gr
distal.overlaps <- findOverlaps(.encode.promoters.distal.gr,loops.gr)
.encode.promoters.distal.gr$enhancer.starts <- NA
.encode.promoters.distal.gr$enhancer.ends <- NA
s <- split(1:length(distal.overlaps),queryHits(distal.overlaps),drop=T)
.encode.promoters.distal.gr$enhancer.starts[as.integer(names(s))] <- sapply(s, function(x) {
  return(loops.gr$enhancer.start[subjectHits(distal.overlaps)[x]])
})
.encode.promoters.distal.gr$enhancer.ends[as.integer(names(s))] <- sapply(s, function(x) {
  return(loops.gr$enhancer.end[subjectHits(distal.overlaps)[x]])
})

.encode.promoters.prox.gr <- ComplexID:::.encode.promoters.prox.gr
prox.overlaps <- findOverlaps(.encode.promoters.prox.gr,loops.gr)
.encode.promoters.prox.gr$enhancer.starts <- NA
.encode.promoters.prox.gr$enhancer.ends <- NA
s <- split(1:length(prox.overlaps),queryHits(prox.overlaps),drop=T)
.encode.promoters.prox.gr$enhancer.starts[as.integer(names(s))] <- sapply(s, function(x) {
  return(loops.gr$enhancer.start[subjectHits(prox.overlaps)[x]])
})
.encode.promoters.prox.gr$enhancer.ends[as.integer(names(s))] <- sapply(s, function(x) {
  return(loops.gr$enhancer.end[subjectHits(prox.overlaps)[x]])
})

# create internal data
.eqtl.gr <- eqtl.gr
names(mcols(.eqtl.gr))[1] <- "genes"

.gene.annotation.gr <- ComplexID:::.gene.annotation.gr
.encode.promoters.distal.gr <- ComplexID:::.encode.promoters.distal.gr
.encode.promoters.distal.gr$enhancer.starts <- NULL
.encode.promoters.distal.gr$enhancer.ends <- NULL
.encode.promoters.prox.gr <- ComplexID:::.encode.promoters.prox.gr
.encode.promoters.prox.gr$enhancer.starts <- NULL
.encode.promoters.prox.gr$enhancer.ends <- NULL
.density <- ComplexID:::.density
.corum.subunits <- ComplexID:::.corum.subunits
.W_norm_t_spars <- ComplexID:::.W_norm_t_spars
.geneToComplex <- ComplexID:::.geneToComplex
.utr.entrez.gr <- ComplexID:::.utr.entrez.gr
.hugoNames <- ComplexID:::.hugoNames
.complexNames <- ComplexID:::.complexNames
.complete.igraph <- ComplexID:::.complete.igraph
.ent.to.hug <- ComplexID:::.ent.to.hug
.eqtl.gr <- ComplexID:::.eqtl.gr
.loops.gr <- ComplexID:::.loops.gr
.all.enhancers.gr <- all.enhancers.gr

use_data(.gene.annotation.gr,.encode.promoters.distal.gr,.encode.promoters.prox.gr,.density,.corum.subunits,.W_norm_t_spars,.geneToComplex,.utr.entrez.gr,.hugoNames,.complexNames,.complete.igraph,.ent.to.hug,.eqtl.gr,.loops.gr,.all.enhancers.gr,internal=T,overwrite=T)

rm(.all.enhancers.gr,.loops.gr,.ent.to.hug,.complete.igraph,.complexNames,.hugoNames,.utr.entrez.gr,.geneToComplex,.W_norm_t_spars,.corum.subunits,
   .density,.encode.promoters.prox.gr,.encode.promoters.distal.gr,.gene.annotation.gr)
# redo loops for 5 human normal cell lines

# rename some of the corum complexes:
corum.annot$ComplexName[corum.annot$ComplexID==629] <- "BLM-TBPL1"
corum.annot$ComplexName[corum.annot$ComplexID==1218] <- "BLM-TERF2"
corum.annot$ComplexName[corum.annot$ComplexID==2318] <- "ITGA6-ITGB4-LAMB1-Laminin10/12 complex"
corum.annot$ComplexName[corum.annot$ComplexID==2319] <- "ITGA6-ITGB4-LAMB2-Laminin10/12 complex"
corum.annot$ComplexName[corum.annot$ComplexID==230] <- "Large Mediator Complex"
corum.annot$ComplexName[corum.annot$ComplexID==5450] <- "Mediator OPA1 Complex"
corum.annot$ComplexName[corum.annot$ComplexID==1193] <- "Rap1-PARP1 complex"
corum.annot$ComplexName[corum.annot$ComplexID==1204] <- "Rap1-POT1 complex"
corum.annot$ComplexName[corum.annot$ComplexID==1132] <- "RFC4-RIalpha complex"
corum.annot$ComplexName[corum.annot$ComplexID==3064] <- "RNA polymerase II TBP complex, chromatin structure modifying"
corum.annot$ComplexName[corum.annot$ComplexID==3065] <- "RNA polymerase II ACTL6A complex, chromatin structure modifying"
corum.annot$ComplexName[corum.annot$ComplexID==3066] <- "RNA polymerase II ERCC3 complex, chromatin structure modifying"

.complexNames <- sapply(ComplexID:::.geneToComplex, function(x) {
  if (x[1] <= length(ComplexID:::.corum.subunits))
    return(paste(unique(corum.annot$ComplexName[match(names(x),corum.annot$ComplexID)]),collapse=";"))
  return("None")
})

use_data(.gene.annotation.gr,.encode.promoters.distal.gr,.encode.promoters.prox.gr,.density,.corum.subunits,.W_norm_t_spars,.geneToComplex,.utr.entrez.gr,.hugoNames,.complexNames,.complete.igraph,.ent.to.hug,.eqtl.gr,.loops.gr,internal=T,overwrite=T)

#####
#####
#####
#####
#####
##### MAKE NEW DATABASE - INSTEAD OF REMOVING GENES THAT DON'T EXIST IN STRING OR CORUM
##### ADD THEM INTO THE DATABASE AS "IMMOBILE" VERTICES
#####
#####
#####
#####
#####

gene.annotations <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\gene_annotations\\ensembl_gene_annotations_protein_coding.txt.gz",stringsAsFactors = F)
#names(gene.annotations) <- c("Ensembl.Gene.ID","Chromosome.Name","Gene.Start..bp.","Gene.End..bp.","Strand","Transcription.Start.Site..TSS.","genes" )
gene.annotations$Strand <- ifelse(gene.annotations$Strand<0,"-","+")
dup <- duplicated(gene.annotations[,c(3:8,11)])
gene.annotations <- gene.annotations[!dup,]

dup <- duplicated(gene.annotations[,4:8]) | duplicated(gene.annotations[,4:8],fromLast = T)
gene.annotations$dup <- dup
gene.annotations$transgene <- paste(gene.annotations$Chromosome.scaffold.name,
                                    gene.annotations$Gene.start..bp.,
                                    gene.annotations$Gene.end..bp.,
                                    gene.annotations$Strand,
                                    gene.annotations$Transcription.start.site..TSS.)

gene.annotations$genecoords <- paste(gene.annotations$Chromosome.scaffold.name,
                                     gene.annotations$Gene.start..bp.,
                                     gene.annotations$Gene.end..bp.,
                                     gene.annotations$Strand)

s <- split(gene.annotations$Gene.stable.ID,gene.annotations$genecoords,drop=T)
s.unique <- sapply(s,function(x) length(unique(x)))

gene.annotations[gene.annotations$genecoords %in% names(s)[s.unique>1],]

gene.annotations$Gene.name[gene.annotations$Gene.name=="APITD1-CORT"] <- "APITD1"
gene.annotations$Gene.name[gene.annotations$Gene.name=="MFRP"] <- "C1QTNF5"
gene.annotations$Gene.name[gene.annotations$Gene.name=="POC1B-GALNT4"] <- "GALNT4"

gene.annotations[gene.annotations$genecoords %in% names(s)[s.unique>1],]

#gene.annotations[gene.annotations$genecoords %in% names(s)[s.unique>1],]







non.protein.annotations <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\gene_annotations\\ensembl_gene_annotations_non_protein_coding.txt.gz",stringsAsFactors = F)
dup <- duplicated(non.protein.annotations[,c(3:8,11)])
non.protein.annotations <- non.protein.annotations[!dup,]

all.protein.annotations <- rbind(gene.annotations[gene.annotations$Protein.stable.ID!="",1:12],
                                 non.protein.annotations[non.protein.annotations$Protein.stable.ID!="",1:12])

non.protein.annotations <- non.protein.annotations[!(non.protein.annotations$Gene.name %in% all.protein.annotations$Gene.name),1:12]

string.ppi <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\string_ppi\\9606.protein.links.v10.5.txt.gz",sep=" ",stringsAsFactors = F)
string.ppi$protein1 <- substring(string.ppi$protein1,first = 6)
string.ppi$protein2 <- substring(string.ppi$protein2,first = 6)

# filter STRING PPI for interactions
string.ppi.filt <- string.ppi[string.ppi$combined_score>700,]
a <- paste(string.ppi.filt$protein1,string.ppi.filt$protein2,sep="_")
b <- paste(string.ppi.filt$protein2,string.ppi.filt$protein1,sep="_")
sum(a %in% b)
which(a %in% b)[1:10]

# annotate the string ppi with ensembl gene names
string.ppi.annot <- merge(string.ppi.filt,all.protein.annotations[,c("Protein.stable.ID","EntrezGene.ID","Gene.name","Source.of.gene.name")],by.x="protein1",by.y="Protein.stable.ID",all.x=T)
string.ppi.annot <- merge(string.ppi.annot,all.protein.annotations[,c("Protein.stable.ID","EntrezGene.ID","Gene.name","Source.of.gene.name")],by.x="protein2",by.y="Protein.stable.ID",all.x=T)
sum(is.na(string.ppi.annot$Gene.name.x) | is.na(string.ppi.annot$Gene.name.y))
string.ppi.annot <- string.ppi.annot[!is.na(string.ppi.annot$Gene.name.x) & !is.na(string.ppi.annot$Gene.name.y),]

string.ordered.name <- pbapply(string.ppi.annot,1,function(x) {
  paste(x[c(4,6)][order(as.character(x[c(4,6)]))],collapse = "_")
})







# load CORUM
corum.annot <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\corum\\coreComplexes.txt",stringsAsFactors = F)
corum.annot <- corum.annot[corum.annot$Organism=="Human",]

corum.annot$subunits.Gene.name.syn. <- gsub("None","",corum.annot$subunits.Gene.name.syn.)
corum.annot$subunits.Gene.name.syn. <- gsub(",","",corum.annot$subunits.Gene.name.syn.)

corum.annot.sep <- do.call(rbind, apply(corum.annot, 1,function(r) {
  orig=trimws(unlist(strsplit(paste0(unlist(r[17]),";"), ";")))
  alt=trimws(unlist(strsplit(paste0(unlist(r[18]),";"), ";")))
  if (length(orig)!=length(alt))
    ret <- data.frame(orig=orig,alt="MISTAKE_IN_CORUM")
  else
    ret <- data.frame(orig=orig,alt=alt)
  #print(as.data.frame(r)[rep(1,nrow(ret)),])
  ret <- cbind(as.data.frame(as.list(r))[rep(1,nrow(ret)),],ret)
  return(ret)
})
)

corum.annot.sep <- do.call(rbind,
                           apply(corum.annot.sep, 1,
                                 function(r) do.call(expand.grid,
                                                     c(unlist(r[-22]),
                                                       strsplit(trimws(as.character(r[22])), " ")))))

sum(temp$orig %in% string.ppi.annot$Gene.name.x & temp$Var22 %in% string.ppi.annot$Gene.name.y)
temp <- unique(corum.annot.sep[,21:22])
temp[temp$Var22 %in% temp$orig,]



corum.annot.sep <- do.call(rbind, apply(corum.annot, 1,function(r) {
  orig=trimws(unlist(strsplit(paste0(unlist(r[6]),";"), ";")))
  alt=trimws(unlist(strsplit(paste0(unlist(r[17]),";"), ";")))
  if (length(orig)!=length(alt))
    ret <- data.frame(orig=orig,alt="MISTAKE_IN_CORUM")
  else
    ret <- data.frame(orig=orig,alt=alt)
  #print(as.data.frame(r)[rep(1,nrow(ret)),])
  ret <- cbind(as.data.frame(as.list(r))[rep(1,nrow(ret)),],ret)
  return(ret)
})
)

corum.annot.sep <- do.call(rbind,
                           apply(corum.annot.sep, 1,
                                 function(r) do.call(expand.grid,
                                                     c(unlist(r[-22]),
                                                       strsplit(trimws(as.character(r[22])), " ")))))

write.table(unique(corum.annot.sep[,"orig"]),"D:\\temp\\corum_uniprot.txt",quote = F,row.names = F,col.names = F)
edited.uni.prot <- sapply(strsplit(as.character(unique(corum.annot.sep[,"orig"])),"-"),"[[",1)

uni.to.hgnc <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\uniprot_mapping\\corum_uniprot_to_hgnc.txt",stringsAsFactors=F)
length(unique(uni.to.hgnc$From))
sum(duplicated(uni.to.hgnc$From))
uni.to.hgnc[duplicated(uni.to.hgnc$From) | duplicated(uni.to.hgnc$From,fromLast=T),]
unique(uni.to.hgnc$From[duplicated(uni.to.hgnc$From) | duplicated(uni.to.hgnc$From,fromLast=T)])

corum.annot.sep$orig[!(corum.annot.sep$orig %in% uni.to.hgnc$From)]
unique(corum.annot.sep$orig[!(corum.annot.sep$orig %in% uni.to.hgnc$From)])
sum(grepl("-$",corum.annot.sep$orig))

edited.uni.prot <- sapply(strsplit(as.character(unique(corum.annot.sep[,"orig"])),"-"),"[[",1)
write.table(matrix(edited.uni.prot,nrow=length(edited.uni.prot)),"D:\\temp\\corum_uniprot.txt",quote = F,row.names = F,col.names = F)

uni.to.hgnc <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\uniprot_mapping\\corum_uniprot_to_hgnc.txt",stringsAsFactors=F)
length(unique(uni.to.hgnc$From))
sum(duplicated(uni.to.hgnc$From))
uni.to.hgnc[duplicated(uni.to.hgnc$From) | duplicated(uni.to.hgnc$From,fromLast=T),]
unique(uni.to.hgnc$From[duplicated(uni.to.hgnc$From) | duplicated(uni.to.hgnc$From,fromLast=T)])

hugo.syn <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\hugo_synonyms\\hugo_synonyms.txt",stringsAsFactors = F)
all.synonyms <- unique(c(unlist(strsplit(hugo.syn$Previous.Symbols,", ")),unlist(strsplit(hugo.syn$Synonyms,", "))))
all.synonyms <- all.synonyms[all.synonyms!=""]
sum(all.synonyms %in% hugo.syn$Approved.Symbol)

uni.to.hugo <- merge(uni.to.hgnc,hugo.syn,by.x=2,by.y=1,all.x=T)
uni.to.hugo[duplicated(uni.to.hugo$From) | duplicated(uni.to.hugo$From,fromLast=T),]
corum.synonyms <- unique(c(unlist(strsplit(uni.to.hugo$Previous.Symbols,", ")),unlist(strsplit(uni.to.hugo$Synonyms,", "))))
corum.synonyms <- corum.synonyms[corum.synonyms!=""]
sum(corum.synonyms %in% hugo.syn$Approved.Symbol)

s <- split(1:nrow(uni.to.hugo),uni.to.hugo$From,drop=T)
keep.all.names <- c("Q16637","P62158")
uni.to.real.hugo <- lapply(1:length(s),function(x) {
  if (names(s)[x] %in% keep.all.names)
    return(uni.to.hugo[unlist(s[x]),3])
  return(uni.to.hugo[s[[x]],3])
})
names(uni.to.real.hugo) <- names(s)

corum.annot$real.hugo <- apply(corum.annot,1,function(x) {
  prot.names <- unlist(strsplit(x[6],";"))
  prot.names <- sapply(strsplit(prot.names,"-"),"[[",1)
  ret <- uni.to.real.hugo[prot.names]
  ret[sapply(ret,length)==0] <- "None"
  return(paste(sapply(ret,paste,collapse=","),collapse=";"))
})

# DONE WITH CHANGING CORUM
# Now update the string PPI, create a dictionary linking a symbol to its primary HUGO ID
# Account for the multi mapping proteins in CORUM

# create a mapping from a synonym to a real hugo name
hugo.syn.edited <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\hugo_synonyms\\hugo_synonyms.txt",stringsAsFactors = F)
hugo.syn.edited$all.synonyms <- paste(hugo.syn.edited$Previous.Symbols,hugo.syn.edited$Synonyms,sep=", ")
hugo.syn.edited$all.synonyms <- gsub("(^, )|(, $)","",hugo.syn.edited$all.synonyms)
# edit specific proteins for CORUM multimapping
multi.s <- s[sapply(s,length)>1 & !(names(s) %in% keep.all.names)]
approved.n <- lapply(multi.s,function(x) {
  uni.to.hugo$Approved.Symbol[x]
})

for (i in approved.n) {
  app <- paste(i[2:length(i)],collapse=", ")
  syn <- paste(hugo.syn.edited$all.synonyms[hugo.syn.edited$Approved.Symbol %in% i],collapse=", ")
  hugo.syn.edited$all.synonyms[hugo.syn.edited$Approved.Symbol==i[1]] <- gsub("(^, )|(, $)","",paste(app,syn,sep=", "))
  remove <- which(hugo.syn.edited$Approved.Symbol %in% i[2:length(i)])
  hugo.syn.edited <- hugo.syn.edited[-remove,]
}

# create dataframe mapping synonyms to approved names:
library(data.table)
syn.to.app <- apply(hugo.syn.edited, 1, function(x) {
  spl <- unlist(strsplit(x[6],split = ", "))
  return(data.frame("Approved.Symbol"=rep(x[2],length(spl)+1),
                    "Synonym"=c(x[2],spl)))
})
syn.to.app <- as.data.frame(rbindlist(syn.to.app))
syn.to.app <- syn.to.app[!grepl("~",syn.to.app$Approved.Symbol),]
syn.to.app <- unique(syn.to.app)
table(table(syn.to.app$Synonym))
synonyms <- hugo.syn.edited$all.synonyms[hugo.syn.edited$all.synonyms!=""]
synonyms <- unique(unlist(strsplit(synonyms,", ")))
dup.syn<-names(table(syn.to.app$Synonym[syn.to.app$Synonym %in% synonyms]))[table(syn.to.app$Synonym[syn.to.app$Synonym %in% synonyms])>1]
dup.syn.all<-names(table(syn.to.app$Synonym))[table(syn.to.app$Synonym)>1]
all.equal(dup.syn,dup.syn.all)
sum(dup.syn %in% syn.to.app$Approved.Symbol)
syn.to.app[syn.to.app$Synonym %in% dup.syn,][1:10,]

# check how many approved symbols are being mapped to other approved symbols:
app.to.app <- table(syn.to.app$Synonym)[names(table(syn.to.app$Synonym)) %in% syn.to.app$Approved.Symbol & table(syn.to.app$Synonym) > 1]
length(app.to.app)

# remove the links between approved symbols except for those between the rolled up proteins:
approved.names.to.keep <- uni.to.hugo$Approved.Symbol[unlist(multi.s)]
to.remove <- which(syn.to.app$Synonym %in% syn.to.app$Approved.Symbol & !(syn.to.app$Synonym %in% approved.names.to.keep) & syn.to.app$Synonym %in% names(app.to.app))
syn.to.app <- syn.to.app[-to.remove,]
dup.syn<-names(table(syn.to.app$Synonym[syn.to.app$Synonym %in% synonyms]))[table(syn.to.app$Synonym[syn.to.app$Synonym %in% synonyms])>1]
sum(dup.syn %in% syn.to.app$Approved.Symbol)

# for the duplicated synonyms, check if the Entrez Gene ID helps?
# deconvolute HUGO IDs for genes in annotation using Entrez Gene ID:
# For each duplicated synonym, look at its entrez gene ID in Ensembl, hopefully its only 1 Entrez ID

entrez.of.dups <- lapply(dup.syn, function(x) {
  ret <- unique(all.protein.annotations$EntrezGene.ID[all.protein.annotations$Gene.name==x])
  #ret <- ret[!is.na(ret)]
  return(ret)
})
names(entrez.of.dups) <- dup.syn
table(sapply(entrez.of.dups,length))
entrez.of.dups[sapply(entrez.of.dups,length)==1]

# There are only 8 multi-mapping HUGO synonyms that exist within the Ensembl database, 2 of which don't have Entrez IDs, therefore there are still 1,287 synonyms that are multimapping
sum(dup.syn %in% all.protein.annotations$Gene.name)

# check if the duplicated synonyms are in the protein annotations:
dup.syn[dup.syn %in% all.protein.annotations$Gene.name]

# check how many protein annotations are in duplicated synonyms:
all.protein.annotations[all.protein.annotations$Gene.name %in% dup.syn & all.protein.annotations$Source.of.gene.name == "HGNC Symbol",]
unique(all.protein.annotations$Gene.name[all.protein.annotations$Gene.name %in% dup.syn & all.protein.annotations$Source.of.gene.name == "HGNC Symbol"])

all.protein.annotations[all.protein.annotations$Gene.name %in% dup.syn,]

# check if the ambiguous proteins are in the PPI:

string.ppi.annot[string.ppi.annot$Gene.name.x=="ATP6C" | string.ppi.annot$Gene.name.y=="ATP6C",]
string.ppi.annot[string.ppi.annot$Gene.name.x=="CSRP2BP" | string.ppi.annot$Gene.name.y=="CSRP2BP",]

# Edit synonyms to remove ambiguity, see report.
dim(syn.to.app)
syn.to.app <- syn.to.app[!(syn.to.app$Approved.Symbol=="B3GNT2" & syn.to.app$Synonym=="B3GNT1"),]
syn.to.app <- syn.to.app[!(syn.to.app$Approved.Symbol=="C11orf98" & syn.to.app$Synonym=="C11orf48"),]
syn.to.app <- syn.to.app[!(syn.to.app$Approved.Symbol=="ATP6V1C1" & syn.to.app$Synonym=="ATP6C"),]
syn.to.app <- syn.to.app[!(syn.to.app$Approved.Symbol=="BHLHE40" & syn.to.app$Synonym=="STRA13"),]
syn.to.app <- syn.to.app[!(syn.to.app$Approved.Symbol %in% c("CEACAM7","PSG2") & syn.to.app$Synonym=="CEA"),]
syn.to.app <- syn.to.app[!(syn.to.app$Approved.Symbol=="FBXL19" & syn.to.app$Synonym=="CXXC11"),]
syn.to.app <- syn.to.app[!(syn.to.app$Approved.Symbol=="LPCAT1" & syn.to.app$Synonym=="AGPAT9"),]
dim(syn.to.app)



dim(string.ppi.annot)
string.ppi.real.hugo <- merge(string.ppi.annot,syn.to.app,by.x="Gene.name.x",by.y="Synonym",all.x=T)
string.ppi.real.hugo <- merge(string.ppi.real.hugo,syn.to.app,by.x="Gene.name.y",by.y="Synonym",all.x=T)
string.ppi.real.hugo <- string.ppi.real.hugo[,c(10,11)]
string.ppi.real.hugo <- unique(string.ppi.real.hugo)
string.ppi.real.hugo <- string.ppi.real.hugo[!is.na(string.ppi.real.hugo$Approved.Symbol.x) & !is.na(string.ppi.real.hugo$Approved.Symbol.y),]
string.ppi.real.hugo <- string.ppi.real.hugo[string.ppi.real.hugo$Approved.Symbol.x != string.ppi.real.hugo$Approved.Symbol.y,]
library(pbapply)
string.ppi.real.hugo$order_name <- pbapply(string.ppi.real.hugo,1,function(x) {
  paste(x[1:2][order(as.character(x[1:2]))],collapse = "_")
})
table(table(string.ppi.real.hugo$order_name))
sum(duplicated(string.ppi.real.hugo$order_name))
sum(duplicated(string.ppi.real.hugo[1:2,]))
# all interactions are in the PPI from both directions
dim(string.ppi.real.hugo)

# edit protein annotation:
sum(all.protein.annotations$Gene.name %in% syn.to.app$Synonym)
dim(all.protein.annotations)
all.protein.annotations.final <- merge(all.protein.annotations,syn.to.app,by.x=11,by.y=2,all.x=T)
dim(all.protein.annotations.final)
all.protein.annotations.final$final_genes <- ifelse(is.na(all.protein.annotations.final$Approved.Symbol),as.character(all.protein.annotations.final$Gene.name),as.character(all.protein.annotations.final$Approved.Symbol))
dup <- duplicated(all.protein.annotations.final[,c(5:9,14)])
all.protein.annotations.final <- all.protein.annotations.final[!dup,]

# edit UTR annotation:
utr <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\utr_annotation\\ensembl_utr.txt.gz",stringsAsFactors = F)

utr <- data.frame(Gene.name=c(utr$Gene.name[!is.na(utr$X5..UTR.start)],utr$Gene.name[!is.na(utr$X3..UTR.start)]),
                  CHR=c(utr$Chromosome.scaffold.name[!is.na(utr$X5..UTR.start)],utr$Chromosome.scaffold.name[!is.na(utr$X3..UTR.start)]),
                  Start=c(utr$X5..UTR.start[!is.na(utr$X5..UTR.start)],utr$X3..UTR.start[!is.na(utr$X3..UTR.start)]),
                  End=c(utr$X5..UTR.end[!is.na(utr$X5..UTR.start)],utr$X3..UTR.end[!is.na(utr$X3..UTR.start)]))
utr <- unique(utr)
utr <- merge(utr,syn.to.app,by.x=1,by.y=2)

# Output final CORUM complexes, PPI, protein annotations, and non-protein annotations:
write.table(corum.annot,"X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\final_edited_annotations\\corum_annotation.txt",sep = "\t",quote = F,row.names = F,col.names = T)
write.table(string.ppi.real.hugo,"X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\final_edited_annotations\\string_ppi_700.txt",sep = "\t",quote = F,row.names = F,col.names = T)
write.table(all.protein.annotations.final,"X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\final_edited_annotations\\protein_coding_annotations.txt",sep = "\t",quote = F,row.names = F,col.names = T)
write.table(non.protein.annotations,"X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\final_edited_annotations\\non_protein_coding_annotations.txt",sep = "\t",quote = F,row.names = F,col.names = T)
write.table(utr,"X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\final_edited_annotations\\utr_annotations.txt",sep = "\t",quote = F,row.names = F,col.names = T)







######
######
######
######
######
######
###### DONE WITH UPDATED CORUM AND STRING PPI WITH LATEST HUGO NAMES
######
######
######
######
######
######
######
######
######
######

# convert CORUM to list connected each complex to all proteins in it
corum.subunits <- strsplit(corum.annot$real.hugo,"(;|,)")
corum.subunits <- lapply(corum.subunits,function(x) x[which(x != "None")])
names(corum.subunits) <- corum.annot$ComplexID
corum.subunits <- corum.subunits[sapply(corum.subunits,length)>1]

library(ggplot2)
ggplot(mapping=aes(sapply(corum.subunits,length)))+geom_histogram()+labs(x="Size of Complex",y="Count of Complexes",title="Distribution of the Sizes of Complexes in CORUM Core Complexes")
corum.subunits[sapply(corum.subunits,length)>70]


# remove duplicate complexes by alphabetizing the subunit lists, and finding the duplicates
corum.subunits <- lapply(corum.subunits, function(x) {
  return(sort(x,decreasing = F))
  #paste(sort(as.integer(x),decreasing = F),collapse=";")
})
corum.subunits <- corum.subunits[!duplicated(corum.subunits)]

# remove very large complexes that are outliers (>70 subunits)
corum.subunits <- corum.subunits[sapply(corum.subunits,length)<=70]

# examine fully-overlapping complexes
find.large.complex <- sapply(corum.subunits, function(x) {
  sapply(corum.subunits, function(y) {
    if (length(x) >= length(y))
      if (sum(y %in% x) == length(y))
        return(T)
    return(F)
    if (sum(x %in% y) == length(x))
      return(T)
    return(F)
  })
})
num.smaller.complexes<-apply(find.large.complex,2,sum)
num.larger.complexes<-apply(find.large.complex,1,sum)

ggplot(mapping=aes(num.smaller.complexes-1))+geom_bar()+scale_x_continuous("Number of CORUM sub-complexes in complex", 0:9, 0:9)+labs(title="Distribution of the Number of Sub-Complexes in a Complex")
sum(num.smaller.complexes>1)

ggplot(mapping=aes(num.larger.complexes-1))+geom_bar()+scale_x_continuous("Number of CORUM complexes in another Complex", 0:25, 0:25)+labs(title="Distribution of the Number of Complexes that are part of another Complex")
sum(num.larger.complexes>1)
# remove all complexes that have larger complex that completely incorporate them
corum.subunits <- corum.subunits[num.larger.complexes==1]

# convert STRING's ENSP-ENSP interactions to Entrez-Entrez gene interactions and remove any duplicates
# string.ppi.entrez <- string.ppi.annot[,c("EntrezGene.ID.x","EntrezGene.ID.y")]
# string.ppi.entrez <- unique(string.ppi.entrez)
# string.ppi.entrez <- rbind(string.ppi.entrez,data.frame(EntrezGene.ID.x=string.ppi.entrez$EntrezGene.ID.y,
#                                                         EntrezGene.ID.y=string.ppi.entrez$EntrezGene.ID.x))
# string.ppi.entrez <- unique(string.ppi.entrez)

# remove interactions between the same gene and itself
# string.ppi.entrez <- string.ppi.entrez[string.ppi.entrez$EntrezGene.ID.x!=string.ppi.entrez$EntrezGene.ID.y,]

# create dictionary, linking each gene to all of its neighbors.
geneToGene <- split(string.ppi.real.hugo$Approved.Symbol.x,string.ppi.real.hugo$Approved.Symbol.y,drop=T)

# check if there are any complexes that do not have a single protein in the PPI
which(sapply(corum.subunits, function(x) { sum(x %in% names(geneToGene)) })==0)
corum.subunits[which(sapply(corum.subunits, function(x) { sum(x %in% names(geneToGene)) })==0)]

# create network of complexes
C2C.mat <- matrix(0,nrow = length(corum.subunits),ncol=length(corum.subunits))
rownames(C2C.mat) <- colnames(C2C.mat) <- names(corum.subunits)

for (i in 1:nrow(C2C.mat)) {
  row.subunits <- corum.subunits[[i]]
  for (j in 1:i) {
    if (i==j)
      next
    col.subunits <- corum.subunits[[j]]
    intersection <- row.subunits %in% col.subunits
    col.intersection <- col.subunits %in% row.subunits
    denom <- (length(row.subunits) - sum(intersection))*(length(col.subunits) - sum(col.intersection))
    col.subunits <- col.subunits[!col.intersection]
    row.subunits.filt <- row.subunits[!intersection]

    row.interactions <- unlist(geneToGene[as.character(row.subunits.filt)])
    #if(sum(as.integer(col.subunits) %in% unlist(row.interactions)))
    C2C.mat[i,j] <- sum(sapply(col.subunits,function(x) length(which(x==row.interactions))))/denom
    #print(j)
  }
  if (i %% 100 == 0)
    print(i)
}
C2C.mat[upper.tri(C2C.mat)] <- t(C2C.mat)[upper.tri(C2C.mat)]

# create list of individual protein complexes (proteins that are in the STRING DB but not in CORUM)
indivToGene <- geneToGene[!(names(geneToGene) %in% unlist(corum.subunits))]

# create I2C and C2I matrices to connect proteins not in CORUM to complexes in CORUM
I2C.mat <- matrix(0,nrow = length(indivToGene),ncol=length(corum.subunits))
rownames(I2C.mat) <- names(indivToGene)
colnames(I2C.mat) <- names(corum.subunits)

indivGeneDegrees <- sapply(indivToGene,length)
cAproteins <- lapply(corum.subunits,function(x) unlist(geneToGene[as.character(x)]) )

for (i in 1:nrow(I2C.mat)) {
  geneI <- names(indivToGene[i])
  for (j in 1:ncol(I2C.mat)) {
    cA <- cAproteins[[j]]
    I2C.mat[i,j] <- sum(geneI == cA)
  }
  if (i %% 100 == 0)
    print(i)
}

C2I.mat <- t(I2C.mat)

I2C.mat <- I2C.mat/indivGeneDegrees
C2I.mat <- C2I.mat/sapply(corum.subunits,length)

# Create I2I matrix to connect proteins not in CORUM with other proteins not in CORUM
I2I.mat <- matrix(0,nrow = length(indivToGene),ncol=length(indivToGene))
rownames(I2I.mat) <- names(indivToGene)
colnames(I2I.mat) <- names(indivToGene)

for (i in 1:nrow(I2I.mat)) {
  I2I.mat[i,which(names(indivToGene) %in% indivToGene[[i]])] <- 1
  if (i %% 100 == 0)
    print(i)
}

I2I.mat <- I2I.mat/indivGeneDegrees

W <- rbind(cbind(C2C.mat,C2I.mat),cbind(I2C.mat,I2I.mat))

W_norm_t <- t(W)
W_norm_t <- sweep(W_norm_t,2,colSums(W_norm_t),'/')
sum(is.na(W_norm_t))
W_norm_t_spars <- Matrix(W_norm_t, sparse = TRUE)


# Calculate density(Ca) for all complexes
density <- sapply(corum.subunits, function(x) {
  all.interactions <- unlist(geneToGene[as.character(x)])
  edges <- sum(sapply(x,function(y) sum(y == all.interactions)))/2
  return((2*edges)/(length(x)*(length(x)-1)))
})

ggplot(mapping=aes(density))+geom_histogram()+labs(x="Density of Complex",title="Distribution of Density of Multi-Protein Complexes")
table(sapply(corum.subunits,length)[density==1])
table(sapply(corum.subunits,length)[density==0])
min(density[density!=0])

density <- density+min(density[density!=0])
density <- ifelse(density>1,1,density)
ggplot(mapping=aes(density))+geom_histogram()+labs(x="Density of Complex",title="Distribution of Density of Multi-Protein Complexes\nAfter Adjustment")

string.ppi.real.hugo$Approved.Symbol.x <- as.character(string.ppi.real.hugo$Approved.Symbol.x)
string.ppi.real.hugo$Approved.Symbol.y <- as.character(string.ppi.real.hugo$Approved.Symbol.y)

all.genes.in.network <- unique(c(string.ppi.real.hugo$Approved.Symbol.x,string.ppi.real.hugo$Approved.Symbol.y,unlist(corum.subunits)))

geneToComplex <- pblapply(all.genes.in.network, function(x) {
  corum <- which(sapply(corum.subunits, function(y) x %in% y))
  indiv.names <- rownames(W)[(length(corum.subunits)+1):nrow(W)]
  return(c(corum,which(as.character(x) == indiv.names)+length(corum.subunits)))
})
names(geneToComplex) <- all.genes.in.network

# rename some of the corum complexes:
corum.annot$ComplexName[corum.annot$ComplexID==629] <- "BLM-TBPL1"
corum.annot$ComplexName[corum.annot$ComplexID==1218] <- "BLM-TERF2"
corum.annot$ComplexName[corum.annot$ComplexID==2318] <- "ITGA6-ITGB4-LAMB1-Laminin10/12 complex"
corum.annot$ComplexName[corum.annot$ComplexID==2319] <- "ITGA6-ITGB4-LAMB2-Laminin10/12 complex"
corum.annot$ComplexName[corum.annot$ComplexID==230] <- "Large Mediator Complex"
corum.annot$ComplexName[corum.annot$ComplexID==5450] <- "Mediator OPA1 Complex"
corum.annot$ComplexName[corum.annot$ComplexID==1193] <- "Rap1-PARP1 complex"
corum.annot$ComplexName[corum.annot$ComplexID==1204] <- "Rap1-POT1 complex"
corum.annot$ComplexName[corum.annot$ComplexID==1132] <- "RFC4-RIalpha complex"
corum.annot$ComplexName[corum.annot$ComplexID==3064] <- "RNA polymerase II TBP complex, chromatin structure modifying"
corum.annot$ComplexName[corum.annot$ComplexID==3065] <- "RNA polymerase II ACTL6A complex, chromatin structure modifying"
corum.annot$ComplexName[corum.annot$ComplexID==3066] <- "RNA polymerase II ERCC3 complex, chromatin structure modifying"

complexNames <- pbsapply(geneToComplex, function(x) {
  if (x[1] <= length(geneToComplex))
    return(paste(unique(corum.annot$ComplexName[match(names(x),corum.annot$ComplexID)]),collapse=";"))
  return("None")
})
complexNames[complexNames=="NA"] <- "None"

# Construct igraph object of PPI, add edges to it of genes not in PPI but in CORUM
library(igraph)
complex.size <- sapply(corum.subunits,length)
colA <- unlist(mapply(rep,as.integer(names(corum.subunits)),complex.size))
corum.edge.list <- matrix(data = c(colA,unlist(corum.subunits)),ncol = 2)
corum.edge.list <- lapply(corum.subunits,function(x) as.data.frame(t(combn(x,2))))
library(data.table)
corum.edge.list <- as.data.frame(rbindlist(corum.edge.list))

string.ppi.real.hugo$color <- ifelse(paste(string.ppi.real.hugo$Approved.Symbol.x,string.ppi.real.hugo$Approved.Symbol.y,sep=";") %in% c(paste(corum.edge.list[,1],corum.edge.list[,2],sep=";"),paste(corum.edge.list[,2],corum.edge.list[,1],sep=";")),"red","black")
corum.edge.list <- cbind(corum.edge.list,
                         matrix(ifelse(paste(corum.edge.list[,1],corum.edge.list[,2],sep=";") %in% c(paste(string.ppi.real.hugo$Approved.Symbol.x,string.ppi.real.hugo$Approved.Symbol.y,sep=";"),paste(string.ppi.real.hugo$Approved.Symbol.y,string.ppi.real.hugo$Approved.Symbol.x,sep=";")),"red","green"),ncol = 1))

colnames(corum.edge.list) <- names(string.ppi.real.hugo)[c(1,2,4)]
complete.edge.list <- unique(rbind(string.ppi.real.hugo[,c(1,2,4)],corum.edge.list))
complete.edge.list$Approved.Symbol.x <- as.character(complete.edge.list$Approved.Symbol.x)
complete.edge.list$Approved.Symbol.y <- as.character(complete.edge.list$Approved.Symbol.y)

ordered <- ifelse(apply(complete.edge.list[,1:2],1,function(x) order(x)[1]) == 1,paste(complete.edge.list$Approved.Symbol.x,complete.edge.list$Approved.Symbol.y,sep=";"),paste(complete.edge.list$Approved.Symbol.y,complete.edge.list$Approved.Symbol.x,sep=";"))
complete.edge.list <- complete.edge.list[!duplicated(ordered),]

#complete.edge.list$color <-

complete.igraph <- graph_from_edgelist(as.matrix(complete.edge.list[,1:2]),directed=F)
E(complete.igraph)$color <- complete.edge.list$color

# create database of eqtls:
eqtl.files <- list.files("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\gtex\\GTEx_Analysis_v6p_eQTL\\",full.names=T)
eqtl.files <- eqtl.files[grepl("signif_snpgene_pairs",eqtl.files)]

eqtl.list <- lapply(eqtl.files, function(x) {
  print(x)
  infile <- read.delim(x,stringsAsFactors = F)
  s <- strsplit(infile[,1],"_")
  infile$snp_chr <- sapply(s,"[[",1)
  infile$snp_pos <- sapply(s,"[[",2)
  return(infile[,c("gene_id","snp_chr","snp_pos")])
})
library(data.table)
eqtl.df <- as.data.frame(rbindlist(eqtl.list))
eqtl.df <- unique(eqtl.df)
eqtl.df$gene_id <- substr(eqtl.df$gene_id,1,15)
eqtl.df <- unique(eqtl.df)

eqtl.df <- merge(eqtl.df,all.protein.annotations.final[,c("Gene.stable.ID","final_genes")],by=1)

#eqtl.df <- eqtl.df[eqtl.df$EntrezGene.ID %in% as.character(all.entrez.genes.in.network),]
eqtl.df <- eqtl.df[,2:4]
eqtl.df <- unique(eqtl.df)
eqtl.gr <- makeGRangesFromDataFrame(eqtl.df,
                                    keep.extra.columns = T,ignore.strand = T,seqnames.field = "snp_chr",
                                    start.field = "snp_pos",end.field="snp_pos")
eqtl.gr$Feature="eQTL"

# rename objects for use in package

all.protein.annotations.final$Strand <- ifelse(all.protein.annotations.final$Strand %in% c("-","+"),all.protein.annotations.final$Strand,
                                               ifelse(all.protein.annotations.final$Strand == "-1","-","+"))
gene.annotation.gr <- makeGRangesFromDataFrame(all.protein.annotations.final[,c(5:8,14,9)],
                                               keep.extra.columns = T,strand.field = "Strand",seqnames.field = "Chromosome.scaffold.name",
                                               start.field = "Gene.Start..bp.",end.field="Gene.End..bp.")

gene.annotation.gr$Feature <- "Gene_Body"
names(mcols(gene.annotation.gr))[1] <- "genes"

.gene.annotation.gr <- ComplexID:::.gene.annotation.gr
.encode.promoters.distal.gr <- ComplexID:::.encode.promoters.distal.gr
.encode.promoters.distal.gr$enhancer.starts <- NULL
.encode.promoters.distal.gr$enhancer.ends <- NULL
.encode.promoters.prox.gr <- ComplexID:::.encode.promoters.prox.gr
.encode.promoters.prox.gr$enhancer.starts <- NULL
.encode.promoters.prox.gr$enhancer.ends <- NULL
.density <- density
.corum.subunits <- corum.subunits
.W_norm_t_spars <- W_norm_t_spars
.geneToComplex <- geneToComplex

utr.gr <- makeGRangesFromDataFrame(utr[,2:5],
                                          keep.extra.columns = T,ignore.strand = T,seqnames.field = "CHR",
                                          start.field = "Start",end.field="End")

utr.gr$Feature <- "UTR"
names(mcols(utr.gr))[1] <- "genes"

.utr.gr <- utr.gr
.complexNames <- complexNames
.complete.igraph <- complete.igraph
.eqtl.gr <- eqtl.gr
.loops.gr <- ComplexID:::.loops.gr
.all.enhancers.gr <- ComplexID:::.all.enhancers.gr
.non.ppi.genes <- unique(gene.annotation.gr$genes[!(gene.annotation.gr$genes %in% all.genes.in.network)])

use_data(.gene.annotation.gr,.encode.promoters.distal.gr,.encode.promoters.prox.gr,.density,.corum.subunits,.W_norm_t_spars,.geneToComplex,.utr.gr,.complexNames,.complete.igraph,.eqtl.gr,.loops.gr,.all.enhancers.gr,.non.ppi.genes,internal=T,overwrite=T)

rm(.non.ppi.genes,.all.enhancers.gr,.loops.gr,.complete.igraph,.complexNames,.utr.gr,.geneToComplex,.W_norm_t_spars,.corum.subunits,
   .density,.encode.promoters.prox.gr,.encode.promoters.distal.gr,.gene.annotation.gr)

# messed up writing package:

all.protein.annotations.final <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\final_edited_annotations\\protein_coding_annotations.txt",stringsAsFactors = F)
all.protein.annotations.final$Strand <- ifelse(all.protein.annotations.final$Strand %in% c("-","+"),all.protein.annotations.final$Strand,
                                               ifelse(all.protein.annotations.final$Strand == "-1","-","+"))
gene.annotation.gr <- makeGRangesFromDataFrame(all.protein.annotations.final[,c(5:8,14,9)],
                                               keep.extra.columns = T,strand.field = "Strand",seqnames.field = "Chromosome.scaffold.name",
                                               start.field = "Gene.Start..bp.",end.field="Gene.End..bp.")

gene.annotation.gr$Feature <- "Gene_Body"
names(mcols(gene.annotation.gr))[1] <- "genes"

.gene.annotation.gr <- gene.annotation.gr
.encode.promoters.distal.gr <- ComplexID:::.encode.promoters.distal.gr
.encode.promoters.prox.gr <- ComplexID:::.encode.promoters.prox.gr
.density <- density
.corum.subunits <- ComplexID:::.corum.subunits
.W_norm_t_spars <- W_norm_t_spars
.geneToComplex <- geneToComplex
.utr.gr <- ComplexID:::.utr.gr
.utr.gr$genes <- as.character(.utr.gr$genes)
.complexNames <- complexNames
.complete.igraph <- complete.igraph
.eqtl.gr <- ComplexID:::.eqtl.gr
names(mcols(.eqtl.gr))[1] <- "genes"
.loops.gr <- ComplexID:::.loops.gr
.all.enhancers.gr <- ComplexID:::.all.enhancers.gr
.non.ppi.genes <- ComplexID:::.non.ppi.genes

use_data(.gene.annotation.gr,.encode.promoters.distal.gr,.encode.promoters.prox.gr,.density,.corum.subunits,.W_norm_t_spars,.geneToComplex,.utr.gr,.complexNames,.complete.igraph,.eqtl.gr,.loops.gr,.all.enhancers.gr,.non.ppi.genes,internal=T,overwrite=T)

rm(.non.ppi.genes,.all.enhancers.gr,.loops.gr,.complete.igraph,.complexNames,.utr.gr,.geneToComplex,.W_norm_t_spars,.corum.subunits,
   .density,.encode.promoters.prox.gr,.encode.promoters.distal.gr,.gene.annotation.gr)


# data set: PGRS univariate test P-VALUE 0.0001

pgrs.uni <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\data\\results_PGRS_SNPs_pleio_Univariate_4outcomes_n656_ANNOTATIONS.txt")

pgrs.uni.sig.wm <- pgrs.uni[pgrs.uni$WMM.Pval <0.0001,c("SNP","CHR","BP")]
pgrs.uni.sig.wm$pheno <- "WM"
pgrs.uni.sig.aro <- pgrs.uni[pgrs.uni$AROUSAL.Pval <0.0001,c("SNP","CHR","BP")]
pgrs.uni.sig.aro$pheno <- "Arousal"
pgrs.uni.sig.hyper <- pgrs.uni[pgrs.uni$PARENT_HYPER.Pval <0.0001,c("SNP","CHR","BP")]
pgrs.uni.sig.hyper$pheno <- "Hyperactivity"
pgrs.uni.sig.inat <- pgrs.uni[pgrs.uni$PARENT_INATT.Pval <0.0001,c("SNP","CHR","BP")]
pgrs.uni.sig.inat$pheno <- "Inattention"

pgrs.uni.sig <- rbind(pgrs.uni.sig.wm,pgrs.uni.sig.aro,pgrs.uni.sig.hyper,pgrs.uni.sig.inat)

pgrs.uni.sig.gr <- makeGRangesFromDataFrame(df = pgrs.uni.sig,keep.extra.columns = T,ignore.strand = T,seqnames.field = "CHR",start.field = "BP",end.field = "BP")

phenosim.mat <- matrix(c("WM","Arousal","Hyperactivity","Inattention",1,1,1,1),ncol=2)

system.time(test <- runComplexID(pgrs.uni.sig.gr,phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T))

test <- annotateHits(Hits = pgrs.uni.sig.gr,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T)

generatePlot(test$scores$HUGO.Gene.Name[test$scores$Gene.in.Network=="Yes"][1:100],order=0,vertex.size=2,vertex.label=NA)

# try to put in non-protein coding annotation:
non.protein.annotations <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\final_edited_annotations\\non_protein_coding_annotations.txt",stringsAsFactors = F)
non.protein.annotations$Gene.name[non.protein.annotations$Gene.name %in% gene.annotation.gr$genes]
non.protein.annotations <- non.protein.annotations[!(non.protein.annotations$Gene.name %in% gene.annotation.gr$genes | non.protein.annotations$Gene.name %in% unlist(corum.subunits)),]
non.protein.annotations$Strand <- ifelse(non.protein.annotations$Strand <0 ,"-","+")

.non.protein.annotations.gr <- makeGRangesFromDataFrame(non.protein.annotations[,c(4:6,11)],
                                               keep.extra.columns = T,ignore.strand = T,seqnames.field = "Chromosome.scaffold.name",
                                               start.field = "Gene.Start..bp.",end.field="Gene.End..bp.")

.non.protein.annotations.gr$Feature <- "Non_protein_region"
names(mcols(.non.protein.annotations.gr))[1] <- "genes"
.non.protein.regions <- unique(.non.protein.annotations.gr$genes)


.gene.annotation.gr <- ComplexID:::.gene.annotation.gr
.encode.promoters.distal.gr <- ComplexID:::.encode.promoters.distal.gr
.encode.promoters.prox.gr <- ComplexID:::.encode.promoters.prox.gr
.density <- ComplexID:::.density
.corum.subunits <- ComplexID:::.corum.subunits
.W_norm_t_spars <- ComplexID:::.W_norm_t_spars
.geneToComplex <- ComplexID:::.geneToComplex
.utr.gr <- ComplexID:::.utr.gr
.complexNames <- ComplexID:::.complexNames
.complete.igraph <- ComplexID:::.complete.igraph
.eqtl.gr <- ComplexID:::.eqtl.gr
.loops.gr <- ComplexID:::.loops.gr
.all.enhancers.gr <- ComplexID:::.all.enhancers.gr
.non.ppi.genes <- ComplexID:::.non.ppi.genes


use_data(.gene.annotation.gr,.encode.promoters.distal.gr,.encode.promoters.prox.gr,.density,.corum.subunits,.W_norm_t_spars,.geneToComplex,.utr.gr,.complexNames,.complete.igraph,.eqtl.gr,.loops.gr,.all.enhancers.gr,.non.ppi.genes,.non.protein.annotations.gr,.non.protein.regions,internal=T,overwrite=T)

rm(.non.protein.annotations.gr,.non.protein.regions,.non.ppi.genes,.all.enhancers.gr,.loops.gr,.complete.igraph,.complexNames,.utr.gr,.geneToComplex,.W_norm_t_spars,.corum.subunits,
   .density,.encode.promoters.prox.gr,.encode.promoters.distal.gr,.gene.annotation.gr)

system.time(test <- runComplexID(pgrs.uni.sig.gr,phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T,non_proteins=T))
system.time(test <- runComplexID(pgrs.uni.sig.gr,phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T))

system.time(test <- runComplexID(pgrs.uni.sig.gr,phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T,non_proteins=T))
generatePlot(test$scores$HUGO.Gene.Name[test$scores$Gene.in.Network=="Yes"][1:100],order=0,vertex.size=2,vertex.label=NA)

test <- annotateHits(Hits = pgrs.uni.sig.gr,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T,non_proteins=T)
test <- annotateHits(Hits = pgrs.uni.sig.gr,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T)

test <- annotateHits(pgc.top.gr[681],promoterRange = 20000,upstream = 0,downstream = 0,promoters = T,gene.body = F,utr = F,eqtl = F,enhancers = T,loopDist = 100000)


# Starting working again on ComplexID on 1/26/2018
# Start by editing eQTL annotation - add tissue annotation so user can filter the eqtl by tissue

eqtl.files <- list.files("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\gtex\\GTEx_Analysis_v6p_eQTL\\",full.names=T)
eqtl.files <- eqtl.files[grepl("signif_snpgene_pairs",eqtl.files)]

eqtl.list <- lapply(eqtl.files, function(x) {
  print(x)
  file.name <- strsplit(x,"\\\\")[[1]]
  file.name <- file.name[length(file.name)]
  tissue <- strsplit(file.name,"_")[[1]][1]
  tissue <- ifelse(tissue=="Whole","Blood",tissue)
  tissue <- ifelse(tissue=="Small","Intestine",tissue)
  infile <- read.delim(x,stringsAsFactors = F)
  infile$tissue <- tissue
  s <- strsplit(infile[,1],"_")
  infile$snp_chr <- sapply(s,"[[",1)
  infile$snp_pos <- sapply(s,"[[",2)
  return(infile[,c("gene_id","snp_chr","snp_pos","tissue")])
})
library(data.table)
eqtl.df <- as.data.frame(rbindlist(eqtl.list))
eqtl.df <- unique(eqtl.df)
eqtl.df$gene_id <- substr(eqtl.df$gene_id,1,15)
eqtl.df <- unique(eqtl.df)

eqtl.df <- merge(eqtl.df,all.protein.annotations.final[,c("Gene.stable.ID","final_genes")],by=1)

#eqtl.df <- eqtl.df[eqtl.df$EntrezGene.ID %in% as.character(all.entrez.genes.in.network),]
eqtl.df <- eqtl.df[,c(2,3,5,4)]
eqtl.df <- unique(eqtl.df)
eqtl.gr <- makeGRangesFromDataFrame(eqtl.df,
                                    keep.extra.columns = T,ignore.strand = T,seqnames.field = "snp_chr",
                                    start.field = "snp_pos",end.field="snp_pos")
eqtl.gr$Feature="eQTL"
eqtl.gr.list <- split(eqtl.gr,eqtl.gr$tissue,drop=T)
eqtl.gr.list <- lapply(eqtl.gr.list,unique)
eqtl.gr.list <- GRangesList(eqtl.gr.list)
eqtl.gr <- unlist(eqtl.gr.list)
names(mcols(eqtl.gr))[1] <- "genes"
rm(eqtl.gr.list)

.gene.annotation.gr <- ComplexID:::.gene.annotation.gr
.encode.promoters.distal.gr <- ComplexID:::.encode.promoters.distal.gr
.encode.promoters.prox.gr <- ComplexID:::.encode.promoters.prox.gr
.density <- ComplexID:::.density
.corum.subunits <- ComplexID:::.corum.subunits
.W_norm_t_spars <- ComplexID:::.W_norm_t_spars
.geneToComplex <- ComplexID:::.geneToComplex
.utr.gr <- ComplexID:::.utr.gr
.complexNames <- ComplexID:::.complexNames
.complete.igraph <- ComplexID:::.complete.igraph
.eqtl.gr <- eqtl.gr
.loops.gr <- ComplexID:::.loops.gr
.all.enhancers.gr <- ComplexID:::.all.enhancers.gr
.non.ppi.genes <- ComplexID:::.non.ppi.genes
.non.protein.annotations.gr <- ComplexID:::.non.protein.annotations.gr
.non.protein.regions <- ComplexID:::.non.protein.regions

use_data(.gene.annotation.gr,.encode.promoters.distal.gr,.encode.promoters.prox.gr,.density,.corum.subunits,.W_norm_t_spars,.geneToComplex,.utr.gr,.complexNames,.complete.igraph,.eqtl.gr,.loops.gr,.all.enhancers.gr,.non.ppi.genes,.non.protein.annotations.gr,.non.protein.regions,internal=T,overwrite=T)

rm(.non.protein.annotations.gr,.non.protein.regions,.non.ppi.genes,.all.enhancers.gr,.loops.gr,.complete.igraph,.complexNames,.utr.gr,.geneToComplex,.W_norm_t_spars,.corum.subunits,
   .density,.encode.promoters.prox.gr,.encode.promoters.distal.gr,.gene.annotation.gr)

# Now edit the enhancer annotation to account for tissues
encode.enhancers.tissues <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\enhancers\\ENCODE_Predictions\\tissue_mapping\\ENCODE_Enhancer_Tissues.txt",stringsAsFactors = F)
encode.enhancers.tissues$Name <- paste0(encode.enhancers.tissues$Name,"_predictions.bed.gz")

encode.list <- list.files("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\enhancers\\ENCODE_Predictions\\",full.names = T,pattern="_prediction")
encode.list.file <- list.files("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\enhancers\\ENCODE_Predictions\\",full.names = F,pattern="_prediction")

encode.enhancers.list <- lapply(encode.list,read.delim,header=F)
encode.enhancers.list <- lapply(1:length(encode.list),function(x) {
  tissue <- encode.enhancers.tissues$tissue[encode.enhancers.tissues$Name==encode.list.file[x]]
  temp <- encode.enhancers.list[[x]]
  temp$V1 <- sapply(temp$V1,substring,first=4)
  temp$tissue <- tissue
  return(temp)
})
library(data.table)
encode.enhancers <- as.data.frame(rbindlist(encode.enhancers.list))
# examine blood, see if merging overlapping regions increases the width of the resulting regions greatly:
blood.tissue <- encode.enhancers[encode.enhancers$tissue=="Blood",]
blood.tissue$V1 <- sapply(blood.tissue$V1,substring,first=4)
blood.tissue.gr <- makeGRangesFromDataFrame(blood.tissue,
                                                keep.extra.columns = F,ignore.strand = T,seqnames.field = "V1",
                                                start.field = "V2",end.field="V3")
library(ggplot2)
dat <- data.frame(width=width(blood.tissue.gr),data="unmerged")
dat <- rbind(dat,data.frame(width=width(reduce(blood.tissue.gr)),data="merged"))
ggplot(dat,mapping=aes(x=width,fill=data))+geom_histogram(alpha=0.75)

ggplot(mapping=aes(x=width(blood.tissue.gr)))+geom_histogram()
ggplot(mapping=aes(x=width(reduce(blood.tissue.gr))))+geom_histogram()

# appears as though merging the enhancers for each tissue is beneficial.
# Merge all ranges for each tissue
encode.enhancers.list <- lapply(encode.enhancers.list,function(x) {
  ret<-makeGRangesFromDataFrame(x,
                           keep.extra.columns = T,ignore.strand = T,seqnames.field = "V1",
                           start.field = "V2",end.field="V3")
  tiss <- ret$tissue[1]
  ret <- reduce(ret)
  ret$tissue <- tiss
  return(ret)
})
encode.enhancers.list <- GRangesList(encode.enhancers.list)
encode.enhancers <- unlist(encode.enhancers.list)

vista.fantom.enhancers <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\EWAS\\Analysis_Peter\\Enhancer_annotation_data\\ensembl_hg19_other_regulatory_regions_vista_fantom_enhancers.txt.gz",stringsAsFactors = F)
vista.gr <- makeGRangesFromDataFrame(vista.fantom.enhancers[vista.fantom.enhancers$Feature.Type=="VISTA Enhancers",],
                                     keep.extra.columns = F,ignore.strand = T,seqnames.field = "Chromosome.Name",
                                     start.field = "Start..bp.",end.field="End..bp.")
vista.gr <- reduce(vista.gr)
fantom.robust.gr <- makeGRangesFromDataFrame(vista.fantom.enhancers[vista.fantom.enhancers$Feature.Type=="FANTOM predictions"
                                                                    & vista.fantom.enhancers$Feature.Type.Description=="FANTOM enhancers, robust",],
                                             keep.extra.columns = F,ignore.strand = T,seqnames.field = "Chromosome.Name",
                                             start.field = "Start..bp.",end.field="End..bp.")
fantom.robust.gr <- reduce(fantom.robust.gr)
vista.gr$tissue <- "VISTA"
fantom.robust.gr$tissue <- "FANTOM"

all.enhancers.gr <- c(encode.enhancers,vista.gr,fantom.robust.gr)

# output to new sysdata file

.gene.annotation.gr <- ComplexID:::.gene.annotation.gr
.encode.promoters.distal.gr <- ComplexID:::.encode.promoters.distal.gr
.encode.promoters.prox.gr <- ComplexID:::.encode.promoters.prox.gr
.density <- ComplexID:::.density
.corum.subunits <- ComplexID:::.corum.subunits
.W_norm_t_spars <- ComplexID:::.W_norm_t_spars
.geneToComplex <- ComplexID:::.geneToComplex
.utr.gr <- ComplexID:::.utr.gr
.complexNames <- ComplexID:::.complexNames
.complete.igraph <- ComplexID:::.complete.igraph
.eqtl.gr <- ComplexID:::.eqtl.gr
.loops.gr <- ComplexID:::.loops.gr
.all.enhancers.gr <- all.enhancers.gr
.non.ppi.genes <- ComplexID:::.non.ppi.genes
.non.protein.annotations.gr <- ComplexID:::.non.protein.annotations.gr
.non.protein.regions <- ComplexID:::.non.protein.regions

use_data(.gene.annotation.gr,.encode.promoters.distal.gr,.encode.promoters.prox.gr,.density,.corum.subunits,.W_norm_t_spars,.geneToComplex,.utr.gr,.complexNames,.complete.igraph,.eqtl.gr,.loops.gr,.all.enhancers.gr,.non.ppi.genes,.non.protein.annotations.gr,.non.protein.regions,internal=T,overwrite=T)

rm(.non.protein.annotations.gr,.non.protein.regions,.non.ppi.genes,.all.enhancers.gr,.loops.gr,.complete.igraph,.complexNames,.utr.gr,.geneToComplex,.W_norm_t_spars,.corum.subunits,
   .density,.encode.promoters.prox.gr,.encode.promoters.distal.gr,.gene.annotation.gr,.eqtl.gr)

system.time(test <- runComplexID(pgrs.uni.sig.gr,phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T))
system.time(test <- runComplexID(pgrs.uni.sig.gr,phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,enhancerTissues = "Blood",utr = T))

# edit promoter annotation to include tissue options:
encode.promoters.tissues <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\Promoter_annotation_data\\tissue_mapping\\ENCODE_Promoter_Tissues.txt",stringsAsFactors = F)
encode.promoters.tissues$Name <- paste0(encode.promoters.tissues$Name,"_predictions.bed.gz")

encode.list <- list.files("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\Promoter_annotation_data\\",full.names = T,pattern="_prediction")
encode.list.file <- list.files("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\Promoter_annotation_data\\",full.names = F,pattern="_prediction")

encode.promoters.list <- lapply(encode.list,read.delim,header=F)
encode.promoters.list <- lapply(1:length(encode.list),function(x) {
  tissue <- encode.promoters.tissues$tissue[encode.promoters.tissues$Name==encode.list.file[x]]
  temp <- encode.promoters.list[[x]]
  temp$V1 <- sapply(temp$V1,substring,first=4)
  temp$tissue <- tissue
  return(temp)
})
library(data.table)
encode.promoters <- as.data.frame(rbindlist(encode.promoters.list))

encode.promoters.list <- lapply(encode.promoters.list,function(x) {
  ret<-makeGRangesFromDataFrame(x,
                                keep.extra.columns = T,ignore.strand = T,seqnames.field = "V1",
                                start.field = "V2",end.field="V3")
  tiss <- ret$tissue[1]
  ret <- reduce(ret)
  ret$tissue <- tiss
  return(ret)
})
encode.promoters.list <- GRangesList(encode.promoters.list)
encode.promoters <- unlist(encode.promoters.list)

.gene.annotation.gr <- ComplexID:::.gene.annotation.gr
.encode.promoters.gr <- encode.promoters
.density <- ComplexID:::.density
.corum.subunits <- ComplexID:::.corum.subunits
.W_norm_t_spars <- ComplexID:::.W_norm_t_spars
.geneToComplex <- ComplexID:::.geneToComplex
.utr.gr <- ComplexID:::.utr.gr
.complexNames <- ComplexID:::.complexNames
.complete.igraph <- ComplexID:::.complete.igraph
.eqtl.gr <- ComplexID:::.eqtl.gr
.loops.gr <- ComplexID:::.loops.gr
.all.enhancers.gr <- ComplexID:::.all.enhancers.gr
.non.ppi.genes <- ComplexID:::.non.ppi.genes
.non.protein.annotations.gr <- ComplexID:::.non.protein.annotations.gr
.non.protein.regions <- ComplexID:::.non.protein.regions

use_data(.gene.annotation.gr,.encode.promoters.gr,.density,.corum.subunits,.W_norm_t_spars,.geneToComplex,.utr.gr,.complexNames,.complete.igraph,.eqtl.gr,.loops.gr,.all.enhancers.gr,.non.ppi.genes,.non.protein.annotations.gr,.non.protein.regions,internal=T,overwrite=T)

rm(.non.protein.annotations.gr,.non.protein.regions,.non.ppi.genes,.all.enhancers.gr,.loops.gr,.complete.igraph,.complexNames,.utr.gr,.geneToComplex,.W_norm_t_spars,.corum.subunits,
   .density,.encode.promoters.gr,.gene.annotation.gr,.eqtl.gr)


system.time(test <- runComplexID(pgrs.uni.sig.gr,phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T))






