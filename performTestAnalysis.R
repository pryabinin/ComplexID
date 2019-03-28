# data set: PGRS univariate test

pgrs.uni <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\GWAS_analysis_Wilmot\\Analysis1_EF\\analysis_1_Pleio\\Results\\SNPS_PGRS\\results_PGRS_SNPs_pleio_Univariate_4outcomes_n656_ANNOTATIONS.txt")

pgrs.uni.sig.wm <- pgrs.uni[pgrs.uni$WMM.Pval <0.05,c("SNP","CHR","BP")]
pgrs.uni.sig.wm$pheno <- "WM"
pgrs.uni.sig.aro <- pgrs.uni[pgrs.uni$AROUSAL.Pval <0.05,c("SNP","CHR","BP")]
pgrs.uni.sig.aro$pheno <- "Arousal"
pgrs.uni.sig.hyper <- pgrs.uni[pgrs.uni$PARENT_HYPER.Pval <0.05,c("SNP","CHR","BP")]
pgrs.uni.sig.hyper$pheno <- "Hyperactivity"
pgrs.uni.sig.inat <- pgrs.uni[pgrs.uni$PARENT_INATT.Pval <0.05,c("SNP","CHR","BP")]
pgrs.uni.sig.inat$pheno <- "Inattention"

pgrs.uni.sig <- rbind(pgrs.uni.sig.wm,pgrs.uni.sig.aro,pgrs.uni.sig.hyper,pgrs.uni.sig.inat)

pgrs.uni.sig.gr <- makeGRangesFromDataFrame(df = pgrs.uni.sig,keep.extra.columns = T,ignore.strand = T,seqnames.field = "CHR",start.field = "BP",end.field = "BP")

phenosim.mat <- matrix(c("WM","Arousal","Hyperactivity","Inattention",1,1,1,1),ncol=2)

system.time(test <- runComplexID(pgrs.uni.sig.gr,phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T))
sapply(gregexpr(";", test$scores$Complexes[1:10]),length)
test$scores$HUGO.Gene.Name[1:10]
sum(test$scores$Num.Hits!=0)

write.table(test$scores,file = "X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\results\\PGRS\\univariate\\PGRS_univariate_pvalue_05.txt",quote = F,sep="\t",row.names = F)
sum(test$scores$Num.Hits!=0)

library(ggplot2)
sum(test$scores$Complexes=="None")
ggplot(mapping = aes(x=sapply(gregexpr(";", test$scores$Complexes[test$scores$Complexes!="None"]),length))) + geom_histogram()+labs(x="Number of Complexes a Gene is Part of",y="Number of Genes",title="Distribution of the Number of Complexes a Gene is part of (must be in at least one complex)")
ggplot(mapping=aes(x=test$scores$score)) +geom_histogram() +labs(x="Score of Gene",y="Number of Genes",title="Distribution of the Number of Complexes a Gene is part of (must be in at least one complex)\nonly showing scores > 1e-2") + xlim(1e-2, 2)


# re-test, using hits with very high p-value:

pgrs.uni.sig.wm.test <- pgrs.uni[pgrs.uni$WMM.Pval >.95,c("SNP","CHR","BP")]
pgrs.uni.sig.wm.test$pheno <- "WM"
pgrs.uni.sig.aro.test <- pgrs.uni[pgrs.uni$AROUSAL.Pval >.95,c("SNP","CHR","BP")]
pgrs.uni.sig.aro.test$pheno <- "Arousal"
pgrs.uni.sig.hyper.test <- pgrs.uni[pgrs.uni$PARENT_HYPER.Pval >.95,c("SNP","CHR","BP")]
pgrs.uni.sig.hyper.test$pheno <- "Hyperactivity"
pgrs.uni.sig.inat.test <- pgrs.uni[pgrs.uni$PARENT_INATT.Pval >.95,c("SNP","CHR","BP")]
pgrs.uni.sig.inat.test$pheno <- "Inattention"

pgrs.uni.sig.test <- rbind(pgrs.uni.sig.wm.test,pgrs.uni.sig.aro.test,pgrs.uni.sig.hyper.test,pgrs.uni.sig.inat.test)

pgrs.uni.sig.test.gr <- makeGRangesFromDataFrame(df = pgrs.uni.sig.test,keep.extra.columns = T,ignore.strand = T,seqnames.field = "CHR",start.field = "BP",end.field = "BP")

system.time(test <- runComplexID(pgrs.uni.sig.test.gr,phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T))
sapply(gregexpr(";", test$scores$Complexes[1:10]),length)
test$scores$HUGO.Gene.Name[1:10]
sum(test$scores$Num.Hits!=0)

#system.time(test <- runComplexID(pgrs.uni.sig.test.gr[4],phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T))
write.table(test$scores,file = "X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\results\\PGRS\\univariate\\PGRS_univariate_pvalue_95.txt",quote = F,sep="\t",row.names = F)


# data set: PGRS univariate test P-VALUE 0.0001

pgrs.uni <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\GWAS_analysis_Wilmot\\Analysis1_EF\\analysis_1_Pleio\\Results\\SNPS_PGRS\\results_PGRS_SNPs_pleio_Univariate_4outcomes_n656_ANNOTATIONS.txt")

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
sum(test$scores$Num.Hits!=0)

ggplot(mapping=aes(x=test$scores$score)) +geom_histogram() +labs(x="Score of Gene",y="Number of Genes",title="Distribution of Scores of Genes from PGRS Analysis, Hit p-values threshold = 0.0001")

sapply(gregexpr(";", test$scores$Complexes[1:10]),length)
test$scores$HUGO.Gene.Name[1:10]
write.table(test$scores,file = "X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\results\\PGRS\\univariate\\PGRS_univariate_pvalue_0001.txt",quote = F,sep="\t",row.names = F)


# adjust phenotype similarity to latent factors from Joel's paper:

phenosim.mat <- matrix(c("WM","Arousal","Hyperactivity","Inattention",.52,.57,.05,.05),ncol=2)

system.time(test <- runComplexID(pgrs.uni.sig.gr,phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T))
test$scores[1:10,]

write.table(test$scores,file = "X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\results\\PGRS\\univariate\\PGRS_univariate_pvalue_0001_phenoSim_latent_corr.txt",quote = F,sep="\t",row.names = F)

# rerun using 77 sites with >0.95 p-values
pgrs.uni.sig.wm.test <- pgrs.uni[pgrs.uni$WMM.Pval >.95,c("SNP","CHR","BP")]
pgrs.uni.sig.wm.test$pheno <- "WM"
pgrs.uni.sig.wm.test <- pgrs.uni.sig.wm.test[1:24,]
pgrs.uni.sig.aro.test <- pgrs.uni[pgrs.uni$AROUSAL.Pval >.95,c("SNP","CHR","BP")]
pgrs.uni.sig.aro.test$pheno <- "Arousal"
pgrs.uni.sig.aro.test <- pgrs.uni.sig.aro.test[1:23,]
pgrs.uni.sig.hyper.test <- pgrs.uni[pgrs.uni$PARENT_HYPER.Pval >.95,c("SNP","CHR","BP")]
pgrs.uni.sig.hyper.test$pheno <- "Hyperactivity"
pgrs.uni.sig.hyper.test <- pgrs.uni.sig.hyper.test[1:20,]
pgrs.uni.sig.inat.test <- pgrs.uni[pgrs.uni$PARENT_INATT.Pval >.95,c("SNP","CHR","BP")]
pgrs.uni.sig.inat.test$pheno <- "Inattention"
pgrs.uni.sig.inat.test <- pgrs.uni.sig.inat.test[1:10,]

pgrs.uni.sig <- rbind(pgrs.uni.sig.wm.test,pgrs.uni.sig.aro.test,pgrs.uni.sig.hyper.test,pgrs.uni.sig.inat.test)

pgrs.uni.sig.gr <- makeGRangesFromDataFrame(df = pgrs.uni.sig,keep.extra.columns = T,ignore.strand = T,seqnames.field = "CHR",start.field = "BP",end.field = "BP")

phenosim.mat <- matrix(c("WM","Arousal","Hyperactivity","Inattention",1,1,1,1),ncol=2)

system.time(test <- runComplexID(pgrs.uni.sig.gr,phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T))
sum(test$scores$Num.Hits!=0)

ggplot(mapping=aes(x=test$scores$score)) +geom_histogram() +labs(x="Score of Gene",y="Number of Genes",title="Distribution of Scores of Genes from PGRS Analysis, 77 Random hits from P-value > 0.95")
ggplot(mapping=aes(x=test$scores$score)) +geom_histogram() +labs(x="Score of Gene",y="Number of Genes",title="Distribution of Scores of Genes from PGRS Analysis, 77 Random hits from P-value > 0.95\nonly showing scores > 1e-2") + xlim(1e-2, 2)


write.table(test$scores,file = "X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\results\\PGRS\\univariate\\PGRS_univariate_pvalue_095_random_77_hits.txt",quote = F,sep="\t",row.names = F)



# Use pleiotropic analysis:
pgrs.pleio <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\GWAS_analysis_Wilmot\\Analysis1_EF\\analysis_1_Pleio\\Results\\SNPS_PGRS\\results_PGRS_SNPs_pleio_combined_4outcomes_n656_ANNOTATIONS.txt",stringsAsFactors = F)
pgrs.pleio <- pgrs.pleio[pgrs.pleio$index.beta!="0",]
pgrs.pleio.hits <- apply(pgrs.pleio,1,function(x) {
  #print(x)
  x <- as.data.frame(t(as.data.frame(x,stringsAsFactors=F)),stringsAsFactors=F)
  #return(x)
  ret <- x[rep(1,length(strsplit(x$index.beta,",")[[1]])),]
  ret$index.beta <- strsplit(x$index.beta,",")[[1]]
  return(ret)
})
library(data.table)
pgrs.pleio.hits <- as.data.frame(rbindlist(pgrs.pleio.hits))
pgrs.pleio.hits$BP <- as.integer(pgrs.pleio.hits$BP)
pgrs.pleio.hits$CHR <- as.integer(pgrs.pleio.hits$CHR)

pgrs.pleio.hits.gr <- makeGRangesFromDataFrame(pgrs.pleio.hits[,c("SNP","index.beta","CHR","BP")],keep.extra.columns = T,ignore.strand = T,seqnames.field = "CHR",
                                        start.field = "BP",end.field="BP")

phenosim.mat <- matrix(c("1","2","3","4",1,1,1,1),ncol=2)
system.time(test <- runComplexID(pgrs.pleio.hits.gr,phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T))
sum(test$scores$Num.Hits!=0)
write.table(test$scores,file = "X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\results\\PGRS\\pleiotropic\\PGRS_pleiotropic.txt",quote = F,sep="\t",row.names = F)

# Use univariate, one endophenotype at a time:

pgrs.uni.sig.wm <- pgrs.uni[pgrs.uni$WMM.Pval <0.0001,c("SNP","CHR","BP")]
pgrs.uni.sig.wm$pheno <- "WM"
pgrs.uni.sig.aro <- pgrs.uni[pgrs.uni$AROUSAL.Pval <0.0001,c("SNP","CHR","BP")]
pgrs.uni.sig.aro$pheno <- "Arousal"
pgrs.uni.sig.hyper <- pgrs.uni[pgrs.uni$PARENT_HYPER.Pval <0.0001,c("SNP","CHR","BP")]
pgrs.uni.sig.hyper$pheno <- "Hyperactivity"
pgrs.uni.sig.inat <- pgrs.uni[pgrs.uni$PARENT_INATT.Pval <0.0001,c("SNP","CHR","BP")]
pgrs.uni.sig.inat$pheno <- "Inattention"

pgrs.uni.sig.gr.wm <- makeGRangesFromDataFrame(df = pgrs.uni.sig.wm,keep.extra.columns = T,ignore.strand = T,seqnames.field = "CHR",start.field = "BP",end.field = "BP")

phenosim.mat <- matrix(c("WM",1),ncol=2)

system.time(test <- runComplexID(pgrs.uni.sig.gr.wm,phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T))
sum(test$scores$Num.Hits!=0)
write.table(test$scores,file = "X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\results\\PGRS\\univariate\\PGRS_univariate_pvalue_0001_working_memory.txt",quote = F,sep="\t",row.names = F)


pgrs.uni.sig.gr.aro <- makeGRangesFromDataFrame(df = pgrs.uni.sig.aro,keep.extra.columns = T,ignore.strand = T,seqnames.field = "CHR",start.field = "BP",end.field = "BP")
phenosim.mat <- matrix(c("Arousal",1),ncol=2)
system.time(test <- runComplexID(pgrs.uni.sig.gr.aro,phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T))
write.table(test$scores,file = "X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\results\\PGRS\\univariate\\PGRS_univariate_pvalue_0001_arousal.txt",quote = F,sep="\t",row.names = F)

pgrs.uni.sig.gr.hyper <- makeGRangesFromDataFrame(df = pgrs.uni.sig.hyper,keep.extra.columns = T,ignore.strand = T,seqnames.field = "CHR",start.field = "BP",end.field = "BP")
phenosim.mat <- matrix(c("Hyperactivity",1),ncol=2)
system.time(test <- runComplexID(pgrs.uni.sig.gr.hyper,phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T))
write.table(test$scores,file = "X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\results\\PGRS\\univariate\\PGRS_univariate_pvalue_0001_hyperactivity.txt",quote = F,sep="\t",row.names = F)

pgrs.uni.sig.gr.inat <- makeGRangesFromDataFrame(df = pgrs.uni.sig.inat,keep.extra.columns = T,ignore.strand = T,seqnames.field = "CHR",start.field = "BP",end.field = "BP")
phenosim.mat <- matrix(c("Inattention",1),ncol=2)
system.time(test <- runComplexID(pgrs.uni.sig.gr.inat,phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T))
apply(test$scores[1:10,],1,function(x) length(gregexpr(";",x[3])[[1]]))
write.table(test$scores,file = "X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\results\\PGRS\\univariate\\PGRS_univariate_pvalue_0001_inattention.txt",quote = F,sep="\t",row.names = F)


# use pgc 2850 data set
pgc2850.uni <- read.csv("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\GWAS_analysis_Wilmot\\Analysis1_EF\\analysis_1_Pleio\\Archive\\Results\\SNPs_top_2850\\results_Univariate_4outcomes.csv")
pgc2850.pleio <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\GWAS_analysis_Wilmot\\Analysis1_EF\\analysis_1_Pleio\\Archive\\Results\\SNPs_top_2850\\results_2850SNPs_4outcomes_n656_pvalpoint1_ANNOTATIONS.txt")

pgc2850.uni <- merge(pgc2850.uni,pgc2850.pleio[,c(1,4,5)],by=1,all.x=T)

chr.bp <- sapply(pgc2850.uni$SNPName[1:19], function(x) {
  s <- unlist(strsplit(as.character(x),"\\."))
  return(c(substring(s[1],4),s[2]))
})
chr.bp <- t(chr.bp)
chr.bp[13,2] <- 103914139

pgc2850.uni$CHR[1:19] <- chr.bp[,1]
pgc2850.uni$BP[1:19] <- chr.bp[,2]

pgc2850.uni.wm <- pgc2850.uni[pgc2850.uni$WMM.Pval <0.05,c("SNPName","CHR","BP")]
pgc2850.uni.wm$pheno <- "WM"
pgc2850.uni.aro <- pgc2850.uni[pgc2850.uni$AROUSAL.Pval <0.05,c("SNPName","CHR","BP")]
pgc2850.uni.aro$pheno <- "Arousal"
pgc2850.uni.hyper <- pgc2850.uni[pgc2850.uni$PARENT_HYPER.Pval <0.05,c("SNPName","CHR","BP")]
pgc2850.uni.hyper$pheno <- "Hyperactivity"
pgc2850.uni.inat <- pgc2850.uni[pgc2850.uni$PARENT_INATT.Pval <0.05,c("SNPName","CHR","BP")]
pgc2850.uni.inat$pheno <- "Inattention"

pgc2850.uni.sig <- rbind(pgc2850.uni.wm,pgc2850.uni.aro,pgc2850.uni.hyper,pgc2850.uni.inat)

pgc2850.uni.sig.gr <- makeGRangesFromDataFrame(df = pgc2850.uni.sig,keep.extra.columns = T,ignore.strand = T,seqnames.field = "CHR",start.field = "BP",end.field = "BP")

phenosim.mat <- matrix(c("WM","Arousal","Hyperactivity","Inattention",1,1,1,1),ncol=2)

system.time(test <- runComplexID(pgc2850.uni.sig.gr,phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T))
sum(test$scores$Num.Hits!=0)
ggplot(mapping=aes(x=test$scores$score)) +geom_histogram() +labs(x="Score of Gene",y="Number of Genes",title="Distribution of Scores of Genes using Hits from PGC meta analysis (p-value < 0.05)")
ggplot(mapping=aes(x=test$scores$score)) +geom_histogram() +labs(x="Score of Gene",y="Number of Genes",title="Distribution of Scores of Genes using Hits from PGC meta analysis (p-value < 0.05)\nonly showing scores > 1e-2") + xlim(1e-2, 2)

write.table(test$scores,file = "X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\results\\PGC2850\\univariate\\PGC2850_univariate_pvalue_05.txt",quote = F,sep="\t",row.names = F)

# decrease p-value cutoff to 0.001

pgc2850.uni.wm <- pgc2850.uni[pgc2850.uni$WMM.Pval <0.001,c("SNPName","CHR","BP")]
pgc2850.uni.wm$pheno <- "WM"
pgc2850.uni.aro <- pgc2850.uni[pgc2850.uni$AROUSAL.Pval <0.001,c("SNPName","CHR","BP")]
pgc2850.uni.aro$pheno <- "Arousal"
pgc2850.uni.hyper <- pgc2850.uni[pgc2850.uni$PARENT_HYPER.Pval <0.001,c("SNPName","CHR","BP")]
pgc2850.uni.hyper$pheno <- "Hyperactivity"
pgc2850.uni.inat <- pgc2850.uni[pgc2850.uni$PARENT_INATT.Pval <0.001,c("SNPName","CHR","BP")]
pgc2850.uni.inat$pheno <- "Inattention"

pgc2850.uni.sig <- pgc2850.uni.inat

pgc2850.uni.sig.gr <- makeGRangesFromDataFrame(df = pgc2850.uni.sig,keep.extra.columns = T,ignore.strand = T,seqnames.field = "CHR",start.field = "BP",end.field = "BP")

phenosim.mat <- matrix(c("Inattention",1),ncol=2)

system.time(test <- runComplexID(pgc2850.uni.sig.gr,phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T))
sum(test$scores$Num.Hits!=0)
write.table(test$scores,file = "X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\results\\PGC2850\\univariate\\PGC2850_univariate_pvalue_001.txt",quote = F,sep="\t",row.names = F)


# p-value cutoff 0.01

pgc2850.uni.wm <- pgc2850.uni[pgc2850.uni$WMM.Pval <0.01,c("SNPName","CHR","BP")]
pgc2850.uni.wm$pheno <- "WM"
pgc2850.uni.aro <- pgc2850.uni[pgc2850.uni$AROUSAL.Pval <0.01,c("SNPName","CHR","BP")]
pgc2850.uni.aro$pheno <- "Arousal"
pgc2850.uni.hyper <- pgc2850.uni[pgc2850.uni$PARENT_HYPER.Pval <0.01,c("SNPName","CHR","BP")]
pgc2850.uni.hyper$pheno <- "Hyperactivity"
pgc2850.uni.inat <- pgc2850.uni[pgc2850.uni$PARENT_INATT.Pval <0.01,c("SNPName","CHR","BP")]
pgc2850.uni.inat$pheno <- "Inattention"

pgc2850.uni.sig <- rbind(pgc2850.uni.wm,pgc2850.uni.aro,pgc2850.uni.hyper,pgc2850.uni.inat)

pgc2850.uni.sig.gr <- makeGRangesFromDataFrame(df = pgc2850.uni.sig,keep.extra.columns = T,ignore.strand = T,seqnames.field = "CHR",start.field = "BP",end.field = "BP")

phenosim.mat <- matrix(c("WM","Arousal","Hyperactivity","Inattention",1,1,1,1),ncol=2)

system.time(test <- runComplexID(pgc2850.uni.sig.gr,phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T))
sum(test$scores$Num.Hits!=0)
ggplot(mapping=aes(x=test$scores$score)) +geom_histogram() +labs(x="Score of Gene",y="Number of Genes",title="Distribution of Scores of Genes using Hits from PGC meta analysis (p-value < 0.05)") + xlim(1e-2, 2)
write.table(test$scores,file = "X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\results\\PGC2850\\univariate\\PGC2850_univariate_pvalue_01.txt",quote = F,sep="\t",row.names = F)

# calculate pearson correlation between endophenotypes and ADHD

endo <- read.csv("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\GWAS_analysis_Wilmot\\Analysis1_EF\\Clinical Outcomes\\Factor Scores\\Factor_Scores_for_PGS_Analyses_n786.csv")
sample.annot <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\GWAS_data\\final_genotypes\\ADHD_GWAS_Sample_Annotation_Cleaned.txt")

m <- merge(endo,sample.annot[,c(1,6)],by=1)
wm.cor <- cor(m$WMM,m$status,use="complete.obs")
aro.cor <- cor(m$AROUSAL,m$status,use="complete.obs")
hyper.cor <- cor(m$PARENT_HYPER,m$status,use="complete.obs")
inatt.cor <- cor(m$PARENT_INATT,m$status,use="complete.obs")

pgrs.uni <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\GWAS_analysis_Wilmot\\Analysis1_EF\\analysis_1_Pleio\\Results\\SNPS_PGRS\\results_PGRS_SNPs_pleio_Univariate_4outcomes_n656_ANNOTATIONS.txt")

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

phenosim.mat <- matrix(c("WM","Arousal","Hyperactivity","Inattention",abs(wm.cor),abs(aro.cor),abs(hyper.cor),abs(inatt.cor)),ncol=2)

system.time(test <- runComplexID(pgrs.uni.sig.gr,phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T))
test$scores[1:10,]

write.table(test$scores,file = "X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\results\\PGRS\\univariate\\PGRS_univariate_pvalue_0001_phenoSim_raw_correlations.txt",quote = F,sep="\t",row.names = F)

# Performing analysis on PGC's top 12 regions, also include genes from gtex, bio_eqtls, hic_gz, and hi_cp:

top12.pgc <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\data\\Supplementary_eTable_2.txt",stringsAsFactors = F)
top12.pgc$start <- as.integer(sapply(strsplit(top12.pgc$bp,"-"), "[[", 1))
top12.pgc$end <- as.integer(sapply(strsplit(top12.pgc$bp,"-"), "[[", 2))
top12.pgc <- top12.pgc[,c(1,2,13,14)]
additional.regions <- as.data.frame(matrix(c("MED8",1,43849588,43855479,
                               "CCDC24",1,44457031,44462200,
                               "ATP6V0B",1,44440159,44443967,
                               "HYI",1,43916824,43919660,
                               "HYI-AS1",1,43919598,43922666,
                               "SZT2",1,43855553,43918321,
                               "SZT2-AS1",1,43913447,43914315,
                               "B4GALT2",1,44444615,44456840,
                               "ATP6V1E1P1",1,43368903,43369583,
                               "EBNA1BP2",1,43629846,43736607,
                               "WDR65",1,43637820,43720029 ,
                               "RP5-898J17.1",1,96719625,96839681,
                               "HPR",16,72088522,72111145,
                               "HP",16,72088491,72094954,
                               "HCCAT5",16,73126248,73127673,
                               "AC107622.1",3,20383960,20392420,
                               "ZNF385D",3,21459915,22414812,
                               "TMEM161B",5,87485450,87565293),
                             byrow = T,ncol=4))
names(additional.regions) <- names(top12.pgc)
top12.pgc <- rbind(top12.pgc,additional.regions)

top12.pgc.gr <- makeGRangesFromDataFrame(df = top12.pgc,keep.extra.columns = T,ignore.strand = T,seqnames.field = "chr",start.field = "start",end.field = "end")
top12.pgc.gr$pheno <- "ADHD"
phenosim.mat <- matrix(c("ADHD",1),ncol=2)
system.time(test <- runComplexID(top12.pgc.gr,phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T))
test$scores[1:10,]
write.table(test$scores,file = "X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\results\\PGC_12_GenomeWide_SNPs\\PGC_12_GenomineWide_SNPs_With_Regions_EQTLS.txt.",quote = F,sep="\t",row.names = F)


# Perform analysis using only 12 SNPs:
top12snps <- matrix(c("rs11420276",1,44184192,44184193,
                    "rs11591402",10,106747354,106747354,
                    "rs1222063",1,96602440,96602440,
                    "rs1427829",12,89760744,89760744,
                    "rs212178",16,72578131,72578131,
                    "rs281324",15,47754018,47754018,
                    "rs28411770",4,31151456,31151456,
                    "rs4858241",3,20669071,20669071,
                    "rs4916723",5,87854395,87854395,
                    "rs5886709",7,114086134,114086135,
                    "rs74760947",8,34352610,34352610,
                    "rs9677504",2,215181889 ,215181889 ),
                    byrow = T,ncol=4)

top12snps <- cbind(top12snps,"ADHD")


phenosim.mat <- matrix(c("ADHD",1),ncol=2)
system.time(test <- runComplexID(top12snps,phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T))
write.table(test$scores,file = "X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\results\\PGC_12_GenomeWide_SNPs\\PGC_12_GenomineWide_SNPs.txt",quote = F,sep="\t",row.names = F)





# see if there are more hugo or entrez gene names:
genenamesmapping <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\ensembl_gene_to_protein\\Gene_names_8_18_17.gz")
length(unique(genenamesmapping$Gene.stable.ID))
length(unique(genenamesmapping$Gene.stable.ID[genenamesmapping$HGNC.symbol==""]))
length(unique(genenamesmapping$Gene.stable.ID[is.na(genenamesmapping$EntrezGene.ID)]))
sum(genenamesmapping$HGNC.symbol=="" & !is.na(genenamesmapping$EntrezGene.ID))
corum.annot$subunits.Gene.name.[1:5]

# see if the protein IDs in the STRING PPI have more HUGO or Entrez gene names:
string.ppi.annot.hugo <- merge(string.ppi.filt,genenamesmapping,by.x="protein1",by.y="Protein.stable.ID",all.x=T)
string.ppi.annot.hugo <- merge(string.ppi.annot.hugo,genenamesmapping,by.x="protein2",by.y="Protein.stable.ID",all.x=T)
string.ppi.annot.hugo <- string.ppi.annot.hugo[!is.na(string.ppi.annot.hugo$HGNC.symbol.x) & !is.na(string.ppi.annot.hugo$HGNC.symbol.y),]
string.ppi.annot.hugo <- unique(string.ppi.annot.hugo[,c(6,10)])
string.ppi.annot.hugo <- unique(rbind(string.ppi.annot.hugo,data.frame(HGNC.symbol.x=string.ppi.annot.hugo$HGNC.symbol.y,
                                                                HGNC.symbol.y=string.ppi.annot.hugo$HGNC.symbol.x)))
dim(string.ppi.annot.hugo)
dim(string.ppi.entrez)

# test loops:
pgrs.uni <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\GWAS_analysis_Wilmot\\Analysis1_EF\\analysis_1_Pleio\\Archive\\Results\\SNPS_PGRS\\new\\results_PGRS_SNPs_pleio_Univariate_4outcomes_n656_ANNOTATIONS.txt")

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

# test loops using 12 genome-wide PGC SNPs, but don't include the eqtl and hic loops:
top12.pgc <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\data\\Supplementary_eTable_2.txt",stringsAsFactors = F)
top12.pgc$start <- as.integer(sapply(strsplit(top12.pgc$bp,"-"), "[[", 1))
top12.pgc$end <- as.integer(sapply(strsplit(top12.pgc$bp,"-"), "[[", 2))
top12.pgc <- top12.pgc[,c(1,2,13,14)]

top12.pgc.gr <- makeGRangesFromDataFrame(df = top12.pgc,keep.extra.columns = T,ignore.strand = T,seqnames.field = "chr",start.field = "start",end.field = "end")
top12.pgc.gr$pheno <- "ADHD"
phenosim.mat <- matrix(c("ADHD",1),ncol=2)
system.time(test <- runComplexID(top12.pgc.gr,phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T))
test$scores[1:10,]
generatePlot(centralGenes = as.character(test$scores$Entrez.Gene.ID[1:10]),colorGenes = as.character(test$scores$Entrez.Gene.ID[1:10][test$scores$Num.Hits[1:10]>0]))
generatePlot(centralGenes = as.character(test$scores$Entrez.Gene.ID[1:50]),colorGenes = as.character(test$scores$Entrez.Gene.ID[1:50][test$scores$Num.Hits[1:50]>0]),vertex.size=4,vertex.label.cex=.75)
generatePlot(centralGenes = as.character(test$scores$Entrez.Gene.ID[1:100]),colorGenes = as.character(test$scores$Entrez.Gene.ID[1:100][test$scores$Num.Hits[1:100]>0]),vertex.size=4,vertex.label.cex=.75)

generatePlot(centralGenes = as.character(test$scores$Entrez.Gene.ID[1:10]),order=1,vertex.size=4,useHugoNames=F,vertex.label=NA)

# Test annotation function:
# Test for top 12 SNP regions:

top12.pgc <- read.delim("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\data\\Supplementary_eTable_2.txt",stringsAsFactors = F)
top12.pgc$start <- as.integer(sapply(strsplit(top12.pgc$bp,"-"), "[[", 1))
top12.pgc$end <- as.integer(sapply(strsplit(top12.pgc$bp,"-"), "[[", 2))
top12.pgc <- top12.pgc[,c(1,2,13,14)]

top12.pgc.gr <- makeGRangesFromDataFrame(df = top12.pgc,keep.extra.columns = T,ignore.strand = T,seqnames.field = "chr",start.field = "start",end.field = "end")
top12.pgc.gr$pheno <- "ADHD"

test <- annotateHits(top12.pgc.gr,promoterRange = 10000,upstream = 0,downstream = 0,promoters = F,gene.body = F,utr = F,eqtl = T,enhancers = F)

# Test annotation function for previously annotated loops:

test <- annotateHits(top12.pgc.gr,promoterRange = 20000,upstream = 0,downstream = 0,promoters = F,gene.body = F,utr = F,eqtl = F,enhancers = T,loopDist = 100000)


pgc.top <- read.csv("X:\\OHSU Shared\\Restricted\\OCTRI\\Informatics\\TransBio\\ADHD\\SNP_Networks\\loops\\mike_pgc_annotation\\adhd_pthresh_1e-5_genes_100kb_ldist_2kb_pdist_updated.csv")
pgc.top <- pgc.top[,c(1,12,13,13)]
pgc.top.gr <- makeGRangesFromDataFrame(df = pgc.top,keep.extra.columns = T,ignore.strand = T,seqnames.field = "chrom",start.field = "bp",end.field = "bp.1")
pgc.top.gr$pheno <- "ADHD"

test <- annotateHits(pgc.top.gr,promoterRange = 20000,upstream = 0,downstream = 0,promoters = F,gene.body = F,utr = F,eqtl = F,enhancers = T,loopDist = 100000)



test <- annotateHits(a,promoterRange = 2000,upstream = 0,downstream = 0,promoters = F,gene.body = F,utr = F,eqtl = F,enhancers = T,loopDist = 1)

a <- findOverlaps(pgc.top.gr[681],ComplexID:::.all.enhancers.gr)
aa <- findOverlaps(ComplexID:::.loops.gr,ComplexID:::.all.enhancers.gr[subjectHits(a)],maxgap = 100000)
aaa <- GRanges(seqnames=as.character(seqnames(ComplexID:::.loops.gr)[queryHits(aa)]),
                        ranges=IRanges(start=as.integer(ComplexID:::.loops.gr$y1)[queryHits(aa)]-100000,
                                       end=as.integer(ComplexID:::.loops.gr$y2)[queryHits(aa)]+100000))
aaa$enhancer.start <- start(ComplexID:::.all.enhancers.gr[subjectHits(a)])[subjectHits(aa)]
aaa$enhancer.end <- end(ComplexID:::.all.enhancers.gr[subjectHits(a)])[subjectHits(aa)]


.gene.annotation.gr <- ComplexID:::.gene.annotation.gr
.encode.promoters.distal.gr <- ComplexID:::.encode.promoters.distal.gr
.encode.promoters.prox.gr <- ComplexID:::.encode.promoters.prox.gr






tss.regions.gr <- GRanges(seqnames=seqnames(.gene.annotation.gr),
                          ranges=IRanges(start=ifelse(strand(.gene.annotation.gr)=="-",.gene.annotation.gr$Transcription.start.site..TSS.,.gene.annotation.gr$Transcription.start.site..TSS.-20000),
                                         end=ifelse(strand(.gene.annotation.gr)=="-",.gene.annotation.gr$Transcription.start.site..TSS.+20000,.gene.annotation.gr$Transcription.start.site..TSS.)),
                          strand=strand(.gene.annotation.gr),mcols(.gene.annotation.gr))

promoter.distal.tss.hits <- findOverlaps(.encode.promoters.distal.gr,tss.regions.gr)
split.hits.idx <- split(1:length(promoter.distal.tss.hits),queryHits(promoter.distal.tss.hits),drop=T)
.encode.promoters.distal.gr$genes <- ""
.encode.promoters.distal.gr$Feature <- "Distal_Promoter"
.encode.promoters.distal.gr$genes[as.integer(names(split.hits.idx))] <- sapply(split.hits.idx, function(x) {
  return(tss.regions.gr$genes[subjectHits(promoter.distal.tss.hits)[x]])
})

promoter.prox.tss.hits <- findOverlaps(.encode.promoters.prox.gr,tss.regions.gr[tss.regions.gr$genes=="CRIM1"])
split.hits.idx <- split(1:length(promoter.distal.tss.hits),queryHits(promoter.prox.tss.hits),drop=T)
.encode.promoters.prox.gr$genes <- ""
.encode.promoters.prox.gr$Feature <- "Distal_Promoter"
.encode.promoters.prox.gr$genes[as.integer(names(split.hits.idx))] <- sapply(split.hits.idx, function(x) {
  return(tss.regions.gr$genes[subjectHits(promoter.prox.tss.hits)[x]])
})

"CRIM1" %in% unlist(.encode.promoters.distal.gr$genes)
"CRIM1" %in% unlist(.encode.promoters.prox.gr$genes)
p <- .encode.promoters.prox.gr[sapply(.encode.promoters.prox.gr$genes, function(x) "CRIM1" %in% x)]
promoters.with.genes <- c(.encode.promoters.distal.gr[.encode.promoters.distal.gr$genes != ""],.encode.promoters.prox.gr[.encode.promoters.prox.gr$genes != ""])
"CRIM1" %in% unlist(promoters.with.genes$genes)



system.time(test <- runComplexID(pgrs.uni.sig.gr,phenoSim=phenosim.mat,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T))
test$scores[1:20,]



# Figure out why large annotations take a long time ~600,000 hits

hits.big <- rep(hits,4000)
system.time(test1 <- annotateHits(hits.big[1:100000],promoterRange=0,gene.body = F,promoters = T,utr=F,eqtl=F,enhancers = F))
system.time(test1 <- annotateHits(hits.big[1:200000],promoterRange=0,gene.body = F,promoters = T,utr=F,eqtl=F,enhancers = F))
system.time(test1 <- annotateHits(hits.big,promoterRange=0,gene.body = F,promoters = T,utr=F,eqtl=F,enhancers = F))

