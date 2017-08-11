#'  ComplexID: Integrated Hit Annotations and Random Walk of Complex Network
#'
#'  A package with a self-contained annotation database and protein complex network for finding causitive genes for a phenotype of interest.
#'
#' @section ComplexID functions:
#' runComplexID
#'
#' @docType package
#' @name ComplexID
#' @author Peter Ryabinin
#' @import Matrix
#' @import IRanges
#' @import S4Vectors
#' @import GenomicRanges
NULL

#' Annotates Hits, Performs Random Walk, and Scores Genes
#'
#' Annotates hits to genes, performs random walk with restarts on a network of protein complexes, and then scores each gene in the network for its association with the pheontype of interest
#'
#' @param Hits Granges object with two meta data columns, or a matrix or data frame with at least 4 columns. \cr If it is a Granges object, then the first meta data column is the site's name. The second meta data columns is a phenotype that the site is causitive for.\cr
#' If it is a matrix or data frame, then the first column must be the Hit's name, the second column must be chromosome designation, the third column must the base pair position, and the fourth column must a phenotype that the site is causitive for. \cr
#' For both Grange objects and matrices/dataframes, each entry/row corresponds to one site that is causitive to one phenotype. If a site is causitive in multiple phenotypes then there would be multiple entries for the same site but all with different values in the phenotype column
#' @param phenoSim matrix or data frame with two columns. The first column are names of phenotypes that match the same phenotypes found in Hits. The second column are phenotype similarity values between the phenotype in that row and the phenotype of interest (values between 0 and 1), with higher values denoting higher similarity
#' @param promoterRange single integer greater than zero. How many bases to look upstream of a TSS of a gene in order to find a promoter region.
#' @param eps single numeric, must be greater than zero. L1 norm threshold between current and previous interations of random walk at which to terminate the random walk
#' @param alpha single numeric in the range of (0,1]. The weight given to the vector of initialized values for the random walk, higher value of alpha means more weight for the initialized values
#' @details
#' Annotates Hits to genes using a built-in annotation database. Gene annotations come from ENSEMBL genes that have Entrez gene IDs and are in the STRING PPI with threshold >700. Promoter regions are
#' from ENCODE annotation, a hit in Hits is in a promoter region for a gene if it lies within a promoter region that is a number of bases upstream equal to promoterRange.\cr
#' \cr
#' After the causitive genes for each endophenotype are identified, it performs a Random Walk with Restarts on a pre-constructed protein complex network as in the RWPCN method.
#' The protein complex network was constructed in a similar way is in the RWPCN method. For a PPI we used STRING with a threshold cutoff of 700. Protein IDs in STRING were then mapped to Entrez gene ids.
#' Protein complexes were retrieved from CORUM. Any complex with no genes in the PPI was removed along with 5 of the largest complexes (more than 70 subunits) \cr
#' \cr
#' A random walk with restarts is initialized and performed as in RWPCN then all genes in the PPI and complexes are scored according to the weights in the complex network.
#' @return data frame of two columns. The first column are Entrez Gene IDs and the second column are the scores of each gene for the phenotype of interest. Sorted from largest to smallest score.
#' @examples
#' data("hits")
#' data("hits.pheno")
#' test <- runComplexID(Hits = hits,phenoSim=hits.pheno,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T)
#' @export
runComplexID <- function(Hits,phenoSim,promoterRange=100000,eps=1e-10,alpha=0.8,upstream=NULL,downstream=NULL,utr=T) {
  # Check for errors in input
  if (promoterRange < 0)
    stop("promoterRange must be greater than zero")
  if (eps <= 0)
    stop("eps must be greater than zero")
  if (alpha <= 0 | alpha > 1)
    stop("alpha must be in the range of (0,1]")
  if (length(upstream)>0 & upstream < 0)
    stop("upstream must be greater than or equal to zero")
  if (length(downstream)>0 & downstream < 0)
    stop("downstream must be greater than or equal to zero")
  if (utr != T & utr != F)
    stop("utr must be either TRUE or FALSE")
  if (ncol(mcols(Hits)) < 2)
    stop("Hits must have at least two meta data columns")
  if (ncol(phenoSim) < 2)
    stop("phenoSim must have at least two columns")
  # Convert matrix to genomic ranges object
  if (is.matrix(Hits) | is.data.frame(Hits)) {
    Hits <- as.data.frame(Hits)
    Hits <- makeGRangesFromDataFrame(df = Hits,keep.extra.columns = T,ignore.strand = T,seqnames.field = names(Hits)[2],start.field = names(Hits)[3],end.field = names(Hits)[3])
  }
  if (sum(mcols(Hits)[,2] %in% phenoSim[,1]) < length(Hits))
    stop("some hits in Hits have a phenotype that is nonexistant in phenoSim")

  # create annotations according to the user's promoter threshold
  annotations <- .createAnnotationDB(promoterRange,upstream,downstream,utr)
  # get set of seed genes linked to a phenotype
  seedgenes <- .getSeedGenes(Hits,annotations)
  # initialize network
  F0 <- .initNetwork(seedgenes$seedgenes,phenoSim)
  # perform random walk with restarts
  Ffinal <- .RWPCN(F0,eps,alpha)
  # calculate final scores for each gene
  geneScores <- .calcGeneScores(Ffinal)
  # Determine if a gene has a hit in it
  geneScores <- merge(geneScores,seedgenes$hitsPerGene,by=1,all.x=T)
  names(geneScores)[ncol(geneScores)] <- "Num.Hits"
  geneScores$Num.Hits <- ifelse(is.na(geneScores$Num.Hits),0,geneScores$Num.Hits)
  geneScores <- geneScores[order(geneScores$score,decreasing = T),]
  row.names(geneScores) <-NULL
  return(list("scores"=geneScores,"missingHits"=seedgenes$missingHits))
}

#' @keywords internal
.createAnnotationDB <- function(promoterRange,upstream,downstream,utr) {
  tss.regions.gr <- GRanges(seqnames=seqnames(.gene.annotation.gr),
                            ranges=IRanges(start=ifelse(strand(.gene.annotation.gr)=="-",.gene.annotation.gr$Transcription.Start.Site..TSS.,.gene.annotation.gr$Transcription.Start.Site..TSS.-promoterRange),
                                           end=ifelse(strand(.gene.annotation.gr)=="-",.gene.annotation.gr$Transcription.Start.Site..TSS.+promoterRange,.gene.annotation.gr$Transcription.Start.Site..TSS.)),
                            strand=strand(.gene.annotation.gr),mcols(.gene.annotation.gr))

  promoter.distal.tss.hits <- findOverlaps(.encode.promoters.distal.gr,tss.regions.gr)
  split.hits.idx <- split(1:length(promoter.distal.tss.hits),queryHits(promoter.distal.tss.hits),drop=T)
  .encode.promoters.distal.gr$genes <- ""
  .encode.promoters.distal.gr$Feature <- "Distal_Promoter"
  .encode.promoters.distal.gr$genes[as.integer(names(split.hits.idx))] <- sapply(split.hits.idx, function(x) {
    return(tss.regions.gr$genes[subjectHits(promoter.distal.tss.hits)[x]])
  })

  promoter.prox.tss.hits <- findOverlaps(.encode.promoters.prox.gr,tss.regions.gr)
  split.hits.idx <- split(1:length(promoter.distal.tss.hits),queryHits(promoter.distal.tss.hits),drop=T)
  .encode.promoters.prox.gr$genes <- ""
  .encode.promoters.prox.gr$Feature <- "Distal_Promoter"
  .encode.promoters.prox.gr$genes[as.integer(names(split.hits.idx))] <- sapply(split.hits.idx, function(x) {
    return(tss.regions.gr$genes[subjectHits(promoter.prox.tss.hits)[x]])
  })

  temp.gr <- .gene.annotation.gr
  temp.gr$Transcription.Start.Site..TSS. <- NULL

  ret <- c(temp.gr,.encode.promoters.distal.gr[.encode.promoters.distal.gr$genes != ""],.encode.promoters.prox.gr[.encode.promoters.prox.gr$genes != ""])

  if (utr)
    ret <- c(ret,.utr.entrez.gr)

  if (length(upstream)>0) {
    upstream.gr <- GRanges(seqnames=seqnames(.gene.annotation.gr),
                              ranges=IRanges(start=ifelse(strand(.gene.annotation.gr)=="-",.gene.annotation.gr$Transcription.Start.Site..TSS.,.gene.annotation.gr$Transcription.Start.Site..TSS.-upstream),
                                             end=ifelse(strand(.gene.annotation.gr)=="-",.gene.annotation.gr$Transcription.Start.Site..TSS.+upstream,.gene.annotation.gr$Transcription.Start.Site..TSS.)),
                              strand=strand(.gene.annotation.gr),genes=mcols(.gene.annotation.gr)[,"genes"])
    upstream.gr$Feature = "Upstream"
    ret <- c(ret,upstream.gr)
  }

  if (length(downstream)>0) {
    downstream.gr <- GRanges(seqnames=seqnames(.gene.annotation.gr),
                           ranges=IRanges(start=ifelse(strand(.gene.annotation.gr)=="-",end(.gene.annotation.gr)-downstream,end(.gene.annotation.gr)),
                                          end=ifelse(strand(.gene.annotation.gr)=="-",end(.gene.annotation.gr),end(.gene.annotation.gr)+downstream)),
                           strand=strand(.gene.annotation.gr),genes=mcols(.gene.annotation.gr)[,"genes"])
    downstream.gr$Feature = "Downstream"
    ret <- c(ret,downstream.gr)
  }
  return(ret)
}

#' @keywords internal
.getSeedGenes <- function(Hits,annotations) {
  snpOverlaps <- findOverlaps(Hits,annotations,ignore.strand=T)
  missingHits <- Hits[!(1:length(Hits) %in% queryHits(snpOverlaps))]
  pheno <- mcols(Hits)[queryHits(snpOverlaps),2]
  time.to.repeat <- sapply(annotations$genes[subjectHits(snpOverlaps)],length)
  snpPheno <- unlist(mapply(rep, pheno, time.to.repeat))

  geneToPheno <- unique(matrix(c(unlist(annotations$genes[subjectHits(snpOverlaps)]),snpPheno),ncol = 2))

  snpToGene <- unique(matrix(c(unlist(mapply(rep, queryHits(snpOverlaps), time.to.repeat)),unlist(annotations$genes[subjectHits(snpOverlaps)])),ncol=2))
  hitsPerGene <- as.data.frame(table(snpToGene[,2]),stringsAsFactors=F)

  return(list("seedgenes"=unique(geneToPheno),"missingHits"=missingHits,"hitsPerGene"=hitsPerGene))
}

#' @keywords internal
.initNetwork <- function(seedgenes,phenoSim) {
  geneToPheno <- split(seedgenes[,2],seedgenes[,1],drop=T)
  seedScore <- sapply(geneToPheno, function(x) {
    sum(as.numeric(phenoSim[as.character(phenoSim[,1]) %in% as.character(x),2]))
  })
  names(seedScore) <- names(geneToPheno)

  corum.f0 <- sapply(.corum.subunits, function(x) sum(seedScore[as.character(x)],na.rm=T))
  corum.f0 <- corum.f0 * .density
  indiv.f0 <- sapply(rownames(.W_norm_t_spars)[(length(.corum.subunits)+1):nrow(.W_norm_t_spars)],function(x) sum(seedScore[as.character(x)],na.rm=T))

  return(matrix(c(corum.f0,indiv.f0)))
}

#' @keywords internal
.RWPCN <- function(F0,eps,alpha) {
  F.diff <- eps+1
  FOld <- F0
  while(F.diff>eps) {
    Fnew <- (1-alpha) * .W_norm_t_spars %*% FOld + alpha * F0
    F.diff <- sum(abs(Fnew-FOld))
    FOld <- Fnew
  }
  return(Fnew)
}

#' @keywords internal
.calcGeneScores <- function(Ffinal) {
  scores <- sapply(.geneToComplex, function(x) {
    sum(Ffinal[x])
  })
  ret <- data.frame("Entrez.Gene.ID"=names(.geneToComplex),"HUGO Gene Name"=.hugoNames,"Complexes"=.complexNames,"score"=scores)
  ret <- ret[order(ret[,2],decreasing = T),]
  return(ret)
}
