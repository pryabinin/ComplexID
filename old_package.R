#'  ComplexID: Integrated Hit Annotations and Random Walk of Complex Network
#'
#'  A package with a self-contained annotation database and protein complex network for finding causitive genes for a phenotype of interest.
#'
#' @section ComplexID functions:
#' runComplexID
#' generatePlot
#' annotateHits
#'
#'
#' @docType package
#' @name ComplexID
#' @author Peter Ryabinin
#' @import Matrix
#' @import IRanges
#' @import S4Vectors
#' @import GenomicRanges
#' @import igraph
NULL

#' Annotates Hits, Performs Random Walk, and Scores Genes
#'
#' Annotates hits to genes, performs random walk with restarts on a network of protein complexes, and then scores each gene in the network for its association with the pheontype of interest
#'
#' @param Hits Granges object with two meta data columns, or a matrix or data frame with at least 4 columns. \cr If it is a Granges object, then the first meta data column is the site's name. The second meta data columns is a phenotype that the site is associated with.\cr
#' If it is a matrix or data frame, then the first column must be the Hit's name, the second column must be chromosome designation, the third column must the base pair position, and the fourth column must a phenotype that the site is causitive for. \cr
#' For both Grange objects and matrices/dataframes, each entry/row corresponds to one site that is causitive to one phenotype. If a site is causitive in multiple phenotypes then there would be multiple entries for the same site but all with different values in the phenotype column
#' @param phenoSim matrix or data frame with two columns. The first column are names of phenotypes that match the same phenotypes found in Hits. The second column are phenotype similarity values between the phenotype in that row and the phenotype of interest (values between 0 and 1), with higher values denoting higher similarity
#' @param promoterRange single integer greater than or equal to zero. How many bases to look upstream of a TSS of a gene in order to find a promoter region for a gene.
#' @param eps single numeric, must be greater than zero. L1 norm threshold between current and previous interations of random walk at which to terminate the random walk
#' @param alpha single numeric in the range of (0,1]. The weight given to the vector of initialized values for the random walk, higher value of alpha means more weight for the initialized values
#' @param upstream single integer or NULL. How far upstream of a transcription start site a hit can be for it to be annotated to that gene. A NULL value is equivalent to a value of zero (no upstream sites will be annotated to a gene unless they lie in a promoter region, see promoterRange parameter).
#' @param downstream single integer or NULL. How far downstream of a transcription start site a hit can be for it to be annotated to that gene. A NULL value is equivalent to a value of zero (no downstream sites will be annotated to a gene).
#' @param utr TRUE or FALSE. If TRUE then it will look for hits in the 3' and 5' UTRs of genes, otherwise it will not.
#' @param eqtl TRUE or FALSE. By default TRUE. If TRUE, then hits may be mapped to eQTL loci, and therefore genes effected by those eQTLs be designated as causitive
#' @param enhancers TRUE or FALSE. By default TRUE. If TRUE, then hits may be mapped to enhancer loci and linked to genes via looping structures and promoters
#' @param loopDist single integer. By default 0. The maximum allowable distance that an enhancer or promoter can be from a looping region to be annotated to it.
#' @details
#' Annotates Hits to genes using a built-in annotation database. Gene annotations come from ENSEMBL genes that have Entrez gene IDs and are in the STRING PPI with threshold >700. Promoter regions are
#' from ENCODE annotation, a hit in Hits is in a promoter region for a gene if it lies within a promoter region that is a number of bases upstream equal to promoterRange.\cr
#' \cr
#' After the causitive genes for each endophenotype are identified, it performs a Random Walk with Restarts on a pre-constructed protein complex network as in the RWPCN method.
#' The protein complex network was constructed in a similar way is in the RWPCN method. For a PPI we used STRING with a threshold cutoff of 700. Protein IDs in STRING were then mapped to Entrez gene ids.
#' Protein complexes were retrieved from CORUM. Any complex with no genes in the PPI was removed along with 5 of the largest complexes (more than 70 subunits) \cr
#' \cr
#' A random walk with restarts is initialized and performed as in RWPCN then all genes in the PPI and complexes are scored according to the weights in the complex network.
#' @return A list with two objects: a data frame called "scores" and a GRanges object "missingHits"
#' The data frame "scores" has six columns showing the scores of each gene, related to how much that gene is important to the query phenotype, as well as other information about the gene. It is ordered with the highest scoring genes first.\cr
#' The first column are Entrez Gene IDs, the second column are HUGO gene names, the third column are the names of the complexes that gene is part of, the fourth columns is the score for the gene, the fifth column are the features of that gene that have a hit in them, and the sixth column is the number of hits that were annotated to that gene.
#' The GRanges object "missingHits" lists all of the input hits that were not mapped to any gene.
#' @examples
#' data("hits")
#' data("hits.pheno")
#' test <- runComplexID(Hits = hits,phenoSim=hits.pheno,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T)
#' @export
runComplexID <- function(Hits,phenoSim,promoterRange=100000,eps=1e-10,alpha=0.8,upstream=0,downstream=0,gene.body=T,promoters=T,utr=T,eqtl=T,enhancers=T,loopDist=0,geneScoring=sum) {
  # Check for errors in input
  if (promoterRange < 0)
    stop("promoterRange must be greater than zero")
  if (eps <= 0)
    stop("eps must be greater than zero")
  if (alpha <= 0 | alpha > 1)
    stop("alpha must be in the range of (0,1]")
  if (length(upstream)>0)
    if (upstream<0)
      stop("upstream must be greater than or equal to zero")
  if (length(downstream)>0)
    if(downstream < 0)
      stop("downstream must be greater than or equal to zero")
  if (utr != T & utr != F)
    stop("utr must be either TRUE or FALSE")
  if (ncol(phenoSim) < 2)
    stop("phenoSim must have at least two columns")
  # Convert matrix to genomic ranges object
  if (is.matrix(Hits) | is.data.frame(Hits)) {
    if (ncol(Hits)<5)
      stop("Hits must have at least 5 columns")
    Hits <- as.data.frame(Hits)
    Hits <- makeGRangesFromDataFrame(df = Hits,keep.extra.columns = T,ignore.strand = T,seqnames.field = names(Hits)[2],start.field = names(Hits)[3],end.field = names(Hits)[4])
  }
  if (ncol(mcols(Hits)) < 2)
    stop("Hits must have at least two meta data columns")
  if (sum(mcols(Hits)[,2] %in% phenoSim[,1]) < length(Hits))
    stop("some hits in Hits have a phenotype that is nonexistant in phenoSim")

  # create annotations according to the user's promoter threshold
  annotations <- .createAnnotationDB(promoterRange,upstream,downstream,gene.body,promoters,utr,eqtl,enhancers,loopDist)
  # get set of seed genes linked to a phenotype
  seedgenes <- .getSeedGenes(Hits,annotations)
  # initialize network
  F0 <- .initNetwork(seedgenes$seedgenes,phenoSim)
  # perform random walk with restarts
  Ffinal <- .RWPCN(F0[[1]],eps,alpha)
  # calculate final scores for each gene
  geneScores <- .calcGeneScores(Ffinal,geneScoring)
  # add scores for genes not in PPI or CORUM
  non.ppi.df <- data.frame("Entrez.Gene.ID"=.non.ppi.genes,
                           "HUGO Gene Name"=names(.non.ppi.genes),
                           "Complexes"="None",
                           "score"=F0[[2]],
                           "Gene.in.Network"="No")
  geneScores <- rbind(geneScores,non.ppi.df)
  # Determine if a gene has a hit in it
  geneScores <- merge(geneScores,seedgenes$hitsPerGene,by=1,all.x=T)
  names(geneScores)[(ncol(geneScores)-1):ncol(geneScores)] <- c("Feature.of.Hits","Num.Hits")
  geneScores$Num.Hits <- ifelse(is.na(geneScores$Num.Hits),0,geneScores$Num.Hits)
  geneScores <- geneScores[order(geneScores$score,decreasing = T),]
  row.names(geneScores) <-NULL
  return(list("scores"=geneScores,"missingHits"=seedgenes$missingHits))
}

#' Produce Network Plot
#'
#' Produces a network plot of the specified genes and their neighbors (if desired)
#'
#' @param centralGenes character vector. The Entrez Gene IDs of the genes that around which the network plot will be centered
#' @param order single integer, by default 0, The degree of neighboring genes that will be included. A value of zero means no neighboring genes are included
#' @param useHugoNames single binary TRUE or FALSE, defaults to TRUE which means the vertices will be labeled by their Hugo names. If FALSE then the vertcies will be labeled by the their EntrezGene IDs. If no labels are desired then set useHugoNames to FALSE and pass "vertex.label=NA" in the ... parameter.
#' @param colorGenes character vector. By default, is equal to centralGenes. The Entrez Gene IDs of the genes that will be colored red. The rest of the genes will be colored orange.
#' @param ... arguments passed to plot.igraph function
#' @details
#' Plots a graph of the network of genes centered around the centralGenes, including neighbors out to the degree of "order". Any additional parameters for the plot.igraph function may be passed as well. \cr\cr
#' Vertices are colored red if they are in centralGenes and orange otherwise. \cr
#' \cr
#' Edges are colored red if the interaction exists in the PPI and the two vertices share a complex. Edges are colored black if the interaction exists in the PPI but the two vertices do not share a compelx. Edges are colored green if the two vertices share a complex but do not interact in the PPI. \cr\cr
#' Complexes are plotted as ellipses circling the vertices. Vertices that are not circled are not part of any complex. Vertices which are circled individually are part of a complex but none of the other genes are included in this graph.
#' @return A plot of the subnetwork
#' @examples
#' data("hits")
#' data("hits.pheno")
#' test <- runComplexID(Hits = hits,phenoSim=hits.pheno,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T)
#' generatePlot(test$scores$Entrez.Gene.ID[1:10])
#' @export
generatePlot <- function(centralGenes,order=0,useHugoNames=T,colorGenes=centralGenes,...) {
  graph.to.plot <- induced_subgraph(.complete.igraph,unlist(neighborhood(.complete.igraph,order = order,nodes = as.character(centralGenes))))
  V(graph.to.plot)$color <- ifelse(names(V(graph.to.plot)) %in% as.character(colorGenes),"red","orange")
  complexes.to.plot <- .corum.subunits[sapply(.corum.subunits,function(x) { sum(as.integer(names(V(graph.to.plot))) %in% x)>0 })]
  complexes.to.plot <- lapply(complexes.to.plot,function(x) { as.character(x)[as.character(x) %in% names(V(graph.to.plot))] })
  if (useHugoNames)
    return(plot.igraph(x=graph.to.plot,mark.groups = complexes.to.plot,vertex.label=.ent.to.hug[names(V(graph.to.plot))],...))
  else
    return(plot.igraph(x=graph.to.plot,mark.groups = complexes.to.plot,...))
}

#' Annotates Hits by Genomic Features
#' @inheritParams runComplexID
#' @export
annotateHits <- function(Hits,promoterRange=100000,upstream=0,downstream=0,gene.body=T,promoters=T,utr=T,eqtl=T,enhancers=T,loopDist=0) {
  # Check for errors in input
  if (promoterRange < 0)
    stop("promoterRange must be greater than zero")
  if (length(upstream)>0)
    if (upstream<0)
      stop("upstream must be greater than or equal to zero")
  if (length(downstream)>0)
    if(downstream < 0)
      stop("downstream must be greater than or equal to zero")
  if (utr != T & utr != F)
    stop("utr must be either TRUE or FALSE")
  # Convert matrix to genomic ranges object
  if (is.matrix(Hits) | is.data.frame(Hits)) {
    if (ncol(Hits)<5)
      stop("Hits must have at least 5 columns")
    Hits <- as.data.frame(Hits)
    Hits <- makeGRangesFromDataFrame(df = Hits,keep.extra.columns = T,ignore.strand = T,seqnames.field = names(Hits)[2],start.field = names(Hits)[3],end.field = names(Hits)[4])
  }
  if (ncol(mcols(Hits)) < 2)
    stop("Hits must have at least two meta data columns")

  # create annotations according to the user's promoter threshold
  annotations <- .createAnnotationDB(promoterRange,upstream,downstream,gene.body,promoters,utr,eqtl,enhancers,loopDist)

  # annotate input hits
  mcols(Hits) <- cbind(mcols(Hits),data.frame(1:length(Hits)))
  snpOverlaps <- findOverlaps(Hits,annotations,ignore.strand=T)
  if (length(snpOverlaps)<1)
    stop("No input hits were mapped to any gene using these annotation parameters")
  missingHits <- Hits[!(1:length(Hits) %in% queryHits(snpOverlaps))]
  time.to.repeat <- sapply(annotations$genes[subjectHits(snpOverlaps)],length)

  entrezToHugo <- as.character(.hugoNames)
  names(entrezToHugo) <- as.character(names(.geneToComplex))
  out.df <- data.frame(snpName=as.character(unlist(mapply(rep,mcols(Hits)[queryHits(snpOverlaps),1],time.to.repeat))),
                       entrez.genes=unlist(annotations$genes[subjectHits(snpOverlaps)]),
                       hugo.names=entrezToHugo[as.character(unlist(annotations$genes[subjectHits(snpOverlaps)]))],
                       features=unlist(mapply(rep,annotations$Feature[subjectHits(snpOverlaps)],time.to.repeat)),
                       order=unlist(mapply(rep,mcols(Hits)[queryHits(snpOverlaps),ncol(mcols(Hits))],time.to.repeat)),
                       stringsAsFactors = F)

  out.df <- aggregate(cbind(entrez.genes,hugo.names,features,snpName)~order, data = unique(out.df), paste, collapse = ";")
  out.df$snpName <- as.character(sapply(strsplit(out.df$snpName,";"),"[[",1))
  out.df <- out.df[,c(5,2,3,4,1)]
  temp.df <- data.frame(snpName=mcols(missingHits)[,1],
                        entrez.genes=NA,
                        hugo.names=NA,
                        features=NA,
                        order=mcols(missingHits)[,ncol(mcols(missingHits))],
                        stringsAsFactors = F)

  out.df <- rbind(out.df,temp.df)
  out.df <- out.df[order(out.df$order,decreasing = F),]
  out.df$order <- NULL
  return(out.df)
}

#' @keywords internal
.createAnnotationDB <- function(promoterRange,upstream,downstream,gene.body,promoters,utr,eqtl,enhancers,loopDist) {
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

  if (gene.body) {
    ret <- .gene.annotation.gr
    ret$Transcription.Start.Site..TSS. <- NULL
  }
  else
    ret <- GRanges()

  if(enhancers) {
    new.loops.gr <- GRanges(seqnames=seqnames(.loops.gr),
                            ranges=IRanges(start=start(.loops.gr)-loopDist,
                                           end=end(.loops.gr)+loopDist))
    mcols(new.loops.gr) <- mcols(.loops.gr)
    enhancer.overlaps <- findOverlaps(new.loops.gr,.all.enhancers.gr,ignore.strand=T)

    new.loops.gr <- GRanges(seqnames=as.character(seqnames(new.loops.gr)[queryHits(enhancer.overlaps)]),
                            ranges=IRanges(start=as.integer(new.loops.gr$y1)[queryHits(enhancer.overlaps)]-loopDist,
                                           end=as.integer(new.loops.gr$y2)[queryHits(enhancer.overlaps)]+loopDist))
    new.loops.gr$enhancer.start <- start(.all.enhancers.gr)[subjectHits(enhancer.overlaps)]
    new.loops.gr$enhancer.end <- end(.all.enhancers.gr)[subjectHits(enhancer.overlaps)]

    promoters.with.genes <- c(.encode.promoters.distal.gr[.encode.promoters.distal.gr$genes != ""],.encode.promoters.prox.gr[.encode.promoters.prox.gr$genes != ""])

    promoter.overlaps <- findOverlaps(new.loops.gr,promoters.with.genes,ignore.strand=T)
    time.to.repeat <- sapply(promoters.with.genes$genes[subjectHits(promoter.overlaps)],length)
    enhancer.annot <- GRanges(seqnames=unlist(mapply(rep,as.character(seqnames(new.loops.gr))[queryHits(promoter.overlaps)],time.to.repeat)),
                              ranges=IRanges(start=unlist(mapply(rep,as.integer(new.loops.gr$enhancer.start)[queryHits(promoter.overlaps)],time.to.repeat)),
                                             end=unlist(mapply(rep,as.integer(new.loops.gr$enhancer.end)[queryHits(promoter.overlaps)],time.to.repeat))))
    enhancer.annot$genes <- unlist(promoters.with.genes$genes[subjectHits(promoter.overlaps)])

    enhancer.annot$Feature <- "Enhancer"
    ret <- c(ret,enhancer.annot)

    .encode.promoters.distal.gr$enhancer.starts <- NULL
    .encode.promoters.distal.gr$enhancer.ends <- NULL

    .encode.promoters.prox.gr$enhancer.starts <- NULL
    .encode.promoters.prox.gr$enhancer.ends <- NULL
  }

  if (promoters)
    ret <- c(ret,.encode.promoters.distal.gr[.encode.promoters.distal.gr$genes != ""],.encode.promoters.prox.gr[.encode.promoters.prox.gr$genes != ""])

  if (utr)
    ret <- c(ret,.utr.entrez.gr)

  if (upstream>0) {
    upstream.gr <- GRanges(seqnames=seqnames(.gene.annotation.gr),
                           ranges=IRanges(start=ifelse(strand(.gene.annotation.gr)=="-",.gene.annotation.gr$Transcription.Start.Site..TSS.,.gene.annotation.gr$Transcription.Start.Site..TSS.-upstream),
                                          end=ifelse(strand(.gene.annotation.gr)=="-",.gene.annotation.gr$Transcription.Start.Site..TSS.+upstream,.gene.annotation.gr$Transcription.Start.Site..TSS.)),
                           strand=strand(.gene.annotation.gr),genes=mcols(.gene.annotation.gr)[,"genes"])
    upstream.gr$Feature = "Upstream"
    ret <- c(ret,upstream.gr)
  }

  if (downstream>0) {
    downstream.gr <- GRanges(seqnames=seqnames(.gene.annotation.gr),
                             ranges=IRanges(start=ifelse(strand(.gene.annotation.gr)=="-",end(.gene.annotation.gr)-downstream,end(.gene.annotation.gr)),
                                            end=ifelse(strand(.gene.annotation.gr)=="-",end(.gene.annotation.gr),end(.gene.annotation.gr)+downstream)),
                             strand=strand(.gene.annotation.gr),genes=mcols(.gene.annotation.gr)[,"genes"])
    downstream.gr$Feature = "Downstream"
    ret <- c(ret,downstream.gr)
  }

  if (eqtl)
    ret <- c(ret,.eqtl.gr)
  return(ret)
}

#' @keywords internal
.getSeedGenes <- function(Hits,annotations) {
  snpOverlaps <- findOverlaps(Hits,annotations,ignore.strand=T)
  if (length(snpOverlaps)<1)
    stop("No input hits were mapped to any gene using these annotation parameters")
  missingHits <- Hits[!(1:length(Hits) %in% queryHits(snpOverlaps))]
  pheno <- mcols(Hits)[queryHits(snpOverlaps),2]
  time.to.repeat <- sapply(annotations$genes[subjectHits(snpOverlaps)],length)
  snpPheno <- as.character(unlist(mapply(rep, pheno, time.to.repeat)))

  geneToPheno <- unique(matrix(c(unlist(annotations$genes[subjectHits(snpOverlaps)]),snpPheno),ncol = 2))

  snpToGene <- unique(data.frame(snp=unlist(mapply(rep, queryHits(snpOverlaps), time.to.repeat)),
                                 gene=unlist(annotations$genes[subjectHits(snpOverlaps)]),
                                 feature=unlist(mapply(rep, annotations$Feature[subjectHits(snpOverlaps)], time.to.repeat))))

  hitsPerGene <- aggregate(feature ~ gene, data = snpToGene, paste, collapse = ",")
  hitsPerGene$Freq <- ifelse(grepl(",",hitsPerGene$feature)==F,1,sapply(gregexpr(",",hitsPerGene$feature),length)+1)
  #hitsPerGene <- as.data.frame(table(snpToGene[,2]),stringsAsFactors=F)

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

  non.ppi.f0 <- sapply(.non.ppi.genes, function(x) sum(seedScore[as.character(x)],na.rm=T))
  return(list(matrix(c(corum.f0,indiv.f0)),non.ppi.f0))
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
.calcGeneScores <- function(Ffinal,geneScoring) {
  scores <- sapply(.geneToComplex, function(x) {
    geneScoring(Ffinal[x])
  })
  ret <- data.frame("Entrez.Gene.ID"=names(.geneToComplex),"HUGO Gene Name"=.hugoNames,"Complexes"=.complexNames,"score"=scores,"Gene.in.Network"="Yes",stringsAsFactors = F)
  ret <- ret[order(ret[,2],decreasing = T),]
  return(ret)
}
