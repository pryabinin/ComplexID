#'  ComplexID: Integrated Hit Annotations and Random Walk of Complex Network
#'
#'  A package with a self-contained annotation database and protein complex network for finding causitive genes for a phenotype of interest.
#'
#' @section ComplexID functions:
#' runComplexID
#' generatePlot
#' annotateHits
#' getPromoterTissues
#' geteQTLTissues
#' getEnhancerTissues
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
#' @import qgraph
NULL

#' Annotates Hits, Performs Random Walk, and Scores Genes
#'
#' Annotates hits to genes, performs random walk with restarts on a network of protein complexes, and then scores each gene in the network for its association with the phenotype of interest
#'
#' @param Hits Granges object with two meta data columns, or a matrix or data frame with at least 5 columns. \cr If it is a Granges object, then the first meta data column is the site's name. The second meta data columns is a phenotype that the site is associated with.\cr
#' If it is a matrix or data frame, then the first column must be the Hit's name, the second column must be chromosome designation, the third column must the starting base pair position, the fourth column must be the ending base pair position (equal to starting bp position for a standard SNP) and the fifth column must a phenotype that the site is associated with. \cr
#' For both Grange objects and matrices/dataframes, each entry/row corresponds to one site that is associated to one phenotype. If a site is associated in multiple phenotypes then there would be multiple entries for the same site but all with different values in the phenotype column
#' @param phenoSim matrix or data frame with two columns. The first column are names of phenotypes that match the same phenotypes found in Hits. The second column are phenotype similarity values between the phenotype in that row and the phenotype of interest (values between 0 and 1), with higher values denoting higher similarity
#' @param promoterRange single integer greater than or equal to zero. How many bases to look upstream of a TSS of a gene in order to find a promoter region for a gene.
#' @param eps single numeric, must be greater than zero. L1 norm threshold between current and previous interations of random walk at which to terminate the random walk
#' @param alpha single numeric in the range of (0,1]. The weight given to the vector of initialized values for the random walk, higher value of alpha means more weight for the initialized values
#' @param upstream single integer. By default 0. How far upstream of a transcription start site a hit can be for it to be annotated to that gene. A NULL value is equivalent to a value of zero (no upstream sites will be annotated to a gene unless they lie in a promoter region, see promoterRange parameter).
#' @param downstream single integer. By default 0. How far downstream of a transcription start site a hit can be for it to be annotated to that gene. A NULL value is equivalent to a value of zero (no downstream sites will be annotated to a gene).
#' @param geneBody TRUE or FALSE, by default TRUE. If TRUE, then hits will be annotated to the bodies (exons and introns) of protein coding genes. If FALSE, hits will not be annotated to those regions.
#' @param promoters TRUE or FALSE, by default TRUE. If TRUE, then hits will be annotated to promoter regions. If FALSE, hits will not be annotated to promoter regions.
#' @param promoterTissues character vector, by default is "all". If "all", then all promoters from all tissues will be included in the annotation, otherwise, only promoter regions from tissues specified by promoterTissues will be used for annotation.
#' @param utr TRUE or FALSE. If TRUE then it will look for hits in the 3' and 5' UTRs of genes, otherwise it will not.
#' @param eqtl TRUE or FALSE. By default TRUE. If TRUE, then hits may be mapped to eQTL loci, and therefore genes effected by those eQTLs be designated as associated to those hits.
#' @param eqtlTissues character vector, by default is "all". If "all", then all eQTLs from all tissues will be included in the annotation, otherwise, only eQTL sites from tissues specified by promoterTissues will be used for annotation.
#' @param enhancers TRUE or FALSE. By default TRUE. If TRUE, then hits may be mapped to enhancer loci and linked to genes via looping structures and promoters
#' @param enhancerTissues character vector, by default is "all". If "all", then all enhancers from all tissues will be included in the annotation, otherwise, only enhancers regions from tissues specified by promoterTissues will be used for annotation.
#' @param loopDist single integer. By default 0. The maximum allowable distance that an enhancer or promoter can be from a looping region to be annotated to it.
#' @param non_proteins TRUE or FALSE. By default FALSE. If TRUE then hits may be mapped to non-protein regions, if FALSE then that annotation will not be used.
#' @param geneScoring a function that takes a vector and outputs a single number. By default the "sum" function. This is the function that will determine the score of a gene based on the scores of the complexes that it belongs to. The input of the function is a vector of numerical values that represent the scores of the complexes that a gene belongs to. Scores are determined by the RWPCN algorithm. The output of the function should be a single numerical value.
#' @param useAllTSS TRUE or FALSE. By default TRUE. If TRUE, then all unique transcription start sites will be considered when looking at upstream regions of a gene (for promoters and upstream regions). If FALSE, it will a single start site for a gene, namely the start of the gene.
#' @details
#' Annotates Hits to genes using a built-in annotation database. Protein coding genes, non-protein coding genes, and UTR annotations come from the ENSEMBL version 89 annotation of GRCH37. Promoter and Enhancer regions are
#' from ENCODE annotation version 3, eQTL are from the gtexportal version 6.\cr
#' \cr
#' After the associated genes for each endophenotype are identified, it performs a Random Walk with Restarts on a pre-constructed protein complex network as in the RWPCN method.
#' The protein complex network was constructed in a similar way is in the RWPCN method. For a PPI we used STRING with a threshold cutoff of 700. Protein IDs in STRING were mapped to approved HUGO names using ENSEMBL and HGNC.
#' Protein complexes were retrieved from CORUM. Any complex with no genes in the PPI was removed along with 5 of the largest complexes (more than 70 subunits) \cr
#' \cr
#' A random walk with restarts is initialized and performed as in RWPCN then all genes in the PPI and complexes are scored according to the weights in the complex network.
#' @return A list with two objects: a data frame called "scores" and a GRanges object "missingHits"
#' The data frame "scores" has seven columns showing the scores of each gene, related to how much that gene is important to the query phenotype, as well as other information about the gene. It is ordered with the highest scoring genes first.\cr
#' The first columns is the HUGO gene names, the second column are the names of the complexes that gene is part of, the third columns is the score for the gene, the fourth column says whether or not the gene was in the PPI and/or a complex, the fifth column says whether or not the gene is a protein coding gene, the sixth column are the features of that gene that have a hit in them, and the seventh column is the number of hits that were annotated to that gene.
#' The GRanges object "missingHits" lists all of the input hits that were not mapped to any gene.
#' @examples
#' data("hits")
#' data("hits.pheno")
#' test <- runComplexID(Hits = hits,phenoSim=hits.pheno,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T)
#' @export
runComplexID <- function(Hits,phenoSim,promoterRange=100000,eps=1e-10,alpha=0.8,upstream=0,downstream=0,geneBody=T,promoters=T,promoterTissues="all",utr=T,eqtl=T,eqtlTissues="all",enhancers=T,enhancerTissues="all",loopDist=0,non_proteins=F,geneScoring=sum,useAllTSS=T) {
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
  annotations <- .createAnnotationDB(promoterRange,upstream,downstream,geneBody,promoters,promoterTissues,utr,eqtl,eqtlTissues,enhancers,enhancerTissues,non_proteins,loopDist,useAllTSS)
  # get set of seed genes linked to a phenotype
  seedgenes <- .getSeedGenes(Hits,annotations)
  # initialize network
  F0 <- .initNetwork(seedgenes$seedgenes,phenoSim,non_proteins)
  # perform random walk with restarts
  Ffinal <- .RWPCN(F0[[1]],eps,alpha)
  # calculate final scores for each gene
  geneScores <- .calcGeneScores(Ffinal,geneScoring)
  # add scores for genes not in PPI or CORUM
  genes.to.add <- .non.ppi.genes
  protein.coding <- rep("Yes",length(genes.to.add))
  if (non_proteins) {
    genes.to.add <- c(genes.to.add,.non.protein.regions)
    protein.coding <- c(protein.coding,rep("No",length(.non.protein.regions)))
  }
  non.ppi.df <- data.frame("HUGO Gene Name"=genes.to.add,
                            "Complexes"="None",
                            "score"=F0[[2]],
                           "Gene.in.Network"="No",
                           "Protein.coding"=protein.coding)
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
#' @param centralGenes character vector. The HUGO names of the genes that around which the network plot will be centered
#' @param order single integer, by default 0, The degree of neighboring genes that will be included. A value of zero means no neighboring genes are included
#' @param colorGenes character vector. By default, is equal to centralGenes. The HUGO names of the genes that will be colored red. The rest of the genes will be colored orange.
#' @param makeTkplot TRUE or FALSE. By default is FALSE. If TRUE, create a tkplot instead of a standard igraph so you can edit it with the GUI.
#' @param spreadNodes TRUE or FALSE. By default is FALSE. If TRUE, uses Fruchterman Reingold algorithm to try to spread the vertices of the graph out further.
#' @param ... arguments passed to plot.igraph function
#' @details
#' Plots a graph of the network of genes centered around the centralGenes, including neighbors out to the degree of "order". Any additional parameters for the plot.igraph function may be passed as well. \cr\cr
#' Vertices are colored red if they are in centralGenes and orange otherwise. \cr
#' \cr
#' Edges are colored red if the interaction exists in the PPI and the two vertices share a complex. Edges are colored black if the interaction exists in the PPI but the two vertices do not share a complex. Edges are colored green if the two vertices share a complex but do not interact in the PPI. \cr\cr
#' Complexes are plotted as ellipses circling the vertices. Vertices that are not circled are not part of any complex. Vertices which are circled individually are part of a complex but none of the other genes are included in this graph.
#' @return A plot of the subnetwork
#' @examples
#' data("hits")
#' data("hits.pheno")
#' test <- runComplexID(Hits = hits,phenoSim=hits.pheno,promoterRange = 10000,upstream = 1000,downstream = 1000,utr = T)
#' generatePlot(test$scores$HUGO.Gene.Name[1:10])
#' @export
generatePlot <- function(centralGenes,order=0,colorGenes=centralGenes,makeTkplot=F,spreadNodes=F,...) {
  graph.to.plot <- induced_subgraph(.complete.igraph,unlist(neighborhood(.complete.igraph,order = order,nodes = as.character(centralGenes))))
  V(graph.to.plot)$color <- ifelse(names(V(graph.to.plot)) %in% as.character(colorGenes),"red","orange")
  complexes.to.plot <- .corum.subunits[sapply(.corum.subunits,function(x) { sum(names(V(graph.to.plot)) %in% x)>0 })]
  complexes.to.plot <- lapply(complexes.to.plot,function(x) { as.character(x)[as.character(x) %in% names(V(graph.to.plot))] })
  if (makeTkplot) {
    if (spreadNodes)
      return(tkplot(graph=graph.to.plot,mark.groups = complexes.to.plot,l=.callqgraphfruchterman,...))
    return(tkplot(graph=graph.to.plot,mark.groups = complexes.to.plot,...))
  }
  if (spreadNodes)
    return(plot.igraph(x=graph.to.plot,mark.groups = complexes.to.plot,l=.callqgraphfruchterman,...))
  return(plot.igraph(x=graph.to.plot,mark.groups = complexes.to.plot,...))
}

#' Annotates Hits by Genomic Features
#' @inheritParams runComplexID
#' @export
annotateHits <- function(Hits,promoterRange=100000,upstream=0,downstream=0,geneBody=T,promoters=T,promoterTissues="all",utr=T,eqtl=T,eqtlTissues="all",enhancers=T,enhancerTissues="all",non_proteins=F,loopDist=0,useAllTSS=T) {
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
  annotations <- .createAnnotationDB(promoterRange,upstream,downstream,geneBody,promoters,promoterTissues,utr,eqtl,eqtlTissues,enhancers,enhancerTissues,non_proteins,loopDist,useAllTSS)

  # annotate input hits
  mcols(Hits) <- cbind(mcols(Hits),data.frame(1:length(Hits)))
  snpOverlaps <- findOverlaps(Hits,annotations,ignore.strand=T)
  if (length(snpOverlaps)<1)
    stop("No input hits were mapped to any gene using these annotation parameters")
  missingHits <- Hits[!(1:length(Hits) %in% queryHits(snpOverlaps))]
  time.to.repeat <- sapply(annotations$genes[subjectHits(snpOverlaps)],length)

  out.df <- data.frame(snpName=as.character(unlist(mapply(rep,mcols(Hits)[queryHits(snpOverlaps),1],time.to.repeat))),
                       hugo.names=unlist(annotations$genes[subjectHits(snpOverlaps)]),
                       features=unlist(mapply(rep,annotations$Feature[subjectHits(snpOverlaps)],time.to.repeat)),
                       order=unlist(mapply(rep,mcols(Hits)[queryHits(snpOverlaps),ncol(mcols(Hits))],time.to.repeat)),
                       stringsAsFactors = F)

  out.df <- aggregate(x = out.df[,1:3], by = list(out.df$order), FUN = paste, collapse = ";")
  names(out.df)[1] <- "order"
  #out.df <- aggregate(cbind(hugo.names,features,snpName)~order, data = unique(out.df), paste, collapse = ";")
  out.df$snpName <- as.character(sapply(strsplit(out.df$snpName,";"),"[[",1))
  out.df <- out.df[,c("snpName","hugo.names","features","order")]
  if (length(missingHits) > 0) {
    temp.df <- data.frame(snpName=mcols(missingHits)[,1],
                          hugo.names=NA,
                          features=NA,
                          order=mcols(missingHits)[,ncol(mcols(missingHits))],
                          stringsAsFactors = F)
    out.df <- rbind(out.df,temp.df)
  }
  out.df <- out.df[order(out.df$order,decreasing = F),]
  out.df$order <- NULL
  return(out.df)
}

#' Lists tissues available for promoter regions
#'
#' Lists the tissues that can be used for the promoterTissues argument in the runComplexID and annotateHits functions
#'
#' @details
#' Lists which tissues can be used for annotating promoter hits.
#' @return A character vector
#' @examples
#' getPromoterTissues()
#' @export
getPromoterTissues <- function() {
  return(unique(.encode.promoters.gr$tissue))
}

#' Lists tissues available for eQTL sites
#'
#' Lists the tissues that can be used for the eqtlTissues argument in the runComplexID and annotateHits functions
#'
#' @details
#' Lists which tissues can be used for annotating eQTL hits.
#' @return A character vector
#' @examples
#' geteQTLTissues()
#' @export
geteQTLTissues <- function() {
  return(unique(.eqtl.gr$tissue))
}

#' Lists tissues available for enhancer regions
#'
#' Lists the tissues that can be used for the enhancerTissues argument in the runComplexID and annotateHits functions
#'
#' @details
#' Lists which tissues can be used for annotating enhancer hits.
#' @return A character vector
#' @examples
#' getEnhancerTissues()
#' @export
getEnhancerTissues <- function() {
  return(unique(.all.enhancers.gr$tissue))
}


#' @keywords internal
.createAnnotationDB <- function(promoterRange,upstream,downstream,geneBody,promoters,promoterTissues,utr,eqtl,eqtlTissues,enhancers,enhancerTissues,non_proteins,loopDist,useAllTSS) {
  if (useAllTSS)
    tss.regions.gr <- GRanges(seqnames=seqnames(.gene.annotation.gr),
                            ranges=IRanges(start=ifelse(strand(.gene.annotation.gr)=="-",.gene.annotation.gr$Transcription.start.site..TSS.,.gene.annotation.gr$Transcription.start.site..TSS.-promoterRange),
                                           end=ifelse(strand(.gene.annotation.gr)=="-",.gene.annotation.gr$Transcription.start.site..TSS.+promoterRange,.gene.annotation.gr$Transcription.start.site..TSS.)),
                            strand=strand(.gene.annotation.gr),mcols(.gene.annotation.gr))
  else {
    tss.regions.gr <- GRanges(seqnames=seqnames(.gene.annotation.gr),
                              ranges=IRanges(start=ifelse(strand(.gene.annotation.gr)=="-",start(.gene.annotation.gr),start(.gene.annotation.gr)-promoterRange),
                                             end=ifelse(strand(.gene.annotation.gr)=="-",start(.gene.annotation.gr)+promoterRange,start(.gene.annotation.gr))),
                              strand=strand(.gene.annotation.gr),mcols(.gene.annotation.gr))
    tss.regions.gr <- unique(tss.regions.gr)
  }

  if (enhancers | promoters) {
    if (promoterTissues == "all")
      current.promoters <- .encode.promoters.gr
    else
      current.promoters <- .encode.promoters.gr[.encode.promoters.gr$tissue %in% promoterTissues]
    current.promoters <- reduce(current.promoters)
    #current.promoters$tissue <- NULL
    promoter.tss.hits <- findOverlaps(current.promoters,tss.regions.gr)
    current.promoters$genes <- ""
    current.promoters$Feature <- "Promoter"
    if (length(promoter.tss.hits) > 0) {
      split.hits.idx <- split(1:length(promoter.tss.hits),queryHits(promoter.tss.hits),drop=T)
      current.promoters$genes[as.integer(names(split.hits.idx))] <- sapply(split.hits.idx, function(x) {
        return(tss.regions.gr$genes[subjectHits(promoter.tss.hits)[x]])
      })
    }
  }

  # promoter.distal.tss.hits <- findOverlaps(.encode.promoters.distal.gr,tss.regions.gr)
  # .encode.promoters.distal.gr$genes <- ""
  # .encode.promoters.distal.gr$Feature <- "Distal_Promoter"
  # if (length(promoter.distal.tss.hits) > 0) {
  #   split.hits.idx <- split(1:length(promoter.distal.tss.hits),queryHits(promoter.distal.tss.hits),drop=T)
  #   .encode.promoters.distal.gr$genes[as.integer(names(split.hits.idx))] <- sapply(split.hits.idx, function(x) {
  #     return(tss.regions.gr$genes[subjectHits(promoter.distal.tss.hits)[x]])
  #   })
  # }
  #
  # promoter.prox.tss.hits <- findOverlaps(.encode.promoters.prox.gr,tss.regions.gr)
  # .encode.promoters.prox.gr$genes <- ""
  # .encode.promoters.prox.gr$Feature <- "Distal_Promoter"
  # if (length(promoter.prox.tss.hits) > 0) {
  #   split.hits.idx <- split(1:length(promoter.prox.tss.hits),queryHits(promoter.prox.tss.hits),drop=T)
  #   .encode.promoters.prox.gr$genes[as.integer(names(split.hits.idx))] <- sapply(split.hits.idx, function(x) {
  #     return(tss.regions.gr$genes[subjectHits(promoter.prox.tss.hits)[x]])
  #   })
  # }

  if (geneBody) {
    ret <- unique(.gene.annotation.gr)
    ret$Transcription.start.site..TSS. <- NULL
  }
  else
    ret <- GRanges()

  if(enhancers) {
    if (enhancerTissues=="all")
      current.enhancers <- .all.enhancers.gr
    else {
      current.enhancers <- .all.enhancers.gr[.all.enhancers.gr$tissue %in% enhancerTissues]
    }
    current.enhancers <- reduce(current.enhancers)
    #current.enhancers$tissue <- NULL
    new.loops.gr <- GRanges(seqnames=seqnames(.loops.gr),
                            ranges=IRanges(start=start(.loops.gr)-loopDist,
                                           end=end(.loops.gr)+loopDist))
    mcols(new.loops.gr) <- mcols(.loops.gr)
    enhancer.overlaps <- findOverlaps(new.loops.gr,current.enhancers,ignore.strand=T)

    new.loops.gr <- GRanges(seqnames=as.character(seqnames(new.loops.gr)[queryHits(enhancer.overlaps)]),
                            ranges=IRanges(start=as.integer(new.loops.gr$y1)[queryHits(enhancer.overlaps)]-loopDist,
                                           end=as.integer(new.loops.gr$y2)[queryHits(enhancer.overlaps)]+loopDist))
    new.loops.gr$enhancer.start <- start(current.enhancers)[subjectHits(enhancer.overlaps)]
    new.loops.gr$enhancer.end <- end(current.enhancers)[subjectHits(enhancer.overlaps)]

    promoters.with.genes <- current.promoters[current.promoters$genes != ""]

    promoter.overlaps <- findOverlaps(new.loops.gr,promoters.with.genes,ignore.strand=T)
    time.to.repeat <- sapply(promoters.with.genes$genes[subjectHits(promoter.overlaps)],length)
    enhancer.annot <- GRanges(seqnames=unlist(mapply(rep,as.character(seqnames(new.loops.gr))[queryHits(promoter.overlaps)],time.to.repeat)),
                              ranges=IRanges(start=unlist(mapply(rep,as.integer(new.loops.gr$enhancer.start)[queryHits(promoter.overlaps)],time.to.repeat)),
                                             end=unlist(mapply(rep,as.integer(new.loops.gr$enhancer.end)[queryHits(promoter.overlaps)],time.to.repeat))))
    enhancer.annot$genes <- unlist(promoters.with.genes$genes[subjectHits(promoter.overlaps)])

    enhancer.annot$Feature <- "Enhancer"
    ret <- c(ret,enhancer.annot)

  }

  if (promoters)
    ret <- c(ret,current.promoters[current.promoters$genes != ""])

  if (utr)
    ret <- c(ret,.utr.gr)

  if (upstream>0) {
    if (useAllTSS)
      upstream.gr <- GRanges(seqnames=seqnames(.gene.annotation.gr),
                             ranges=IRanges(start=ifelse(strand(.gene.annotation.gr)=="-",.gene.annotation.gr$Transcription.start.site..TSS.,.gene.annotation.gr$Transcription.start.site..TSS.-upstream),
                                            end=ifelse(strand(.gene.annotation.gr)=="-",.gene.annotation.gr$Transcription.start.site..TSS.+upstream,.gene.annotation.gr$Transcription.start.site..TSS.)),
                             strand=strand(.gene.annotation.gr),genes=mcols(.gene.annotation.gr)[,"genes"])
    else
      upstream.gr <- GRanges(seqnames=seqnames(.gene.annotation.gr),
                             ranges=IRanges(start=ifelse(strand(.gene.annotation.gr)=="-",start(.gene.annotation.gr),start(.gene.annotation.gr)-upstream),
                                            end=ifelse(strand(.gene.annotation.gr)=="-",start(.gene.annotation.gr)+upstream,start(.gene.annotation.gr))),
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

  if (eqtl) {
    if (eqtlTissues=="all")
      temp <- .eqtl.gr
    else
      temp <- .eqtl.gr[.eqtl.gr$tissue %in% eqtlTissues]
    temp <- unique(temp)
    temp$tissue <- NULL
    ret <- c(ret,temp)
  }
  if (non_proteins)
    ret <- c(ret,.non.protein.annotations.gr)
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
.initNetwork <- function(seedgenes,phenoSim,non_proteins) {
  geneToPheno <- split(seedgenes[,2],seedgenes[,1],drop=T)
  seedScore <- sapply(geneToPheno, function(x) {
    sum(as.numeric(phenoSim[as.character(phenoSim[,1]) %in% as.character(x),2]))
  })
  names(seedScore) <- names(geneToPheno)

  corum.f0 <- sapply(.corum.subunits, function(x) sum(seedScore[as.character(x)],na.rm=T))
  corum.f0 <- corum.f0 * .density
  indiv.f0 <- sapply(rownames(.W_norm_t_spars)[(length(.corum.subunits)+1):nrow(.W_norm_t_spars)],function(x) sum(seedScore[as.character(x)],na.rm=T))

  if (non_proteins)
    non.ppi.f0 <- sapply(c(.non.ppi.genes,.non.protein.regions), function(x) sum(seedScore[as.character(x)],na.rm=T))
  else
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
  ret <- data.frame("HUGO Gene Name"=names(.geneToComplex),"Complexes"=.complexNames,"score"=scores,"Gene.in.Network"="Yes","Protein.coding"="Yes",stringsAsFactors = F)
  ret <- ret[order(ret[,2],decreasing = T),]
  return(ret)
}

#' @keywords internal
.callqgraphfruchterman <- function(g) {
  return(qgraph.layout.fruchtermanreingold(edgelist = get.edgelist(g),vcount=vcount(g),area=8*(vcount(g)^2),repulse.rad=(vcount(g)^3.1)))
}
