### Fucntions which I find useful in mulitple different projects/analyses ###


#' Function to make a cummulative x axis for chromosomes across the genome
#'
#'
#'
#'
#' @export
make.axis.adjust.cummulative <- function(input.with.cummulative=NA, chrlength=NA){
  
  
  if(!all(is.na(input.with.cummulative))){
  chrname <- names(input.with.cummulative)[names(input.with.cummulative) %in% c("Chromosome","chromosome","Chr","chr","chromosome_name","seqnames")][1]

  #if using X chromosome relapse with orderable numberic
  chromosomes <- as.character(input.with.cummulative[, chrname])
  chromosomes[chromosomes == "X"] <- 23
  chromosomes <- as.numeric(chromosomes)
  
  means <-sapply(unique(chromosomes), function(chr) ((max(input.with.cummulative[chromosomes==chr,"End.cummulitive"]) - min(input.with.cummulative[chromosomes==chr,"Start.cummulitive"])) / 2) + min(input.with.cummulative[chromosomes==chr,"Start.cummulitive"]) )
  out <- data.frame(chr = unique(input.with.cummulative[, chrname]),
                    center = means)
}
  if(!all(is.na(chrlength))){
    chrlengths.cum <- sapply(1:length(chrlength), function(i) sum(chrlength[1:i]))
    means <- sapply(1:length(chrlength), function(i){
      if(i==1) lastchr <- 1 else lastchr <- chrlengths.cum[i-1]
      return(((chrlengths.cum[i] - lastchr)/2) + lastchr)
    })
    out <- data.frame(chr = names(chrlength),
                      center = means)
  }
  return(out)
}

#' add cummulaative start / end columns too order the genome chr1 -> chrX
#'
#'
#'
#'
#' @export
add.cummulative.start.end <- function(input, add.GrCh37.scale=TRUE){
  
  chrname <- names(input)[names(input) %in% c("Chromosome","chromosome","Chr","chr","chromosome_name","seqnames")][1]
  startname <- names(input)[names(input) %in% c("start","Start","startpos","Startpos","start_position")][1]
  endname <- names(input)[names(input) %in% c("end","End","endpos","Endpos","end_position")][1]
  
  #if using X chromosome relapse with orderable numberic
  chromosomes <- as.character(input[, chrname])
  chromosomes[chromosomes == "X"] <- 23
  chromosomes <- as.numeric(chromosomes)
  
  input <- input[order(as.numeric(chromosomes)),]
  
  #add bin.index
  input$bin.index <- 1:nrow(input)
  

  #get rid of factors
  input <- as.data.frame(apply(input,2,as.character),stringsAsFactors=F)
  input[,2:3] <- as.data.frame(apply(input[,2:3],2,as.numeric),stringsAsFactors=F)
  
  if(add.GrCh37.scale){
    load(file=paste0(home,"/proj-tracerx-lung/tctProjects/frankella/Tx_exome/Inputs/GrCh37.scale.rda"))
  } 
  if(!add.GrCh37.scale){
    scale <- lapply(unique(chromosomes), function(chr) max(input[chromosomes==chr,endname]))
    names(scale) <- unique(chromosomes)
  }
  
  nCHR <- length(scale)
  input$Start.cummulitive <- NA
  s <- 0
  nbp <- c()
  for (chr in as.numeric(names(scale))){
    nbp[chr] <- scale[[which(names(scale)==chr)]]
    input[chromosomes == chr,"Start.cummulitive"] <- input[chromosomes == chr,startname] + s
    s <- s + nbp[chr]
  }
  input$End.cummulitive <- NA
  s <- 0
  nbp <- c()
  for (chr in as.numeric(names(scale))){
    nbp[chr] <- scale[[which(names(scale)==chr)]]
    input[chromosomes == chr,"End.cummulitive"] <- input[chromosomes == chr,endname] + s
    s <- s + nbp[chr]
  }
  
  return(input)
  
}


#' function to round to a specific base 
#'
#'
#'
#'
#' @export
mround <- function( x, base ) base * round( x / base )

#' function to roun down to a specific base 
#'
#'
#'
#'
#' @export
mfloor <- function( x, base ) base * floor( x / base )

#' function to round up to a specific base 
#'
#'
#'
#'
#' @export
mceiling <- function( x, base ) base * ceiling( x / base )


#' function to overlay tree information onto a TRACERx mutTable
#'
#'
#'
#'
#' @export
Overlay.tree.parents <- function(mutTable, trees, clone.field, patient.field){
  mutTable$Tree.parent <- NA
  pats <- unique(mutTable[,patient.field])
  pat.clone.table <- do.call(rbind,lapply(pats, function(pat){
    print(pat)
    if( !any(names(trees) == pat)) return(NULL)
    pat.tree <- trees[[which(names(trees)==pat)]]
    pat.clones <- unique(mutTable[mutTable[,patient.field]==pat,clone.field])
    pat.tree <- remove.clones.on.tree(pat.tree,clones.to.keep = pat.clones)
    pat.tree <- logically.order.tree(pat.tree)
    
    if(class(pat.tree)=="character" | class(pat.tree)=="numeric") pat.tree <- matrix(pat.tree,ncol = 2,byrow = TRUE)
    root <-pat.tree[1,1]
    out <- data.frame(SampleID = pat, parent = c("root",pat.tree[,1]), daughter = c(root,pat.tree[,2]),stringsAsFactors = F)
    
    if(all(as.character(pat.tree)==root)){
      out <- data.frame(SampleID = pat, parent = "root", daughter = root,stringsAsFactors = F)
    }
    return(out)
    
  }))
  
  order.change <- order(paste(mutTable[,patient.field],mutTable[,clone.field]))
  mutTable <- mutTable[order.change,]
  
  #overlay
  mutTable$Tree.parent <- unlist(lapply(unique(paste(mutTable[,patient.field],mutTable[,clone.field])), function(pat.clone){
    mut.no <- sum( paste(mutTable[,patient.field],mutTable[,clone.field]) ==  pat.clone)
    parent <- pat.clone.table[paste(pat.clone.table$SampleID,pat.clone.table$daughter) == pat.clone,"parent"]
    if(length(parent)==0) parent <- NA
    out <- rep(parent,mut.no)
    return(out)
  }))
  
  #reverse order change
  mutTable <- mutTable[match(1:nrow(mutTable),order.change),]
  
  return(mutTable)
  
}

#' function to exact tree that hase been overlaid onto a TRACERx mutTable
#'
#'
#'
#'
#' @export
extract.tree.mutTable <- function(mutTable, clone.field, patient.field=NA, parent.field){
  if(!is.na(patient.field)){
    pats <- unique(mutTable[,patient.field])
    
    out <- lapply(pats, function(pat){
      tree <- as.matrix(mutTable[mutTable[,patient.field]==pat ,c(parent.field,clone.field)])
      
      tree <- tree[!duplicated(tree[,2]),]
      if(class(tree)=="character" | class(tree)=="numeric") tree <- matrix(tree,ncol = 2,byrow = TRUE)
      
      tree <- tree[!(is.na(tree[,2])| is.na(tree[,1])),]
      if(class(tree)=="character" | class(tree)=="numeric") tree <- matrix(tree,ncol = 2,byrow = TRUE)
      
      if(any(!tree[,1]=="root")){
        tree <- tree[!tree[,1]=="root",]
      } else {
        tree[tree[,1]=="root",1] <- tree[tree[,1]=="root",2]
      }
      if(class(tree)=="character" | class(tree)=="numeric") tree <- matrix(tree,ncol = 2,byrow = TRUE)
      
      tree <- apply(tree,2,as.numeric)
      if(class(tree)=="character" | class(tree)=="numeric") tree <- matrix(tree,ncol = 2,byrow = TRUE)
      
      tree <- apply(tree,2,as.character)
      if(class(tree)=="character" | class(tree)=="numeric") tree <- matrix(tree,ncol = 2,byrow = TRUE)
      
      colnames(tree) <- NULL
      return(tree)
    })
    names(out) <- pats
    
  } else {
    
    tree <- as.matrix(mutTable[ ,c(parent.field,clone.field)])
    
    tree <- tree[!duplicated(tree[,2]),]
    if(class(tree)=="character" | class(tree)=="numeric") tree <- matrix(tree,ncol = 2,byrow = TRUE)
    
    tree <- tree[!(is.na(tree[,2])| is.na(tree[,1])),]
    if(class(tree)=="character" | class(tree)=="numeric") tree <- matrix(tree,ncol = 2,byrow = TRUE)
    
    if(any(!tree[,1]=="root")){
      tree <- tree[!tree[,1]=="root",]
    } else {
      tree[tree[,1]=="root",1] <- tree[tree[,1]=="root",2]
    }
    if(class(tree)=="character" | class(tree)=="numeric") tree <- matrix(tree,ncol = 2,byrow = TRUE)
    
    tree <- apply(tree,2,as.numeric)
    if(class(tree)=="character" | class(tree)=="numeric") tree <- matrix(tree,ncol = 2,byrow = TRUE)
    
    tree <- apply(tree,2,as.character)
    if(class(tree)=="character" | class(tree)=="numeric") tree <- matrix(tree,ncol = 2,byrow = TRUE)
    
    colnames(tree) <- NULL
    out <- tree
  }
  
  
  return(out)
}


#mutTableAllpriority will keep all its mutations and those mutations in mutTableAllextra that are duplicated in mutTableAllpriority will be removed

#' function to merge MutTables ensuring the sampleID formats etc are consistent 
#'
#' 
#'
#'
#' @export
integrate.mutTables <- function(mutTableAllpriority,mutTableAllextra){
  
  #match IDs
  mutTableAllpriority$SampleID <- trx_rename.fn(mutTableAllpriority$SampleID)
  mutTableAllextra$SampleID <- trx_rename.fn(mutTableAllextra$SampleID)
  
  mutTableAllpriority$mutation_id <- rename.TxIDs.in.mutid(mutTableAllpriority$mutation_id)
  mutTableAllextra$mutation_id <- rename.TxIDs.in.mutid(mutTableAllextra$mutation_id)
  
  #remove dups in extra
  mutTableAllextra <- mutTableAllextra[!mutTableAllextra$mutation_id %in% unique(mutTableAllpriority$mutation_id),]
  
  #full join (missing cols in present in the other will be NA)
  suppressPackageStartupMessages(library(dplyr))
  
  mutTableAllout <- suppressMessages(full_join(mutTableAllpriority,mutTableAllextra))

  return(mutTableAllout)  
  
}

#' Function to rename any TRACERx patient ids to the latest version
#'
#'
#'
#'
#' @export
trx_rename.fn <- function(trxid, trialID='LTX'){
  new_trxid <- as.character(trxid)
  #Includes hospital id? 
  hospital<-grepl('\\w_', substr(trxid[1],1,2))
  if(hospital){
    new_trxid <- substring(new_trxid,3)
  }
  new_trxid <- as.numeric(sub(trialID, '', new_trxid))
  new_trxid <- sprintf("%04d", new_trxid)
  new_trxid <- paste(trialID, new_trxid,sep='')
  return(new_trxid)
}


#'  function to rename TRACERx mutation ids with correct most up to datte patient ids
#'
#'
#'
#'
#' @export
rename.TxIDs.in.mutid <- function(mutids) {
  seperatedrownames <- strsplit(mutids,split=":")
  return(
    paste( trx_rename.fn( sapply( seperatedrownames,function(x) x[1])),
           sapply( seperatedrownames,function(x){ paste(x[2:4], collapse  = ":") }),
           sep = ":")
  )
}


#' Function to adjust any quantification (subject) for cellularity (cell) and then average accross regions to get 1 quantification per tumour
#'
#'
#'
#'
#' @export
cal.Subject.accross.regions <- function(SubjectbyRegion,CellbyRegion){
  if(is.na(SubjectbyRegion)){
    
    return(NA)
    
  }
  
  #Extract subject in each region
  Subject <- data.frame(Subject = sapply(strsplit(strsplit(SubjectbyRegion,split=";")[[1]],split = ":"), function (x) x[2]),
                        Region = as.character(sapply(strsplit(strsplit(SubjectbyRegion,split=";")[[1]],split = ":"), function (x) x[1])),stringsAsFactors = F)
  #Extract cellularity in each region
  Cell <- data.frame(Cellularity = sapply(strsplit(strsplit(CellbyRegion,split=";")[[1]],split = ":"), function (x) x[2]),
                     Region = as.character(sapply(strsplit(strsplit(CellbyRegion,split=";")[[1]],split = ":"), function (x) x[1])),stringsAsFactors = F)
  
  #someimtes a "-" is used instrad of a "." seperater in region names - make consistant either way
  Subject$Region <- gsub("\\-","\\.",Subject$Region)
  
  #remove NAs & where no cellularity is availible
  nodata.region <- Cell[Cell$Cellularity == "Region.Absent.Cellularity.Data", "Region"]
  if(length(nodata.region)>0){
    Cell <- Cell[ !Cell$Region == nodata.region, ]
    Subject <- Subject[ !Subject$Region == nodata.region, ]
  }
  
  #if no data left after removal of NAs simpy return NA
  if(nrow(Cell) == 0 | nrow(Subject) == 0){
    
    return(NA)
    
  } 
  
  #make remaining data numeric
  Cell$Cellularity <- as.numeric(Cell$Cellularity)
  Subject$Subject <- as.numeric(Subject$Subject) 
  
  #if we have extra regions in ASCAT output but not in Muttable (failed mutect eg)
  Cell <- Cell[Cell$Region %in% Subject$Region,]
  #also remove for now if not in ASCAT (perhasp failed for now) but in muttable
  Subject <-Subject[Subject$Region %in% Cell$Region,]
  
  #reorder 
  Cell <- Cell[match( Cell$Region, Subject$Region ),]
  
  #Calcuate the fraction of total tumour cells in each sample, adjusted by if tumour cells spread evenly accross all regions (ie 1.1 = if 10% more tumour cells here than would be expected in even spread)
  adjCell <- Cell$Cellularity / mean(Cell$Cellularity)
  
  #estimate meanSubject accross whole tumour
  return ( mean( Subject$Subject * adjCell ) ) 
  
}

#' Function to FDR correct a set of p-values even when duplicates exist and only a subset is being considered for correction 
#'
#'
#'
#'
#' @export
correct.specific.ps <- function(table,field.to.deduplicate=NA,field.to.correct, indices.to.correct = "all",method ="BH"){
  
  if(all(indices.to.correct =="all")){
    indices.to.correct <- rep(TRUE,nrow(table))
  }
  
  #check classes are correct
  table[,field.to.deduplicate] <- as.character(table[,field.to.deduplicate])
  table[,field.to.correct] <- as.numeric(table[,field.to.correct])
  
  #need to ensure all dedup fields ordered together (makes code a lot faster) then can restore order at the end
  order.change <- order(paste(table[,field.to.deduplicate]))
  table <- table[order.change,]
  indices.to.correct <- indices.to.correct[order.change]
  
  if(!is.na(field.to.deduplicate)){
    table.filtered <- table[indices.to.correct,]
    ps.to.correct <- table.filtered[!duplicated(table.filtered[,field.to.deduplicate]),field.to.correct]
    unique.dedups <- table.filtered[!duplicated(table.filtered[,field.to.deduplicate]),field.to.deduplicate]
  } else {
    ps.to.correct <- table[indices.to.correct,field.to.correct]
  }
  
  if(grepl("\\.p$|_p$",field.to.correct)) new.field.name <- gsub("\\.p$|_p$",".q",field.to.correct) else new.field.name <- paste0(field.to.correct,".q")
  
  table[,new.field.name] <- NA
  ps.corrected <- p.adjust(ps.to.correct, n = length(ps.to.correct), method = method)
  
  if(!is.na(field.to.deduplicate)){
    table[,new.field.name][indices.to.correct] <- round(unlist(lapply(1:length(ps.corrected), function(i) rep(ps.corrected[i],sum(table[,field.to.deduplicate][indices.to.correct] == unique.dedups[i])))),digits=4)
  } else {
    table[,new.field.name][indices.to.correct] <- ps.corrected
  }
  
  #restore order
  table <- table[match(1:nrow(table),order.change),]
  
  #make the q column next to the p column
  field.to.correct.col <- which(names(table)==field.to.correct)
  table <- table[,c(1:field.to.correct.col,ncol(table),(field.to.correct.col+1):(ncol(table)-1))]
  
  return(table)
  
}

#' Function to extract CCF tables from a TRACERx mutTable containing many samples
#'
#'
#'
#'
#' @export
extract_ccf_table <- function( muttable_df, ccf_col = "PhyloCCF_MPhase", sampleid_col = "tumour_id", cluster_col = "PyCloneCluster_MPhase", driver_col = "DriverMut" ){
  
  ## malnipulate as data/table bt return to original class (tibble, data.frame) at the end ##
  
  orig_class <- class(muttable_df)
  muttable_df <- data.table::as.data.table(muttable_df)
  
  ## which samples in the table have data in CCf column ##
  
  samples <- muttable_df[ !is.na( get(ccf_col) ), unique( get(sampleid_col) ) ] 
  
  ## output list of CCF tables for ech sample ##
  
  out <- lapply(samples, function(sample){
    
    muttable_dfsamp <- muttable_df[ get(sampleid_col) == sample ]
    
    CCFs <- muttable_dfsamp[ !is.na( get(cluster_col) ), .(get(cluster_col), get(ccf_col)) ]
    setnames(CCFs, c("cluster", "CCFs")) 
    
    ccf_by_region <- data.table::as.data.table( do.call(rbind, lapply( strsplit(CCFs[, CCFs],split=";"), function(x) sapply(strsplit(x,split = ":"), function (x) x[2]))) )
    data.table::setnames(ccf_by_region, sapply(strsplit( strsplit(CCFs[, CCFs],split=";")[[1]],split = ":"), function (x) x[1]) )
    
    ccf_by_region <- apply( ccf_by_region, 2, as.numeric)
    
    CCFs <- cbind( CCFs, ccf_by_region)
    
    region_cols <- names(CCFs)[3:ncol(CCFs)]
    
    meanCCF <- rowMeans( CCFs[, ..region_cols ] )
    
    CCFs[, meanCCF := meanCCF ]
    
    mean_region_cols <- paste0("mean", region_cols)
    
    ## couldn't find a good data table wy of doing all regiono cols at once
    for( i in 1:length(region_cols) )  CCFs[, (mean_region_cols[i]) := mean( get(region_cols[i])), by = cluster, ]
    
    ### now summarise by cluster ###
    
    CCF.table <- cbind( CCFs[, .(meanCCF = mean(meanCCF)), by = cluster, ], CCFs[ !duplicated(cluster), ..mean_region_cols] )
    
    ## add on other info (no of muts and whther there are drivers) ##
    
    CCF.table[, driver := sapply(cluster, function(cluster_id) any( muttable_dfsamp[ get(cluster_col) == cluster_id, get(driver_col) ] == 'TRUE' ))]
    CCF.table[, no_muts := sapply(cluster, function(cluster_id) muttable_dfsamp[ get(cluster_col) == cluster_id, .N ])  ]
    
    setcolorder(CCF.table, c("cluster", "no_muts", "driver", "meanCCF"))
    
    
    class(CCF.table) <- orig_class
    
    return(CCF.table)
    
  })
  
  names(out) <- samples
  
  return(out)
  
}


###################################################
### Function to caulcated wgII from seg CN data ###
###################################################

## code written by Nicolai ##
## Input is an ASCAT output matrix, with nA in column 7, nB in column 8 #

#' Function to caulcated wgII from seg CN data
#'
#'
#'
#'
#' @export
calc.wgii <- function(seg, threshold = 0.6, check.names=FALSE, include.sexchrom=FALSE){
  # Edit 20140826: gii & wgii: exclude sex chromosomes
  # Edit 20150316: Use wMajor as ploidy
  seg[,1] <- as.character(seg[,1])
  if(check.names){
    seg <- check.names.fn(seg)
  }
  output.mat <- matrix(NA, nrow=length(unique(seg[,1])), ncol=4)
  colnames(output.mat) <- c('GII', 'wGII', 'FLOH', 'wFLOH')
  rownames(output.mat) <- unique(seg[,1])
  chrom.length <- setNames(rep(NA, length(unique(seg[,2]))), unique(seg[,2]))
  if(! all(seg[,8] <= seg[,7]) ){
    cat("Warning!! nB  not always <= nA!!  -- Correcting for internal use (only!)\n") # In case ASCAT people change the algorithm
    tmp <- seg
    seg[tmp[,8] > tmp[,7],7]  <- tmp[tmp[,8] > tmp[,7],8]
    seg[tmp[,8] > tmp[,7],8]  <- tmp[tmp[,8] > tmp[,7],7]
  }
  cat("Setting sample ploidy according to major proportion copy number. Ploidy < 2 re-set at 2\n")
  
  message( "Calculating pliody for each sample..\n" )
  
  ploidy <- calc.ploidy(seg, check.names=check.names)
  ploidy <- setNames(ploidy[,2], rownames(ploidy))
  ploidy[ploidy < 2] <- 2
  for(i in names(chrom.length)){chrom.length[i] <- max(seg[seg[,2] %in% i,4]) - min(seg[seg[,2] %in% i,3])} # restrict to the part that can be measured
  chrom.length <- setNames(as.numeric(chrom.length), names(chrom.length))
  
  message( "Calculating metrics for each sample..\n" )
  
  pb <- txtProgressBar(1, length(unique(seg[,1])), style=3, width = 30)
  
  for(i in unique(seg[,1])){
    sample.seg <- seg[seg[,1] %in% i,] # Restrict to sample
    sample.seg.cn <- sample.seg[sample.seg[,6] < (ploidy[i] - threshold) | sample.seg[,6] > (ploidy[i] + threshold),] # Restrict to aberrant segments
    sample.seg.loh <- sample.seg[sample.seg[,8] == 0 & sample.seg[,7] != 0,] # Restrict to LOH segments
    output.mat[i,1] <- sum(as.numeric(sample.seg.cn[,4]) - as.numeric(sample.seg.cn[,3]))/sum(chrom.length) # GII
    output.mat[i,3] <- sum(as.numeric(sample.seg.loh[,4] - sample.seg.loh[,3]))/sum(chrom.length) # FLOH
    wgii <- vector()
    wfloh <- vector()
    if(!include.sexchrom){
      chrom.length <- chrom.length[!names(chrom.length) %in% c('X','Y','x','y',23,24)]
    }
    for(j in names(chrom.length)){
      sample.seg.chr <- sample.seg.cn[sample.seg.cn[,2] %in% j,]
      sample.seg.chr.loh <- sample.seg.loh[sample.seg.loh[,2] %in% j,]
      wgii <- c(wgii, sum(sample.seg.chr[,4] - sample.seg.chr[,3])/chrom.length[j])
      wfloh <- c(wfloh, sum(sample.seg.chr.loh[,4] - sample.seg.chr.loh[,3])/chrom.length[j])
    }
    output.mat[i,2] <- sum(wgii)/length(chrom.length)
    output.mat[i,4] <- sum(wfloh)/length(chrom.length)
    
    setTxtProgressBar(pb, which( unique(seg[,1]) == i ))
    
  }
  return(output.mat)
}


#' Function to calculate ploidy distribution and main ploidy per sample from seg input (Nicolai)
#'
#'
#'
#'
#' @export
calc.ploidy <- function(seg, cnCol=6,check.names=FALSE){
  #20140311 : added wMajor, weighted ploidy by chromosome
  #20151209 : Added the option to call mean raw copy if it exists
  #seg = ASCAT segmented output (or any segmented output), with total CN in column 6
  seg[,1] <- as.character(seg[,1])
  if(check.names){
    seg <- check.names.fn(seg)
  }
  samples <- unique(seg[,1])
  out.ploidy <- matrix(nrow=length(samples), ncol=14)
  rownames(out.ploidy) <- samples
  colnames(out.ploidy) <- c('Major','wMajor','MeanRaw',0:10) #wMajor is weighted by chromosome, MeanRaw is mean of the raw CN values from ASCAT 2.3
  
  pb <- txtProgressBar(1, length(samples), style=3, width = 30)
  
  for(j in samples){
    sample.seg <- seg[seg[,1] %in% j,]
    ploidy <- vector()
    for(k in unique(sample.seg[,cnCol])){
      tmp <- sample.seg[sample.seg[,cnCol] %in% k,]
      ploidy <- c(ploidy, setNames(sum(tmp[,4]-tmp[,3]), k))
    }
    ploidy <- (ploidy/sum(ploidy))[order(ploidy/sum(ploidy),decreasing=T)]
    out.ploidy[j,1] <- as.numeric(names(ploidy)[1])
    out.ploidy[j,3:13] <- ploidy[colnames(out.ploidy)[3:13]]
    wMajor <- c()
    for(k in unique(sample.seg[,2])){
      tmp <- sample.seg[sample.seg[,2] %in% k,]
      tmp <- setNames(tmp[,4]-tmp[,3], tmp[,cnCol])
      wMajor <- c(wMajor,setNames(as.numeric(names(sort(tapply(tmp, names(tmp),sum), decreasing=T)[1])), k))
    }
    out.ploidy[j,2] <- as.numeric(names(sort(table(wMajor), decreasing=T)[1]))
    if(any(grepl('nAraw', colnames(seg)))){
      a <- rowSums(sample.seg[,c('nAraw','nBraw')])
      b <- sample.seg[,4]-sample.seg[,3]
      out.ploidy[j,'MeanRaw'] <- sum(a * b/sum(b))
    }
    setTxtProgressBar(pb, which( samples == j ))
  }
  if(!any(grepl('nAraw', colnames(seg)))){
    out.ploidy <- out.ploidy[,-c(3)]
  }
  out.ploidy
}




#===================================#
# interpret_following_Args function #
#===================================#

###  function designed to assign objects to following arguemnts extracted by commandArgs(). ###
###  takes commardArgs output (args) and a table (arg.table) specifying a automatic order   ###
###  of object names to assign to each argument and a list of flags and matching object     ###
###  names so that flags can be used in any order to specify what objects following         ###
###  arguments should be assigned                                                           ###

# # example data:
# 
# arg.table <- data.table( arguments  = c("mut_table_path", "scna_table_path", "trees_path", "clin_path", "outputs.folder"),
#                          flags = c("--muts", "--scnas", "--trees", "--clin", "--out") )
# 
# args <- c("--scnas", "path/to/scna.table.txt", "--muts", "path/to/mut.table.txt")
# # ( in terminal = 'Rscript script.path --scnas path/to/scna.table.txt --muts path/to/mut.table.txt'
# #  with in script args <- commandArgs() )
# 
# # or
# 
# args <- c("path/to/mut.table.txt", "path/to/scna.table.txt")
# # ( in terminal = 'Rscript script.path path/to/mut.table.txt path/to/scna.table.txt'
# #   with in script args <- commandArgs() )
# 
# # both egs above of args will be intrepretted the same. The first uses flags to assign arguments and the second uses the order
# # of the arguuemnts field in the arg.table
# 
# # in these cases the function will assign "../../scna.table.txt" to an object called scna_table_path
# # and "../../mut.table.txt" to object called mut_table_path

#' Function to allow R script to interpret following arguments specified in terminal with flags
#'
#'
#'
#'
#' @export
interpret_following_Args <- function(args, arg.table){
  
  if(length(args) > 0){
    
    using_flags <- any( grepl("--", args) )  ## have flags been used in the following args - if not just use the default order 
    ## default order specified by order of args field in args.table 
    ## this means input arguments cannot contain "--" symbols which are reserved for flags
    
    ## if no flags just use the default order to specify ##
    
    if( using_flags == FALSE ){
      
      for( i in 1:length(args) ) assign( arg.table[ i, arguments ], args[i], envir = parent.frame() ) 
      
    }
    
    if( using_flags == TRUE ){
      
      args.input.table <- data.table( input.arguments = args[ seq(2, length(args), by = 2 ) ],    ## odd args will be input arguments
                                      input.flags = args[ seq(1, (length(args) - 1), by = 2 ) ] ) ## even args will be flags
      
      # limit arg table to those present in flags #
      arg.table <- arg.table[ flags %in% sapply( strsplit( args, split = " " ), "[[", 1 ) ]
      
      arg.table[, arg.input := args.input.table[ match( flags, input.flags ), input.arguments] ]
      
      for( i in 1:nrow(arg.table) ) assign( arg.table[ i, arguments ], arg.table[ i, arg.input ], envir = parent.frame() ) 
      
    }
    
  }
  
}

#' Function to extract the tree information from all patients in a tree output directory
#'
#'
#'
#'
#' @export
read_trees <- function( tree_path, return_all = FALSE, remove_missing = TRUE ){
  
  tumours <- system( paste('ls', tree_path ), intern = TRUE )
  tumours <- tumours[ grepl('LTX', tumours) ]
  
  tree_data <- lapply(tumours, function(tumour){
     
    tumour_no_cluster <- gsub( '_Cluster.{1}', '', tumour)
    tree_path <- paste0(tree_path, '/', tumour, '/', tumour_no_cluster, '.tree.RDS' )
    
    if( file.exists(tree_path) ){
    tree <- readRDS( tree_path ) 
    } else return( NULL )
    
  })
  names(tree_data) <- tumours
 
  if( remove_missing ) tree_data <- tree_data[ !sapply( tree_data, function(tree) all(is_null( tree )) ) ]
  
  if( return_all ) return( tree_data )
   
  tree_data <- lapply(tree_data, function(tree) tree$graph_pyclone$Corrected_tree )
  
  return(tree_data)
  
}


#=====#
# END #
#=====#



# suppressPackageStartupMessages(library(sf))
# suppressPackageStartupMessages(library(units))
# suppressPackageStartupMessages(library(rgeos))
# suppressPackageStartupMessages(library(qlcMatrix))
# suppressPackageStartupMessages(library(RColorBrewer))
# suppressPackageStartupMessages(library(smoothr))


