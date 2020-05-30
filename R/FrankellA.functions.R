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



#' function to order tree root -> branches -> leaves
#'
#'
#'
#'
#' @export
logically.order.tree <- function(tree){
  
  # if only 1 clone on tree then just output it as it is #
  if(all(tree==unique(tree)[1])) return(tree)
  
  ### assign levels to each parent clone depending on how near the trunk ###
  
  # work out the root (the only clone tht's never a daughter ) #
  root <- tree[,1] [! tree[,1] %in% tree[,2] ]
  
  # make empty list for levels of ecah clone #
  levels <- rep(NA,nrow(tree))
  
  # aasign the root to level 1 #
  levels[tree[,1] %in% root ] <- 1
  
  # list aall the daughter's of root to work out the level #
  daughters <- tree[ tree[,1] %in% root,2]
  
  # go down the tree in levels and asssign clenns correct levels #
  l <- 1
  repeat{
    l <- l + 1
    parents <- daughters
    levels[tree[,1] %in% parents ] <- l
    daughters <- tree[ tree[,1] %in% parents,2]
    if(all(!is.na(levels))) break
  }
  tree <- tree[order(levels),]
  return(tree)
}



#' function to remove clones from a phenogenetic tree matrix
#'
#' This maintains parent -> daughter relationships,  even if an intermediate
#' (branch) clone is being removed
#'
#'
#' @export
remove.clones.on.tree <- function(tree, clones.to.remove = NA, clones.to.keep = NA){
  
  # check that info have been provided on what clones to remove
  if(all(is.na(clones.to.remove)) & all(is.na(clones.to.keep))){
    stop("you have not provided info on which clones to remove")
  }
  
  # check classes are correct #
  clones.to.keep <- as.character( clones.to.keep )
  tree <- as.matrix( tree )
  
  # define list of all clones in the tree #
  all.clones <- unique(as.character(tree))
  
  # if clones to keep specific rather than clones         ##
  # to remove then work out from this which to be removed ##
  if(all(is.na(clones.to.remove))){
    clones.to.remove <- setdiff(all.clones,clones.to.keep)
  }
  
  # ensure class is correct and all clones specified to keep are in the tree matrix #
  clones.to.remove <- as.character(clones.to.remove)
  clones.to.keep <- clones.to.keep[clones.to.keep %in% all.clones]
  
  # if nothing to remove that is ctually ono the tree then just return the tree with a warning #
  if(length(clones.to.remove)==0){
    warning( "no clones specified to remove that are on the tree\n")
    return(tree)
  }
  
  # order the tree trunk -> branches -> leaves #
  tree <- logically.order.tree(tree)
  
  # get the root clone #
  root <- tree[ 1, 1 ]
  
  # if only the root left then return root as parent and daughter to maintaian the same structure #
  if(all(clones.to.keep==root)){
    return(matrix(c(root,root),ncol = 2, byrow = TRUE))
  }
  
  # now loop round ecah clone to remove, get rid of all relationships its involved in and  #
  # then reassign aany daughter(s) to its parent                                           #
  for(clone in clones.to.remove){
    parent <- tree[tree[,2]==clone,1]
    if(any(tree[,1]==clone)){
      daughters <- tree[tree[,1]==clone,2]
      tree <- rbind(tree, matrix(c(rep(parent,length(daughters)),daughters),ncol = 2))
    }
    tree <- tree[!(tree[,1]==clone | tree[,2]==clone),]
    if(class(tree)=="character" | class(tree)=="numeric") tree <- matrix(tree,ncol = 2,byrow = TRUE)
  }
  
  # return pruned tree #
  return(tree)
  
  # END #
  
}


#' Function to correct CCFs when sum of daughters CCF > parent CCF 
#'
#'
#'
#'
#' @export
make.CCFs.tree.consistant <- function( tree.mat, CCF.data, warning.limit = 1 , parent.adjust = 1,
                                       decrease.daughters = TRUE, increase.parents = FALSE ){
  
  # one of decrease daughters or parents must be true
  if( all( !decrease.daughters & !increase.parents ) ) cat( "Please set either decrease.daughters or decrease.parents arguments to TRUE or cannot correct tree\n" )
  
  # order tree trunk -> leaf #
  tree.mat <- logically.order.tree( tree.mat )
  
  # limit to clones in CCF table #
  tree.mat <- remove.clones.on.tree( tree.mat, clones.to.keep = CCF.data$clones )
  
  # get root #
  root <- tree.mat[ 1, 1 ]
  
  # get ordered list of parents #
  parent.order <- rev(unique(tree.mat[,1]))
  
  # detect if fractions or percentages #
  is_frac <- CCF.data[ CCF.data$clones ==  root, "CCF" ] < 2
  if( is_frac ) clonal_CCF <- 1 else clonal_CCF <- 100
  
  # for each parent check if it has CCF < sum of daughters if so correct it #
  for(parent in parent.order){
    
    # get names of daughteer clones #
    daughters <- tree.mat[ tree.mat[, 1] == parent, 2 ]
    
    # get row indices for daughters and parent #
    parentrow <- which( CCF.data$clones == parent )
    daughterrows <- which( CCF.data$clones %in% daughters )
    
    # get CCF values #
    parent.CCF <- CCF.data[ parentrow, "CCF" ]
    daughter.total.CCF <- sum( CCF.data[ daughterrows, "CCF" ] )
    
    # if decrease parent == TRUE & daughter CCF > parent CCF then increase parent CCF to match daughters, then if #
    # parent CCF > 1 then adjustt all CCFs on the tree to allow parent CCF < 1 #
    # if decrease daughter == TRUE & daughter CCF > parent CCF then decrease daughter CCFs to match paarent, then if #
    # parent adjust allows soome wiggle room for niose in CCF data - default = 0  howeever #
    if( parent.adjust * parent.CCF < daughter.total.CCF ){
      if( daughter.total.CCF / parent.CCF > warning.limit ){
        if( decrease.daughters  )  type <- "Decreasing daughter CCFs proportionally" else type <- "Increasing parent CCF"
        cat( paste0("        ", "clone ", parent, " has daughters with ", signif( (daughter.total.CCF * clonal_CCF) / parent.CCF, 3 ), "% CCF of itself. ", type, "\n") )
      }
      if( increase.parents  ) CCF.data[ parentrow, "CCF" ] <- daughter.total.CCF * parent.adjust
      if( decrease.daughters  ) CCF.data[ daughterrows, "CCF" ] <- sapply(daughterrows, function(rowi) (CCF.data[ rowi, "CCF" ] / daughter.total.CCF) * parent.CCF )
    }
  }
  
  # normalise to root #
  # if the clonal cluster CCF needed to be increased then adjust this back down too 1 and adjust all other clones by similar margin #
  if( any( CCF.data$CCF > clonal_CCF ) ) cat( paste0("        Clonal CCF needed increasing to accommodate daughters. Therefore decreasing all CCFs proportionally so clonal CCF == 1 again\n") )
  
  # ensure clonal cluster = 1 #
  CCF.data$CCF <- ( (CCF.data$CCF * clonal_CCF) / CCF.data[ CCF.data$clones == root, "CCF"] ) 
  
  return( CCF.data )
}



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
  muttable_df <- as.data.table(muttable_df)
  
  ## which samples in the table have data in CCf column ##
  
  samples <- muttable_df[ !is.na( get(ccf_col) ), unique( get(sampleid_col) ) ] 
  
  ## output list of CCF tables for ech sample ##
  
  out <- lapply(samples, function(sample){
    
    muttable_dfsamp <- muttable_df[ get(sampleid_col) == sample ]
    
    CCFs <- muttable_dfsamp[ !is.na( get(cluster_col) ), .(get(cluster_col), get(ccf_col)) ]
    setnames(CCFs, c("cluster", "CCFs")) 
    
    ccf_by_region <- as.data.table( do.call(rbind, lapply( strsplit(CCFs[, CCFs],split=";"), function(x) sapply(strsplit(x,split = ":"), function (x) x[2]))) )
    setnames(ccf_by_region, sapply(strsplit( strsplit(CCFs[, CCFs],split=";")[[1]],split = ":"), function (x) x[1]) )
    
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

#=====#
# END #
#=====#



# suppressPackageStartupMessages(library(sf))
# suppressPackageStartupMessages(library(units))
# suppressPackageStartupMessages(library(rgeos))
# suppressPackageStartupMessages(library(qlcMatrix))
# suppressPackageStartupMessages(library(RColorBrewer))
# suppressPackageStartupMessages(library(smoothr))


