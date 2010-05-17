############################################
# analysis flow 
############################################
#
# 1. find attractors
# 2. remove flat genes 
# 3. find synexpression groups 
# 4. functional enrichment within synexpression groups
# 5. correlated sets 
#
############################################
# Step 1 
############################################
## function: findAttractors
############################################
# input arguments: 
#	dat.fr			data matrix, rows = all genes (even unannotated ones), columns = all samples to be considered in GSEA model
#					note: data.fr rownames must correspond to gene names! 
#	class.vector	character vector of sample classes
#	log2			logical indicator, true if data in dat.fr has been log2-transformed, if log2 = FALSE, data is transformed within the function
#	numPerm 		number of permutations to be performed for the GSEAlm permutation P-value step
#	min.pwaysize	defaults to 5, minimum size of pathway to consider (>= 5)
#	...				other arguments 
############################################

findAttractors <- function(myEset, cellTypeTag, nperm=100, min.pwaysize=5, direction="upper", annotation="illuminaHumanv2BeadID.db", ...){
	
	#require(GSEAlm)
	#require(KEGG.db)
	require(annotation, character.only=TRUE)

	ann <- strsplit(annotation, ".db")[[1]]
	loadNamespace(annotation)
	pathEnv <- get(paste(ann, "PATH", sep=""), as.environment(match(paste("package:", annotation, sep=""), search())))
	dat.fr <- exprs(myEset) 
	all.probes <- rownames(dat.fr)
	dat.fr <- as.matrix(dat.fr)

	class.vector <- as.factor(pData(myEset)[,colnames(pData(myEset)) %in% cellTypeTag])
	
	# checking KEGG annotations
	path.hits <- mget(intersect(all.probes, ls(pathEnv)), pathEnv) 

	list.wpway <- sort(names(path.hits)[unlist(sapply(path.hits, flagPwayExists))])
	all.pways <- unlist(mget(list.wpway, pathEnv))
	all.pways <- unique(all.pways)													 

	# making expression data object for kegg annotated probes only 
	dat.detect.wkegg <- dat.fr[rownames(dat.fr) %in% list.wpway,]
	dat.detect.wkegg <- as.matrix(dat.detect.wkegg)
	
	# make geneset incidience matrix 
	kegg.incidence.matrix <- buildKeggIncidenceMatrix(all.pways, rownames(dat.detect.wkegg), annotation)

	keep.pways <- apply(kegg.incidence.matrix, 1, sum) >= min.pwaysize 
	kegg.incidence.matrix <- kegg.incidence.matrix[keep.pways,]

	# making the eset with the current data 
	eset <- new("ExpressionSet")
	eset@assayData <- new.env()
	assign("exprs", dat.detect.wkegg, eset@assayData)

	pheno.dat <- data.frame(colnames(dat.detect.wkegg), class.vector)
	colnames(pheno.dat) <- c("ChipID", cellTypeTag)
	p.eset <- new("AnnotatedDataFrame", data=pheno.dat)
	eset@phenoData <- p.eset
	p.kegg <- gsealmPerm(eset, as.formula(paste("~ ", cellTypeTag, sep="")), kegg.incidence.matrix, nperm = nperm, na.rm=TRUE)

	rankpmat <- formatrankpathways(p.kegg, kegg.incidence.matrix, direction=direction) 
	
	out <- new("AttractorModuleSet") 
	out@eSet <- eset
	out@incidenceMatrix <- kegg.incidence.matrix
	out@rankedPathways <- rankpmat
	out@cellTypeTag <- cellTypeTag
	return(out)	
}

# eventually convert this to a summary method for the find.attractors output 

formatrankpathways <- function(gseaResMat, incidenceMat, direction="upper"){
	#require(KEGG.db)
	if( !direction %in% c("upper", "lower") ){
		stop("direction must be either upper or lower")
	}
	if( direction == "upper" ){ pvals <- gseaResMat[,2] }
	else{ pvals <- gseaResMat[,1] }
	pways <- rownames(gseaResMat) 
	pways.names <- unlist(mget(pways, KEGGPATHID2NAME))
	p.size <- apply(incidenceMat, 1, sum)	
	mat <- data.frame(pways[order(pvals)], pways.names[order(pvals)], pvals[order(pvals)], p.size[order(pvals)])
	colnames(mat) <- c("KEGGID", "PATHWAY", "Pvalue", "Size")
	# sort the pathway hits with a Pvalue = 0, by size 
	# maybe expand this secondary ranking to all other pathways with the same Pvalue?
	mat[mat$Pvalue == 0,] <- mat[order(mat[mat$Pvalue == 0,]$Size, decreasing=TRUE),]
	return(mat)
}
		
## internal functions

flagPwayExists <- function(x){ 
	if( length(x) == 1 ){
		if( is.na(x) ){ flag <- FALSE }
		else{ flag <- TRUE }
	}
	else{ flag <- TRUE }
}

buildKeggIncidenceMatrix <- function(kegg.ids, gene.ids, annotation){
	
	require(annotation, character.only=TRUE)
	#require(KEGG.db)
	ann <- strsplit(annotation, ".db")[[1]]
	path2probeEnv <- get(paste(ann, "PATH2PROBE", sep=""))
	pway.genes <- mget(kegg.ids, path2probeEnv)
	
	convert.to.row <- function(x.genes, row.genes){
		res <- integer(length(row.genes)) ; res[row.genes %in% x.genes] <- 1 
		return(res)
	}
	
	xmat <- t(sapply(lapply(pway.genes, convert.to.row, gene.ids), cbind))
	rownames(xmat) <- kegg.ids ; colnames(xmat) <- gene.ids
	return(xmat)
}

buildCorMatrix <- function(dat.fr, module.genes, cor.cutoff){
		second.list <- NULL 
		non.module.genes <- setdiff(rownames(dat.fr), module.genes) 		
		calc.corr.onegene <- function(onegene, dat.fr, non.module.genes, cor.cutoff){
			x <- as.numeric(dat.fr[rownames(dat.fr) %in% onegene,])
			cvals <- cor(x, t(dat.fr[rownames(dat.fr) %in% non.module.genes,]))
			if( sum( cvals >= cor.cutoff ) > 0 ){
				return(non.module.genes[cvals >= cor.cutoff])
			}
		}
		second.list <- unlist(lapply(as.list(module.genes), calc.corr.onegene, dat.fr, non.module.genes, cor.cutoff))
		return(unique(second.list))
}
