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

findAttractors <- function(myEset, cellTypeTag, min.pwaysize=5, annotation="illuminaHumanv2BeadID.db", ...){
	
	require(KEGG.db)
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

	new.order <- order(class.vector, colnames(dat.detect.wkegg))
	dat.detect.wkegg <- dat.detect.wkegg[,new.order] 
	class.vector <- class.vector[new.order]
	
	fstat <- apply(dat.detect.wkegg, 1, function(y,x){ anova(lm(y ~ x))[[4]][1] }, x=class.vector)
	fstat <- log(fstat, 2)
			
	evalPway <- function(index, global){
		pway.vals <- global[index==1]
		t.test(pway.vals, global)$p.value 
	}
	
	t.pvals <- apply(kegg.incidence.matrix, 1, evalPway, global=fstat)
	t.pvals <- p.adjust(t.pvals, "BH") 

	size <- apply(kegg.incidence.matrix, 1, sum)
	
	tab <- data.frame(KEGGID = rownames(kegg.incidence.matrix), KEGGNAME = unlist(mget(rownames(kegg.incidence.matrix), KEGGPATHID2NAME)), AdjustedPvalues = t.pvals, NumberDetectedGenes = size)
	tab <- tab[order(t.pvals),]
	
	# making the eset with the current data 
	eset <- new("ExpressionSet")
	eset@assayData <- new.env()
	assign("exprs", dat.detect.wkegg, eset@assayData)

	pheno.dat <- data.frame(colnames(dat.detect.wkegg), class.vector)
	colnames(pheno.dat) <- c("ChipID", cellTypeTag)
	p.eset <- new("AnnotatedDataFrame", data=pheno.dat)
	eset@phenoData <- p.eset
	
	out <- new("AttractorModuleSet") 
	out@eSet <- eset
	out@incidenceMatrix <- kegg.incidence.matrix
	out@rankedPathways <- tab
	out@cellTypeTag <- cellTypeTag
	return(out)	
}

# eventually convert this to a summary method for the find.attractors output 

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
	require(KEGG.db)
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
