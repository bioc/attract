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

#require("AnnotationDbi")
#require(org.Hs.eg.db)
#require("KEGGREST")
#library(reactome.db)

buildCustomIncidenceMatrix <- function(geneSetFrame, geneNames) {
    
    pwayToGeneList <- apply(geneSetFrame, 1, function(x) x[x!=""]) #pathway to gene List
    
    convert.to.row <- function(x.genes, row.genes){ # x.genes are list of pathways that has each gene. row.genes are all genes that have a kegg annotation
        res <- integer(length(row.genes)) ; res[row.genes %in% x.genes] <- 1 #see which probes are in the pathway
        return(res)
    }
    
    xmat <- t(sapply(lapply(pwayToGeneList, convert.to.row, geneNames), cbind))
    colnames(xmat) <- geneNames
    return(xmat)
}

filterDataSet <- function(data,filterPerc=0.75){ #only necessary for RNAseq data
    keep.rows <- apply(data, 1, function(x) sum(x==0)/length(x)) < filterPerc
    data <- data[keep.rows,]
    data <- log(data +1,2)
    return(data)
}
findAttractors <- function(myEset, cellTypeTag, min.pwaysize=5, annotation="illuminaHumanv2.db", database="KEGG", analysis="microarray", databaseGeneFormat=NULL, ...){

    require(annotation, character.only=TRUE)
	ann <- strsplit(annotation, ".db")[[1]]
	loadNamespace(annotation)
    envPos <- match(paste("package:", annotation, sep=""), search())
	dat.fr <- exprs(myEset)
    #keep.rows <- apply(dat.fr, 1, sum) >0
    #dat.fr <- dat.fr[keep.rows,] # remove rows where sum of expression data is 0
	all.probes <- rownames(dat.fr)
	dat.fr <- as.matrix(dat.fr)

	class.vector <- as.factor(pData(myEset)[,colnames(pData(myEset)) %in% cellTypeTag])
	
    #do if using a custom database
    if(database!="KEGG" & database!= "reactome") {
        # create data frame for custom geneSet and vector of all the genes
        maxColumns <- max(count.fields(database, sep="\t"))
        geneSetFrame <- read.table(database,header=FALSE,fill=TRUE,col.names=1:maxColumns,sep="\t",quote="",colClasses = "character")
        rownames(geneSetFrame) <- geneSetFrame[,1]
        geneSetFrame <- geneSetFrame[,c(-1,-2)]
        geneSetFrame <- as.matrix(geneSetFrame)
        geneNames <- unique(as.character(geneSetFrame))
        quoteIndex <- match("",geneNames)
        geneNames <- geneNames[-quoteIndex] # vector of all genes
        custom.incidence.matrix <- buildCustomIncidenceMatrix(geneSetFrame, geneNames) # create incidence matrix
        keep.pways <- apply(custom.incidence.matrix, 1, sum) >= min.pwaysize
        custom.incidence.matrix <- custom.incidence.matrix[keep.pways,] # filter incidence matrix
        #create list.wpway: the genes from your expression data actually in the pathways
        if(analysis=="RNAseq") {
            #convert genes to ENSEMBL IDs
            genestoENSEMBL <- select(get(annotation,envPos, as.environment(envPos)), keys=geneNames, keytype=databaseGeneFormat, columns=c(databaseGeneFormat,"ENSEMBL") )
            geneSettoENSEMBL <- unique(genestoENSEMBL$ENSEMBL) #some genes have same ensembl IDS
            list.wpway <- sapply(rownames(dat.fr),function(x) x%in%geneSettoENSEMBL)
            list.wpway <- names(which(list.wpway))
        } else { # convert genes to probes
            genestoProbe <- select(annotation, keys=geneNames, keytype=databaseGeneFormat, columns=c(databaseGeneFormat,"PROBEID")) # annotation must be package name. not str
            geneSettoProbe <- genestoProbe$PROBEID
            list.wpway <- sapply(rownames(dat.fr),function(x) x%in%geneSettoProbe)
            list.wpway <- names(which(list.wpway))
        }
        
        dat.detect.wkegg <- dat.fr[rownames(dat.fr) %in% list.wpway,] # only keep genes actually in a pathway
        
        new.order <- order(class.vector, colnames(dat.detect.wkegg))
        dat.detect.wkegg <- dat.detect.wkegg[,new.order]
        class.vector <- class.vector[new.order]
        
        fstat <- apply(dat.detect.wkegg, 1, function(y,x){ anova(lm(y ~ x))[[4]][1] }, x=class.vector)
        fstat <- log(fstat, 2)
        
        evalPway <- function(index, global){
            pway.vals <- global[index==1] # there are NAs where there are ensembl genes in the incidence matrix not in fstat
            #print("NEXT")
            #print(pway.vals)
            t.test(pway.vals, global)$p.value
        }
        
        t.pvals <- apply(custom.incidence.matrix, 1, evalPway, global=fstat)
        t.pvals <- p.adjust(t.pvals, "BH")
        
        size <- apply(custom.incidence.matrix, 1, sum)
        tab <- data.frame(pathID = rownames(custom.incidence.matrix), pathNAME = rownames(custom.incidence.matrix), AdjustedPvalues = t.pvals, NumberDetectedGenes = size)
        tab <- tab[order(t.pvals),]
        
        eset <- new("ExpressionSet")
        eset@assayData <- new.env()
        assign("exprs", dat.detect.wkegg, eset@assayData)
        
        pheno.dat <- data.frame(colnames(dat.detect.wkegg), class.vector)
        colnames(pheno.dat) <- c("ChipID", cellTypeTag)
        p.eset <- new("AnnotatedDataFrame", data=pheno.dat)
        eset@phenoData <- p.eset
        
        out <- new("AttractorModuleSet")
        out@eSet <- eset
        out@incidenceMatrix <- custom.incidence.matrix
        out@rankedPathways <- tab
        out@cellTypeTag <- cellTypeTag
        return(out)

    }
    
    #checking reactome annotations
    if( database=="reactome") {
        if( analysis=="microarray") {
            myEnv <- get(paste(ann, "ENTREZID", sep=""), envPos, as.environment(envPos)) # probe to gene env
            gene.hits <- mget(intersect(all.probes, ls(myEnv)), myEnv)
        } else if (analysis == "RNAseq") { #RNAseq option only supports ensembl genes as of now
            #myEnv <- get(paste(ann, "ENSEMBL", sep=""), envPos, as.environment(envPos)) # ensembl to gene env
            #gene.hits <- mget(intersect(all.probes, ls(org.Hs.egENSEMBL)), org.Hs.egENSEMBL)
            ensemblToEntrez <- select(get(annotation,envPos, as.environment(envPos)), keys=all.probes, keytype="ENSEMBL", columns=c("ENSEMBL","ENTREZID") )
            #combinedRows <- aggregate(ensmblToEntrez[,2],FUN=identity, by=list(ensmblToEntrez[,1]))
            #names(combinedRows$x) <- combinedRows$Group.1
            gene.hits <- ensemblToEntrez$ENTREZID
            #gene.hits <- as.list(combinedRows$x)
        }
        gene.hits <- unique(as.numeric(gene.hits)[!is.na(as.numeric(gene.hits))])
        path.hits <- mget(intersect(gene.hits, ls(reactomeEXTID2PATHID)),reactomeEXTID2PATHID) #get all pathways
        list.wpway <- sort(names(path.hits)[unlist(sapply(path.hits, flagPwayExists))]) # only get entrez genes that are in pathways
        all.pways <- unlist(mget(list.wpway, reactomeEXTID2PATHID))
        all.pways <- unique(all.pways)
    } else if( database=="KEGG") { # checking KEGG annotations
        if( analysis=="microarray") {
            pathEnv <- get(paste(ann, "PATH", sep=""), envPos, as.environment(envPos)) #####
            path.hits <- mget(intersect(all.probes, ls(pathEnv)), pathEnv) #get probes to pathways list of our probes
            list.wpway <- sort(names(path.hits)[unlist(sapply(path.hits, flagPwayExists))])
            all.pways <- unlist(mget(list.wpway, pathEnv))
            all.pways <- unique(all.pways)
        } else if (analysis=="RNAseq") {
            #ensembl=useMart("ensembl")
            #ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
            #ensmblToEntrez <- getBM(filters= "ensembl_gene_id",attributes=c('ensembl_gene_id','entrezgene'),values=all.probes,mart=ensembl)
            #ensemblToEntrez <- select(get(annotation,envPos, as.environment(envPos)), keys=all.probes, keytype="ENSEMBL", columns=c("ENSEMBL","ENTREZID") )
            #ENSEMBLEnv <- get(paste(ann, "ENSEMBL", sep=""), envPos, as.environment(envPos))
            #ensemblToEntrez <- select(org.Hs.eg.db, keys=all.probes, keytype="ENSEMBL", columns=c("ENSEMBL","ENTREZID") )
            #ensemblToEntrez <- mget(intersect(all.probes, ls(revmap(ENSEMBLEnv))), revmap(ENSEMBLEnv))
            #list.wGenes <- sort(names(ensemblToEntrez)[unlist(sapply(ensemblToEntrez, flagPwayExists))])
            #all.Genes <- unlist(mget(list.wGenes, revmap(ENSEMBLEnv)))
            ensemblToEntrez <- select(get(annotation,envPos, as.environment(envPos)), keys=all.probes, keytype="ENSEMBL", columns=c("ENSEMBL","ENTREZID") )
            gene.hits <- ensemblToEntrez$ENTREZID
            gene.hits <- unique(as.numeric(gene.hits)[!is.na(as.numeric(gene.hits))]) # vector of unique entrez IDs with NAs taken out
            pathEnv <- get(paste(ann, "PATH", sep=""), envPos, as.environment(envPos))
            path.hits <- mget(intersect(gene.hits, ls(pathEnv)), pathEnv)
            list.wpway <- sort(names(path.hits)[unlist(sapply(path.hits, flagPwayExists))])
            all.pways <- unlist(mget(list.wpway, pathEnv))
            all.pways <- unique(all.pways)
            #combinedRows <- aggregate(ensmblToEntrez[,2],FUN=identity, by=list(ensmblToEntrez[,1]))
            #names(combinedRows$x) <- combinedRows$Group.1
            #gene.hits <- as.list(combinedRows$x)
            #gene.hits <- sapply(all.probes, function(x) ensmblToEntrez[grep(x, ensmblToEntrez[,1]),2], USE.NAMES=TRUE)
            #gene.hits <- gene.hits[lapply(gene.hits,length)>0]
            #gene.hits <- sapply(gene.hits, function(x) x[!is.na(x)])
            #myEnv <- get(paste(ann, "ENSEMBL", sep=""), envPos, as.environment(envPos))
            #ensmblID  <- unique(as.character(org.Hs.egENSEMBL)) # get all ensembl IDs that map to an entrez gene
            #gene.hits <- mget(intersect(all.probes, ls(org.Hs.egENSEMBL)), org.Hs.egENSEMBL) # ensemble to genes

#            entrezIDsToKegg <- keggLink("pathway", "hsa")
#            path.hits <- sapply(paste("hsa:",gene.hits,sep=""), function(x) entrezIDsToKegg[names(entrezIDsToKegg) == x]) # get all kegg pathways that contains above entrez gene
#            list.wpway <- sort(names(path.hits)[unlist(sapply(path.hits, flagPwayExists))]) # remove entrez genes that do not have any pathways
            #list.wpway <- sapply(strsplit(list.wpway,":"), function(x) x[2])
#            all.pways <- sapply(list.wpway, function(x) entrezIDsToKegg[names(entrezIDsToKegg) == x])
#            all.pways <- unlist(all.pways)
#            all.pways <- unique(all.pways)
#            list.wpway <- sapply(strsplit(list.wpway,":"), function(x) x[2]) #remove hsa from entrez gene list
            #remove prefix from the pathway path:hsa
#            all.pways <- sapply(strsplit(all.pways,":"), function(x) x[2])
#            all.pways <- substring(all.pways,4)
            #gene.hits <- unique(as.numeric(gene.hits)[!is.na(as.numeric(gene.hits))])
            #path.hits <- mget(intersect(gene.hits, ls(reactomeEXTID2PATHID)),reactomeEXTID2PATHID) # genes to path
            #list.wpway <- sort(names(path.hits)[unlist(sapply(path.hits, flagPwayExists))]) #sorted entrez IDs
            #all.pways <- select(reactome.db, keys=list.wpway, keytype="ENTREZID", columns=c("ENTREZID","PATHID") )
            #all.pways <- unique(all.pways$PATHID)
            #all.pways <- unlist(mget(list.wpway, reactomeEXTID2PATHID))
            #all.pways <- unique(all.pways)
        }
    }
	# making expression data object for pathway database annotated probes only
    if(analysis=="RNAseq" | database=="reactome") {
        allEntrezGenes <- ensemblToEntrez$ENTREZID
        names(allEntrezGenes) <- ensemblToEntrez$ENSEMBL # you can have many pathways with a single entrez gene.
        dat.detect.wkegg <- dat.fr[unique(names(allEntrezGenes[allEntrezGenes %in% list.wpway])),] # find which entrez IDs are in list.wpway and then get ensembl names of them and expression data
    } else {
        dat.detect.wkegg <- dat.fr[rownames(dat.fr) %in% list.wpway,] # get all probes that are actually in a kegg pathway that are in your data set
    }
    dat.detect.wkegg <- as.matrix(dat.detect.wkegg)

	# make geneset incidience matrix
	kegg.incidence.matrix <- buildKeggIncidenceMatrix(all.pways, rownames(dat.detect.wkegg), annotation, database, analysis, envPos)

	keep.pways <- apply(kegg.incidence.matrix, 1, sum) >= min.pwaysize 
	kegg.incidence.matrix <- kegg.incidence.matrix[keep.pways,]

	new.order <- order(class.vector, colnames(dat.detect.wkegg))
	dat.detect.wkegg <- dat.detect.wkegg[,new.order]
    #keep.rows <- apply(dat.detect.wkegg, 1, sum) >0
    #dat.detect.wkegg <- dat.detect.wkegg[keep.rows,] # remove rows where sum of expression data is 0
	class.vector <- class.vector[new.order]
	
	fstat <- apply(dat.detect.wkegg, 1, function(y,x){ anova(lm(y ~ x))[[4]][1] }, x=class.vector)
	fstat <- log(fstat, 2)
			
	evalPway <- function(index, global){
        pway.vals <- global[index==1] # there are NAs where there are ensembl genes in the incidence matrix not in fstat
        #print("NEXT")
        #print(pway.vals)
		t.test(pway.vals, global)$p.value
	}
	
	t.pvals <- apply(kegg.incidence.matrix, 1, evalPway, global=fstat)
	t.pvals <- p.adjust(t.pvals, "BH") 

	size <- apply(kegg.incidence.matrix, 1, sum)
    #######
    
    ##Put all info with Pvalues in table
    if( database=="KEGG") {
    #Get KEGG pathway IDs to NAMES
    KEGGPATHNAMES <- keggList("pathway")
    names(KEGGPATHNAMES) <- substring(names(KEGGPATHNAMES),9)
    ##get Kegg pathway names that are in the matrix
    namesInMatrix <- KEGGPATHNAMES[sapply(names(KEGGPATHNAMES), function(x) x %in% rownames(kegg.incidence.matrix))]
    correctOrder <- sapply(rownames(kegg.incidence.matrix), function(x) match(x, names(namesInMatrix))) #indexes of where the pathway names are in which rows of incidence matrix
    namesInMatrix <- namesInMatrix[correctOrder] # sort the names of the pathways by the rownames of the kegg incidence matrix
	tab <- data.frame(KEGGID = rownames(kegg.incidence.matrix), KEGGNAME = namesInMatrix, AdjustedPvalues = t.pvals, NumberDetectedGenes = size)
    } else if(database=="reactome") {
        #tab <- data.frame(KEGGID = rownames(kegg.incidence.matrix), KEGGNAME = unlist(mget(rownames(kegg.incidence.matrix), reactomePATHID2NAME)), AdjustedPvalues = t.pvals, NumberDetectedGenes = size)
        tab <- data.frame(reactomeID = rownames(kegg.incidence.matrix), reactomeNAME = select(reactome.db, keys=rownames(kegg.incidence.matrix), keytype="PATHID", columns=c("PATHID","PATHNAME") )$PATHNAME, AdjustedPvalues = t.pvals, NumberDetectedGenes = size)

    }
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
    if( length(x) == 0 ) {
        flag <- FALSE
    }
	else if( length(x) == 1 ){
		if( is.na(x) ){ flag <- FALSE }
		else{ flag <- TRUE }
	}
	else{ flag <- TRUE }
}

buildKeggIncidenceMatrix <- function(kegg.ids, gene.ids, annotation, database, analysis, envPos){
	if( analysis=="RNAseq") {
        if( database=="reactome") {
            pway.genes <- mget(kegg.ids, reactomePATHID2EXTID) #converts reactome pathways back to entrez IDs
        } else {
            require(annotation, character.only=TRUE)
            ann <- strsplit(annotation, ".db")[[1]]
            path2probe <- select(get(annotation,envPos, as.environment(envPos)), keys=kegg.ids, keytype="PATH", columns=c("PATH","ENTREZID") )
            pway.genes <- sapply(kegg.ids, function(x) path2probe[path2probe$PATH==x,2])
            #pway.genes <- path2probe$ENTREZID
            #pway.genes <- unique(pway.genes) # vector of unique entrez IDs with NAs taken out
            #path2probeEnv <- get(paste(ann, "PATH2EG", sep=""))
            #pway.genes <- mget(kegg.ids, path2probeEnv) #a list of pathways that has each probe in the pathway
#            keggEntrezToPway <- keggLink("pathway", "hsa") # get vector where names are entrez IDs and elements are pways
#            pway.genes <- sapply(paste("path:hsa",kegg.ids,sep=""), function(x) names(keggEntrezToPway[keggEntrezToPway == x])) # convert kegg pways I have to entrez gene IDs
            #pway.genes <- unlist(pway.genes)
            #pway.genes <- unique(pway.genes) # all entrez IDs in my kegg pathways
            #pway.genes <- lapply(strsplit(pway.genes,":"), function(x) x[2]) # remove hsa from entrez IDs
#            pwayNames <- sapply(strsplit(names(pway.genes),":"), function(x) x[2])
#            pwayNames <- substring(pwayNames,4)
#            names(pway.genes) <- pwayNames
        }
        # convert.to.row <- function(x.genes, row.genes){
        #    x.genes <- sapply(strsplit(x.genes,":"), function(x) x[2])
        #    entrez2Ensmbl <- select(org.Hs.eg.db, keys=x.genes, keytype="ENTREZID", columns=c("ENTREZID","ENSEMBL") ) #there are entrez IDs that match to several ensembl IDS
        #    x.genes <- entrez2Ensmbl$ENSEMBL
        #    x.genes <- unique(x.genes) # convert entrez genes to ensembl
        #    res <- integer(length(row.genes)) ; res[row.genes %in% x.genes] <- 1
        #    return(res)
        #}
    }
    else {
        require(annotation, character.only=TRUE)
        ann <- strsplit(annotation, ".db")[[1]]
        path2probeEnv <- get(paste(ann, "PATH2PROBE", sep=""))
        pway.genes <- mget(kegg.ids, path2probeEnv) #a list of pathways that has each probe in the pathway
		
    }
    convert.to.row <- function(x.genes, row.genes, envPos){ # x.genes are list of pathways that has each gene. row.genes are all genes that have a kegg annotation
        if (analysis == "RNAseq") {
            #if( database == "KEGG") {
                #x.genes <- sapply(strsplit(x.genes,":"), function(x) x[2])
            #}
            entrez2Ensmbl <- select(get(annotation,envPos, as.environment(envPos)), keys=x.genes, keytype="ENTREZID", columns=c("ENTREZID","ENSEMBL") )
            x.genes <- entrez2Ensmbl$ENSEMBL
            x.genes <- unique(x.genes)
            if(length(which(is.na(x.genes))) > 0) {
                x.genes <- x.genes[-which(is.na(x.genes))]
            }
        }
        res <- integer(length(row.genes)) ; res[row.genes %in% x.genes] <- 1 #see which probes are in the pathway
        return(res)
    }

    xmat <- t(sapply(lapply(pway.genes, convert.to.row, gene.ids, envPos), cbind))
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
