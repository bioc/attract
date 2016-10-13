
###########################################################################
# Step 4 - functional enrichment for the synexpression groups
###########################################################################

calcFuncSynexprs <- function(mySynExpressionSet, myAttractorModuleSet, ontology="BP", min.pvalue=.05, min.pwaysize=5, annotation="illuminaHumanv2.db", analysis="microarray", expressionSetGeneFormat=NULL, ...){
    #require(GOstats)
    if (analysis=="RNAseq") {
        require(annotation, character.only=TRUE)
        loadNamespace(annotation)
        envPos <- match(paste("package:", annotation, sep=""), search())
        if(expressionSetGeneFormat == "ENTREZID") {
            gene.uni <- rownames(exprs(myAttractorModuleSet@eSet))
            entrez.uni <- gene.uni
        } else {
            gene.uni <- select(get(annotation,envPos, as.environment(envPos)), keys=rownames(exprs(myAttractorModuleSet@eSet)), keytype=expressionSetGeneFormat, columns=c(expressionSetGeneFormat,"ENTREZID") )
            entrez.uni <- gene.uni$ENTREZID
            entrez.uni <- unique(entrez.uni)
            entrez.uni <- entrez.uni[!is.na(entrez.uni)]
            entrezEnv <- NULL
        }
    } else {
        ann <- strsplit(annotation, ".db")[[1]]
        entrezEnv <- get(paste(ann, "ENTREZID", sep=""))
        gene.uni <- rownames(exprs(myAttractorModuleSet@eSet))
        entrez.uni <- unique(unlist(mget(intersect(gene.uni, ls(entrezEnv)), entrezEnv)))
        entrez.uni <- entrez.uni[!is.na(entrez.uni)]
    }
    calc.func.onesyngroup <- function(onesyngroup, entrezEnv, entrez.uni, annotation, ontology, analysis,expressionSetGeneFormat){
        if (analysis=="RNAseq") {
            require(annotation, character.only=TRUE)
            loadNamespace(annotation)
            envPos <- match(paste("package:", annotation, sep=""), search())
            if(expressionSetGeneFormat == "ENTREZID") {
                ensemblToEntrez <- onesyngroup
                entrez.mgenes <- ensemblToEntrez
            } else {
                ensemblToEntrez <- select(get(annotation,envPos, as.environment(envPos)), keys=onesyngroup, keytype=expressionSetGeneFormat, columns=c(expressionSetGeneFormat,"ENTREZID") )
                entrez.mgenes <- unique(ensemblToEntrez$ENTREZID)
                entrez.mgenes <- entrez.mgenes[!is.na(entrez.mgenes)]
            }
        } else {
            entrez.mgenes <- unique(unlist(mget(intersect(onesyngroup, ls(entrezEnv)), entrezEnv)))
            entrez.mgenes <- entrez.mgenes[!is.na(entrez.mgenes)]
        }
        if( ontology != "KEGG" ){
            hyperObj <- new("GOHyperGParams", geneIds=entrez.mgenes, universeGeneIds=entrez.uni, annotation=annotation,
            ontology=ontology, pvalueCutoff=1)
            xtab <- summary(hyperGTest(hyperObj))
            apval <- p.adjust(xtab[,2], "BH")
            keep.rows <- apval <= min.pvalue & xtab[,5] >= min.pwaysize
            if( length(keep.rows) > 0 ){
                xtab <- data.frame(xtab[keep.rows,], AdjPval=apval[keep.rows])
            }
            else{
                xtab <- NA
            }
        }
        
        return(xtab)
    }
    return(lapply(mySynExpressionSet@groups, calc.func.onesyngroup, entrezEnv, entrez.uni, annotation, ontology, analysis,expressionSetGeneFormat))
}

findCorrPartners <- function(mySynExpressionSet, myEset, removeGenes=NULL, cor.cutoff=0.85, ...){ 
	dat.fr <- exprs(myEset)
	dat.fr <- dat.fr[!rownames(dat.fr) %in% removeGenes,]
	calc.onegroup <- function(onegroup, removeGenes, dat.fr, cor.cutoff){
		onegroup <- onegroup[!onegroup %in% removeGenes]
		cor.genes <- buildCorMatrix(dat.fr, onegroup, cor.cutoff=cor.cutoff)
		return(cor.genes)
	}
	corr.synexprs <- lapply(mySynExpressionSet@groups, calc.onegroup, removeGenes, dat.fr, cor.cutoff)
	avgmat <- lapply(corr.synexprs, function(xgenes, exprsdat){ apply(exprsdat[rownames(exprsdat) %in% xgenes,], 2, mean, na.rm=TRUE) },
						exprsdat=dat.fr)

	res <- new("SynExpressionSet")
	res@groups <- corr.synexprs
	res@profiles <- t(sapply(avgmat, rbind))
	return(res)
}

