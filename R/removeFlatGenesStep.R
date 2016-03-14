############################################
# Step 2 - remove flat genes
############################################

removeFlatGenes <- function(eSet, cellTypeTag, contrasts=NULL, limma.cutoff=0.05, ...){
	dat.fr <- exprs(eSet)
	class.vector <- as.factor(pData(eSet)[,colnames(pData(eSet)) %in% cellTypeTag])

#require(limma)
	design <- model.matrix(~0+class.vector)
	lev <- levels(class.vector)
	colnames(design) <- lev

	fit <- lmFit(dat.fr, design)
	if( is.null(contrasts) ){
		my.contrasts <- NULL 
		
		for( i in 1:(length(lev)-1) ){
			my.contrasts <- c(my.contrasts, paste(lev[i], lev[i+1], sep=" - "))
		}
		cont.diff <- makeContrasts(contrasts=my.contrasts, levels=design)
	}
	else{
		my.contrasts <- contrasts
		# should probably check that if contrasts is supplied, it's a valid contrast vector
		cont.diff <- makeContrasts(contrasts=contrasts, levels=lev)
	}
	
	fit2 <- contrasts.fit(fit, cont.diff)
	fit2 <- eBayes(fit2)
	fit2$genes <- rownames(dat.fr)
	tt.all <- topTable(fit2, coef=(1:(length(my.contrasts))), adjust="fdr", n=nrow(dat.fr))

	remove.genes <- tt.all[tt.all$adj.P.Val > limma.cutoff,1]
	return(remove.genes) 
}

