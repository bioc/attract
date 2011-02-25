# attract-specific classes

setClass("AttractorModuleSet", representation(eSet="ExpressionSet", 
			cellTypeTag="character", incidenceMatrix="matrix", rankedPathways="data.frame"))

setClass("SynExpressionSet", representation(groups="list", profiles="matrix"))
