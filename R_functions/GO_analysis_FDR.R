GO_analysis_FDR <- function(q, b, c) {
    
    background.genes <- q$gene_name
    geneUniverse <- background.genes
    
    if(b=="upregulated"){
        test.genes <- filter(q, FDR < 0.05, logFC > 0)$gene_name
    } else if (b=="downregulated"){
        test.genes <- filter(q, FDR < 0.05, logFC < 0)$gene_name
    } else {
        print("Incorect fold change parameter")
        stop()
    }
    
    genesOfInterest <- test.genes
    geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
    names(geneList) <- geneUniverse
    myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot=annFUN.org,    mapping="org.Mm.eg.db", ID = "alias", nodeSize=20)
    print(myGOdata)
    
    resultFisher <- runTest(myGOdata, algorithm = "weight01", statistic = "fisher")
    #resultKS <- runTest(myGOdata, algorithm = "classic", statistic = "ks")
    #resultKS.elim <- runTest(myGOdata, algorithm = "elim", statistic = "ks")
    #classicKS = resultKS, elimKS = resultKS.elim, - add later
    
    allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "classicFisher", topNodes = c)
    
    #Calculate enrichment
    padding = 5
    allRes$Enrichment <- log2((allRes$Significant + padding)/(allRes$Expected + padding))
    
    #showSigOfNodes(myGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')
    #nodes_plot <- recordPlot(showSigOfNodes(myGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all'))
    
    #Building a df of DE genes belonging to top 20 GO BP caegories
    DE_genes_in_top_GO_cat <- function(r){
        fisher.go <- allRes[r,1]
        #print(allRes[x,c(1,2)])
        fisher.ann.genes <- genesInTerm(myGOdata, whichGO=fisher.go)
        df <- data.frame(GO.ID = allRes[r,c(1)], Term = allRes[r,c(2)], gene_name=intersect(as.character(fisher.ann.genes[[1]]), q$gene_name))
        df <- filter(df, gene_name %in% test.genes)
        df
    }
    
    list(allRes, lapply(1:20, DE_genes_in_top_GO_cat))
}
