.libPaths("C:/Users/diala/Documents/R/win-library/3.2")


root_path = "C:/Firelite/Projects/PCA_Metformin_1"
setwd(root_path)

in_path   = "./Input_Data/"
code_path = "./Code/R/"
out_path  = "./Results/norm_9_no_pval_cutoff/"

library(biomaRt)
mart.snp = useMart(biomart="ENSEMBL_MART_SNP", host="mar2016.archive.ensembl.org", path="/biomart/martservice",
                   dataset="hsapiens_snp")


getENSG = function(rs="rs3043732", mart=mart.snp) {
    results = getBM(attributes=c("refsnp_id","ensembl_gene_stable_id","associated_gene"),
                    filters="snp_filter", values=rs, mart=mart)
    return(results)
}

overlap.results = data.frame(disease="", gwas.studies.no=0, associated.snps.no=0, hs.mapped.genes.no=0, homolog.orfs.no=0,
                             tested.orfs.no.in.pca=0, common.orfs.no=0, hyper.test.pval=0, stringsAsFactors=F)

DF           = read.csv(paste(out_path,"CTRL and MTF Total.csv",sep=""), sep=",", stringsAsFactors=F)
mtf.hits     = unique(subset(DF, pVal<0.01)$ORF)
gwas.studies = read.table(paste(in_path,"gwas_catalog_cancer.txt",sep=""), header=T, sep="\t", quote="\"",stringsAsFactors=F)
diseases     = sort(unique(gwas.studies$disease.trait))

for ( i in 1:length(diseases) ) {
    disease = diseases[i]
    overlap.results[i,"disease"] = disease
    gwas.study = subset(gwas.studies, disease.trait %in% disease)
    overlap.results[i,"gwas.studies.no"]    = length(unique(gwas.study$pubmedid))
    overlap.results[i,"associated.snps.no"] = length(unique(gwas.study$snps))

    gwas.mapped.genes = getENSG(rs=gwas.study$snps)
    gwas.mapped.ids = gwas.mapped.genes$ensembl_gene_stable_id
    gwas.mapped.genes = gwas.mapped.genes$associated_gene
    gwas.mapped.ids = tryCatch(unlist(sapply(gwas.mapped.ids, function(x) unlist(strsplit(x,",")))),
                                 error=function(cond) { return(c()) } )
    gwas.mapped.genes = tryCatch(unlist(sapply(gwas.mapped.genes, function(x) unlist(strsplit(x,",")))),
                                 error=function(cond) { return(c()) } )
    
    if ( !is.null(gwas.mapped.genes) ){
        gwas.mapped.genes = setdiff(unique(gwas.mapped.genes),c("NR","NA"))
    }
    
    if ( !is.null(gwas.mapped.ids) ){
      gwas.mapped.ids = setdiff(unique(gwas.mapped.ids),c("NR","NA"))
    }
    
    overlap.results[i,"hs.mapped.genes.no"] = length(gwas.mapped.genes)
    overlap.results[i,"hs.mapped.ids.no"] = length(gwas.mapped.ids)

    human = useMart(host='mar2016.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
    yeast = useMart(host='mar2016.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='scerevisiae_gene_ensembl')
    yeast.homolog = tryCatch(getLDS(attributes=c("hgnc_symbol","chromosome_name", "start_position"), filters="hgnc_symbol",
                           values=gwas.mapped.genes, mart=human, attributesL=c("ensembl_peptide_id", "chromosome_name",
                           "start_position"), martL=yeast), error=function(cond) { return(c()) } )
    
    yeast.homolog.by.id = tryCatch(getLDS(attributes=c("hgnc_symbol","chromosome_name", "start_position"), filters="ensembl_gene_id",
                                    values=gwas.mapped.ids, mart=human, attributesL=c("ensembl_peptide_id", "chromosome_name",
                                                                                        "start_position"), martL=yeast), error=function(cond) { return(c()) } )

    orf.gwas = unique(yeast.homolog$Ensembl.Protein.ID)
    orf.gwas.by.id = unique(yeast.homolog.by.id$Ensembl.Protein.ID)
    overlap.results[i,"homolog.orfs.no"] = length(orf.gwas)
    overlap.results[i,"homolog.orfs.by.id.no"] = length(orf.gwas.by.id)
    orf.gwas = intersect(orf.gwas, DF$ORF)
    orf.gwas.by.id = intersect(orf.gwas.by.id, DF$ORF)
    overlap.results[i,"tested.orfs.no.in.pca"] = ifelse(length(orf.gwas)==0, NA, length(orf.gwas))
    overlap.results[i,"tested.orfs.by.id.no.in.pca"] = ifelse(length(orf.gwas.by.id)==0, NA, length(orf.gwas.by.id))

    univ = union(orf.gwas, DF$ORF)
    q = length(intersect(orf.gwas, mtf.hits))
    m = length(mtf.hits)
    n = length(univ)-m
    k = length(orf.gwas)
    pval = phyper(q-1,m,n,k,lower.tail=F)

    overlap.results[i,"common.orfs.no"]  = ifelse(length(orf.gwas)==0, NA, q)
    overlap.results[i,"hyper.test.pval"] = ifelse(length(orf.gwas)==0, NA, pval)
}


i = nrow(overlap.results)+1
gwas.studies = read.table(paste(in_path,"gwas_catalog_t2d.txt",sep=""), header=T, sep="\t", quote="\"",stringsAsFactors=F)
diseases     = sort(unique(gwas.studies$disease.trait))[c(1,4)]

gwas.study = subset(gwas.studies, disease.trait %in% diseases)
overlap.results[i,"disease"] = paste(diseases, collapse=" | ")
overlap.results[i,"gwas.studies.no"]    = length(unique(gwas.study$pubmedid))
overlap.results[i,"associated.snps.no"] = length(unique(gwas.study$snps))

gwas.mapped.genes = getENSG(rs=gwas.study$snps)
gwas.mapped.genes = gwas.mapped.genes$associated_gene
gwas.mapped.genes = tryCatch(unlist(sapply(gwas.mapped.genes, function(x) unlist(strsplit(x,",")))),
                                 error=function(cond) { return(c()) } )
if ( !is.null(gwas.mapped.genes) )
    gwas.mapped.genes = setdiff(unique(gwas.mapped.genes),c("NR","NA"))
overlap.results[i,"hs.mapped.genes.no"] = length(gwas.mapped.genes)

human = useMart(host='mar2016.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
yeast = useMart(host='mar2016.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='scerevisiae_gene_ensembl')
yeast.homolog = tryCatch(getLDS(attributes=c("hgnc_symbol","chromosome_name", "start_position"), filters="hgnc_symbol",
                         values=gwas.mapped.genes, mart=human, attributesL=c("ensembl_peptide_id", "chromosome_name",
                         "start_position"), martL=yeast), error=function(cond) { return(c()) } )

orf.gwas = unique(yeast.homolog$Ensembl.Protein.ID)
overlap.results[i,"homolog.orfs.no"] = length(orf.gwas)
orf.gwas = intersect(orf.gwas, DF$ORF)
overlap.results[i,"tested.orfs.no.in.pca"] = ifelse(length(orf.gwas)==0, NA, length(orf.gwas))

univ = union(orf.gwas, DF$ORF)
q = length(intersect(orf.gwas, mtf.hits))
m = length(mtf.hits)
n = length(univ)-m
k = length(orf.gwas)
pval = phyper(q-1,m,n,k,lower.tail=F)

overlap.results[i,"common.orfs.no"]  = ifelse(length(orf.gwas)==0, NA, q)
overlap.results[i,"hyper.test.pval"] = ifelse(length(orf.gwas)==0, NA, pval)

write.table(overlap.results, file=paste(out_path,"gwas.pca.enrichment.txt",sep=""), sep="\t", quote=F, row.names=F)


