library(biomaRt)
library(reshape2)
library(stringr)
library(readr)

##### function to retrieve genes from SNP #####
getENSG = function(rs="rs3043732", mart=mart.snp) {
  results = getBM(attributes=c("refsnp_id","ensembl_gene_stable_id","associated_gene"),
                  filters="snp_filter", values=rs, mart=mart)
  return(results)
}


##### initialize ensembl SNP database #####
mart.snp = useMart(biomart="ENSEMBL_MART_SNP", path="/biomart/martservice",dataset="hsapiens_snp")

##### read unique yeast-human homologies #####
unique_homologies <- read_csv("Documents/PROJECTS/GWAS-enrich/references/ensembl_homologies.csv")

##### read associations #####
associations.all <- read_delim("Documents/PROJECTS/GWAS-enrich/data/raw/traits/prostate_carcinoma.tsv",
                    "\t", escape_double = FALSE, trim_ws = TRUE)

##### drop intergenic variants #####
associations.no.intergenic = associations.all[associations.all$CONTEXT!="intergenic_variant",]

##### initialize results dataframe #####
overlap.results = data.frame(trait="", gwas.studies.total=0, associated.snps.total=0, mapped.human.genes.total=0, mapped.human.genes.in.universe=0,
                             mapped.yeast.genes.in.universe=0, homologies.drawn.by.pca=0, ORFs.drawn.by.pca=0, common.hits=0, hypergeom.test.pval=0, stringsAsFactors=F)

##### read PCA data #####
pca.hits <- read_delim("Documents/PROJECTS/GWAS-enrich/references/hdPCA/hdPCA_metformin_hits.tsv",
                              "\t", escape_double = FALSE, trim_ws = TRUE)

valid_draws = merge(pca.hits,unique_homologies,by.x = "ORF",by.y = "yeast.id")
unique_draws = as.data.frame(unique(valid_draws$ORF))
colnames(unique_draws)="ORF"
Universe = nrow(unique_homologies)

gwas.traits = sort(unique(associations.all$`DISEASE/TRAIT`))

for ( i in 1:length(gwas.traits) ) {
  current.trait = gwas.traits[i]
  overlap.results[i,"trait"] = current.trait
  current.associations = subset(associations.all, `DISEASE/TRAIT` %in% current.trait)
  
  overlap.results[i,"gwas.studies.total"]    = length(unique(current.associations$PUBMEDID))
  overlap.results[i,"associated.snps.total"] = length(unique(current.associations$SNPS))
  
  if(length(unique(current.associations$SNPS))<10){
       next
  }
  
  gwas.mapped.genes = getENSG(rs=current.associations$SNPS)
  
  if(is.null(gwas.mapped.genes)){
    next
  }
  
  gwas.mapped.ids.alt = unlist(sapply(current.associations$SNP_GENE_IDS, function(x) unlist(strsplit(x,","))))
  if ( !is.null(gwas.mapped.ids.alt) ){
    gwas.mapped.ids.alt = setdiff(unique(gwas.mapped.ids.alt),c("NR","NA"))
  }
  
  gwas.mapped.ids.all = gwas.mapped.genes$ensembl_gene_stable_id
  gwas.mapped.ids.all = unlist(sapply(gwas.mapped.ids.all, function(x) unlist(strsplit(x,","))))
  
  if ( !is.null(gwas.mapped.ids.all) ){
      gwas.mapped.ids.all = setdiff(unique(gwas.mapped.ids.all),c("NR","NA"))
  } else {
    next
  }
  
  overlap.results[i,"mapped.human.genes.total"] = length(gwas.mapped.ids.all)
  
  mapped.ensembl.ids <- data.frame(sort(unique(gwas.mapped.ids.all)))
  colnames(mapped.ensembl.ids)="id"
  
  gwas.mapped.ids.in.universe = merge.data.frame(x = mapped.ensembl.ids, y = unique_homologies,by.x = "id",by.y = "human.id")
  
  colnames(gwas.mapped.ids.in.universe)[1]="human.id"
  
  overlap.results[i,"mapped.human.genes.in.universe"] = length(unique(gwas.mapped.ids.in.universe$human.id))
  
  overlap.results[i,"mapped.yeast.genes.in.universe"] = length(unique(gwas.mapped.ids.in.universe$yeast.id))
  
  common.hits = merge.data.frame(x = gwas.mapped.ids.in.universe,y = pca.hits,by.x = "yeast.id",by.y = "ORF")
  
  nMarked = nrow(gwas.mapped.ids.in.universe)
  nDrawn = nrow(valid_draws)
  
  overlap.results[i,"homologies.drawn.by.pca"] = nDrawn
  
  overlap.results[i,"ORFs.drawn.by.pca"] = nrow(unique_draws)
  
  nFound = nrow(common.hits)
  
  pval = phyper(nFound,nMarked,Universe-nMarked,nDrawn,lower.tail = FALSE)
  
  overlap.results[i,"common.hits"]  = ifelse(nFound==0, NA, nFound)
  overlap.results[i,"hypergeom.test.pval"] = ifelse(nFound==0, NA, pval)
}

write_csv(intersection,path = "Documents/PROJECTS/GWAS-enrich/data/processed/traits/PR_interval_homologies.csv",col_names = TRUE,quote_escape = FALSE)

P = function(N){
  nFound <- 15
  nDrawn <- 370
  nMarked <- 72
  nNotMarked <- 2696 - nMarked
  phyper(nFound,nMarked,nNotMarked,nDrawn,lower.tail = FALSE)
}