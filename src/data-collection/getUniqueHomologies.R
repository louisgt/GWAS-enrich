library(biomaRt)

mart.hs <- useMart("ensembl", "hsapiens_gene_ensembl")
mart.sc <- useMart("ensembl", "scerevisiae_gene_ensembl")
mart.snp = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", path="/biomart/martservice",
                   dataset="hsapiens_snp")

human.homolog = tryCatch(getLDS(attributes=c("chromosome_name","ensembl_gene_id"), filters="ensembl_gene_id",
values=hdPCA_list$ID, mart=mart.sc, attributesL=c("chromosome_name","hgnc_symbol","ensembl_gene_id","scerevisiae_homolog_orthology_confidence"), martL=mart.hs), error=function(cond) { return(c()) } )

orf.gwas = as.data.frame(unique(human.homolog$Gene.stable.ID.1))
homologies = human.homolog[,c(2,4,5)]
colnames(homologies)=c("yeast_id","human_name","human_id")
unique_homologies = homologies[!duplicated(t(apply(homologies,1,sort))),]
df = merge(unique_homologies,hdPCA_list,by.x = "yeast_id",by.y = "X3")
df$X1=NULL
df2 = df[,c(4,1,2,3)]
colnames(df2)[1]="yeast_name"
write_csv(df2,path = "Documents/PROJECTS/GWAS-enrich/references/unique_homologies.csv",col_names = TRUE,quote_escape = FALSE)
