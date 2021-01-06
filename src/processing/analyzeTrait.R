library(biomaRt)
library(reshape2)
library(ggplot2)
library(scales)
library(stringr)

unique_homologies <- read_csv("Documents/PROJECTS/GWAS-enrich/references/unique_homologies.csv")

trait <- read_delim("Documents/PROJECTS/GWAS-enrich/data/raw/traits/PR_interval.tsv",
"\t", escape_double = FALSE, col_types = cols(CHR_ID = col_skip(),
CHR_POS = col_skip(), CNV = col_skip(),
DATE = col_skip(), DOWNSTREAM_GENE_DISTANCE = col_skip(),
`FIRST AUTHOR` = col_skip(), `GENOTYPING TECHNOLOGY` = col_skip(),
`INITIAL SAMPLE SIZE` = col_skip(),
JOURNAL = col_skip(), LINK = col_skip(),
MAPPED_TRAIT_URI = col_skip(), MERGED = col_skip(),
`P-VALUE (TEXT)` = col_skip(), `PLATFORM [SNPS PASSING QC]` = col_skip(),
PUBMEDID = col_skip(), PVALUE_MLOG = col_skip(),
REGION = col_skip(), `REPLICATION SAMPLE SIZE` = col_skip(),
STUDY = col_skip(), `STUDY ACCESSION` = col_skip(),
UPSTREAM_GENE_DISTANCE = col_skip()),
trim_ws = TRUE)

# drop intergenic variants
trait = trait[trait$CONTEXT!="intergenic_variant",]

# select only matching trait
trait = trait[trait$`DISEASE/TRAIT`=="PR interval",]

trait = trait[rowSums(is.na(trait)) != ncol(trait), ]

ntotal_associations = nrow(trait)
nunique_associations = length(unique(trait$SNPS))

hit_genes = trait$SNP_GENE_IDS

hit_genes = data.frame(hit_genes,do.call(rbind,str_split(hit_genes,", ")))

hit_genes = hit_genes[rowSums(is.na(hit_genes)) != ncol(hit_genes), ]

df = melt(hit_genes,id.vars = "hit_genes")

GWAS_genes = as.data.frame(unique(df$value))
colnames(GWAS_genes)="ID"

GWAS_set = merge(GWAS_genes,unique_homologies,by.x = "ID",by.y = "human_id")

nMarked = nrow(GWAS_set)

#size = nrow(unique_homologies)

######

#intersection with a dataset

hdPCA_metformin <- read_delim("Documents/PROJECTS/GWAS-enrich/references/hdPCA/hdPCA_metformin_hits.tsv",
"\t", escape_double = FALSE, trim_ws = TRUE)

valid_draws = merge(hdPCA_metformin,unique_homologies,by.x = "ORF",by.y = "yeast_id")

unique_draws = as.data.frame(unique(valid_draws$ORF))
colnames(unique_draws)="ORF"
#nDrawn = length(unique_draws$ORF)

nDrawn = nrow(valid_draws)

intersection = merge(GWAS_set,unique_draws,by.x = "yeast_id","ORF")

write_csv(intersection,path = "Documents/PROJECTS/GWAS-enrich/data/processed/traits/PR_interval_homologies.csv",col_names = TRUE,quote_escape = FALSE)

nFound = nrow(intersection)

#size = 2696
#size = 4454

phyper(nFound,nMarked,2696-nMarked,nDrawn,lower.tail = FALSE)

phyper(nFound,nMarked,4454-nMarked,nDrawn,lower.tail = FALSE)


P = function(N){
  nFound <- 15
  nDrawn <- 370
  nMarked <- 72
  nNotMarked <- 2696 - nMarked
  phyper(nFound,nMarked,nNotMarked,nDrawn,lower.tail = FALSE)
}

