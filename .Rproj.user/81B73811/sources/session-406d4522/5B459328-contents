install.packages("enrichR")
library(enrichR)
library(dplyr)
library(tibble)
library(ggplot2)
dbs <- listEnrichrDbs()
drug_dbs <- dbs[grep("Drug|LINCS|pert", dbs$libraryName, ignore.case = T), ]
selected_dlibs <- c(
  "LINCS_L1000_Chem_Pert_up",
  "LINCS_L1000_Chem_Pert_down",
  "Drug_Perturbations_from_GEO_up",
  "Drug_Perturbations_from_GEO_down",
  "NIBR_DRUGseq_2025_up",
  "NIBR_DRUGseq_2025_down",
  "DrugMatrix",
  "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO",
  "DGIdb_Drug_Targets_2024"
)
up_genes <- c("MMP9", "KDM6B", "JUN", "FOS", "HIF1A", "NFKBIA") 
down_genes <- c("TP53", "CDKN1A", "BCL2", "GADD45A", "ATM", "CDK1")

up_enrich <- enrichr(up_genes, selected_dlibs)
down_enrich <- enrichr(down_genes, selected_dlibs)

head(up_enrich[["LINCS_L1000_Chem_Pert_down"]], 10)
head(down_enrich[["LINCS_L1000_Chem_Pert_up"]], 10)


head(up_enrich[["Drug_Perturbations_from_GEO_down"]], 10)
head(down_enrich[["Drug_Perturbations_from_GEO_up"]], 10)

top_drugs <- up_enrich[["LINCS_L1000_Chem_Pert_down"]]
top_drugs <- top_drugs[order(top_drugs$Adjusted.P.value), ]
head(top_drugs, 5)

clean_enrichr <- function(
    enrichr_df, label, top_n = 15) {
  enrichr_df %>% 
    as_tibble() %>%
    select(Term, Overlap, Adjusted.P.value, Combined.Score, Genes) %>%
    rename_with(~paste0(label, "_", .), -Term) %>%
    arrange(!!sym(paste0(label, "_Adjusted.P.value"))) %>%
  slice_head(n = top_n)
}

top_up_drugs <- clean_enrichr(up_enrich[["Drug_Perturbations_from_GEO_down"]], "Up")
top_up_drugs_lincs <- clean_enrichr(up_enrich[["LINCS_L1000_Chem_Pert_down"]], "Up")
top_down_drugs_lincs <- clean_enrichr(down_enrich[["LINCS_L1000_Chem_Pert_up"]], "Down")

top_down <- down_enrich[["LINCS_L1000_Chem_Pert_up"]][1:10, ]
ggplot(top_down, aes(x=reorder(Term, Combined.Score), y=Combined.Score)) +
  geom_col(fill="tomato") +
  coord_flip() +
  labs(title="Top Drugs Upregulating Downregulated Genes", x="Drug", y="Combined Score") +
  theme_minimal()
