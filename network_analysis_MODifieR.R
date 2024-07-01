#####################################################################
# Network Analysis using Module Identification with MODifieR        #
#####################################################################

# Input: DMC lists for condition/s of interest
# Output: MODifieR module and consensus gene lists

# Author: David Martínez-Enguita (2023)

# 1. Importation and formatting of DMCs
# 2. Generation and storage of MODifieR input
# 3. Module identification with MODifieR
# 4. Consensus of MODifieR modules (frequency)
# 5. Pathway enrichment of modules
# 6. Graph visualization of consensus modules
# 7. Venn diagram of consensus module genes
# 8. Mapping of hyper/hypo methylated probes in consensus modules

# Set library path to update and load required packages
pack_R <- c("MODifieR", "doParallel", "clusterProfiler", "VennDiagram", "riverplot",
            "dplyr", "RColorBrewer", "xlsx", "readxl", "purrr", 
            "IlluminaHumanMethylationEPICanno.ilm10b2.hg19", "ReactomePA",
            "org.Hs.eg.db", "DOSE")
for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

set.seed(777)
main_dir <- "./"
input_dir <- "Data/"
output_dir <- "Results/"

#####################################################################
# 1. Importation and formatting of DMP table                        #
#####################################################################

# Import main methylation results
DMC_Exp_vs_HC <- read.csv(paste0(main_dir, input_dir, "DMC_Exp_vs_HC_padj0.05_0.2.csv"), sep = ";")
DMC_Pat_vs_Exp <- read.csv(paste0(main_dir, input_dir, "DMC_Pat_vs_Exp_padj0.05_0.2.csv"), sep = ";")
DMC_Pat_vs_HC <- read.csv(paste0(main_dir, input_dir, "DMC_Pat_vs_HC_padj0.05_0.2.csv"))

# Annotation auxiliary functions
symbol.to.entrez <- function(x) {
  x <- mget(x, org.Hs.egSYMBOL2EG, ifnotfound = NA)
  if (length(x) == 0) { return("NA") 
  } else {unlist(lapply(x, function(i) { i[1]}))}
}

entrez.to.symbol <- function(x) {
  x <- mget(x, org.Hs.egSYMBOL, ifnotfound = NA)
  if (length(x) == 0) { return("NA") 
  } else {unlist(lapply(x, function(i) { i[1]}))}
}

# Retrieve methylation array annotation
data(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
list_annot_EPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)@listData
annot_EPIC <- as.data.frame(cbind(list_annot_EPIC$Name, list_annot_EPIC$Relation_to_Island, 
                                  list_annot_EPIC$UCSC_RefGene_Name, list_annot_EPIC$UCSC_RefGene_Accession,
                                  list_annot_EPIC$UCSC_RefGene_Group))
colnames(annot_EPIC) <- c("CpG_ID", "Relation_to_Island", "Gene_Symbol", "Gene_Accession", "Gene_Group")

# Format rows with multiple annotations per probe
annot_EPIC <- tidyr::separate_rows(annot_EPIC, Gene_Symbol)
row_index <- c(TRUE, rowSums(tail(annot_EPIC, -1) == head(annot_EPIC, -1)) != ncol(annot_EPIC))
annot_EPIC <- annot_EPIC[row_index, ]

annot_EPIC <- annot_EPIC %>% 
  mutate(key = CpG_ID) %>%
  filter(key != lag(key, default="0")) %>% 
  select(-key)

# Correct deprecated annotation terms
annot_correct <- read.csv(paste0(main_dir, "Data/annotation_correction_table.csv"))

annot_EPIC[annot_EPIC == ""] <- "NA"
annot_EPIC <- annot_EPIC %>% mutate(Gene_Symbol = recode(Gene_Symbol,
                                                         !!!setNames(annot_correct$updated_ID, annot_correct$deprecated_ID)))

# Map to Entrez IDs
annot_EPIC[, "Entrez_ID"] <- symbol.to.entrez(annot_EPIC$Gene_Symbol)

# Format to true NA
annot_EPIC <- annot_EPIC %>%
  mutate(across(everything(), ~ifelse(. == "NA", NA, as.character(.))))
annot_EPIC <- as.data.frame(annot_EPIC)

# Check final output
length(unique(annot_EPIC$CpG_ID))
length(unique(annot_EPIC$Gene_Symbol))
length(unique(annot_EPIC$Entrez_ID))

# Store or load methylation array annotation
write.csv(annot_EPIC, file = paste0(main_dir, input_dir, "/annotation_EPIC.csv"), 
          quote = FALSE, row.names = FALSE)

annot_EPIC <- read.csv(file = paste0(main_dir, input_dir, "/annotation_EPIC.csv"))
annot_EPIC[] <- lapply(annot_EPIC, as.character)

# Annotate and store corrected DMG sets
DMG_Exp_vs_HC <- annot_EPIC[annot_EPIC$CpG_ID %in% DMC_Exp_vs_HC$Name, ]
length(na.omit(unique(DMG_Exp_vs_HC$Gene_Symbol)))
entrez_Exp_vs_HC <- na.omit(unique(annot_EPIC[annot_EPIC$CpG_ID %in% DMC_Exp_vs_HC$Name, "Entrez_ID"]))
length(entrez_Exp_vs_HC)

DMG_Pat_vs_Exp <- annot_EPIC[annot_EPIC$CpG_ID %in% DMC_Pat_vs_Exp$Name, ]
length(na.omit(unique(DMG_Pat_vs_Exp$Gene_Symbol)))
entrez_Pat_vs_Exp <- na.omit(unique(annot_EPIC[annot_EPIC$CpG_ID %in% DMC_Pat_vs_Exp$Name, "Entrez_ID"]))
length(entrez_Pat_vs_Exp)

DMG_Pat_vs_HC <- annot_EPIC[annot_EPIC$CpG_ID %in% DMC_Pat_vs_HC$Name, ]
length(na.omit(unique(DMG_Pat_vs_HC$Gene_Symbol)))
entrez_Pat_vs_HC <- na.omit(unique(annot_EPIC[annot_EPIC$CpG_ID %in% DMC_Pat_vs_HC$Name, "Entrez_ID"]))
length(entrez_Pat_vs_HC)

write.csv(DMG_Exp_vs_HC, file = paste0(main_dir, output_dir, "DMG_Exp_vs_HC.csv"), 
          row.names = FALSE)
write.csv(DMG_Pat_vs_Exp, file = paste0(main_dir, output_dir, "DMG_Pat_vs_Exp.csv"), 
          row.names = FALSE)
write.csv(DMG_Pat_vs_HC, file = paste0(main_dir, output_dir, "DMG_Pat_vs_HC.csv"), 
          row.names = FALSE)

#####################################################################
# 2. Generation and storage of MODifieR input                       #
#####################################################################

# Create and store MODifieR input
input_Exp_vs_HC <- create_custom_microarray_input_object(diff_genes = 
                                                           as.data.frame(matrix(data = c(entrez_Exp_vs_HC, rep(0, length(entrez_Exp_vs_HC))), ncol = 2)))
input_Exp_vs_HC$diff_genes$pvalue <- 0
input_Pat_vs_Exp <- create_custom_microarray_input_object(diff_genes = 
                                                            as.data.frame(matrix(data = c(entrez_Pat_vs_Exp, rep(0, length(entrez_Pat_vs_Exp))), ncol = 2)))
input_Pat_vs_Exp$diff_genes$pvalue <- 0
input_Pat_vs_HC <- create_custom_microarray_input_object(diff_genes = 
                                                           as.data.frame(matrix(data = c(entrez_Pat_vs_HC, rep(0, length(entrez_Pat_vs_HC))), ncol = 2)))
input_Pat_vs_HC$diff_genes$pvalue <- 0

saveRDS(object = input_Exp_vs_HC, file = paste0(main_dir, "Results/RDS_files/input_Exp_vs_HC.RDS"))
saveRDS(object = input_Pat_vs_Exp, file = paste0(main_dir, "Results/RDS_files/input_Pat_vs_Exp.RDS"))
saveRDS(object = input_Pat_vs_HC, file = paste0(main_dir, "Results/RDS_files/input_Pat_vs_HC.RDS"))

input_Exp_vs_HC <- readRDS(paste0(main_dir, "Results/RDS_files/input_Exp_vs_HC.RDS"))
input_Pat_vs_Exp <- readRDS(paste0(main_dir, "Results/RDS_files/input_Pat_vs_Exp.RDS"))
input_Pat_vs_HC <- readRDS(paste0(main_dir, "Results/RDS_files/input_Pat_vs_HC.RDS"))

#####################################################################
# 3. Module identification with MODifieR                            #
#####################################################################

## Import PPI network
network_score <- "700"
network_category <- "combined_score"
ppi_network <- read.csv(paste0(main_dir, input_dir, "PPI_network/", network_category, "/entrez_", 
                               network_category, "_", network_score, ".txt"), stringsAsFactors = FALSE)

# Module identification with MODifieR
module_list <- list()
mod_input <- input_Pat_vs_HC
mod_name <- "Pat_vs_HC"

## MCODE
mcode_out <- mcode(MODifieR_input = mod_input, 
                   ppi_network    = ppi_network, 
                   hierarchy      = 1, 
                   vwp            = 0.5,     
                   haircut        = FALSE, 
                   fluff          = FALSE, 
                   fdt            = 0.8, 
                   loops          = TRUE, 
                   deg_cutoff     = 0.05, 
                   module_cutoff  = 3.5, 
                   dataset_name   = mod_name)
if (length(mcode_out$module_genes) >= 1) {
  module_list[[length(module_list) + 1]] <- c(mcode_out$module_genes)
  names(module_list)[[length(module_list)]] <- paste0("mcode_size_", length(c(mcode_out$module_genes)))
} else { 
  module_list[[length(module_list) + 1]] <- NA
  names(module_list)[[length(module_list)]] <- paste0("mcode_size_", length(c(mcode_out$module_genes)))
}

## Clique SuM
# SQLite database construction
#build_clique_db(ppi_network = ppi_network, 
#db_folder   = paste0(main_dir, "Results/RDS_files/"), 
#db_name     = "SQL_db")

cliquesum_out <- clique_sum_permutation(MODifieR_input      = mod_input, 
                                        db                  = paste0(main_dir, "Results/RDS_files/SQL_db.sqlite"), 
                                        n_iterations        = 10000, 
                                        clique_significance = 0.01, 
                                        min_clique_size     = 5, 
                                        multiple_cores      = FALSE, 
                                        n_cores             = 1, 
                                        dataset_name        = NULL)
if (length(cliquesum_out$module_genes) >= 1 ) {
  module_list[[length(module_list) + 1]] <- c(cliquesum_out$module_genes)
  names(module_list)[[length(module_list)]] <- paste0("cliquesum_size_", length(cliquesum_out$module_genes))
} else { 
  module_list[[length(module_list) + 1]] <- "NA"
  names(module_list)[[length(module_list)]] <- paste0("cliquesum_size_", length(cliquesum_out$module_genes))
}

## DIAMOnD
diamond_out <- diamond(MODifieR_input = mod_input, 
                       ppi_network    = ppi_network, 
                       deg_cutoff     = 0.05, 
                       n_output_genes = 200, 
                       seed_weight    = 10, 
                       include_seed   = TRUE, 
                       dataset_name   = mod_name)
if (length(diamond_out$module_genes) >= 1 ) {
  module_list[[length(module_list) + 1]] <- c(diamond_out$module_genes)
  names(module_list)[[length(module_list)]] <- paste0("diamond_size_", length(diamond_out$module_genes))
} else { 
  module_list[[length(module_list) + 1]] <- "NA" 
  names(module_list)[[length(module_list)]] <- paste0("diamond_size_", length(diamond_out$module_genes))
}

# Merge the module genes to store
module_length <- lapply(module_list, length)
module_longest <- which.max(module_length)
genes_longest <- unlist(module_list[module_longest])
module_bind <- as.data.frame(matrix(data = NA, ncol = length(genes_longest), nrow = 0))

for (i in 1:length(module_list)) {
  mod_row <- c(NA, module_list[[i]], rep(NA, length(genes_longest) - length(module_list[[i]])))
  module_bind <- rbind(module_bind, mod_row, stringsAsFactors = FALSE)
  if (is.null(names(module_list)[i])) {
    rownames(module_bind)[i] <- paste0("NULL", i)
  } else if (is.na(names(module_list)[i])) {
    rownames(module_bind)[i] <- paste0("NA", i)
  } else { rownames(module_bind)[i] <- names(module_list)[i] }
}
module_out <- t(module_bind[, -1])

write.csv(module_out, file = paste0(main_dir, "Results/Modules/", mod_name, "_modules.csv"), 
          row.names = FALSE, quote = FALSE)

#####################################################################
# 4. Consensus of MODifieR modules (frequency)                      #
#####################################################################

## Across dataset (all methods)
dmg_list <- list(na.omit(DMG_Pat_vs_HC$Entrez_ID))
module_list <- module_list_Pat_vs_HC
mod_name <- "Pat_vs_HC"

names(dmg_list) <- paste0("DMGs_size_", length(dmg_list[[1]]))

# Generation of consensus modules by frequency
mod_genes_freq <- as.data.frame(table(unlist(module_list, recursive = TRUE, use.names = FALSE)))
mod_genes_freq$Var1 <- as.character(mod_genes_freq$Var1)
consensus_list <- list()

for (i in 1:length(module_list)) {
  con_mod <- mod_genes_freq[mod_genes_freq$Freq >= i, 1]
  consensus_list[[length(consensus_list) + 1]] <- con_mod
  names(consensus_list)[length(consensus_list)] <- paste0("consensus_", i, "_3_size_", length(con_mod))
}

# Merge the consensus module genes to store
combined_list <- c(dmg_list, module_list)
combined_list <- c(combined_list, consensus_list)
module_length <- lapply(combined_list, length)
module_longest <- which.max(module_length)
genes_longest <- unlist(combined_list[module_longest])
module_bind <- as.data.frame(matrix(data = NA, ncol = length(genes_longest), nrow = 0))

for (i in 1:length(combined_list)) {
  mod_row <- c(NA, combined_list[[i]], rep(NA, length(genes_longest) - length(combined_list[[i]])))
  module_bind <- rbind(module_bind, mod_row, stringsAsFactors = FALSE)
  if (is.null(names(combined_list)[i])) {
    rownames(module_bind)[i] <- paste0("NULL", i)
  } else if (is.na(names(combined_list)[i])) {
    rownames(module_bind)[i] <- paste0("NA", i)
  } else { rownames(module_bind)[i] <- names(combined_list)[i] }
}
module_out <- t(module_bind[, -1])

write.csv(module_out, file = paste0(main_dir, "Results/Modules/", mod_name, "_consensus.csv"),
          row.names = FALSE, quote = FALSE)

# Store module genes as Gene Symbols
module_out_symbol <- module_out
for (i in 1:ncol(module_out_symbol)) {
  entrez_sel <- na.omit(unlist(module_out_symbol[, i]))
  symbol_sel <- entrez.to.symbol(entrez_sel)
  module_out_symbol[1:length(symbol_sel), i] <- symbol_sel
}

write.csv(module_out_symbol, file = paste0(main_dir, "Results/Modules/", mod_name, "_consensus_GeneSymbol.csv"), 
          row.names = FALSE)

#####################################################################
# 5. Pathway enrichment of modules                                  #
#####################################################################

## Across dataset (all methods)
mod_name <- "Pat_vs_HC"
mod_out <- read.csv(file = paste0(main_dir, "Results/Modules/", mod_name, "_consensus.csv"))

top_cuts <- colnames(mod_out)

## Pathway analysis (clusterProfiler - KEGG)
enrich_list <- list()
for (top_sel in top_cuts) {
  sel_cpgs <- as.character(unlist(na.omit(mod_out[top_sel])))
  cpgs_kegg <- enrichKEGG(gene = sel_cpgs, organism = "hsa", keyType = "kegg", 
                          pvalueCutoff = 0.05, pAdjustMethod = "BH")
  sig_cpgs_kegg <- cpgs_kegg@result[cpgs_kegg@result$p.adjust < 0.05, ]
  if (nrow(sig_cpgs_kegg) == 0) {
    sig_cpgs_kegg <- cpgs_kegg@result[1, ]
    print(paste0("No enrichment for ", top_sel))
  }
  
  enrich_list[[length(enrich_list) + 1]] <- as.data.frame(cpgs_kegg@result)
  names(enrich_list)[[length(enrich_list)]] <- top_sel
  
  if (top_sel == top_cuts[1]) {
    write.xlsx(cpgs_kegg@result, file = paste0(main_dir, output_dir, "kegg_", mod_name, "_all.xlsx"), 
               sheetName = as.character(top_sel), row.names = FALSE, append = FALSE)
    
    write.xlsx(sig_cpgs_kegg, file = paste0(main_dir, output_dir, "kegg_", mod_name, "_signif.xlsx"), 
               sheetName = as.character(top_sel), row.names = FALSE, append = FALSE)
    
  } else {
    write.xlsx(cpgs_kegg@result, file = paste0(main_dir, output_dir, "kegg_", mod_name, "_all.xlsx"), 
               sheetName = as.character(top_sel), row.names = FALSE, append = TRUE)
    
    write.xlsx(sig_cpgs_kegg, file = paste0(main_dir, output_dir, "kegg_", mod_name, "_signif.xlsx"), 
               sheetName = as.character(top_sel), row.names = FALSE, append = TRUE)
  }
  assign("kegg_enrich", enrich_list)
}

# Apply styles to output file
gc()
wb <- loadWorkbook(file = paste0(main_dir, output_dir, "kegg_", mod_name, "_all.xlsx"))
sheets <- getSheets(wb)
for (sheet in sheets) {
  rows <- getRows(sheet)
  cells <- getCells(rows, colIndex = 9)
  values <- lapply(cells, getCellValue)
  
  # Format column width and column names
  cells <- getCells(rows, colIndex = 1:9)
  setColumnWidth(sheet, colIndex = 1, colWidth = 10)
  setColumnWidth(sheet, colIndex = c(2:9), colWidth = 8)
  setColumnWidth(sheet, colIndex = 2, colWidth = 50)
  setColumnWidth(sheet, colIndex = 8, colWidth = 20)
  
  # Save workbook sheet
  saveWorkbook(wb, file = paste0(main_dir, output_dir, "kegg_", mod_name, "_all.xlsx"))
}

# Apply styles to output file
gc()
wb <- loadWorkbook(file = paste0(main_dir, output_dir, "kegg_", mod_name, "_signif.xlsx"))
sheets <- getSheets(wb)
for (sheet in sheets) {
  rows <- getRows(sheet)
  cells <- getCells(rows, colIndex = 9)
  values <- lapply(cells, getCellValue)
  
  # Format column width and column names
  cells <- getCells(rows, colIndex = 1:9)
  setColumnWidth(sheet, colIndex = 1, colWidth = 10)
  setColumnWidth(sheet, colIndex = c(2:9), colWidth = 8)
  setColumnWidth(sheet, colIndex = 2, colWidth = 50)
  setColumnWidth(sheet, colIndex = 8, colWidth = 20)
  
  # Save workbook sheet
  saveWorkbook(wb, file = paste0(main_dir, output_dir, "kegg_", mod_name, "_signif.xlsx"))
}

# Visualize the top 20 significantly enriched pathways of consensus as a dot plot
pdf(file = paste0(main_dir, "Results/consensus_", mod_name, "_kegg_dotplot.pdf"), width = 10, height = 14)
dotplot(cpgs_kegg, showCategory = 20, title = paste0("Top enriched KEGG pathways (module consensus ", mod_name, ", p.adj < 0.05)"))
dev.off()

## Pathway analysis (ReactomePA - Reactome)
enrich_list <- list()
for (top_sel in top_cuts) {
  sel_cpgs <- as.character(unlist(na.omit(mod_out[top_sel])))
  cpgs_reactome <- enrichPathway(gene = sel_cpgs, organism = "human", 
                                 pvalueCutoff = 0.05, pAdjustMethod = "BH")
  sig_cpgs_reactome <- cpgs_reactome@result[cpgs_reactome@result$p.adjust < 0.05, ]
  if (nrow(sig_cpgs_reactome) == 0) {
    sig_cpgs_reactome <- cpgs_reactome@result[1, ]
    print(paste0("No enrichment for ", top_sel))
  }
  
  enrich_list[[length(enrich_list) + 1]] <- as.data.frame(cpgs_reactome@result)
  names(enrich_list)[[length(enrich_list)]] <- top_sel
  
  if (top_sel == top_cuts[1]) {
    write.xlsx(cpgs_reactome@result, file = paste0(main_dir, output_dir, "reactome_", mod_name, "_all.xlsx"), 
               sheetName = as.character(top_sel), row.names = FALSE, append = FALSE)
    
    write.xlsx(sig_cpgs_reactome, file = paste0(main_dir, output_dir, "reactome_", mod_name, "_signif.xlsx"),  
               sheetName = as.character(top_sel), row.names = FALSE, append = FALSE)
    
  } else {
    write.xlsx(cpgs_reactome@result, file =paste0(main_dir, output_dir, "reactome_", mod_name, "_all.xlsx"),  
               sheetName = as.character(top_sel), row.names = FALSE, append = TRUE)
    
    write.xlsx(sig_cpgs_reactome, file = paste0(main_dir, output_dir, "reactome_", mod_name, "_signif.xlsx"), 
               sheetName = as.character(top_sel), row.names = FALSE, append = TRUE)
  }
  assign("reactome_enrich", enrich_list)
}

# Apply styles to output file
gc()
wb <- loadWorkbook(file = paste0(main_dir, output_dir, "reactome_", mod_name, "_all.xlsx"))
sheets <- getSheets(wb)
for (sheet in sheets) {
  rows <- getRows(sheet)
  cells <- getCells(rows, colIndex = 9)
  values <- lapply(cells, getCellValue)
  
  # Format column width and column names
  cells <- getCells(rows, colIndex = 1:9)
  setColumnWidth(sheet, colIndex = 1, colWidth = 10)
  setColumnWidth(sheet, colIndex = c(2:9), colWidth = 8)
  setColumnWidth(sheet, colIndex = 2, colWidth = 50)
  setColumnWidth(sheet, colIndex = 8, colWidth = 20)
  
  # Save workbook sheet
  saveWorkbook(wb, file = paste0(main_dir, output_dir, "reactome_", mod_name, "_all.xlsx"))
}

# Apply styles to output file
gc()
wb <- loadWorkbook(file = paste0(main_dir, output_dir, "reactome_", mod_name, "_signif.xlsx"))
sheets <- getSheets(wb)
for (sheet in sheets) {
  rows <- getRows(sheet)
  cells <- getCells(rows, colIndex = 9)
  values <- lapply(cells, getCellValue)
  
  # Format column width and column names
  cells <- getCells(rows, colIndex = 1:9)
  setColumnWidth(sheet, colIndex = 1, colWidth = 10)
  setColumnWidth(sheet, colIndex = c(2:9), colWidth = 8)
  setColumnWidth(sheet, colIndex = 2, colWidth = 50)
  setColumnWidth(sheet, colIndex = 8, colWidth = 20)
  
  # Save workbook sheet
  saveWorkbook(wb, file = paste0(main_dir, output_dir, "reactome_", mod_name, "_signif.xlsx"))
}

# Visualize the top 20 significantly enriched pathways of consensus as a dot plot
pdf(file = paste0(main_dir, "Results/consensus_", mod_name, "_reactome_dotplot.pdf"), width = 10, height = 14)
dotplot(cpgs_reactome, showCategory = 20, title = paste0("Top enriched Reactome pathways (module consensus ", mod_name, ", p.adj < 0.05)"))
dev.off()

## DO enrichment (DOSE - DO categories)
enrich_list <- list()
for (top_sel in top_cuts) {
  sel_cpgs <- as.character(unlist(na.omit(mod_out[top_sel])))
  cpgs_do <- enrichDO(gene = sel_cpgs, ont = "DO",
                      pvalueCutoff = 0.05, pAdjustMethod = "BH")
  sig_cpgs_do <- cpgs_do@result[cpgs_do@result$p.adjust < 0.05, ]
  if (nrow(sig_cpgs_do) == 0) {
    sig_cpgs_do <- cpgs_do@result[1, ]
    print(paste0("No enrichment for ", top_sel))
  }
  
  enrich_list[[length(enrich_list) + 1]] <- as.data.frame(cpgs_do@result)
  names(enrich_list)[[length(enrich_list)]] <- top_sel
  
  if (top_sel == top_cuts[1]) {
    write.xlsx(cpgs_do@result, file = paste0(main_dir, output_dir, "do_", mod_name, "_all.xlsx"), 
               sheetName = as.character(top_sel), row.names = FALSE, append = FALSE)
    
    write.xlsx(sig_cpgs_do, file = paste0(main_dir, output_dir, "do_", mod_name, "_signif.xlsx"), 
               sheetName = as.character(top_sel), row.names = FALSE, append = FALSE)
    
  } else {
    write.xlsx(cpgs_do@result, file = paste0(main_dir, output_dir, "do_", mod_name, "_all.xlsx"), 
               sheetName = as.character(top_sel), row.names = FALSE, append = TRUE)
    
    write.xlsx(sig_cpgs_do, file = paste0(main_dir, output_dir, "do_", mod_name, "_signif.xlsx"), 
               sheetName = as.character(top_sel), row.names = FALSE, append = TRUE)
  }
  assign("do_enrich", enrich_list)
}

# Apply styles to output file
gc()
wb <- loadWorkbook(file = paste0(main_dir, output_dir, "do_", mod_name, "_all.xlsx"))
sheets <- getSheets(wb)
for (sheet in sheets) {
  rows <- getRows(sheet)
  cells <- getCells(rows, colIndex = 9)
  values <- lapply(cells, getCellValue)
  
  # Format column width and column names
  cells <- getCells(rows, colIndex = 1:9)
  setColumnWidth(sheet, colIndex = 1, colWidth = 10)
  setColumnWidth(sheet, colIndex = c(2:9), colWidth = 8)
  setColumnWidth(sheet, colIndex = 2, colWidth = 50)
  setColumnWidth(sheet, colIndex = 8, colWidth = 20)
  
  # Save workbook sheet
  saveWorkbook(wb, file = paste0(main_dir, output_dir, "do_", mod_name, "_all.xlsx"))
}

# Apply styles to output file
gc()
wb <- loadWorkbook(file = paste0(main_dir, output_dir, "do_", mod_name, "_signif.xlsx"))
sheets <- getSheets(wb)
for (sheet in sheets) {
  rows <- getRows(sheet)
  cells <- getCells(rows, colIndex = 9)
  values <- lapply(cells, getCellValue)
  
  # Format column width and column names
  cells <- getCells(rows, colIndex = 1:9)
  setColumnWidth(sheet, colIndex = 1, colWidth = 10)
  setColumnWidth(sheet, colIndex = c(2:9), colWidth = 8)
  setColumnWidth(sheet, colIndex = 2, colWidth = 50)
  setColumnWidth(sheet, colIndex = 8, colWidth = 20)
  
  # Save workbook sheet
  saveWorkbook(wb, file = paste0(main_dir, output_dir, "do_", mod_name, "_signif.xlsx"))
}

# Visualize the top 20 significantly enriched terms of consensus as a dot plot
pdf(file = paste0(main_dir, "Results/consensus_", mod_name, "_DO_dotplot.pdf"), width = 10, height = 14)
dotplot(cpgs_do, showCategory = 20, title = paste0("Top enriched DO terms (module consensus ", mod_name, ", p.adj < 0.05)"))
dev.off()

## GO enrichment (clusterProfiler, AnnotationHub - GO terms)
enrich_list <- list()
for (top_sel in top_cuts) {
  gc()
  sel_cpgs <- as.character(unlist(na.omit(mod_out[top_sel])))
  cpgs_go <- enrichGO(gene = sel_cpgs, OrgDb = org.Hs.eg.db, ont = "BP",
                      pvalueCutoff = 0.05, pAdjustMethod = "BH")
  sig_cpgs_go <- cpgs_go@result[cpgs_go@result$p.adjust < 0.05, ]
  gc()
  if (nrow(sig_cpgs_go) == 0) {
    sig_cpgs_go <- cpgs_go@result[1, ]
    print(paste0("No enrichment for ", top_sel))
  }
  
  enrich_list[[length(enrich_list) + 1]] <- as.data.frame(cpgs_go@result)
  names(enrich_list)[[length(enrich_list)]] <- top_sel
  
  if (top_sel == top_cuts[1]) {
    write.xlsx(cpgs_go@result, file = paste0(main_dir, output_dir, "go_", mod_name, "_all.xlsx"),
               sheetName = as.character(top_sel), row.names = FALSE, append = FALSE)
    
    write.xlsx(sig_cpgs_go, file = paste0(main_dir, output_dir, "go_", mod_name, "_signif.xlsx"),
               sheetName = as.character(top_sel), row.names = FALSE, append = FALSE)
    
  } else {
    write.xlsx(cpgs_go@result, file = paste0(main_dir, output_dir, "go_", mod_name, "_all.xlsx"),
               sheetName = as.character(top_sel), row.names = FALSE, append = TRUE)
    
    write.xlsx(sig_cpgs_go, file = paste0(main_dir, output_dir, "go_", mod_name, "_signif.xlsx"),
               sheetName = as.character(top_sel), row.names = FALSE, append = TRUE)
  }
  assign("go_enrich", enrich_list)
}

gc()
wb <- loadWorkbook(file = paste0(main_dir, output_dir, "go_", mod_name, "_all.xlsx"))
sheets <- getSheets(wb)
for (sheet in sheets) {
  rows <- getRows(sheet)
  cells <- getCells(rows, colIndex = 9)
  values <- lapply(cells, getCellValue)
  
  # Format column width and column names
  cells <- getCells(rows, colIndex = 1:9)
  setColumnWidth(sheet, colIndex = 1, colWidth = 10)
  setColumnWidth(sheet, colIndex = c(2:9), colWidth = 8)
  setColumnWidth(sheet, colIndex = 2, colWidth = 50)
  setColumnWidth(sheet, colIndex = 8, colWidth = 20)
  
  # Save workbook sheet
  saveWorkbook(wb, file = paste0(main_dir, output_dir, "go_", mod_name, "_all.xlsx"))
}

# # Apply styles to output file
gc()
wb <- loadWorkbook(file = paste0(main_dir, output_dir, "go_", mod_name, "_signif.xlsx"))
sheets <- getSheets(wb)
for (sheet in sheets) {
  rows <- getRows(sheet)
  cells <- getCells(rows, colIndex = 9)
  values <- lapply(cells, getCellValue)
  
  # Format column width and column names
  cells <- getCells(rows, colIndex = 1:9)
  setColumnWidth(sheet, colIndex = 1, colWidth = 10)
  setColumnWidth(sheet, colIndex = c(2:9), colWidth = 8)
  setColumnWidth(sheet, colIndex = 2, colWidth = 50)
  setColumnWidth(sheet, colIndex = 8, colWidth = 20)
  
  # Save workbook sheet
  saveWorkbook(wb, file = paste0(main_dir, output_dir, "go_", mod_name, "_signif.xlsx"))
}

# Visualize the top 20 significantly enriched terms of consensus as a dot plot
pdf(file = paste0(main_dir, "Results/consensus_", mod_name, "_GO_dotplot.pdf"), width = 10, height = 14)
dotplot(cpgs_go, showCategory = 20, title = paste0("Top enriched GO terms (module consensus ", mod_name, ", p.adj < 0.05)"))
dev.off()

#####################################################################
# 6. Graph visualization of consensus modules                       #
#####################################################################

## Across dataset (all methods)
mod_name <- "Pat_vs_Exp"
mod_genes <- as.character(na.omit(read.csv(file = paste0(main_dir, "Results/Modules/", mod_name, "_consensus_GeneSymbol.csv"))[, 7]))

## Import Gene Symbol PPI network
# Score: 0, 500, 600, 700, 800, 900
# Category: coexpression, combined_score, cooccurence, database, experimental, 
# fusion, neighborhood, textmining
network_score <- "700"
network_category <- "combined_score"
ppi_network <- read.csv(paste0(main_dir, input_dir, "PPI_network/", network_category, "/hgnc_", 
                               network_category, "_", network_score, ".txt"), stringsAsFactors = FALSE)

# Map module genes to PPI network interactions
mod_net <- ppi_network[ppi_network[, 1] %in% mod_genes, ]
mod_net <- mod_net[mod_net[, 2] %in% mod_genes, ]
mod_graph <- graph_from_data_frame(mod_net, directed = FALSE)

# Remove loops and multiple edges from the graph
mod_graph <- igraph::simplify(mod_graph, remove.multiple = TRUE,
                              remove.loops = TRUE)

# Plot igraph representation
pdf(file = paste0(main_dir, "Results/consensus_", mod_name, "_graph.pdf"), width = 14, height = 14)
plot(mod_graph, vertex.size = 4,
     vertex.color = "steelblue2",
     vertex.frame.color = "#555555",
     vertex.label = V(mod_graph)$type,
     vertex.label.color = "black",
     vertex.label.cex = 0.4,
     edge.arrow.size = 0.4, 
     layout = layout_with_fr(mod_graph),
     main = paste0(mod_name))
dev.off()

#####################################################################
# 7. Venn diagram of consensus module genes                         #
#####################################################################

# Build the set of gene lists for the Venn diagram
set1 <- as.character(na.omit(read.csv(file = paste0(main_dir, "Results/Modules/Exp_vs_HC_consensus.csv"))[, 7]))
set2 <- as.character(na.omit(read.csv(file = paste0(main_dir, "Results/Modules/Pat_vs_Exp_consensus.csv"))[, 7]))
set3 <- as.character(na.omit(read.csv(file = paste0(main_dir, "Results/Modules/Pat_vs_HC_consensus.csv"))[, 7]))

pdf(file = paste0(main_dir, "Results/consensus_venndiagram.pdf"), width = 10, height = 10)
draw.triple.venn(area1 = length(set1),
                 area2 = length(set2),
                 area3 = length(set3),
                 n12 = length(intersect(set1, set2)),
                 n23 = length(intersect(set2, set3)),
                 n13 = length(intersect(set1, set3)),
                 n123 = length(intersect(intersect(set1, set2), set3)),
                 category = c("Exp_vs_HC", "Pat_vs_Exp", "Pat_vs_HC"),
                 lty = "blank",
                 fill = c("steelblue1", "yellowgreen", "indianred1") ,
                 cex = 2, cat.cex = 2, cat.fontfamily = rep("serif", 3))
dev.off()

# Extract gene lists from all overlaps between consensus modules
set1 <- as.character(na.omit(read.csv(file = paste0(main_dir, "Results/Modules/Exp_vs_HC_consensus.csv"))[, 7]))
set2 <- as.character(na.omit(read.csv(file = paste0(main_dir, "Results/Modules/Pat_vs_Exp_consensus.csv"))[, 7]))
set3 <- as.character(na.omit(read.csv(file = paste0(main_dir, "Results/Modules/Pat_vs_HC_consensus.csv"))[, 7]))
overlap_list <- list(setdiff(setdiff(set1, set2), set3), setdiff(setdiff(set2, set1), set3), 
                     setdiff(setdiff(set3, set2), set1), setdiff(intersect(set1, set2), set3), 
                     setdiff(intersect(set1, set3), set2), setdiff(intersect(set2, set3), set1), 
                     intersect(intersect(set1, set2), set3))
names(overlap_list) <- c("Exp_vs_HC_unique", "Pat_vs_Exp_unique", 
                         "Pat_vs_HC_unique", "Exp_vs_HC_+_Pat_vs_Exp", 
                         "Exp_vs_HC_+_Pat_vs_HC", "Pat_vs_Exp_+_Pat_vs_HC", 
                         "Exp_vs_HC_+_Pat_vs_Exp_+_Pat_vs_HC")
module_length <- lapply(overlap_list, length)
module_longest <- which.max(module_length)
genes_longest <- unlist(overlap_list[module_longest])
module_bind <- as.data.frame(matrix(data = NA, ncol = length(genes_longest), nrow = 0))

for (i in 1:length(overlap_list)) {
  mod_row <- c(overlap_list[[i]], rep(NA, length(genes_longest) - length(overlap_list[[i]])))
  module_bind <- rbind(module_bind, mod_row, stringsAsFactors = FALSE)
  if (is.null(names(overlap_list)[i])) {
    rownames(module_bind)[i] <- paste0("NULL", i)
  } else if (is.na(names(overlap_list)[i])) {
    rownames(module_bind)[i] <- paste0("NA", i)
  } else { rownames(module_bind)[i] <- names(overlap_list)[i] }
}
module_out <- t(module_bind[, -1])

write.table(module_out, file = paste0(main_dir, "Results/consensus_all_overlaps.csv"), 
            row.names = FALSE, quote = FALSE, sep = ",", dec = ".")

# Store module genes as Gene Symbols
module_out_symbol <- module_out
for (i in 1:ncol(module_out_symbol)) {
  entrez_sel <- na.omit(unlist(module_out_symbol[, i]))
  symbol_sel <- entrez.to.symbol(entrez_sel)
  module_out_symbol[1:length(symbol_sel), i] <- symbol_sel
}

write.csv(module_out_symbol, file = paste0(main_dir, "Results/consensus_all_overlaps_GeneSymbol.csv"), 
          row.names = FALSE)

#####################################################################
# 8. Mapping of hyper/hypo methylated probes in consensus modules   #
#####################################################################

# Annotation auxiliary functions
symbol.to.entrez <- function(x) {
  x <- mget(x, org.Hs.egSYMBOL2EG, ifnotfound = NA)
  if (length(x) == 0) { return("NA") 
  } else {unlist(lapply(x, function(i) { i[1]}))}
}

entrez.to.symbol <- function(x) {
  x <- mget(x, org.Hs.egSYMBOL, ifnotfound = NA)
  if (length(x) == 0) { return("NA") 
  } else {unlist(lapply(x, function(i) { i[1]}))}
}

# Import consensus module genes
mod_name <- "Pat_vs_Exp"

# Import DMC results
DMC_Exp_vs_HC <- read.csv(paste0(main_dir, input_dir, "DMC_Exp_vs_HC_padj0.05_0.2.csv"), sep = ";")
DMC_Pat_vs_Exp <- read.csv(paste0(main_dir, input_dir, "DMC_Pat_vs_Exp_padj0.05_0.2.csv"), sep = ";")
DMC_Pat_vs_HC <- read.csv(paste0(main_dir, input_dir, "DMC_Pat_vs_HC_padj0.05_0.2.csv"))

# Annotate and store corrected DMG sets
DMG_Exp_vs_HC <- read.csv(paste0(main_dir, output_dir, "DMG_Exp_vs_HC.csv"))
DMG_Pat_vs_Exp <- read.csv(paste0(main_dir, output_dir, "DMG_Pat_vs_Exp.csv"))
DMG_Pat_vs_HC <- read.csv(paste0(main_dir, output_dir, "DMG_Pat_vs_HC.csv"))

DMC_table <- DMC_Pat_vs_Exp
DMG_table <- DMG_Pat_vs_Exp
consensus_table <- Pat_vs_Exp_consensus

## Identification of hyper- and hypomethylated genes in consensus module
# All probes for a gene should indicate the same direction of methylation to be considered hyper/hypo
merged_table <- merge(DMG_table, DMC_table[, c("logFC", "Name")], by.x = "CpG_ID", by.y = "Name")

meth_level_consensus <- merge(merged_table, as.data.frame(na.omit(consensus_table[, 7])), by.x = "Entrez_ID", by.y = "na.omit(consensus_table[, 7])")
meth_level_consensus[is.na(meth_level_consensus$logFC), "logFC"] <- 0

meth_level_consensus[, "Direction"] <- "Mixed"
for (i in 1:nrow(meth_level_consensus)) {
  if (meth_level_consensus[i, "logFC"] > 0) {
    meth_level_consensus[i, "Direction"] <- "Hyper"
  } else if (meth_level_consensus[i, "logFC"] < 0) {
    meth_level_consensus[i, "Direction"] <- "Hypo"
  }
}

assign(paste0("meth_level_", mod_name), meth_level_consensus)

write.csv(meth_level_consensus, file = paste0(main_dir, "Results/Modules/consensus_3_3_", mod_name, "_direction.csv"), 
          row.names = FALSE)
          
