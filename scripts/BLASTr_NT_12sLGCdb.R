---
title: "BLASTr  NT + 12sLGCdb"
author: "Hilario, HO; Mendes, GA"
date: "02/2024"
---
  
  ## Carregando bibliotecas ----
{
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggbreak)
  library(phyloseq)
  library(Biostrings)
  library(Matrix)
  library(ShortRead)
  library(dada2)
  library(DECIPHER)
  library(future)
  library(ggh4x)
  library(vegan)
  library(plotly)
  library(ggtext)
  library(BLASTr)
  library(writexl)
  library(readxl)
}
  
## Caminhos ----
{
  prjct_path <- "~/projetos/eDNA_Doce/"
  results_path <- paste0(prjct_path,"results/")
  figs_path <- paste0(results_path,"figures/")
  tbl_path <- paste0(prjct_path,"tables/")
  env_path <- paste0(prjct_path,"environments/")
  blast_path <- paste0(results_path,"blastr/")
}

## Obtencao dos dados ----
{
  # final_results <- read_excel(paste0(tbl_path,"/", "eDNA_Doce-Complete_analysis_results-2024-02-19.xlsx")) %>% tibble()
  final_results <- read.csv(paste0(tbl_path,"/","eDNA_Doce-Complete_analysis_results-2024-02-19_Resultados completos.csv"), sep = ",", check.names = FALSE) 
}

## Alteracoes pos-pipe ----

less_final_results <- final_results %>% 
  select(Unique_File_name,
         "ASV header",
         "ASV (Sequence)",
         "Primer expected length",
         "ASV absolute abundance",
         "Relative abundance on sample",
         "Relative abundance to all samples",
         "Sample total abundance") %>% 
  tibble()

View(less_final_results)

## BLASTn ----

# select ASVs for BLASTn search 
asvs_blast_all <- less_final_results %>%
  unique() %>% 
  pull("ASV (Sequence)") %>% 
  as.character()

asvs_blast_all
length(asvs_blast_all)

# Visualizando os tamanhos das ASVs
asvs_blast_all %>% nchar() %>% table() %>% plot()

# fazendo para apenas 1 ASV

asvs_blast_4434 <- less_final_results %>%
  filter("ASV header" %in% ">ASV_44334_196bp")

asvs_blast_less <- asvs_blast_4434 %>%
  unique() %>% 
  pull("ASV (Sequence)") %>% 
  as.character()

## Usando o BLASTn no R em paralelo com o BLASTr

# BLASTr NT + 12slGCdb
{
tictoc::tic("Parallel - Furrr 2 threads")
Sys.time()

blast_res_1 <- BLASTr::parallel_blast(
  # db_path = '"/data/databases/nt_jun2023/nt /home/gabriel/projetos/peixes-eDNA/databases/12sLGCdb/LGC12Sdb_complete_noGaps_2024-0403.fasta"',
  db_path = '"/data/databases/nt_jun2023/nt /home/gabriel/projetos/peixes-eDNA/databases/12sLGCdb/old_db/LGC12Sdb_complete_noGaps.fasta"',
  # asvs = asvs_blast_all[1:1000],
  asvs = asvs_blast_less,
  out_file = paste0(blast_path, "blast_out_res_1.csv"),
  out_RDS = paste0(blast_path, "blast_out_res_1.RDS"),
  total_cores = 80,
  perc_id = 80,
  num_threads = 2,
  perc_qcov_hsp = 80,
  num_alignments = 3,
  blast_type = "blastn"
)

blast_res_backup <- blast_res_1

tictoc::toc()
Sys.time()
}

{
tictoc::tic("Parallel - Furrr 2 threads")
Sys.time()

blast_res_2 <- BLASTr::parallel_blast(
  db_path = '"/data/databases/nt_jun2023/nt /home/gabriel/projetos/peixes-eDNA/databases/12sLGCdb/LGC12Sdb_complete_noGaps_2024-0403.fasta"' ,
  asvs = asvs_blast_all[1001:1840],
  out_file = paste0(blast_path, "blast_out_res_2.csv"),
  out_RDS = paste0(blast_path, "blast_out_res_2.RDS"),
  total_cores = 80,
  perc_id = 80,
  num_threads = 2,
  perc_qcov_hsp = 80,
  num_alignments = 3,
  blast_type = "blastn"
)

blast_res_backup2 <- blast_res_2



tictoc::toc()
Sys.time()
}

blast_res <- bind_rows(
  blast_res_1
  ,
  blast_res_2
  # blast_res_3,
  # blast_res_4,
  # blast_res_5
  )

View(blast_res)

# Saving environment 
base::save.image(paste0(env_path, "env_", Sys.Date(), "_postBLAST.RData"))  

## Editando o output do BLASTr 

# Editando os hits obtidos com o 12sLGCdb para ficarem semelhantes aos hits obtidos com o NT

blast_res_f <- blast_res %>%
  mutate("1_DB" = ifelse("1_subject header" == "\n", "12sLGCdb", "NT"),
         "1_subject header" = ifelse("1_DB" == "12sLGCdb", str_extract("1_subject", "(?<=_).+"), "1_subject header"),
         "1_subject header" = gsub("_", " ", "1_subject header"),
         "1_subject" = ifelse("1_DB" == "12sLGCdb", str_extract("1_subject", "^[^_]+"), "1_subject"),
         "2_DB" = ifelse("2_subject header" == "\n", "12sLGCdb", "NT"),
         "2_subject header" = ifelse("2_DB" == "12sLGCdb", str_extract("2_subject", "(?<=_).+"), "2_subject header"),
         "2_subject header" = gsub("_", " ", "2_subject header"),
         "2_subject" = ifelse("2_DB" == "12sLGCdb", str_extract("2_subject", "^[^_]+"), "2_subject"),
         "3_DB" = ifelse("3_subject header" == "\n", "12sLGCdb", "NT"),
         "3_subject header" = ifelse("3_DB" == "12sLGCdb", str_extract("3_subject", "(?<=_).+"), "3_subject header"),
         "3_subject header" = gsub("_", " ", "3_subject header"),
         "3_subject" = ifelse("3_DB" == "12sLGCdb", str_extract("3_subject", "^[^_]+"), "3_subject")
  )

View(blast_res_f)

## Adicionando o TaxID aos hits obtidos com o 12sLGCdb

# Função para recuperar o TaxID
retrieve_taxid <- function(organism_name) {
  Sys.sleep(0.1)
  result <- tryCatch({
    myTAI::taxonomy(
      organism = organism_name, 
      db = "ncbi", 
      output = "taxid"
    )
  }, error = function(e) {
    NULL
  })
  
  if (!is.null(result) && !inherits(result, "try-error")) {
    return(result$id[1])
  } else {
    return(NA)
  }
}

# Loop para recuperar o TaxID dos hits 1,2 e 3
for (hit_number in 1:3) {
  db_column <- paste0(hit_number, "_DB")
  subject_header_column <- paste0(hit_number, "_subject header")
  staxid_column <- paste0(hit_number, "_staxid")
  
  spp_db <- unique(blast_res_f[[db_column]][blast_res_f[[db_column]] == "12sLGCdb"])
  
  for (db_value in spp_db) {
    filtered_data <- blast_res_f %>%
      filter(.data[[db_column]] == db_value)
    
    for (subject_header in filtered_data[[subject_header_column]]) {
      organism_name <- word(subject_header, 1)
      
      taxid <- retrieve_taxid(organism_name)
      
      if (!is.na(taxid)) {
        blast_res_f[[staxid_column]][blast_res_f[[db_column]] == db_value & blast_res_f[[subject_header_column]] == subject_header] <- taxid
      } else {
        cat("Failed to retrieve TaxID for:", organism_name, "\n")
      }
    }
  }
  
  # ASVs sem TaxID
  no_taxid <- blast_res_f %>%
    filter(.data[[staxid_column]] == "N/A")
  
  for (subject_header in no_taxid[[subject_header_column]]) {
    organism_name <- word(subject_header, 1)
    
    taxid <- retrieve_taxid(organism_name)
    
    if (!is.na(taxid)) {
      blast_res_f[[staxid_column]][blast_res_f[[subject_header_column]] == subject_header] <- taxid
    } else {
      cat("Failed to retrieve TaxID for:", organism_name, "\n")
    }
  }
}

# Verifique o resultado final para saber se existem ASVs sem TaxID

sapply(blast_res_f[, c("1_staxid", "2_staxid", "3_staxid")], function(x) {
  any(str_detect(x, "\\bN/A\\b"))}) # Se tudo der NA, podemos prosseguir

blast_res_full <- blast_res_f

colnames(blast_res_full)

blast_res_full <- blast_res_full %>% 
  select(-c("OTU")) %>%
  filter(!is.na("1_subject header"))

## Retrieving complete taxonomies for BLAST results ----

#overview the identifications ----
blast_res_full$"1_subject header" %>% unique() %>% sort()

#set hits with poor names to remove from results
bad_1res_IDs <- c(
  "Uncultured organism clone",
  "Uncultured prokaryote",
  "Eukaryotic synthetic construct",
  "16S rRNA amplicon fragment",
  "Uncultured Candidatus",
  "Uncultured bacterium",
  "Uncultured archaeon clone",
  "Complete Metagenome-Assembled"
) %>% 
  paste0(collapse = "|")


blast_res_full <- blast_res_full %>%
  mutate(`blast ID` = `blast ID`,
         "blast ID Origin" = "blast ID Origin",
         "query_taxID" = "query_taxID")


# pick BLASTn res IDs and mark result origin ----

for (asv in 1:nrow(blast_res_full)) {
  if (!is.na(blast_res_full$`1_subject header`[asv]) &&
      stringr::str_detect(string = blast_res_full$`1_subject header`[asv], pattern = bad_1res_IDs)) {
    
    if (!is.na(blast_res_full$`2_subject header`[asv])) {
      blast_res_full$`blast ID`[asv] <- substr(as.character(blast_res_full$`2_subject header`[asv]), 1, 40)
      blast_res_full$`blast ID Origin`[asv] <- "2_"
      blast_res_full$query_taxID[asv] <- blast_res_full$`2_staxid`[asv]
      
      if (!is.na(blast_res_full$`3_subject header`[asv]) &&
          stringr::str_detect(string = blast_res_full$`2_subject header`[asv], pattern = bad_1res_IDs)) {
        
        blast_res_full$`blast ID`[asv] <- substr(as.character(blast_res_full$`3_subject header`[asv]), 1, 40)
        blast_res_full$`blast ID Origin`[asv] <- "3_"
        blast_res_full$query_taxID[asv] <- blast_res_full$`3_staxid`[asv]
        
        if (stringr::str_detect(string = blast_res_full$`3_subject header`[asv], pattern = bad_1res_IDs)) {
          blast_res_full$`blast ID`[asv] <- "Match_not_reliable"
          blast_res_full$`blast ID Origin`[asv] <- NA
          blast_res_full$query_taxID[asv] <- NA
        }
      }
    } else {
      blast_res_full$`blast ID`[asv] <- "Match_not_reliable"
      blast_res_full$`blast ID Origin`[asv] <- NA
      blast_res_full$query_taxID[asv] <- NA
    }
  } else {
    blast_res_full$`blast ID`[asv] <- substr(as.character(blast_res_full$`1_subject header`[asv]), 1, 40)
    blast_res_full$`blast ID Origin`[asv] <- "1_"
    blast_res_full$query_taxID[asv] <- blast_res_full$`1_staxid`[asv]
  }
}


blast_res_full$`blast ID` %>% unique() %>% sort()
blast_res_full$`blast ID Origin` %>% table() 

blast_res_full$`blast ID`<-  blast_res_full$`blast ID` %>% 
  stringr::str_remove(pattern = "^  |^ |Uncultured |uncultured |Candidatus |MAG:|MAG: |MAG TPA_asm: |TPA_asm: |^Cf. |\n|candidate division ") %>% 
  stringr::str_remove(pattern = "^  |^ |Uncultured |uncultured |Candidatus |MAG:|MAG: |MAG TPA_asm: |TPA_asm: |^Cf. |\n|candidate division ") %>% 
  stringr::str_remove_all(pattern = "\\[|\\]") %>% 
  # stringr::str_remove(pattern = "\\:.*.$") %>% 
  stringr::str_replace(pattern = "cf\\. ",replacement = "") %>% 
  stringr::str_replace(pattern = "nr\\. ",replacement = "") %>% 
  stringr::str_replace(pattern = "sp\\. ",replacement = "sp\\.") %>% 
  stringr::str_replace(pattern = "\\,",replacement = "") %>% 
  stringr::str_replace(pattern = "sp\\.",replacement = "sp\\. ")

# blast_res_full_bckp2 <- blast_res_full
# blast_res_full <- blast_res_full_bckp2

blast_res_full$`blast ID` %>% unique() %>% sort()
blast_res_full$`blast ID` %>% unique() %>% sort(decreasing = T)

# selecting just the first 2 names of BLAST result
for (row in 1:nrow(blast_res_full)) {
  
  blast_res_full$`blast ID`[row] <- stringr::str_split_fixed(string = blast_res_full$`blast ID`[row], pattern = " ",n = 3)[1:2] %>% 
    paste0(collapse = " ")
  
}

blast_res_full$`blast ID` %>% unique() %>% sort()


# correct confusing labels to unify identities ----

blast_res_full$`blast ID`[blast_res_full$`blast ID` %in% c("Human DNA",
                                                           "Eukaryotic synthetic",
                                                           "Human chromosome")] <- "Homo sapiens"

blast_res_full$`blast ID` %>% unique() %>% sort()


blast_res_full %>% filter(`blast ID` %in% c(" ")) %>% View()

## Retrieve complete taxonomy ----

# greate a genus colum to be able to join tax results


# blast_res_full_bckp3 <- blast_res_full


blast_res_full$`blast ID` %>% unique() %>% sort()


blast_res_full <- blast_res_full %>%
  mutate(max_tax = case_when(str_detect(`blast ID`,pattern = "PREDICTED:") ~ paste0(str_remove(`blast ID`,pattern = "PREDICTED: ")," (PREDICTED)"),
                             TRUE ~  str_remove(`blast ID`,pattern = " .*$")))

blast_res_full <- blast_res_full %>%
  relocate(`blast ID`,`blast ID Origin`,`query_taxID`,`max_tax`)


blast_res_full$max_tax %>% unique() %>% sort()

blast_res_full %>% 
write_xlsx(path = paste0(results_path, Sys.Date(),"_blast_res_full.xlsx"), 
           col_names = TRUE,
           format_headers = TRUE)


# FUNCTION to retrieve tax ranks using organism taxID ----
#test ----
source("/home/heron/prjcts/ecomol/R/extract_taxonomy_taxID.R")

#caraaaaalhoooooooo ta comendo o primeiro caracterpq?????
extract_taxonomy_taxID("42548")

#buscando as classificações

future::plan(future::multisession(), workers = 78)

taxIDs2search <- blast_res_full$query_taxID %>% unique() %>% na.omit() %>% as.character()

taxIDs2search %>% class()

taxonomy_df <- furrr::future_map_dfr(.x = taxIDs2search,
                                     .f = extract_taxonomy_taxID,
                                     .options = furrr::furrr_options(seed = TRUE))


# (A) chech the ones that have not been retrieved
taxIDs2search[!(taxIDs2search %in% c(unique(taxonomy_df$query_taxID)))]

# (B) search the missing ones

taxonomy_df1 <- furrr::future_map_dfr(.x = taxIDs2search[!(taxIDs2search %in% c(unique(taxonomy_df$query_taxID)))],
                                      .f = extract_taxonomy_taxID,
                                      .options = furrr::furrr_options(seed = TRUE))

# (C) combine the results
taxonomy_df <-  bind_rows(taxonomy_df,taxonomy_df1) %>% unique()

# (D) check if is there any other still missing
taxIDs2search[!(taxIDs2search %in% c(unique(taxonomy_df$query_taxID)))]

# (E) repeat A, B, C and D until no one is missing

taxonomy_df <- taxonomy_df %>% 
  filter(!Rank %in% c("no rank","clade"))


taxonomy_tbl <- taxonomy_df %>% 
  # select(-c("TaxId","Sci_name")) %>%
  dplyr::select(-c("TaxId")) %>%
  unique() %>%
  # dplyr::filter(Rank %in% c("kingdom","phylum","class","order","family")) %>%
  filter(Rank %in% c("superkingdom","kingdom","phylum","subphylum","class","subclass","order","suborder","family","subfamily","genus")) %>% 
  tidyr::pivot_wider(
    id_cols = c(query_taxID,Sci_name),
    names_from = Rank,
    values_from = c(ScientificName)) %>%
  # tidyr::pivot_wider(names_from = Rank,values_from = c(ScientificName,TaxId)) %>%
  # dplyr::select(max_tax,dplyr::starts_with("Scie")) %>% 
  relocate("Sci_name","query_taxID","superkingdom","kingdom","phylum","subphylum","class","subclass","order","suborder","family","subfamily","genus")



saveRDS(object = taxonomy_tbl, file = paste0(results_path,"taxonomy_df_from_taxIDs.rds"))


taxonomy_tbl <- readRDS(file = paste0(results_path,"/taxonomy_df_from_taxIDs.rds"))

# complete taxonomy tbl missing ranks

# taxonomy_tbl_bckp <- taxonomy_tbl
# orgs_tbl_bckp <- orgs_tbl

taxonomy_tbl %>% colnames() %>% paste0(collapse = "\n") %>% cat()
#fill NA tax with combination of max_tax and rank
for (line in 1:nrow(taxonomy_tbl)) {
  # if (taxonomy_tbl$genus[line] %in% c("NA",NA,"")) {
  if (taxonomy_tbl$superkingdom[line] %in% c("NA",NA,"")) {
    taxonomy_tbl$superkingdom[line] <- paste0("superkingdom of ", taxonomy_tbl$kingdom[line]) }
  
  if (taxonomy_tbl$kingdom[line] %in% c("NA",NA,"")) {
    taxonomy_tbl$kingdom[line] <- paste0("kingdom of ", taxonomy_tbl$superkingdom[line]) }
  
  if (taxonomy_tbl$phylum[line] %in% c("NA",NA,"")) {
    taxonomy_tbl$phylum[line] <- paste0("phylum of ", taxonomy_tbl$kingdom[line]) }
  
  if (taxonomy_tbl$subphylum[line] %in% c("NA",NA,"")) {
    taxonomy_tbl$subphylum[line] <- paste0("subphylum of ", taxonomy_tbl$phylum[line]) }
  
  if (taxonomy_tbl$class[line] %in% c("NA",NA,"")) {
    taxonomy_tbl$class[line] <- paste0("class of ", taxonomy_tbl$subphylum[line]) }
  
  if (taxonomy_tbl$subclass[line] %in% c("NA",NA,"")) {
    taxonomy_tbl$subclass[line] <- paste0("subclass of ", taxonomy_tbl$class[line]) }
  
  if (taxonomy_tbl$order[line] %in% c("NA",NA,"")) {
    taxonomy_tbl$order[line] <- paste0("order of ", taxonomy_tbl$subclass[line]) }
  
  if (taxonomy_tbl$suborder[line] %in% c("NA",NA,"")) {
    taxonomy_tbl$suborder[line] <- paste0("suborder of ", taxonomy_tbl$order[line]) }
  
  if (taxonomy_tbl$family[line] %in% c("NA",NA,"")) {
    taxonomy_tbl$family[line] <- paste0("family of ", taxonomy_tbl$suborder[line]) }
  
  if (taxonomy_tbl$subfamily[line] %in% c("NA",NA,"")) {
    taxonomy_tbl$subfamily[line] <- paste0("subfamily of ", taxonomy_tbl$family[line]) }
  
  if (is.na(taxonomy_tbl$genus[line])) {
    taxonomy_tbl$genus[line] <- paste0("genus of ", taxonomy_tbl$subfamily[line]) }
  # 
}

taxonomy_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()

taxonomy_tbl <- taxonomy_tbl %>% 
  dplyr::rename(
    "Superkingdom (BLASTn)" = "superkingdom",
    "Kingdom (BLASTn)" = "kingdom",
    "Phylum (BLASTn)" = "phylum",
    "Subphylum (BLASTn)" = "subphylum",
    "Class (BLASTn)" = "class",
    "Subclass (BLASTn)" = "subclass",
    "Order (BLASTn)" = "order",
    "Suborder (BLASTn)" = "suborder",
    "Family (BLASTn)" = "family",
    "Subfamily (BLASTn)" = "subfamily",
    "Genus (BLASTn)" = "genus")


# taxonomy_tbl_bckp2 <- taxonomy_tbl

#10- bind tax rank cols to DB_tbl ----
blast_res_tax <- left_join(x = blast_res_full, 
                           y = taxonomy_tbl,
                           by = "query_taxID")

blast_res_tax[is.na(blast_res_tax$`Superkingdom (BLASTn)`),] %>% View()

all_ps_tbl %>% unique()

blast_res_tax %>% filter(`Genus (BLASTn)` %in% c(NA)) %>% View()

blast_res_full %>% colnames()

## Fim do codigo do Heron

## Melhorando o output do BLASTr ----

blast_res_tax_fil <- blast_res_tax %>% 
  filter(!"Superkingdom (BLASTn)" %in% NA) %>% ## retirando NA
  mutate("BLASTn pseudo-score" = (`1_indentity`*`1_qcovhsp`/100)) %>% 
  rename("ASV (Sequence)" = "Sequence") %>% 
  left_join(less_final_results,
            by = c("ASV (Sequence)")) %>% 
  select(`Class (BLASTn)`,`blast ID`,`BLASTn pseudo-score`,`ASV absolute abundance`,
         `Relative abundance on sample`,`Relative abundance to all samples`,
         `Sample total abundance`,`Primer expected length`, `Superkingdom (BLASTn)`, 
         `Kingdom (BLASTn)`,`Phylum (BLASTn)`,`Subphylum (BLASTn)`,`Class (BLASTn)`,`Subclass (BLASTn)`,
         `Order (BLASTn)`,`Suborder (BLASTn)`,`Family (BLASTn)`, `Subfamily (BLASTn)`,`Genus (BLASTn)`,
         `Unique_File_name`,`query_taxID`, `max_tax`, `ASV header`, `ASV (Sequence)`, `1_subject header`,
          `1_subject`, `1_indentity`,`1_length`,`1_mismatches`,`1_gaps`,`1_query start`, `1_query end`,
         `1_subject start`,`1_subject end`,`1_e-value`,`1_bitscore`,`1_qcovhsp`,`1_staxid`,`1_DB`,
         `2_subject header`,`2_subject` ,`2_indentity`,`2_length`,`2_mismatches`, `2_gaps`,`2_query start`,`2_query end`,
         `2_subject start`,`2_subject end`,`2_e-value`,`2_bitscore`,`2_qcovhsp`,`2_staxid`,`2_DB`,
         `3_subject header`,`3_subject`,`3_indentity`,`3_length`,`3_mismatches`,`3_gaps`,`3_query start`,
         `3_query end`,`3_subject start`,`3_subject end`,`3_e-value`,`3_bitscore`,`3_qcovhsp`,`3_staxid`,`3_DB`, `OTU`)             
          
# Salvando em excel

blast_res_tax_fil %>% 
write_xlsx(path = paste0(results_path, Sys.Date(),"_blast_res_tax_fil.xlsx"), 
           col_names = TRUE,
           format_headers = TRUE)



