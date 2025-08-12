#Zizhen Zhong Permafrost S25 Project Analysis  
#June 7, 2025 - August 12, 2025 

###----Packages Used----
library(tidyverse)
library(dplyr)

###----Funtion1 ----
#test <- read.csv("T.cylindrosporum_LO42_blastpresultv2", header = FALSE,sep = "#") 
process_blastp_results <- function(file_path){
  
  originalfile <- read.csv(file_path, header = FALSE,sep = "#") 
  
  sep <- grepl("^sp", originalfile$V1)
  
  D1 <- originalfile[!sep,] %>% 
    separate(col = V6, sep = "\t",into = paste0("V6_", 1:12), convert = TRUE)
  
  colnames(D1)[1] <- "Protein_ID"
  colnames(D1)[2] <- "PHI_ID"
  colnames(D1)[3] <- "Gene_name"
  colnames(D1)[5] <- "Pathogen_species"
  colnames(D1)[6] <- "Phenotype"
  
  
  
  D2 <- originalfile[sep,] %>% 
    separate(col = V1, sep = "\t", into = paste0("V6_", 1:12), convert = TRUE) %>%
    separate(V6_1, sep = "\\|", into = c("prefix", "Protein_ID", "Gene_name"),) %>% 
    select(-prefix)
  
  D2<-D2[,-c(14:18)]
  
  combined_data <- bind_rows(D1, D2)
  
  new_names <- c("FunannotateIDs", "Identity", "Total_length", "Mismatch", "Gapsopen","Qstart", "Qend", "Sstart", "Send", "evalue", "bitscore")
  colnames(combined_data)[7:17] <- new_names
  
  combined_data$Gene_name <- dplyr::recode(combined_data$Gene_name,"EFG1_CANAL" = "EFG1","PHR1_CANAL" = "PHR1","PHR2_CANAL" = "PHR2","PLB3_CANAX" = "PLB3","PLB1_CANAL" = "PLB1","CARP6_CANAL" = "SAP6","CARP1_CANAL" = "SAP1","CARP2_CANAL" = "SAP2","CARP3_CANAL" = "SAP3","CARP4_CANAL" = "SAP4","CARP5_CANAL" = "SAP5","CARP7_CANAL" = "SAP7","TEC1_CANAL" = "TEC1")
  
  combined_data <- combined_data[,-c(2,4)]
  
  
  combined_data$Pathogen_species[combined_data$Gene_name %in% c("TEC1","SAP7","SAP5","SAP4","SAP3","SAP2","SAP1","SAP6","PLB1","PLB3","PHR2","PHR1","EFG1")] <- "Candida_albicans"
  
  return(combined_data)
}


###---- Function2 ----

explore_pathogenic_genes <- function(blast_data,output_name, output_name2) {
  
  original_geneids <- unique(blast_data$FunannotateIDs) #original numbers of fungal genes match to known pathogenic genes that passed the 1e-5 paramter in blastp process but many may be redundant because blast tends to generate multiple matches to the same place 
  
  # Primary filtering for significance 
  filtered_data <- blast_data %>% 
    group_by(FunannotateIDs) %>% 
    filter(Identity >= 50) %>% 
    ungroup() # provides putative homolog to the pathogenic genes, but 1 gene ID may have many matches 
  
  filtered_data$Gene_name <- toupper(filtered_data$Gene_name)
  filtered_data$Gene_name <- sub(" \\(.*\\)", "", filtered_data$Gene_name)
  
  potential_sig_numb <- nrow(filtered_data)
  
  #cal <- filtered_data %>% filter(Gene_name == "Calcium-transporting_ATPase_3")
  
  #print(cal)
    
  #more biologically meaningful set of genes for pathogenic function 
  high_identity_genes <- filtered_data %>%
    filter(Identity >= 70) %>% 
    group_by(FunannotateIDs) %>%
    slice_max(order_by = bitscore, n = 1, with_ties = FALSE) %>% 
    ungroup()
  #retain only matches with v. high similarity identity per FunannotateIDs 
  
  assign(output_name, high_identity_genes, envir = .GlobalEnv)
  
  high_conf_numb <- nrow(high_identity_genes)
  
  # Get Potential Conserved Genes that are matched to different protein + species since the protein ID is the uniprot ID is unique 
  #How many different pathogenic proteins are matched to this single gene? 
  match_counts <- filtered_data %>%
    group_by(FunannotateIDs) %>%
    summarise(MatchCount = n_distinct(Protein_ID)) %>%
    arrange(desc(MatchCount))
  
  #Those appearing in ≥5 matched to more than 5 distinct proteins IDs
  #these are likely large conserved gene family 
  top_df <- filtered_data %>%
    group_by(FunannotateIDs) %>%
    summarise(MatchCount = n_distinct(Protein_ID)) %>%
    arrange(desc(MatchCount)) %>% 
    filter(MatchCount >= 5)
  
  view(top_df)
  
  top_IDs <- top_df$FunannotateIDs
  
  top_matches <- filtered_data %>%
    filter(FunannotateIDs %in% top_IDs) %>%
    arrange(desc(Identity))
  
  # Calculate mean identity per gene
  mean_identity <- sapply(top_IDs, function(id) {
    mean(top_matches$Identity[top_matches$FunannotateIDs == id], na.rm = TRUE)
  })
  
  # Find potential gene duplicates
  putative_duplicate_candidates <- filtered_data %>%
    group_by(Gene_name) %>%
    filter(n_distinct(FunannotateIDs) > 1) %>%
    arrange(Gene_name, desc(Identity))
  
  high_confident_duplicates <- filtered_data %>%
    filter(Identity >= 80) %>%
    group_by(Gene_name) %>%
    filter(n_distinct(FunannotateIDs) > 1) %>%
    arrange(Gene_name, desc(Identity)) %>%
    ungroup()
  
  high_confident_duplicates <- high_confident_duplicates %>% 
    select(Protein_ID, Gene_name, FunannotateIDs, Identity, Total_length)
  
  high_dup <- length(unique(high_confident_duplicates$Gene_name))
  
  if (nrow(high_confident_duplicates) == 0) {
    cat("No confident gene duplicates with multiple FunannotateIDs and protein matches found.\n")
  } else {
    cat("Confident gene duplicates found, but need to further confirmation:\n")
    
    print(high_confident_duplicates)
  }
  
  view(high_confident_duplicates)
  
  filtered_joined<- high_identity_genes %>%
    left_join(data1 %>% 
       select(Protein_ID, Function, GO_annotation, Host_descripton),
      by = "Protein_ID"
    ) 
  
  # STEP 1: Count how many Protein_IDs per FunannotateID
  id_counts <- filtered_joined %>%
    group_by(FunannotateIDs) %>%
    mutate(ID_count = n()) %>%
    ungroup()
  
  # STEP 2: Split into duplicates vs. unique
  duplicates <- id_counts %>%
    filter(ID_count > 1) %>%
    group_by(FunannotateIDs) %>%
    slice_max(order_by = bitscore, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  uniques <- id_counts %>%
    filter(ID_count == 1)
  
  # STEP 3: Combine both back
  cleaned_filtered_df <- bind_rows(duplicates, uniques)
  
  
  assign(output_name2, cleaned_filtered_df, envir = .GlobalEnv)
  
    
  # Check for specific virulence genes
  virulence_genes <- c("EFG1", "ALS3", "PLB1","PLB3", "PHR1","PHR2", 
                       "SAP1", "SAP2", "SAP3", "SAP4", "SAP5", "SAP6", "SAP7")
  virulence_hits <- filtered_data %>%
    filter(Gene_name %in% virulence_genes) %>%
    select(Gene_name, FunannotateIDs, Identity)
  
  # Print summary
  cat("Originally, there are", length(original_geneids), "fungal genes found to match during `primary` BLASTp to PHI-databse. These matches may include weak/distant homologs. \n\n")
  
  cat("Out of", potential_sig_numb, "moderately significant matches with >= 50% identiy + >= 50 bitscore filter parameter,", high_conf_numb, "matches show >= 70 identity.\n\n")

  cat("Thus, there are", length(unique(filtered_data$FunannotateIDs)), " putative pathogenic gene homologs from this isolate.\n\n")

  cat("There are", length(unique(high_identity_genes$FunannotateIDs)), "high confidence genes that has conserved functions as query pathogenic genes found in this isolate. \n\n")

  cat("With", length(top_IDs), "genes show match with more than 5 different protein IDs 'suggesting strong signal to be part of a large protein family'. The greater number in matched proteins IDs suggest the more conserved and important. These FunannotateID of:", top_IDs, "\n\n")
  
  cat("There are", length(unique(putative_duplicate_candidates$Gene_name)), "are putative gene duplicates.\n\n")
  
  cat(high_dup, "possible duplicated genes, but needs further investigation for confirmation.\n\n")
  
  
  
  # Return all results in a structured list
  return(list(
    Initial_blast_pathogenic_match= filtered_data,
    Conserved_genes = top_matches,
    Conserved_gene_frequency = top_df,
    Mean_identity = mean_identity,
    Putative_Gene_duplicates = putative_duplicate_candidates,
    True_Gene_duplicates = high_confident_duplicates,
    Similar_key_virulence_C.albicans = virulence_hits
  ))
}



##----T.cylindrosporum Blastp Result Check----
T.cylindrosporum_CBS71870 <- process_blastp_results("T.cylindrosporum_CBS71870_blastpresultv2") 
explore_pathogenic_genes(T.cylindrosporum_CBS71870, "T.cylindrosporum_CBS71870_filtered", "T.cylindrosporum_CBS71870_cleaned")
#There are 235 potential pathogenic gene homologs in CBS71870 isolate, 104/235 are high confident potential pathogenic genes.
#8 potential conserved pathogenic genes, only 1/8 "FUN_002513-T1" showed 17 matched queries with 84.22% mean of identity 
#68 putative gene duplicates, 3/68 are potential gene duplicates:ACL1 FGA2 Ste11, however their loci are too close together, I would say they are likely not duplicated... 
table(T.cylindrosporum_CBS71870_cleaned$Host_descripton)
#68 monocots hosts 
#21 eudicot hosts
#5 moths hosts / 1 grasshopper
#8 rodents 
#1 nematodes 


T.cylindrosporum_CBS71870_cleaned %>% 
  filter(Host_descripton == "nematodes")
#out of all the 104 pathogenics genes found in CBS718.70, 87% of the gene homologs were found to have plant hosts, however,there is PHI-base database bias where it contains alot more plant pathogenic genes in this database  
#1. PHI-base database bias in plant host information
#2. only means it shares more pathogenic machineries with other plant pathogens but likely still an insect pathogen??? 

unique(T.cylindrosporum_CBS71870_cleaned$Function)
T.cylindrosporum_CBS71870_cleaned$Function <- trimws(T.cylindrosporum_CBS71870_cleaned$Function)

#according to names, assign the function into few categories for better investigation
T.cylindrosporum_CBS71870_cleaned <- T.cylindrosporum_CBS71870_cleaned %>%
  mutate(Function_Group = case_when(
    grepl("MAPK|GTPase|G alpha|G beta|kinase", Function, ignore.case = TRUE) ~ "Signal Transduction",
    grepl("transcription|histone|bZip", Function, ignore.case = TRUE) ~ "Transcription Regulation",
    grepl("kinase|isomerase|aminotransferase|dehydrogenase|lyase", Function, ignore.case = TRUE) ~ "Metabolism & Energy",
    grepl("tubulin|actin|motor|clathrin", Function, ignore.case = TRUE) ~ "Cytoskeleton/Transport",
    grepl("pathogenicity|protease|oxidase", Function, ignore.case = TRUE) ~ "Pathogenicity-related",
    grepl("autophagy|BAX|regulator", Function, ignore.case = TRUE) ~ "Stress/Development",
    grepl("hypothetical|uncharacterized", Function, ignore.case = TRUE) ~ "Unknown/Hypothetical",
    TRUE ~ "Other"
  ))

T.cylindrosporum_CBS71870_barplot <- ggplot(T.cylindrosporum_CBS71870_cleaned, aes(x = Function_Group)) +
  geom_bar(fill = "steelblue") +
  coord_flip() +
  labs(title = "Functional Grouping of T.cylindrosporum_CBS718.70 Pathogenic Genes",
       x = "Functional Category",
       y = "Number of Genes") +
  theme_minimal()

T.cylindrosporum_CBS71870_barplot
table(T.cylindrosporum_CBS71870_cleaned$Function_Group)
#Cytoskeleton/Transport      Metabolism & Energy 
#4                       7
#Other    Pathogenicity-related 
#48                        5 
#Signal Transduction       Stress/Development 
#22                        7 
#Transcription Regulation     Unknown/Hypothetical 
#8                       3


#view(T.cylindrosporum_CBS71870_cleaned)
#install.packages("ontologyIndex")
#library(ontologyIndex)
#GO <- get_ontology("go-basic.obo", extract_tags = "everything")
#head(GO$id)
#GO
#molecular function GO:0003674
#biological process GO:0008150
#cellular component GO:0005575 

#T.cylindrosporum_CBS71870_cleaned$GO_annotation
#GO$name
#manually look up Uniprot to find and enter the corresponding GO term into the *cleaned_df 

#T.cylindrosporum_CBS71870_cleaned$GO_annotation[T.cylindrosporum_CBS71870_cleaned$Protein_ID == "A0A0A2JQN4"] <- "GO:0008150"


Conserved_gene_names <- T.cylindrosporum_CBS71870_filtered$Gene_name[T.cylindrosporum_CBS71870_filtered$FunannotateIDs=="FUN_002513-T1"]
#"CMK1"   "PTK1"   "CPMK1"  "MAK2"   "VMK1"   "AMK1"   "Mak2"   "AaFUS3" "Bbmpk1" "MgFus3" "MaMk1"  "ChMK1" 
#data1$GO_annotation[data1$Gene_name %in% Conserved_gene_names]



T.cylindrosporum_IBT41711 <- process_blastp_results("T.cylindrosporum_IBT41711_blastpresultv2") 
explore_pathogenic_genes(T.cylindrosporum_IBT41711, output_name = "T.cylindrosporum_IBT41711_filtered", "T.cylindrosporum_IBT41711_cleaned")
#There are 243 potential pathogenic gene homologs in IBT41711, and 100/391 are high confidence potential pathogenic genes. 
#IBT41711 shows 8 potential conserved genes, only 2/8 genes is high confident potential conserved gene, "FUN_000809-T1" 17 matched query 84.22%, "FUN_009630-T1" 14 matched queries 91.044%
#There are 13 high confident potential gene duplicates / 67 putative gene duplicates. 
#however, only ACL1 seem to be possible duplicates with no short length match 

table(T.cylindrosporum_IBT41711_cleaned$Host_descripton)
#22 eudicots 
#63 monocot hosts 
#1 grasshoppers, 6 moths 
#7 rodents 
#1 nematodes 
#similar to CBS718.70, IBT41711 seem to show more similar pathogenic genes with plant host pathogens 

T.cylindrosporum_IBT41711_cleaned <- T.cylindrosporum_IBT41711_cleaned %>%
  mutate(Function_Group = case_when(
    grepl("MAPK|GTPase|G alpha|G beta|kinase", Function, ignore.case = TRUE) ~ "Signal Transduction",
    grepl("transcription|histone|bZip", Function, ignore.case = TRUE) ~ "Transcription Regulation",
    grepl("kinase|isomerase|aminotransferase|dehydrogenase|lyase", Function, ignore.case = TRUE) ~ "Metabolism & Energy",
    grepl("tubulin|actin|motor|clathrin", Function, ignore.case = TRUE) ~ "Cytoskeleton/Transport",
    grepl("pathogenicity|protease|oxidase", Function, ignore.case = TRUE) ~ "Pathogenicity-related",
    grepl("autophagy|BAX|regulator", Function, ignore.case = TRUE) ~ "Stress/Development",
    grepl("hypothetical|uncharacterized", Function, ignore.case = TRUE) ~ "Unknown/Hypothetical",
    TRUE ~ "Other"
  ))

T.cylindrosporum_IBT41711_barplot<- ggplot(T.cylindrosporum_IBT41711_cleaned, aes(x = Function_Group)) +
  geom_bar(fill = "blue") +
  coord_flip() +
  labs(title = "Functional Grouping of T.cylindrosporum_IBT41711 Pathogenic Genes",
       x = "Functional Category",
       y = "Number of Genes") +
  theme_minimal()

T.cylindrosporum_IBT41711_barplot
table(T.cylindrosporum_IBT41711_cleaned$Function_Group)
#Cytoskeleton/Transport      Metabolism & Energy 
#3                        6 
#Other    Pathogenicity-related 
#12                        3 
#Signal Transduction       Stress/Development 
#19                        4 
#Transcription Regulation     Unknown/Hypothetical 
#3                        2


T.cylindrosporum_IBT41711_filtered$Gene_name[T.cylindrosporum_IBT41711_filtered$FunannotateIDs=="FUN_000809-T1"]
#"CMK1"   "PTK1"   "CPMK1"  "MAK2"   "VMK1"   "AMK1"   "Mak2"   "AaFUS3" "Bbmpk1" "MgFus3" "MaMk1"  "ChMK1" 



T.cylindrosporum_LO42 <- process_blastp_results("T.cylindrosporum_LO42_blastpresultv2")
explore_pathogenic_genes(T.cylindrosporum_LO42, output_name = "T.cylindrosporum_LO42_filtered", "T.cylindrosporum_LO42_cleaned")
#there are 246 potential pathogenic genes in LO42, 101 are high confidence potential pathogenic genes. 
#7 potential conserved pathogenic genes, 1/7 showed high confidence with >= 80 identity "FUN_002651-T1" 84.22% identity
#73 putative gene duplicates, but only 1 gene is true potential gene duplicate ACL1, however, base on the length of the matched sequence, it doesnt seem to be a duplicated gene


table(T.cylindrosporum_IBT41711_cleaned$Host_descripton)
#matched to 63 monocots plant host pathogen's genes with >= 70% identity 
#22 eudicots 
#1 grasshoppers , 6 moths 
#7 rodents 
#1 nematodes 


T.cylindrosporum_LO42_cleaned <- T.cylindrosporum_LO42_cleaned %>%
  mutate(Function_Group = case_when(
    grepl("MAPK|GTPase|G alpha|G beta|kinase", Function, ignore.case = TRUE) ~ "Signal Transduction",
    grepl("transcription|histone|bZip", Function, ignore.case = TRUE) ~ "Transcription Regulation",
    grepl("kinase|isomerase|aminotransferase|dehydrogenase|lyase", Function, ignore.case = TRUE) ~ "Metabolism & Energy",
    grepl("tubulin|actin|motor|clathrin", Function, ignore.case = TRUE) ~ "Cytoskeleton/Transport",
    grepl("pathogenicity|protease|oxidase", Function, ignore.case = TRUE) ~ "Pathogenicity-related",
    grepl("autophagy|BAX|regulator", Function, ignore.case = TRUE) ~ "Stress/Development",
    grepl("hypothetical|uncharacterized", Function, ignore.case = TRUE) ~ "Unknown/Hypothetical",
    TRUE ~ "Other"
  ))

T.cylindrosporum_LO42_barplot<- ggplot(T.cylindrosporum_LO42_cleaned, aes(x = Function_Group)) +
  geom_bar(fill = "blue") +
  coord_flip() +
  labs(title = "Functional Grouping of T.cylindrosporum_LO42 Pathogenic Genes",
       x = "Functional Category",
       y = "Number of Genes") +
  theme_minimal()

T.cylindrosporum_LO42_barplot
table(T.cylindrosporum_LO42_cleaned$Function_Group)
#Cytoskeleton/Transport      Metabolism & Energy 
#3                        5 
#Other    Pathogenicity-related 
#13                        4 
#Signal Transduction       Stress/Development 
#18                        4 
#Transcription Regulation     Unknown/Hypothetical 
#3                        2 


T.cylindrosporum_LO42_filtered$Gene_name[T.cylindrosporum_LO42_filtered$FunannotateIDs =="FUN_002651-T1"]
#"CMK1"   "PTK1"   "CPMK1"  "MAK2"   "VMK1"   "AMK1"   "Mak2"   "AaFUS3"
#"Bbmpk1" "MgFus3" "MaMk1"  "ChMK1"

LO42_conserved_genes <- T.cylindrosporum_LO42_filtered$Gene_name[T.cylindrosporum_LO42_filtered$FunannotateIDs =="FUN_002651-T1"]

data1$Function[data1$Gene_name %in% LO42_conserved_genes]
#these are matched to MAP kinase category 
data1$GO_annotation[data1$Gene_name %in% LO42_conserved_genes]
#GO:0009405; GO:0044409";GO:0044412;GO:0004707;GO:0030437;GO:0009847



###----Funtion3----

process_all_blastp_files <- function(pattern = "_blastpresultv2$") {
  files <- list.files(pattern = pattern)
  
  for (x in files) {
    cat("🔍 Processing file:", x, "\n\n")
    
    obj_name <- gsub("_blastpresultv2$", "", x)
    
    result <- process_blastp_results(x)
    assign(obj_name, result, envir = .GlobalEnv)
    
    filtered_output <- explore_pathogenic_genes(result, output_name = paste0(obj_name, "_filtered"), output_name2=paste0(obj_name, "_cleaned"))
    
    # This will ensure full output is shown:
    cat("\n📢 Summary Output for", obj_name, "\n\n")
    print(filtered_output)  # <- This is key to seeing $Initial_blast_pathogenic_match etc.
    
    cat("\n----------------------------------------\n\n")
  }
}

###----T.album----
process_all_blastp_files(pattern = "_blastpresultv2$")
#CBS86973
  #243 potential pathogenic gene homologs in this isolate, 109 are high confident potential pathogenic gene homologs.
  #8 potential putative conserved pathogenic genes, only 1/8 gene "FUN_007284-T1" showed mean 84.02% identity with 17 different queries 
  #70 putative gene duplicates, but 4/70 true potential gene duplicates at different loci, these gene names are ACL1, Cdc25, FGA2, FgGT2 
  #only FGGT2 showed 3 different loci with long enough total length, when compared
  #showed 2 similar key virulent genes as C.albicans: EFG1, Phr1 both 70% identity 
T.album_CBS86973_filtered$Gene_name[T.album_CBS86973_filtered$FunannotateIDs=="FUN_007284-T1"]

table(T.album_CBS86973_cleaned$Host_descripton)
#22 eudicots 
#1 grasshoppers, 6 moths
#68 monocots 
#11 rodents 
#1 nematodes 

T.album_CBS86973%>% 
  filter(Protein_ID %in% nematode_host_proteinIDs) %>% 
  group_by(Protein_ID) %>% 
  slice_max(order_by = Identity, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  select(Protein_ID,Identity)
#there's 2 similarity with > 70 identity 
#G1XD68 P00445 S8A502


T.album_CBS86973_nematode_genes <- T.album_CBS86973 %>% 
  filter(Protein_ID %in% c("S8A502", "P00445","G1XD68")) %>% 
  group_by(FunannotateIDs)
view(T.album_CBS86973_nematode_genes)
#FUN_004056-T1 78.571% identity to AoRac_(AOL_s00079g171) of Orbilia_oligospora with 140 length alignment 
#FUN_007117-T1 78.125% identity to MhSlt2_(XP_011115715) with 416 alignment 
#FUN_002782-T1 73.377% with Sod1 of Saccharomyces_cerevisiae with 154 alignement in length 


T.album_CBS86973_cleaned <- T.album_CBS86973_cleaned %>%
  mutate(Function_Group = case_when(
    grepl("MAPK|GTPase|G alpha|G beta|kinase", Function, ignore.case = TRUE) ~ "Signal Transduction",
    grepl("transcription|histone|bZip", Function, ignore.case = TRUE) ~ "Transcription Regulation",
    grepl("kinase|isomerase|aminotransferase|dehydrogenase|lyase", Function, ignore.case = TRUE) ~ "Metabolism & Energy",
    grepl("tubulin|actin|motor|clathrin", Function, ignore.case = TRUE) ~ "Cytoskeleton/Transport",
    grepl("pathogenicity|protease|oxidase", Function, ignore.case = TRUE) ~ "Pathogenicity-related",
    grepl("autophagy|BAX|regulator", Function, ignore.case = TRUE) ~ "Stress/Development",
    grepl("hypothetical|uncharacterized", Function, ignore.case = TRUE) ~ "Unknown/Hypothetical",
    TRUE ~ "Other"
  ))

table(T.album_CBS86973_cleaned$Function_Group)
#Cytoskeleton/Transport      Metabolism & Energy 
#3                        6 
#Other    Pathogenicity-related 
#13                        5 
#Signal Transduction       Stress/Development 
#19                        4 
#Transcription Regulation     Unknown/Hypothetical 
#3                        2 


#IQ58 
  #237 potential pathogenic gene homologs in this isolate, 106/237 are high confident potential pathogenic gene homologs. 
  #8 potential conserved pathogenic genes, only 1/8 is high confident "FUN_005476-T1" 17 matched queries with mean 84.02%
  #65 putative gene duplicates, only 2/65 are potential gene duplicates/paralogs, gene names are: ACL1, FGA2
  #however, both are too closely located, may not be true gene duplicates, thus none 

table(T.album_IQ58_cleaned$Host_descripton)
#22 eudicots 
#68 monocots 
#1 grasshopper
#5 moths
#9 rodents 
#1 nematodes 

T.album_IQ58%>% 
  filter(Protein_ID %in% nematode_host_proteinIDs) %>% 
  group_by(Protein_ID) %>% 
  slice_max(order_by = Identity, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  select(Protein_ID,Identity)
#there's 2 similarity with > 70 identity 
#G1XD68 P00445 S8A502


T.album_IQ58_nematode_genes <- T.album_IQ58 %>% 
  filter(Protein_ID %in% c("S8A502", "P00445","G1XD68")) %>% 
  group_by(FunannotateIDs)
view(T.album_IQ58_nematode_genes)
#FUN_008129-T1 78.652% 	AoRac_(AOL_s00079g171 178 alignment #FUN_000468-T1 77.885% MhSlt2_(XP_011115715) 416 alignment
#FUN_006080-T1 71.333% Sod1 150 alignment 

#MS506
  #238 potential pathogenic gene homologs, 96/238 are high confident potential pathogenic genes
  #7 potential putative conserved pathogenic genes, only 1/7 "FUN_004975-T1" with 17 queries matched mean of 84.02% identity 
  #2/63 potential gene duplicates are high confident, ACL1, FGA2; both too closely located thus none duplicates 


table(T.album_MS506_cleaned$Host_descripton)
#20 eudicots 
#61 monocots 
#5 moths, 1 grasshoppers
#7 rodents 
#1 nematode
#1 rabbits and hares 


###----T.amazonense----
process_all_blastp_files(pattern = "_blastpresultv2$")
#LA108
  #233 potential pathogenic gene homologs, 102/233 are high confident potential pathogenic gene homologs
  #8 potential conserved genes, 1/8 high confidence conserved gene "FUN_007426-T1" with 17 matched queries mean 84.02% identity 
  #1/61 potential putative gene duplicates is high confident, ACL1; too closely located + huge total length variation; Thus NONE 
  #similar to C.albicans EFG1 70% identity 

table(T.amazonense_LA108_cleaned$Host_descripton)
#19 eudicots 
#66 monocots 
#1 grashoppers, 6 moths 
#9 rodents 
#1 nematodes 



#MS503 
  #105/238 are high confident potential pathogenic gene homologs 
  #1/7 potential conserved pathogenic gene homolog is high confident FUN_006825-T1, with 17 matched queries mean of 83.989% identity 
  #3/66 putative gene duplicates are high confidence, ACL1,Cdc25, FGA2; CDC25 may be possible, however the variation in length matched between the 2 loci leads to uncertainty 
  #similar to C.albicans EFG1 77.6%, Phr1 50.8%
table(T.amazonense_MS503_cleaned$Host_descripton)
#19 eudicots 
#69 monocots 
#7 moths, 1 grasshoppers
#8 rodents
#1 nematodes 



###----T.endophyticum----
process_all_blastp_files(pattern = "_blastpresultv2$")
#MX560 
  #96/235 high confidence potential pathogenic gene homologs 
  #1/6 potential conserved pathogenic genes maybe true conserved gene, "FUN_003838-T1" 17 matched queries mean of 84.02%
  #2/63 are true potential gene duplicates, ACL1, FGA2; 
  #NONE duplicates, too closely located 
  #similar to C.albicans, EFG1 60.3%, Phr2 50.5%
table(T.endophyticum_MX560_cleaned$Host_descripton)
#20 eudicots 
#60 monocots 
#1 grasshoper, 6 moths
#8 rodents
#1 nematodes


#MX66
  #100/239 are high confidence potential pathogenic gene homologs 
  #1/8 potential conserved pathogenic genes are confident, FUN_005085-T1 17 queries mean of 84.02% identity 
  #1/65 is true potential gene duplicates, ACL1 
  #NONE, too closely located + variation in total length 
  #simialr to C.albicans, EFG1 70.1%, Phr 50.9% 
table(T.endophyticum_MX66_cleaned$Host_descripton)
#21 eudicots 
#61 monocots 
#1 grasshopers, 6 moths
#9 rodents 
#1 nematode
#1 rabbits and hares 


###----T.inflatum----
process_all_blastp_files(pattern = "_blastpresultv2$")
#CBS56784
  #100/238 high confidence potential pathogenic gene homologs 
  #1/7 potential conserved pathogenic genes is high, "FUN_005733-T1" mean 84.187%identity match to 17 queries 
  #2/65 true potential gene duplicates, ACL1,FGA2
  #NONE dup, too closely located 
  #similar to C.albicans, EFG1 54.5%
table(T.inflatum_CBS56784_cleaned$Host_descripton)
#20 eudicots 
#66 monocots 
#1 grasshoper, 5 moths 
#7 rodents 
#1 nematodes 

#CBS71470
  #104/236 high confidence potential pathogenic gene homologs 
  #1/7 true potential conserved pathogenic genes, "FUN_003952-T1" mean 84.187% 17 queries 
  #2/65 true potential gene duplicates, ACL1,FGA2 
  #NONE dup, too closely located 
  #similar to C.albicans EFG1 59.6% 
table(T.inflatum_CBS71470_cleaned$Host_descripton)
#19 eudicots 
#69 monocots 
#1 grasshoper, 6 moths 
#8 rodents 
#1 nematodes 


#CBS82470
  #105/230 high confidence potential pathogenic gene homologs 
  #1/7 potential conserved pathogenic gene is high "FUN_005354-T1" mean 84.187% 17 queries 
  #2/63 true potential gene duplicates, ACL1,FGA2
  #NONE dup 
  #similar C.albicans EFG1 54.5% 
table(T.inflatum_CBS82470_cleaned$Host_descripton)
#19 eudicots 
#71 monocots 
#1grass 5 moths 
#8 rodents
#1 nematodes 


#NBRC31671
  #102/238 high confidence potential pathogenic gene homologs 
  #1/6 high confident potential conserved pathogenic gene, "FUN_005374-T1" mean 84.187% identity 17 queries 
  #1/68 is true potential gene duplicates, ACL1
  #NONE dup 
  #similar to C.albicans EFG1 59.6% 
table(T.inflatum_NBRC31671_cleaned$Host_descripton)
#18 eudicots 
#69 monocots 
#1 grasshoper, 6 moths 
#7 rodents 
#1 nematodes 


#NBRC31975 
  #103/240 high confidence potential pathogenic gene homologs 
  #1/7 high confidence potential conserved pathogenic genes, "FUN_005679-T1" mean 84.187% identity 17 queries 
  #1/65 is true potential gene duplicates, ACL1 
  #NONE dup 
  #similar to C.albicans EFG1 70.1% 
table(T.inflatum_NBRC31975_cleaned$Host_descripton)
#19 eudicots
#68 monocots 
#1 grass, 5 moths
#9 rodents
#1 nematodes 


#NRRL8044
  #101/235 high confidence potential pathogenic gene homologs 
  #1/7 is true potential conserved pathogenic genes, "FUN_002555-T1" mean 84.187% 17 queries 
  #2/64 are true potential gene duplicates, ACL1, FGA2 
  #NONE dup too closely located 
  #similar C.albican EFG1 52.8%
table(T.inflatum_NRRL8044_cleaned$Host_descripton)
#20 eudicots
#66 monocots 
#1 grasshoper, 5 moths
#8 rodents 
#1 nematodes 

#Q84
  #105/235 high confidence potential pathogenic gene homologs 
  #1/8 high confident potential conserved pathogenic gene, "FUN_008863-T1" mean 86.311% 17 queries 
  #3/68 may be true potential gene duplicates, ACL1,BrlA, FGA2 
  #BRLA may be or may not be dup, however total length variation is large while the loci of the gene ID is still consider to be closely linked 
table(T.inflatum_Q84_cleaned$Host_descripton)
#20 eudicots
#67 monocots
#7 moths 
#10 rodents 
#1 nematodes 


###----Other Tolypocladium isolates----
process_all_blastp_files(pattern = "_blastpresultv2$")
#T.capitatum_CBS113982
  #48/219 high confidence potential pathogenic gene homologs 
  #1/8 high confidence potential conserved pathogenic gene, "FUN_006101-T1" mean 83.722% 17 queries matched 
  #2/54 true potential gene duplicates, ACL1, MgHog1
  #similar to C.albicans EFG1 70.1 
table(T.capitatum_CBS113982_cleaned$Host_descripton)
#20 eudicots
#62 monocots
#5 moths
#7 rodents 
#1 nematodes 

#T.geodes_CBS72370 
  #55/245 high confidence potential pathogenic gene homologs 
  #1/8 high confidence potential conserved pathogenic gene, "FUN_007146-T1" mean 83.735% 17 queries matched 
  #2/67 high confidence true potential gene duplicates, ACL1, FGA2 
  #simialr C.albicans EFG1 70.1% 
table(T.geodes_CBS72370_cleaned$Host_descripton)
#23 eudicots 
#67 monocots 
#1 grasshopers, 6 moths 
#7 rodents 
#1 nematodes 
#1 rabbits and hares 


#T.guangdongense_GD115
  #56/238 high confidence potential pathogenic gene homologs 
  #1/9 high confidence potential conserved pathogenic genes, "FUN_002528-T1" mean 84.04% 17 
  #4/60 high confidence true potential gene duplicates, ACL1,CpCOT1,FGA2,Rac
  #similar EFG1 70.1%
table(T.guangdongense_GD115_cleaned$Host_descripton)
#20 eudicots
#66 monocots 
#1 grass, 5 moths 
#8 rodents 


#T.nubicola_CBS56884
  #56/230 high confidence potential pathogenic gene homologs 
  #1/9 high confidence potential conserved pathogenic gene, "FUN_006831-T1" mean 84.22% 17 
  #2/68 high confidence true potential gene duplicates, ACL1,Cdc25 
  #simialr EFG1 69.2% 
table(T.nubicola_CBS56884_cleaned$Host_descripton)
#18 eudicots
#69 monocots
#1 grass, 7 moths
#8 rodents
#1 nematodes 


#T.pustulatum_CBS110433
  #50/236 high confidence potential pathogenic gene 
  #1/7 potential conserved pathogenic genes "FUN_004664-T1" mean 84.005% 17 queries 
  #2/66 true potential gene duplicates,ACL1,FGA2 
  #similar EFG1 61% 
table(T.pustulatum_CBS110433_cleaned$Host_descripton)
#20 eudicots 
#65 monocots 
#1 grass, 6 moths 
#8 rodents 
#1 nematodes 


#T.sp_sup5_PDA1
  #61/251 high confidence potential pathogenic genes
  #1/8 potential conserved pathogenic genes "FUN_007837-T1" 84.02% 17 queries 
  #5/77 true potential gene duplicates, ACL1,ACL2,FGA2,FgGT2,MoSAS3_(MGG_04615)
  #simialr EFG1 65.1%
table(T.sp_sup5_PDA1_cleaned$Host_descripton)
#21 eudicots 
#74 monocots 
#1 grasshoper, 6 moths
#6 rodents 
#1 nematodes 
#2 rabbits and hares 


#T.tundrense_CBS56984
  #53/239 high confidence potential pathogenic gene
  #1/7 high confidence potential conserved pathogenic gene "FUN_001534-T1" mean 84.22% 17 queries 
  #1/68 true potential gene duplicate, ACL1 
  #similar EFG1 69.2% 
table(T.tundrense_CBS56984_cleaned$Host_descripton)
#19 eudicots 
#67 monocots 
#1 grasshopers, 6 moths 
#7 rodents 


#T.ovalisporum_CBS70092
  #51/234 high confidence potential pathogenic gene
  #1/7 high confidence potential conserved pathogenic genes,"FUN_005166-T1" mean 84.189% 17 
  #2/72 high confidence true potential gene duplicate, ACL1,FGA2 
  #EFG1 70.1% 
table(T.ovalisporum_CBS70092_cleaned$Host_descripton)
#21 eudicots 
#67 monocots
#5 moth
#5 rodents 
#1 rabbit and hares 


#T.paradoxum_NRBC100945 
  #49/223 high confidence potential pathogenic gene 
  #1/7 high confidence potential conserved pathogenic gene "FUN_003651-T1" mean 84.204% 17 
  #1/53 true potential gene duplicate, ACL1 
  #EFG1 70.1% 
table(T.paradoxum_NRBC100945_cleaned$Host_descripton)
#19 eudicots
#63 monocots 
#1 grasshopers, 5 moths
#9 rodents
#1 nematodes 
#1 rabbits and hares 


#T.salcha_MEA2 
  #52/243 high confidence potential pathogenic gene
  #1/8 potential conserved pathogenic gene "FUN_007197-T1" mean 82.614% 17 
  #2/69 true potential gene duplicates, ACL1, FGA2
  #EFG1 70.1% 
table(T.salcha_MEA2_cleaned$Host_descripton)
#20 eudicots
#66 monocots
#5 moths
#9 rodents 
#1 nematodes 


#T.tropicale_MX337 
  #50/236 high confidence potential pathogenic gene
  #1/9 high confidence potential conserved pathogenic gene "FUN_006928-T1" mean 84.00% 17 
  #1/63 true potential gene duplicate, ACL1
  #EFG1 70.1% Phr1 50.5% 
table(T.tropicale_MX337_cleaned$Host_descripton)
#19 eudicots
#60 monocots
#1 grass, 6 moths
#9 rodents 
#1 nematodes 
#1 rabbits and hares 


#T.ophioglossoides_CBS100239
  #53/232 high confidence potential pathogenic gene homologs
  #1/7 high confidence potential conserved pathogenic gene, "FUN_002051-T1" mean 86.075% identity, 17 queries 
  #3/65 high confidence potential gene duplicate, FGA2,MET6,Ste11
  #EFG1 69.2% 
table(T.ophioglossoides_CBS100239_cleaned$Host_descripton)
#23 eudicots
#63 monocots 
#5 moths 
#7rodents 
#1 nematodes 
#1 rabbits and hares 
T.ophioglossoides_CBS100239_cleaned %>% 
  filter(Protein_ID %in% nematode_host_proteinIDs) %>% 
  group_by(Protein_ID) %>% 
  slice_max(order_by = Identity, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  select(Protein_ID,Identity)

view(T.ophioglossoides_CBS100239_cleaned %>% 
  filter(Protein_ID =="P00445"))
#matched to Sod1 with 70% similarity 


#Thus, The high proportion of monocot and eudicot host matches observed across all Tolypocladium strains likely reflects a bias in the PHI-base database toward plant pathogens, rather than true host specificity of the isolates."

table(data1$Host_descripton)
#the original file of host description are majority monocots and eudicots with 540 + 212 compare to other host types, this means PHI base database is slightli bias with the amount of plant pathogen species pathogogenic genes it has. 


###----Nematode Genes finding in T.cylindrosporum----
nematode_host_proteinIDs <- data1$Protein_ID[data1$Host_descripton=="nematodes"]
#"P19880" "P24813" "P00445" "J9VUY6" "G1XP67" "S8A502" "G1XD68" "G1XP67" "G1X8J9"
#these are genes known to have nematodes associations 

table(T.cylindrosporum_LO42_cleaned$Protein_ID %in% nematode_host_proteinIDs)
#there are blast sequence matches to known nematode host pathogenic genes 

T.cylindrosporum_LO42_cleaned %>% 
  filter(Protein_ID %in% nematode_host_proteinIDs) %>% 
  group_by(Protein_ID) %>% 
  slice_max(order_by = Identity, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  select(Protein_ID,Identity)
#show the highest identity match for each ProteinID matches

T.cylindrosporum_LO42_cleaned$FunannotateIDs[T.cylindrosporum_LO42_cleaned$Protein_ID=="P00445"]
view(T.cylindrosporum_LO42_cleaned %>% 
  filter(Protein_ID=="P00445"))
#"FUN_000208-T1" 70.130% identity, gene name: Sod1 154 alignment length function of ROS oxidize 

#identity >= 70 suggests that they share similar aa sequence, with longer alignment region, more likely a similar function 


T.cylindrosporum_CBS71870_cleaned %>% 
  filter(Protein_ID %in% nematode_host_proteinIDs) %>% 
  group_by(Protein_ID) %>% 
  slice_max(order_by = Identity, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  select(Protein_ID,Identity)

CBS71870_nematode_genes <- T.cylindrosporum_CBS71870 %>% 
  filter(Protein_ID == "P00445") %>% 
  group_by(FunannotateIDs)
view(CBS71870_nematode_genes)
#"FUN_008262-T1" 70% identity to Sod1 with 150 alignment length 



T.cylindrosporum_IBT41711_cleaned %>% 
  filter(Protein_ID %in% nematode_host_proteinIDs) %>% 
  group_by(Protein_ID) %>% 
  slice_max(order_by = Identity, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  select(Protein_ID,Identity)
#P00445


IBT41711_nematode_genes <- T.cylindrosporum_IBT41711_cleaned %>% 
  filter(Protein_ID %in% "P00445") %>% 
  group_by(FunannotateIDs)
view(IBT41711_nematode_genes)
#"FUN_007094-T1" 70.130% identity to Sod1, 154 alignment length 

#All 3 T.cylindrosporum strains showed 70% similarity to Sod1 gene from Saccharomyces_cerevisiae suggesting this gene is a functionally conserved in this species, it aligns with the nematode killing phenotype previously observed in LO42, but the fact that other strains also exhibits similar identity may suggest they all are possible in nematode killing??? 




###----Core Pathogenic Genes in Tolypocladium Genus---- 
files <- ls(pattern = "_cleaned$")
files

pathogenic_gene_list <- lapply(files, function(x) {
  
  obj <- get(x)
  obj[[1]]
})
names(pathogenic_gene_list) <- files
names(pathogenic_gene_list) <- gsub("_cleaned$", "", basename(files))
names(pathogenic_gene_list)

shared_genes <- Reduce(intersect, pathogenic_gene_list)
length(unique(shared_genes)) # 64 core pathogenic protein IDs are matched across Tolypoclaidum genus isolates 
shared_genes


core_genes_info <- data1 %>%
  filter(Protein_ID %in% shared_genes) %>%
  distinct(Gene_name, Protein_ID, Pathogen_species, Host_descripton, Function)


view(core_genes_info)
core_genes_info$Function[core_genes_info$Protein_ID =="Q9UWF0"] <- "regulates cell proliferation, programmed cell death and autophagy"
#https://doi.org/10.1016/j.bbamcr.2013.10.021 
core_genes_info$Function[core_genes_info$Protein_ID =="Q51MW4"] <- "Autophagy-related protein 8"
#https://www.uniprot.org/uniprotkb/Q51MW4/entry
core_genes_info$Function[core_genes_info$Protein_ID =="G4MNE6"] <- "calcium sensor"
#https://pubmed.ncbi.nlm.nih.gov/19536897/#:~:text=calcium%2Dindependent%20manner-,Calcineurin%20regulatory%20subunit%20B%20is%20a%20unique%20calcium%20sensor%20that,(3)%3A612%2D23.
#10.1002/prot.22474
core_genes_info$Function[core_genes_info$Protein_ID =="G4MTZ0"] <- "chaperon"
#https://www.uniprot.org/uniprotkb/A0A5K1JSC2/entry

core_genes_info$Function[core_genes_info$Protein_ID =="A0A0B2XG51"] <- "transcription regulation"
#https://www.uniprot.org/uniprotkb/A0A0B2XG51/entry

core_genes_info$Function[core_genes_info$Protein_ID =="E9EMI7"] <- "transcription regulation"
#https://www.uniprot.org/uniprotkb/E9EMI7/entry

table(core_genes_info$Host_descripton)
#12 of these matched proteins IDs are found to have host of eudicots 
#44 monocots 
#5 moths 
#3 rodents 
#Despite the plant-pathogen bias in PHI-base, the host descriptions indicate that these core pathogenic genes may be involved in interactions with a broad range of hosts, this suggest the range of host this genus may impact is plants, insects, mammalians 

#however, base on previous finding in nematode host, majority of the strains in this genus showed a gene matched >=70% identity to Sod1 gene that was previously experimented to be pathogenic towards nematodes
#strains except T.guangdongense, T.ovalisporum, and T.tundrense 

#This suggests that most species in the genus Tolypocladium may have functional capacity to affect or infect nematode hosts, likely through conserved mechanisms such as ROS detoxification, which plays a critical role in evading host immune responses


#USE BLASTKOALA of LO42 CORE GENES TO FIND KEGG PATHWAYS 
core_genes_info$Function 
#core_genes_info$GO <- NA
#view(core_genes_info)
#core_genes_info$GO[core_genes_info$Function == "regulatory protein"] <- "GO:0006417"
#https://www.ebi.ac.uk/QuickGO/term/GO:0006417
#core_genes_info$GO[core_genes_info$Protein_ID=="O42621"] <- "GO:0004424"
#core_genes_info$GO[core_genes_info$Protein_ID=="O42784"] <- "GO:0001664"
#core_genes_info$GO[core_genes_info$Protein_ID=="	Q8J286"] <- "GO:0005507"




###-----Unique Pathogenic Genes Acorss Species----
unique_per_sample <- lapply(names(pathogenic_gene_list), function(x) {
  this_sample <- pathogenic_gene_list[[x]]
  other_samples <- unlist(pathogenic_gene_list[names(pathogenic_gene_list) != x])
  setdiff(this_sample, other_samples)
  
})
names(unique_per_sample) <- names(pathogenic_gene_list)

unique_per_sample
#T.amazonense MS503 "MOMCM1" "GAA1"
#T.capitatum CBS113982 "CAC1"
#T.guangdongense GD115 "CT-PKAC" "FGRAV2"
#T.inflatum Q84 "CMK1" 
#T.nubicola CBS56884 "KIN4" "FAM1" "MOATG12"
#T.ovalisporum CBS70092 "FGLEU4_(FGSG_12952)"
#T.sp_sup5_PDA1 "MOSAS3_(MGG_04615)"
#T.tropicale MX337 "CPK1" 


library(ggplot2)
#install.packages("ggtext")
library(ggtext)

unique_counts <- data.frame(
  Sample = names(unique_per_sample),
  Count = sapply(unique_per_sample, length)
)
unique_counts$Label <- ifelse(
  unique_counts$Sample == "T.cylindrosporum_LO42_cleaned",
  "<b>T.cylindrosporum_LO42</b>",
  unique_counts$Sample
)

barplot1 <- ggplot(unique_counts, aes(x = Label, y = Count)) +
  geom_bar(stat = "identity", fill = "salmon") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_markdown(angle = 90, hjust = 1),  # enables bold text
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Number of Unique Pathogenic Genes per Isolate",
    x = "Tolypocladium Isolates",
    y = "Number of Unique Genes"
  )

barplot1
#not putting this in poster 



###----Compare Between T.cylindrosporum isolates----
TC_files <- ls(pattern = "^T\\.cylindrosporum_.*_cleaned$")
TC_files
TCpathogenic_gene_list <- lapply(files, function(x) {
  
  obj <- get(x)
  obj[[2]]
})
names(pathogenic_gene_list) <- TC_files
TCshared_genes <- Reduce(intersect, TCpathogenic_gene_list)
TCshared_genes # 64 shared core pathogenic genes across T.cylindrosporum isolates 


library(VennDiagram)
#install.packages("gridBase")
library(grid)
# Prepare gene lists
genes_lo42 <- T.cylindrosporum_LO42_cleaned$Protein_ID
genes_cbs <- T.cylindrosporum_CBS71870_cleaned$Protein_ID
genes_ibt <- T.cylindrosporum_IBT41711_cleaned$Protein_ID

# Plot Venn diagram
venn.plot <- venn.diagram(
  x = list(LO42 = genes_lo42, CBS718.70 = genes_cbs, IBT41711 = genes_ibt),
  filename = NULL, #to plot 
  fill = c("skyblue", "pink", "lightgreen"),
  alpha = 0.5,
  cat.cex = 1.2,
  margin = 0.1,
)
grid.draw(venn.plot)
grid.text("Shared and Unique Pathogenic Genes, T.cylindrosporum Isolates", y = 0.95, gp = gpar(fontsize = 14, fontface = "bold"))

TC_genes_list <- list(
  CBS = genes_cbs,
  IBT = genes_ibt,
  LO42 = genes_lo42
)
unique_CBS <- setdiff(TC_genes_list$CBS, union(TC_genes_list$IBT, TC_genes_list$LO42))
unique_CBS
#"PAX1","CST1","FGBNI4_(FGRRES_07218)","MOSNF4" unique to only CBS718.70 when compared against the other 2 T.cylindrosporum strains 

unique_IBT <- setdiff(TC_genes_list$IBT, union(TC_genes_list$CBS, TC_genes_list$LO42))
unique_IBT 
#"VMK1" unique to only IBT41711

unique_LO42 <- setdiff(TC_genes_list$LO42, union(TC_genes_list$CBS, TC_genes_list$IBT))
unique_LO42
#"MOHIK5" unique to only LO42 

view(T.cylindrosporum_LO42_cleaned %>% 
  filter(Gene_name == "MOHIK5"))
#function : histidine kinase; not annotated using KEGG blastkoala 
#FUN_009459-T1 
#Magnaporthe_oryzae 
#host: monocot 
#as though the primary source of this PHI host is monocots, however the host range of LO42 target might be different due to the 73.301% identity with 413 alignment length 

#Interestingly, MOHIK5, a Group V HHK known to regulate virulence and stress response in M. oryzae, was uniquely matched in T. cylindrosporum LO42, while its found absence in the other two T.cylindrosporum strains may partially explain LO42’s distinct nematode-lethality phenotype.

#scaffold_750	funannotate	mRNA	22	1240	.	+	.	ID=FUN_009459-T1;Parent=FUN_009459;product=hypothetical protein;Ontology_term=GO:0046777,GO:0007165,GO:0009927,GO:0000155,GO:0005886;Dbxref=InterPro:IPR000014,InterPro:IPR035965,InterPro:IPR005467,InterPro:IPR003661,InterPro:IPR036890,PFAM:PF00512,InterPro:IPR036097,InterPro:IPR000700;note=EggNog:ENOG503NWSP,COG:T;

#PF00512 → Histidine kinase-like ATPase domain
#IPR000014 → His kinase A (phospho-acceptor) domain
#IPR005467 / IPR003661 / IPR036890 / IPR036097 → HAMP and sensory domains
#IPR000700 → Response regulator receiver domain (REC)
#This domain profile matches Group V hybrid histidine kinases, such as MoHik5, known for regulating fungal development, stress responses, and virulence.



#Since the presence of this histidine kinase (HK) homolog is unique to LO42 and absent in the other two T. cylindrosporum strains, it is reasonable to hypothesize that this gene may contribute to broadening the host infection range. While T. cylindrosporum is classically recognized as an insect pathogen, the presence of this group V HK, may enhance LO42’s ability to sense and respond to alternative hosts, such as nematodes. This aligns with preliminary experimental results showing that LO42 exhibits high nematode-killing activity. Notably, although all three T. cylindrosporum strains possess a Sod1 homolog previously associated with nematode pathogenicity, only LO42 harbors this unique histidine kinase, suggesting it may play an additive role in enhancing virulence or environmental adaptability toward non-insect hosts. 

unique(T.cylindrosporum_LO42_cleaned$Gene_name)

###----Investigate core/conservative pathogenic genes in LO42----

LO42_conserved_genes <- T.cylindrosporum_LO42_cleaned$FunannotateIDs[T.cylindrosporum_LO42_cleaned$Protein_ID %in% core_genes_info $Protein_ID] #66 /100 genes of the high confidence pathogenic genes found in LO42 are core and conservative pathogenic gene homologs to all species within genus of Tolypocladium 
writeLines(LO42_conserved_genes, "LO42_high_conf_fungal_geneIDs.txt")
#then go to computecanada to use seqkit grep -f for the corresponding sequence and plug into BLASTKOALA webpage 

#Major group of KEGG category is enviornmental cue signnaling 21
#cellular process 8 
#genetic information processing 7 
#carbohydrate metabolism 6 
#amino acid metabolism 3
#energy metabolism 3 
#The distribution of KEGG pathway categories among your core pathogenic genes suggests that Tolypocladium species rely heavily on environmental sensing and signaling mechanisms, alongside general metabolic and cellular functions, for their pathogenic lifestyle.


#the absence of those key virulence genes used in C.albicans are absent in LO42 blastp search. Only a homologous gene to EFG1 showed 67% identity to C.albicans, this likely suggest that LO42 and C.albicans used different infection mechanisms. Also, other species in Tolypocladium also lacks the same key virulence genes such as HWP1, ALS2, SAP1-7... suggesting this entire genus uses a rather complete different infection mechanisms. This is also supported by other research papers where Tolypoclaidum genus is know for insect and plant pathogen as well as knwon for the production of various secondary metabolites such as cyclosprin A... 

LO42_strainspecific_genes <- T.cylindrosporum_LO42_cleaned %>% 
  filter(!FunannotateIDs %in% LO42_conserved_genes)

writeLines(LO42_strainspecific_genes$FunannotateIDs, "LO42_strain_specific_geneIDs.txt")


core_KOs_df<- read.csv("Core_KOs.txt", sep = "\t", header = FALSE)
core_KOs <- core_KOs_df$V2
core_KOs<- core_KOs[core_KOs != ""]
core_KOs # 62 are KEGG annotated 
writeLines(core_KOs, "core_KOs_for_KEGG.txt")

LO42_strain_KOs_df <- read.csv("LO42_KOs.txt", sep = "\t", header = FALSE)
LO42_strain_KOs<- LO42_strain_KOs_df$V2
LO42_strain_KOs <- LO42_strain_KOs[LO42_strain_KOs != ""]
LO42_strain_KOs #26 are KEGG annotated 
writeLines(LO42_strain_KOs, "LO42_strain_KOs_for_KEGG.txt")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("clusterProfiler")
#install.packages("ggplot2")
library(ggplot2)
library(clusterProfiler)
enrich_result <- enrichKEGG(gene = LO42_strain_KOs, 
           universe = core_KOs,
           organism = "ko",
           pvalueCutoff = 0.05,
           pAdjustMethod = "BH")

enrich_result@result
#no enrichment.... 




###----Compare Between C.albicans & LO42----
process_all_blastp_files(pattern = "_blastpresultv2$")
#Candida albicans does not have any duplicated pathogenic genes, there are 182 pathogenic associated genes found in C.albicans 

setdiff(Candida_albicans_filtered$Gene_name, T.cylindrosporum_LO42_filtered$Gene_name)
#noticed there are 55 unique pathogenic associated genes found in C.albicans but not in LO42, these include Hwp1, Phr1, PLB1, PLB3, TEC1, ALS3, EFG1 

setdiff(T.cylindrosporum_LO42_filtered$Gene_name,Candida_albicans_filtered$Gene_name)
#there are 65 pathogenic genes unique to LO42 
#Conclusion: Even though they both showed the same nematode killing phenotype, but they do not share similar virulence genes, particularly those key virulent genes in C.albican cannot be found in LO42 with significance. 


###----LO42 vs T.tundrense----
setdiff(T.tundrense_CBS56984_cleaned$Protein_ID, T.cylindrosporum_LO42_cleaned$Protein_ID)
#"G4MNA5","Q01143","Q4WI19","F5HCK8"
view(T.tundrense_CBS56984_cleaned %>% 
  filter(Protein_ID %in% c("G4MNA5","Q01143","Q4WI19","F5HCK8")))
#host of these 4 homologs genes are 3 monocots and Q4WI19 1 rodents 

setdiff(T.cylindrosporum_LO42_cleaned$Protein_ID,T.tundrense_CBS56984_cleaned$Protein_ID)
#"E3Q2F8","P00445","A0A7J6IPU7","C5M7A4", G4MKP6" 
view(T.cylindrosporum_LO42_cleaned %>% 
  filter(Protein_ID %in% c("E3Q2F8","P00445","A0A7J6IPU7","C5M7A4", "G4MKP6")))
#2 monocots, 1 eudicots, P00445 1 nematodes, 	C5M7A4 1 rodents 
#the absence of P00445 and C5M7A4 in T.tundrense suggested it has smaller host range compared to LO42

setdiff(T.nubicola_CBS56884_cleaned$Protein_ID, T.cylindrosporum_LO42_cleaned$Protein_ID)
# "G4MPJ4" "G4MWY5" "Q8NK75" "Q4WI19" "I1RST4" "G4NGW8" "N4VSY8" "A9YDN6" "Q51P78"

view(T.nubicola_CBS56884_cleaned %>% 
  filter(Protein_ID %in% c("G4MPJ4","G4MWY5","Q8NK75","Q4WI19","I1RST4","G4NGW8","N4VSY8","A9YDN6","Q51P78") ))
#all are monocots and eudicot, only 1 rodents in host desciption section 

setdiff(T.cylindrosporum_LO42_cleaned$Protein_ID,T.nubicola_CBS56884_cleaned$Protein_ID)
#"G4NG12"     "Q4P9P1"     "G4MUR7"     "G2X7U7"    "A0A0A2JQN4" "E3Q2F8"     "A0A7J6IPU7" "G2X869"  
view(T.cylindrosporum_LO42_cleaned %>% 
  filter(Protein_ID %in% c("G4NG12","Q4P9P1","G4MUR7","G2X7U7","A0A0A2JQN4","E3Q2F8","A0A7J6IPU7","G2X869")))
#half eudicots half monocots 

#there are more differential protein IDs when comparing between T.nubicola and LO42. This is consistent with the phylogenetic tree generated by concatenating 4 core pathogenic genes, in which although LO42 and T.nubicolas and T.tundrense are placed in the same clade suggesting pathogenic functional wise they are closely similar, however T.nubicola is slightly more different compared to the others within this clade. 
#However, the increased number of non-shared pathogenic homologs between T. nubicola and LO42—relative to the LO42 vs T. tundrense comparison—points to a greater degree of functional divergence within the clade

#Interestingly, phylogenetic placement based on concatenated core pathogenic genes placed T. cylindrosporum strains closer to T. nubicola and T. tundrense, while traditional barcoding loci (28S + ITS + TUB1) grouped them with T. paradoxum and T. sp sup PDA1. This discrepancy possibly reflecting selective pressures on pathogenic gene content associated with environmental or host-specific niches.

###-----LO42 Functions 
view(T.cylindrosporum_LO42_cleaned)
T.cylindrosporum_LO42_cleaned$Protein_ID <- trimws(T.cylindrosporum_LO42_cleaned$Protein_ID) 
T.cylindrosporum_LO42_cleaned$Function[T.cylindrosporum_LO42_cleaned$Protein_ID=="G4MTZ0"] <- "ER chaperone" 
#https://www.uniprot.org/uniprotkb/G4MTZ0/entry
T.cylindrosporum_LO42_cleaned$Function[T.cylindrosporum_LO42_cleaned$Protein_ID=="Q9UWF0"] <- "calmodulin-Ca2+ complex" 
T.cylindrosporum_LO42_cleaned$Function[T.cylindrosporum_LO42_cleaned$Protein_ID=="G4MNE6"] <- "regulate cellular stress response" 
#https://www.nature.com/articles/s41598-022-20507-x
T.cylindrosporum_LO42_cleaned$Function[T.cylindrosporum_LO42_cleaned$Protein_ID=="	
A0A0B2XG51"] <- "C2H2 type master regulator of conidiophore development BrlA" 
#https://www.uniprot.org/uniprotkb/A0A0B2XG51/entry
T.cylindrosporum_LO42_cleaned$Function[T.cylindrosporum_LO42_cleaned$Protein_ID=="	
G4N703"] <- "autophagy; lipid degradation" 
#https://www.uniprot.org/uniprotkb/G4N703/entry
T.cylindrosporum_LO42_cleaned$Function[T.cylindrosporum_LO42_cleaned$Protein_ID=="	
E9EMI7"] <- "Conidiophore development regulator abaA" 
#https://www.uniprot.org/uniprotkb/E9EMI7/entry

T.cylindrosporum_LO42_cleaned$Function
shared_genes
core_genes_info$Gene_name
view(T.inflatum_CBS56784_cleaned %>% 
  filter(Host_descripton=="nematodes"))


#没有用的部分
###----Extracting FunannotateIDs for LO42 & C.albicans ----
LO42_high_geneIDs <- unique(T.cylindrosporum_LO42_filtered$FunannotateIDs)
writeLines(LO42_high_geneIDs, "LO42_high_conf_fungal_geneIDs.txt")

process_all_blastp_files(pattern = "_blastpresultv2$")
CA_high_geneIDs <- unique(Candida_albicans_filtered$FunannotateIDs)
writeLines(CA_high_geneIDs, "C.albicans_high_conf_fungal_geneIDs.txt")

rm(Candida_albicans)
rm(Candida_albicans_filtered)





###---Collection Sites----
diff_gene_list <- lapply(pathogenic_gene_list, function(x){
  x[!x %in% shared_genes]
})


#install.packages("readxl")
library(readxl)
collection_site <- as.data.frame(
  read_excel("Collection_site.xlsx", sheet = "Sheet1"))
class(collection_site)

rownames(collection_site) <- collection_site[[1]]

collection_site <- collection_site[-1]
dim(collection_site) #29 rows match 

library(tibble)

transposed_df <- as.data.frame(t(collection_site)) #flip the dataframe 

transposed_df <- transposed_df %>%
  mutate(Latitudes = as.numeric(Latitudes)) %>%
  mutate(Climate_Zone = case_when(
    abs(Latitudes) <= 23.5 ~ "Tropical",
    #(Feeley & Stroud, 2018) 
    abs(Latitudes) > 23.5 & abs(Latitudes) <= 40 ~ "Subtropical", 
    abs(Latitudes) > 40 & abs(Latitudes) < 50 ~ "Temperate", 
    abs(Latitudes) >= 50 & abs(Latitudes) <= 65 ~ "Boreal", 
    abs(Latitudes) > 65 ~ "Polar"
  ))
  #(A Windy.com Company, 2025)

transposed_df$Climate_Zone[rownames(transposed_df)=="T.nubicola_CBS56884"] <- "Alpine" 

view(transposed_df)
dim(transposed_df) #26 samples, 3 samples had no information about collection site including 2 T.endophyticum strains and T.amazonense MS503

   
   
###---- Phylogenetic tree creation----

shared_genes

files

#df_list <- lapply(files, function(x) {
  
#  obj <- get(x)
#  obj %>% 
#    filter(Gene_name %in% shared_genes & Total_length >= 200 & Total_length <= 300)
  
#  print(filtered_genes_by_length[, c("Gene_name", "FunannotateIDs", "Total_length")])
    
#})

df_list <- lapply(files, function(x) {
  
  obj <- get(x) 
  
  filtered <- obj %>%
    filter(Gene_name %in% shared_genes & Total_length >= 200 & Total_length <= 300)
  
  cat("\n---", x, "---\n")
  print(filtered[, c("Gene_name", "FunannotateIDs", "Total_length")])
  
  return(filtered)  # Return the filtered data
})

#notice majority of the samples show MgRho3, Fgvps60_(FGSG_10503),UvBI-1, and Sec4_(FOIG_08151), since they are part of the core pathogenic genes, they are also in length between 200-300bp, not recognizing them from gene duplicates stated earlier. Meeting the requirements talked about by https://doi.org/10.1038/nmicrobiol.2016.48 for generating concatenated genes of new phylogenetic tree. 

target_genes <- c("FGVPS60_(FGSG_10503)", "GCN5","NOX1", "SEC4_(FOIG_08151)")

geneids_list <- list()

for (gene in target_genes) {
  row_ids <- c()
  
  for (file in files) {
    
    df <- get(file)
    id <- df$FunannotateIDs[df$Gene_name == gene]
    row_ids <- c(row_ids, id)
    
  }
  
  geneids_list[[gene]] <- row_ids
  
}

geneids_list #successful 

# Convert to large data.frame
gene_id_df <- as.data.frame(do.call(rbind, geneids_list))

# clean column names 
colnames(gene_id_df) <- gsub("_cleaned$", "", files)
gene_id_df <- cbind(Gene = rownames(gene_id_df), gene_id_df)

colnames(gene_id_df)

gene_id_df <- gene_id_df %>% 
  select(-Gene)#29x4 dimensions, great! 

write.csv(gene_id_df, "gene_funannotate_matrix.csv", row.names = TRUE)

#write.table(gene_id_df,
#            file = "genes_to_extract.txt",
#            row.names = FALSE,
#            col.names = FALSE,
#            quote = FALSE)

#gene_annotateIDS<- read.csv("gene_funannotate_matrix.csv",stringsAsFactors = FALSE )
#rownames(gene_annotateIDS) <- gene_annotateIDS[,1] #succesfully changed names 

#gene_annotateIDS <- gene_annotateIDS %>% 
#  select(-X) #29x4 dim matched as original 
library(Biostrings)
faafiles <- list.files(pattern = "\\.proteins\\.fa$")
faa_list <- lapply(faafiles, readAAStringSet)
names(faa_list) <- gsub("\\.proteins\\.fa$", "", faafiles)
unique(names(faa_list))
all(colnames(gene_id_df) %in% names(faa_list)) #TRUE, file names match 

#names(faa_list$T.album_CBS86973) <- sapply(strsplit(names(faa_list$T.album_CBS86973), " "), `[`, 1)
#faa_list$T.album_CBS86973["FUN_007852-T1"]

faa_list <- lapply(faa_list, function(x) {
  names(x) <- sapply(strsplit(names(x), " "), `[`, 1)
  return(x)
})

faa_list #perfect!

#extracted_sequences <- as.list()

#for (x in rownames(gene_annotateIDS)) {
  
#  for (y in colnames(gene_annotateIDS)){
    #looping through each column then each row, extract the gene_ids wanted 
#    gene_id <- gene_annotateIDS[x,y]
    
#    if(y %in% names(faa_list)){
#      sequences <- faa_list[[y]]
      
#      if (x %in% names(sequences)){
#        extracted_sequences[[paste(x,y, sep = "\t")]] <- sequences[[gene_id]]
#      }
#    }
    
#  }
#}

#extracted_sequences

extracted_sequences <- list()

# Add debug counter
found_sequences <- 0

for (gene_name in rownames(gene_id_df)) {
  for (isolate in colnames(gene_id_df)) {
    
    # Clean the gene ID
    gene_id <- trimws(gene_id_df[gene_name, isolate])
    
    # Skip if empty
    if (is.na(gene_id) || gene_id == "") next
    
    # Debug print
    cat("Searching for:", gene_id, "in", isolate, "\n")
    
    # Check isolate exists
    if (isolate %in% names(faa_list)) {
      
      # Get sequences for this isolate
      isolate_seqs <- faa_list[[isolate]]
      
      # Check if gene exists (exact match)
      if (gene_id %in% names(isolate_seqs)) {
        extracted_sequences[[paste(gene_name, isolate, sep="|")]] <- isolate_seqs[[gene_id]]
        found_sequences <- found_sequences + 1
        cat("FOUND", gene_id, "\n")
      }
    }
  }
}

# Results summary
cat("\nTOTAL SEQUENCES FOUND:", found_sequences, "\n")

extracted_sequences
Fgvps60_sequences<- extracted_sequences[1:29] 
GCN5_sequences <- extracted_sequences[30:58] 
NOX1_sequences <- extracted_sequences[59:87]
Sec4_sequences <- extracted_sequences[88:116]

writeXStringSet(
  AAStringSet(unlist(Fgvps60_sequences)),
  filepath = "Fgvps60_sequences.faa",
  format = "fasta")

writeXStringSet(
  AAStringSet(unlist(GCN5_sequences)),
  filepath = "GCN5_sequences.faa",
  format = "fasta")

writeXStringSet(
  AAStringSet(unlist(NOX1_sequences)),
  filepath = "NOX1_sequences.faa",
  format = "fasta")

writeXStringSet(
  AAStringSet(unlist(Sec4_sequences)),
  filepath = "Sec4_sequences.faa",
  format = "fasta")


library(muscle)
#v3.46.0

Fgvps60_sequences <- AAStringSet(unlist(Fgvps60_sequences))
Fgvps60_sequences
aligned_Fgvps60 <- muscle::muscle(AAStringSet(unlist(Fgvps60_sequences)))
class(Fgvps60_sequences)
aligned_Fgvps60
writeXStringSet(AAStringSet(aligned_mg_rho3),"MgRho3_aligned.fasta")

Fgvps60_aa <- AAStringSet(unlist(Fgvps60_sequences))
Fgvps60_aa
aligned_Fgvps60_aa <- muscle::muscle(AAStringSet(Fgvps60_aa))
class(aligned_Fgvps60_aa)
aligned_Fgvps60_aa
writeXStringSet(AAStringSet(aligned_Fgvps60_aa),"Fgvps60_aligned.fasta")

aligned_UvBI1<- muscle::muscle(AAStringSet(UvBI1_sequences))
aligned_UvBI1
writeXStringSet(AAStringSet(aligned_UvBI1),"UvBI1_aligned.fasta")


aligned_Sec4<- muscle::muscle(AAStringSet(Sec4_sequences))
aligned_Sec4
writeXStringSet(AAStringSet(aligned_Sec4),"Sec4_aligned.fasta")

###went into compute canada to trim off gaps < 0.05%, and manually removed T.endophyticum MX560 because it distorts the tree from previous trial
aligned4 <- readAAStringSet("Sec4_aligned_trimmed.fasta")
aligned2 <- readAAStringSet("NOX1_aligned_aligned.fasta")
aligned3 <- readAAStringSet("GCN5_aligned_trimmed.fasta")
aligned1 <- readAAStringSet("Fgvps60_aligned_trimmed.fasta")

T.cylindrosporum_LO42_filtered$Pathogen_species[T.cylindrosporum_LO42_filtered$Gene_name=="UvBI-1"]

concatenated_seq <- mapply(function(a, b, c,d) {
  AAString(paste0(as.character(a), as.character(b), as.character(c), as.character(d)))
}, aligned1,aligned2, aligned3, aligned4) #concatenating all 4 alignments end by end
writeXStringSet(AAStringSet(concatenated_seq),"concatenated_4core.fasta")

library(ape)
tree<- read.tree("RAxML_bestTree.result")
class(trees) #df 

library(ggtree)
p <- ggtree(tree) +
  geom_tiplab(size = 3.5) +
  theme_tree2()+
#geom_text(aes(label=node), hjust=-0.02) #label clade 52 
#geom_cladelabel(node=52, label="Polar Zone", color="red", offset=0.02)+
 geom_hilight(node=48, fill="gold", extend = 0.035)+
  xlim(0, 0.13)
  
p

#climate_df <- data.frame(
#  label = rownames(transposed_df),                # tree tip labels
#  val = transposed_df$Climate_Zone                # climate zone values
#)
#climate_df$val_numeric <- as.numeric(factor(climate_df$val))

#p2 <- facet_plot(p, panel="Climate Zone", data=climate_df, geom=geom_point, aes(x=val_numeric), color='red3')
#p2 + theme_tree2()

#p2 <- facet_plot(
#  p,panel = "Climate Zone",data = climate_df,geom = geom_point,mapping = aes(x = val_numeric, color = factor(val_numeric))) +
#scale_color_manual(
#    values = c("1" = "red", "2" = "salmon", "3" = "lightblue", "4" = "blue"  ) +
#  theme_tree2()


tree$tip.label #sample names 
rownames(transposed_df)
head(transposed_df)




aligned6 <- readDNAStringSet("ITS2_aligned_trimmed.fasta")
#setdiff(names(aligned5),names(aligned6))
aligned5 <- readDNAStringSet("ITS1_aligned_trimmed.fasta")
#setdiff(names(aligned7), names(aligned5))
#setdiff(names(aligned5),names(aligned7))
aligned7 <- readDNAStringSet("28S_aligned_trimmed.fasta")
#its1_width <- width(aligned6)[1]
#pad_Ns <- AAString(paste(rep("N", its1_width), collapse = "")) #add Ns for 157bp 

#combined <- AAStringSet()

#read in each sequence, if the name is seen in ITS then move on, or else add 157bp in length of Ns
#for (i in names(aligned4)) {
#  gene_seq <- aligned4[[i]]  
#  if (i %in% names(aligned6)) {
#    its_seq <- aligned6[[i]]
#  } else {
#    its_seq <- pad_Ns
#  }
#  combined[[i]] <- AAString(paste0(as.character(gene_seq), as.character(its_seq)))
}

#combined #successful concatenation of ITS + genes sequences 

#writeXStringSet(AAStringSet(combined),"concatenated_18SITS.fasta")



concatenated_ITS28 <- mapply(function(a, b, c) {
  AAString(paste0(as.character(a), as.character(b), as.character(c)))
}, aligned7, aligned6, aligned5)
#concatenating all 4 alignments end by end
library(Biostrings)
writeXStringSet(AAStringSet(concatenated_ITS28), "concatenated_ITS28.fasta")
names(concatenated_ITS28)

library(ggtree)
tree1<- read.tree("RAxML_bestTree.result")
class(tree1) #df 

p1 <- ggtree(tree1) +
  geom_tiplab(size = 3.5) +
  theme_tree2()+
  #geom_text(aes(label=node), hjust=-0.02) #label clade 52 
  #geom_cladelabel(node=52, label="Polar Zone", color="red", offset=0.02)+
  #geom_hilight(node=48, fill="gold", extend = 0.035)+
  xlim(0, 0.40)
p1




###----Heatmap----

all_genes <- unique(unlist(pathogenic_gene_list)) 
all_genes#137 protein query IDs matched in total for all strains in this genus 

gene_matrix <- matrix(0, nrow = length(pathogenic_gene_list), ncol = length(all_genes),dimnames = list(names(pathogenic_gene_list), all_genes))
view(gene_matrix) #all 0 to begin with to allow me for presence and absence count 

#pathogenic_gene_list[["T.album_CBS86973_cleaned"]]

#gene_matrix

for (x in names(pathogenic_gene_list)){
  
  found <- pathogenic_gene_list[[x]]
  
  gene_matrix[x, found] <-1
  
}

gene_matrix <- as.data.frame(gene_matrix)
rownames(gene_matrix)<- gsub("_cleaned", "", rownames(gene_matrix))
class(gene_matrix)

library(pheatmap)

#pheatmap(gene_matrix,
#         cluster_rows = TRUE,
#         cluster_cols = TRUE,
#         color = c("white", "#2874A6"), 
#         show_colnames = TRUE,
#         fontsize_row = 12,
#         fontsize_col = 5
#         )

row_anno <- transposed_df[rownames(gene_matrix), "Climate_Zone", drop = FALSE]

zone_cols <- c(
  Tropical    = "#1a9850",
  Subtropical = "#66bd63",
  Temperate   = "#8073ac",
  Boreal      = "deepskyblue1",
  Polar       = "deepskyblue4",
  Alpine      = "tomato4"
)

pheatmap(gene_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = c("white", "black"),
         show_colnames = TRUE,
         show_rownames = TRUE,
         annotation_row = row_anno,
         annotation_colors = list(Climate_Zone = zone_cols),
         border_color = NA, 
         fontsize_row = 12,
         fontsize_col = 5
         )

core_genes_info$Protein_ID
#there are in total of 137 protein IDs that are matched to the species, while there is an obvious pattern of core presence of certain protein IDs conserved in all strains.

#Also, base on the heatmap we notice that T.capitatum CBS113982, all strains of T.inflatum, all strains of T.cylindrosporum, and T.pustulatum, T.sp_salcha_MEA2 are cluster together. Indicating that they share more similarity in the presence/absence in the overall pathogenicty genes detected. Which means they may have used similar infection mechanisms or simialrity virulence strategies in host infection (more similar pathogenic potentials). They shared a common enviornmental niche. 
#while phylogenetic tree did not cluster these together, likely suggest species under go convergent selection pressure such as similar enviornmental and/or host selection pressure.
#Both T.cylindrosporum and T.inflatum are known as insect pathogens. 
#Notice T.tundrense, T.ophioglossoides, and T.geodes are placed in the sister clade to the T.inflatum,T.cylindrosporum, T.sp_salcha, T.capitatum clade.   


boreal_species <- rownames(transposed_df)[transposed_df$Climate_Zone == "Boreal", ]
#16 species from boreal zone 

boreal_matrix <- gene_matrix[boreal_species, ]
view(boreal_matrix)

# Step 3: Find columns (PAGs) present (value = 1) in ALL boreal species
conserved_boreal_PAGs <- colnames(boreal_matrix)[colSums(boreal_matrix) == nrow(boreal_matrix)]

# Output
conserved_boreal_PAGs

#可能无用
###---- Comparison Between LO42 and T.tundrense, T.nubicola----
#chose T.ophiogloissoides and T.pustulatum as they are both collected from boreal region however with different known host 


genes_ophio <- unique(T.ophioglossoides_CBS100239_cleaned$Protein_ID)
genes_pustulatum <- unique(T.pustulatum_CBS110433_cleaned$Protein_ID)
genes_lo42 <- unique(T.cylindrosporum_LO42_cleaned$Protein_ID)

# Plot Venn diagram
venn.plot1 <- venn.diagram(
  x = list(LO42 = genes_lo42, T.pustulatum_CBS110433 = genes_pustulatum, T.ophioglossoides_CBS100239 = genes_ophio),
  filename = NULL, #to plot 
  fill = c("lightgreen", "orange", "yellow"),
  alpha = 0.5,
  cat.cex = 1.2,
  margin = 0.1,
)
grid.draw(venn.plot1)
grid.text("Overall PAGs Profile Comparison Between LO42 & 2 Distinct Host Strains", y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))

setdiff(
  unique(T.cylindrosporum_LO42_cleaned$Protein_ID),
  union(
    unique(T.ophioglossoides_CBS100239_cleaned$Protein_ID),
    unique(T.pustulatum_CBS110433_cleaned$Protein_ID)
  )
)

T.cylindrosporum_LO42_cleaned$FunannotateIDs[
  T.cylindrosporum_LO42_cleaned$Protein_ID %in% c("E3Q2F8", "G4MKP6")
]
#"FUN_009183-T1" "FUN_009459-T1"

setdiff(
  unique(T.ophioglossoides_CBS100239_cleaned$Protein_ID),
  union(
    unique(T.cylindrosporum_LO42_cleaned$Protein_ID),
    unique(T.pustulatum_CBS110433_cleaned$Protein_ID)
  )
)
T.ophioglossoides_CBS100239_cleaned$FunannotateIDs[T.ophioglossoides_CBS100239_cleaned$Protein_ID %in% c("H1W5W6", "O42773","Q5EKM7","B9DR51","Q4WI19","G4N5S6")]
# "FUN_000621-T1" "FUN_001510-T1" "FUN_002051-T1" "FUN_002815-T1" "FUN_003968-T1" "FUN_006797-T1"
T.ophioglossoides_CBS100239_cleaned$Function[T.ophioglossoides_CBS100239_cleaned$Protein_ID %in% c("H1W5W6", "O42773","Q5EKM7","B9DR51","Q4WI19","G4N5S6")]                                                          

###----T.cylindrosporum-Genome Statistics----
#install.packages("gt")
library(gt)
table <- data.frame(
  N50 = c(229237, 30235, 41578),
  L50 = c(37, 323, 225),
  Total_Length = c(31220495, 32018899, 31356200),
  GC_Content = c(57.79, 57.89, 58.01)
)
rownames(table) <- c("T.cylindrosporum_LO42", "T.cylindrosporum_IBT41711","T.cylindrosporum_CBS718.70")

table$Strain <- rownames(table)
table <- table[, c("Strain", "N50", "L50", "Total_Length", "GC_Content")]

table %>%
  gt() %>%
  tab_header(
    title = "T.cylindrosporum Genome Statistics"
  ) %>%
  fmt_number(
    columns = c(N50, L50, Total_Length),
    use_seps = TRUE,
    decimals = 0
  ) %>%
  fmt_number(
    columns = GC_Content,
    decimals = 2
  ) %>%
  cols_label(
    Strain = "Strains",
    N50 = "N50 (bp)",
    L50 = "L50 (contigs)",
    Total_Length = "Total Length (bp)",
    GC_Content = "GC Content (%)"
  ) %>%
  tab_options(
    table.font.size = "small"
  ) %>% 
  tab_style(      #color the title 
    style = list(
      cell_text(color = "white", weight = "bold"),
      cell_fill(color = "#1B4F72")
    ),
    locations = cells_title(groups = "title")
  ) %>% 
  tab_style(
    style = list(
      cell_text(color = "white", weight = "bold"),
      cell_fill(color = "#2874A6")
    ),
    locations = cells_column_labels(everything())
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "#ECF0F1"),
      cell_text(weight = "bold")
    ),
    locations = cells_body(columns = "Strain")
  )


###----Checking For Up or Down Pathogenicity----###
table(T.cylindrosporum_LO42_cleaned$Phenotype[T.cylindrosporum_LO42_cleaned$Protein_ID%in% core_genes_info$Protein_ID])
#where they have investigated both lethal and loss of pathogenicity phenotypes in their experiment 
#under one condition the mutation of this gene was lethal, it was required for living; under another condition, it caused loss of pathogenicity
#the dash line -- means it has been tested and seen differences in virulence under different conditions or hosts 
#3/66 core pathogenic genes are with lethal__loss of pathogenicity
#54/66 are with loss of pathogenicity phenotype 
#5/66 are loss of pathogenicity__reduced virulence
#1/66 loss of pathogenicity--unaffected pathogenicity 
#3/66 reduced virulence__unaffected__loss of pathogenicity

table(T.cylindrosporum_LO42_cleaned$Host_descripton[T.cylindrosporum_LO42_cleaned$Phenotype =="loss_of_pathogenicity"])

loss_LO42 <- T.cylindrosporum_LO42_cleaned %>% 
  filter(Protein_ID %in% shared_genes) %>% 
  filter(Phenotype == "loss_of_pathogenicity")
table(loss_LO42$Host_descripton)

#install.packages("webr")
library(webr)

all_LO42_core_genes <- T.cylindrosporum_LO42_cleaned %>% 
  filter(Protein_ID %in% shared_genes)

unique(all_LO42_core_genes$Phenotype) #now only 5 types of phenotype 
all_LO42_core_genes$Phenotype[all_LO42_core_genes$Phenotype == "reduced_virulence__loss_of_pathogenicity"] <- "loss_of_pathogenicity__reduced_virulence"  
all_LO42_core_genes$Phenotype[all_LO42_core_genes$Phenotype == "unaffected_pathogenicity__reduced_virulence__loss_of_pathogenicity"] <- "reduced_virulence__unaffected_pathogenicity__loss_of_pathogenicity"  
table(all_LO42_core_genes$Phenotype)


all_LO42_core_genes$Host_descripton[all_LO42_core_genes$Host_descripton == "monocot"] <- "monocots"
all_LO42_core_genes$Host_descripton[all_LO42_core_genes$Host_descripton == "eudicot"] <- "eudicots"
#unique(all_LO42_core_genes$Host_descripton)

PD <- all_LO42_core_genes %>% 
  group_by(Phenotype, Host_descripton) %>% 
  summarise(n = n(), .groups = "drop")

PieDonut(PD, aes(Phenotype, Host_descripton, count = n ), 
         title = "Conserved Pathogenic Homologs with Host Profiles", showDonutName = FALSE, showPieName = FALSE, explode = 5, explodeDonut=TRUE)

PieDonut(PD, aes(Phenotype, Host_descripton, count = n), 
         title = "Conserved Pathogenic Homologs with Host Profiles") +
  theme(
    axis.title = element_blank(),  # Remove axis titles
    axis.text = element_blank(),    # Remove axis text
    axis.ticks = element_blank()  # Remove axis tick      # Remove legend
  )

T.cylindrosporum_LO42_filtered$Phenotype[T.cylindrosporum_LO42_filtered$Gene_name %in% shared_genes]=="unaffected_pathogenicity__reduced_virulence__loss_of_pathogenicity"

T.cylindrosporum_LO42_filtered$Gene_name[T.cylindrosporum_LO42_filtered$Phenotype=="unaffected_pathogenicity__reduced_virulence__loss_of_pathogenicity"]
#pyruvate_kinase

T.cylindrosporum_LO42_filtered %>% 
  filter(Phenotype=="unaffected_pathogenicity__reduced_virulence__loss_of_pathogenicity")

data1$PHI_MolConn_ID[data1$ProteinID=="G4MXS1"] #PHI:8611, used to check on PHI website for papers; since the main paper suggested knockout of this gene results in loss in pathogenicity thus, the complex phenotype of this gene can be reassigned


identity_values1 <- lapply(files, function(x) {
  
  cat("\n---", x, "---\n")
  
  obj <- get(x) 
  
  filtered <- obj %>%
    filter(Query_ID == "G4MXS1") %>% 
    pull(Identity)
  
  print(filtered)
    # Return the filtered data
})

identity_values_flat1 <- unlist(identity_values1)
mean(identity_values_flat1)
#mean identity of this homologous genes found in this genus is 82.622%


T.cylindrosporum_LO42_filtered$Gene_name[T.cylindrosporum_LO42_filtered$Phenotype=="lethal__loss_of_pathogenicity"]
#TUB1 90.9% identity, GFA1 84.5% identity 

T.cylindrosporum_LO42_filtered %>% 
  filter(Phenotype=="lethal__loss_of_pathogenicity")
#Query_ID:Q4WKG5, E9R5E2

data1$PHI_MolConn_ID[data1$ProteinID=="E9R5E2"] #"PHI:2513"

data1$PHI_MolConn_ID[data1$ProteinID== "Q4WKG5"] #"PHI:2530"
#TUB1 not directly related to pathogenicity, however, the depress expression of this gene lead to growth defects and loss in pathogenicity

T.cylindrosporum_LO42_filtered %>% 
  filter(Phenotype=="loss_of_pathogenicity__reduced_virulence")
#VMK1, Erv14_(MGG_08132), GPI, Gcn5 
#Q5EKM7, G4MYQ8, A0A2H3GW15, A0A2H3G405
#96.3, 89.8, 89.7, 89.7

data1$PHI_MolConn_ID[data1$ProteinID== "Q5EKM7"] #PHI:483
#lost of pathogenicity in some host while reduce pathogenicity in some other hosts
data1$PHI_MolConn_ID[data1$ProteinID== "G4MYQ8"] #"PHI:123965"
#involved in intracellular transport of few receptor proteins, link to pathogenicity indirectly 

data1$PHI_MolConn_ID[data1$ProteinID== "A0A2H3G405"] #PHI:124279 
#Gcn5 found to be involved in vegetaion growth 

data1$PHI_MolConn_ID[data1$ProteinID== "A0A2H3GW15"] #PHI:124811


###----Boreal Zone species PAG Analysis----
#extract only boreal species 
rownames(transposed_df[transposed_df$Climate_Zone == "Boreal", ])
boreal_species_df <- transposed_df %>% 
  filter(!is.na(Climate_Zone)) %>% 
  filter(Climate_Zone == "Boreal") 
boreal_species<- rownames(boreal_species_df)
boreal_species


pathogenic_genename_list <- lapply(files, function(x) {
  
  obj <- get(x)
  obj[[2]]
})
pathogenic_genename_list
names(pathogenic_genename_list) <- gsub("_cleaned$", "", basename(files))
names(pathogenic_genename_list)


conserved_boreal_list <- pathogenic_gene_list[names(pathogenic_gene_list) %in% boreal_species]
unique(conserved_boreal_list)
#there are 100 shared PAGs among these boreal zone species 
boreal_list<- setdiff(conserved_boreal_list[[1]], shared_genes)
#when compared with 66 conserved PAGs, the 35 remaining PAGs are boreal specific 
boreal_list

borealPAG_functions<- unique(data1$Function[data1$Protein_ID %in% boreal_list])

data1$Gene_name[data1$Function=="guards cells against oxidative damage by neutralizing reactive oxygen species"]
#SOD1 
borealPAGs_df <- data1 %>% 
  filter(Function %in% borealPAG_functions) %>% 
  filter(Protein_ID %in% boreal_list) %>% 
  select(Protein_ID,Gene_name, Function)

view(borealPAGs_df)

LO42_specific_genes_df<- T.cylindrosporum_LO42_cleaned %>% 
  filter(!Protein_ID %in% shared_genes)

LO42_specific_genes_df$FunannotateIDs
