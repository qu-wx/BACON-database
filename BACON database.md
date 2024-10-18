
# BACON database 
---------
# 1.[Intra-species communication](https://github.com/qu-wx/BACON-database/tree/main/Intra-species)
BACON database for the intra-species communication is curated from existing databases and literature mining, and contains 43 communication entries of 22 species. Users can easily download the corresponding database for the species of interest.

# 2.Cross-species communication
For the synthesis group, it contains 899 communication entries derived from 663 species. For the reception group, it consists of 701 communication entries from 538 species.Users can select their species of interest based on the following pipeline to construct the dataset required for BACON, then can easily input this data into BACON for use.
## 1.Select species from synthesis group
For example, we aim to study the communication among Bacteroides dorei and Acinetobacter baumannii, Aeromonas caviae, Anaerosalibacter massiliensis, and Anaerostipes rhamnosivorans. 

```r
synthesis_group <- read.csv('microbiome_db_synthesis_group.csv',header = T)
species_i_choose <- c("Bacteroides dorei",'Acinetobacter baumannii','Aeromonas caviae','Anaerosalibacter massiliensis','Anaerostipes rhamnosivorans') 
```

## 2.Select species from reception group

```r
reception_group <- read.csv('microbiome_db_reception_group.csv',header = T)
```


## 3.Combine the database from synthesis group and reception group
```c
process_dbs <- function(synthesis_group, reception_group) {
  
  synthesis_i_choose_db <- synthesis_group[c(which(synthesis_group$ligand_species%in%c(species_i_choose))),]
  reception_i_choose_db <- reception_group[c(which(reception_group$target_species%in%c(species_i_choose))),]
  
  as.data.frame(table(synthesis_i_choose_db$ligand_system)) -> synthesis_system_list
  as.data.frame(table(reception_i_choose_db$target_system)) -> reception_system_list
  
  intersect(synthesis_system_list$Var1,reception_system_list$Var1) -> work_system
  synthesis_i_choose_db <- synthesis_i_choose_db[c(which(synthesis_i_choose_db$ligand_system%in%c(work_system))),]
  reception_i_choose_db <- reception_i_choose_db[c(which(reception_i_choose_db$target_system%in%c(work_system))),]
  
  unique_systems <- unique(c(synthesis_i_choose_db$ligand_system, reception_i_choose_db$target_system))
  
  system_lists <- list()
  
  for (system in unique_systems) {
   
    synthesis_count <- sum(synthesis_i_choose_db$ligand_system == system)
    reception_count <- sum(reception_i_choose_db$target_system == system)
    
    if (synthesis_count == reception_count && synthesis_count > 0) {
     
      system_info1 <- synthesis_i_choose_db[c(which(synthesis_i_choose_db$ligand_system==system)),]
      system_info2 <- reception_i_choose_db[c(which(reception_i_choose_db$target_system==system)),]
      combined_info <- cbind(system_info1, system_info2)
      system_lists[[paste0("system_list_", system)]] <- combined_info
    } else if (synthesis_count > 0 && reception_count > 0) {
      
      system_info1 <- synthesis_i_choose_db[c(which(synthesis_i_choose_db$ligand_system==system)),]
      system_info2 <- reception_i_choose_db[c(which(reception_i_choose_db$target_system==system)),]
      A <- nrow(system_info1)
      B <- nrow(system_info2)
      C <- ncol(system_info1)
      combined_result <- data.frame(matrix(ncol = C + C, nrow = A * B))
      colnames(combined_result) <- c(colnames(system_info1), colnames(system_info2))
      row_index <- 1
      for (i in 1:A) {
        for (j in 1:B) {
          combined_result[row_index, ] <- c(system_info1[i, ], system_info2[j, ])
          row_index <- row_index + 1
        }
      }
      system_lists[[paste0("system_list_", system)]] <- combined_result
    }
  }
  system_db_total <- c()
  for (i in 1:length(system_lists)){
    as.data.frame(system_lists[[i]]) -> system_db_1
    system_db_total <- rbind(system_db_total,system_db_1)
  }
  return(system_db_total)
}

```
`system_db_total <- process_dbs(synthesis_group,reception_group)` #get a data frame containing species, genes related to QS communication information

```r
BACON_db_construction <- function(system_db_total){
  if(length(grep(",",system_db_total$ligand_gene))!=0){system_db_total$ligand_gene <- strsplit(system_db_total$ligand_gene,split = ",")}
  if(length(grep(",",system_db_total$ligand_complex))!=0){system_db_total$ligand_complex <- strsplit(system_db_total$ligand_complex,split = ",")}
  for (i in 1:nrow(system_db_total)){system_db_total$ligand_complex[[i]] <- as.double(system_db_total$ligand_complex[[i]])}
  if(length(grep(",",system_db_total$ligand_complex_number))!=0){system_db_total$ligand_complex_number <- strsplit(system_db_total$ligand_complex_number,split = ",")}
  for (i in 1:nrow(system_db_total)){system_db_total$ligand_complex_number[[i]] <- as.double(system_db_total$ligand_complex_number[[i]])}
  
  if(length(grep(",",system_db_total$target_gene))!=0){system_db_total$target_gene <- strsplit(system_db_total$target_gene,split = ",")}
  if(length(grep(",",system_db_total$target_complex))!=0){system_db_total$target_complex <- strsplit(system_db_total$target_complex,split = ",")}
  
  for (i in 1:nrow(system_db_total)){system_db_total$target_complex[[i]] <- as.double(system_db_total$target_complex[[i]])}
  if(length(grep(",",system_db_total$target_complex_number))!=0){system_db_total$target_complex_number <- strsplit(system_db_total$target_complex_number,split = ",")}
  for (i in 1:nrow(system_db_total)){system_db_total$target_complex_number[[i]] <- as.double(system_db_total$target_complex_number[[i]])}
  system_db_total$interaction_name <- paste(system_db_total$ligand_species,system_db_total$target_species,system_db_total$target_system,sep = "_")
  #single_db_ecoli$interaction_name <- single_db_ecoli$language
  
  
  micro_db_i_select <- list()
  for (i in 1:nrow(system_db_total)){
    list_db_i <- list(#organism = single_db_ecoli$ligand_sp[i],
      
      #To_organism = single_db_ecoli$targer_sp[i],
      
      interaction_name = system_db_total$interaction_name[i],
      
      ligand_gene = system_db_total$ligand_gene[i][[1]],
      
      target_gene = system_db_total$target_gene[i][[1]],
      
      ligand_complex = system_db_total$ligand_complex[i][[1]],
      
      ligand_complex_number = system_db_total$ligand_complex_number[i][[1]],
      
      target_complex = system_db_total$target_complex[i][[1]],
      
      target_complex_number = system_db_total$target_complex_number[i][[1]],
      
      QS_language = system_db_total$target_system[i]
      
    )
    micro_db_i_select[[system_db_total$interaction_name[i]]] <- list_db_i
  } 
  
  return(micro_db_i_select)
}

```
`
micro_db_i_select <- BACON_db_construction(system_db_total)` ## get micro_db_i_select used for running BACON


---


