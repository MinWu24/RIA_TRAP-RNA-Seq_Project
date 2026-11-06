library("dplyr")
library("data.table")
library("varhandle")
library("stringr")
library("stringi")
library("readr")
##CHANGE -> search for this to find the sections where the file names need to be changed
##PARAMETER-> search for this to find all the specific parameters set

#Reading in the pairwise sample comparisons from the Majiq analysis

sample1_vs_sample2<-read.table(file='/op-pa.deltapsi.tsv',sep='\t',na.strings=c("","","NA"),header=TRUE,stringsAsFactors=FALSE)
sample1_vs_sample2<-read.table(file='/op-gacA.deltapsi.tsv',sep='\t',na.strings=c("","","NA"),header=TRUE,stringsAsFactors=FALSE)
sample1_vs_sample2<-read.table(file='/op_pak.deltapsi.tsv',sep='\t',na.strings=c("","","NA"),header=TRUE,stringsAsFactors=FALSE)
sample1_vs_sample2<-read.table(file='PA-gacA.deltapsi.tsv',sep='\t',na.strings=c("","","NA"),header=TRUE,stringsAsFactors=FALSE)
sample1_vs_sample2<-read.table(file='PA-PAK.deltapsi.tsv',sep='\t',na.strings=c("","","NA"),header=TRUE,stringsAsFactors=FALSE)
sample1_vs_sample2<-read.table(file='gacA-PAK.deltapsi.tsv',sep='\t',na.strings=c("","","NA"),header=TRUE,stringsAsFactors=FALSE)

#CHANGE
sample_name<- "gacA_pak"
sample_name1<- str_extract(sample_name, "^(.+)(?=_)")
sample_name2<- str_extract(sample_name, "(?<=_)(.+$)")

colnames(sample1_vs_sample2)<- c("Gene.Name", "Gene.ID", "LSV.ID", "mean_dPSI", "prob_changing", "prob_notchanging", "PSI1", "PSI2", "LSV.type", "Num.junctions", "Num.exons", "denovo", "Chromosome", "strand", "Junc.coords", "Exon.coords", "IR.coords", "UCSC.LSV.Link")


#function taking the sample1_vs_sample2 Majiq dPSI file and removing low confidence entries 
split_values <- function(two_sample){
  #take all the relevant columns with ; split values and make a matrix data frame where each row has the
  #corresponding values from different columns. Not doing enything with the exon junctions column yet because
  #those values don't directly correspond to the other ; split values.
  mean_dPSI <- as.numeric(unlist(strsplit(as.character(two_sample["mean_dPSI"]), ";")))
  prob_changing <- as.numeric(unlist(strsplit(as.character(two_sample["prob_changing"]), ";")))
  prob_notchanging <- as.numeric(unlist(strsplit(as.character(two_sample["prob_notchanging"]), ";")))
  sample_onePSI <- as.numeric(unlist(strsplit(as.character(two_sample["PSI1"]), ";")))
  sample_twoPSI <- as.numeric(unlist(strsplit(as.character(two_sample["PSI2"]), ";")))
  denovo <- as.numeric(unlist(strsplit(as.character(two_sample["denovo"]), ";")))
  junc_coords <- strsplit(as.character(two_sample["Junc.coords"]), ";")
  #save the Num.junctions and the Num.exons
  Num.exons <- as.numeric(two_sample["Num.exons"])
  Num.junctions <- as.numeric(two_sample["Num.junctions"])
  IR_coords <- two_sample["IR.coords"]
  
  #have to have modified IR coords that are the last exon coords to add to list of kept junction coords
  #if the meanPSI is high enough because the IR coords don't match the exon coords
  if(!is.na(IR_coords)){
    IR_coords_add<- unlist(strsplit(as.character(two_sample["IR.coords"]), ";"))
    IR_coords_add <-read.table(text = IR_coords_add, sep = "-", na.strings = c("na"), colClasses = "numeric")
    IR_coords_add$V1<- IR_coords_add$V1 - 1
    IR_coords_add$V2<- IR_coords_add$V2 + 1
    IR_coords_add<- with(IR_coords_add, paste(V1, V2, sep="-"))
    IR_coords_add<- as.data.frame(as.numeric(unlist(strsplit(as.character(IR_coords_add), "-"))))
    colnames(IR_coords_add) <- "coords_keep"
  }
  #data frame with all split values
  newdf <- data.frame(mean_dPSI, prob_changing, prob_notchanging, sample_onePSI, sample_twoPSI, denovo, junc_coords, Num.exons, Num.junctions, IR_coords)
  colnames(newdf) <- c("mean_dPSI", "prob_changing", "prob_notchanging", "sample_onePSI", "sample_twoPSI", "denovo", "junc_coords", "Num.exons", "Num.junctions", "IR_coords")
  
  ##PARAMETER
  #select for all the entries where one or both samples have a PSI of >=0.05
  df_select <- subset(newdf, ((sample_onePSI>=0.05 | sample_twoPSI>=0.05) & (sample_onePSI<=0.95 | sample_twoPSI<=0.95)))
  if (nrow(df_select)<1){
    return()
  }
  
  #save the new number of junctions
  Num.junctions <- length(df_select$junc_coords)
  
  ##Need the junction coordinates of the eliminated junctions to remove the irrelevant exon junctions.
  ##Find the eliminated coordinates -> find which coordinates out of these are not present in the junction coords
  ##that were kept(also want to keep those exons) -> find the exons with that junction and eliminate
  #recall the entries that were eliminated by the cutoff
  df_eliminate <- subset(newdf, ((sample_onePSI<0.05 & sample_twoPSI<0.05) | (sample_onePSI>0.95 & sample_twoPSI>0.95)))
  
  
  #get the junction coordinates of the removed junctions
  if (nrow(df_eliminate)>0){
    coords_eliminate <- as.data.frame(as.numeric(unlist(strsplit(as.character(df_eliminate$junc_coords), "-"))))
    colnames(coords_eliminate)<- "coords_eliminate"
    
    coords_keep <- as.data.frame(as.numeric(unlist(strsplit(as.character(df_select$junc_coords), "-"))))
    colnames(coords_keep) <- "coords_keep"
    
    #In case the kept coords correspond to the IR coords
    if(!is.na(IR_coords)){
      
      IR_coords_split <- unique(as.data.frame(as.numeric(unlist(strsplit(as.character(df_select$IR_coords), "-")))))
      colnames(IR_coords_split) <- "IR_coords_split"
      
      #Check if the coords_keep contain the original IR_coords; if they do, add the modified IR_coords to the list
      # of coords_keep so they will match to the right exon_coord
      is_IR_coords<- subset(coords_keep, (coords_keep %in% IR_coords_split$IR_coords_split))
      
      if (nrow(is_IR_coords)>0){
        coords_keep<- rbind(coords_keep, IR_coords_add)
      }
      else{
        coords_keep <- as.data.frame(as.numeric(unlist(strsplit(as.character(df_select$junc_coords), "-"))))
        colnames(coords_keep) <- "coords_keep"
        IR_coords<- NA
      }
    }
    else{
      coords_keep <- as.data.frame(as.numeric(unlist(strsplit(as.character(df_select$junc_coords), "-"))))
      colnames(coords_keep) <- "coords_keep"
    }
    
    if(nrow(coords_keep)>0){colnames(coords_keep) <- "coords_keep"
    
    #get the junction coordinates to be eliminated that are not kept
    coords_not_in <- subset(coords_eliminate, !(coords_eliminate %in% coords_keep$coords_keep))
    coords_not_in<-as.data.frame(unlist(apply(coords_not_in, 2, function(x)as.numeric(x))))
    
    if (nrow(coords_not_in)>0){
      colnames(coords_not_in) <- "coords_eliminate"
      #split the exon coordinates
      exon_coords <- unlist(strsplit(as.character(two_sample["Exon.coords"]), ";"))
      #need to change the 'na' values to NA to get numeric values for the coordinates
      exon_coords <-read.table(text = exon_coords, sep = "-", na.strings = c("na"), colClasses = "numeric")
      exon_coords_keep<-apply(coords_keep, 1, function(x)subset(exon_coords, between(x, exon_coords[1], exon_coords[2], NAbounds = FALSE)))
      
      exon_coords_keep<-unique(bind_rows(exon_coords_keep))
      
      #need to anti-join the exon coords table with the eliminated junction coords table -> sequentially change
      #the junction start and junction end column names to eliminate rows of exon coordinates that match 
      #coords_not_in. The values in the table after an anti-join are only the values in exon_coords that are
      #not in coords_not_in. Next step is to check if these exon coords are among the coordinates that match the
      #kept junction coordinates
      colnames(exon_coords)[1]<-"coords_eliminate"
      exon_coords <- anti_join(exon_coords, coords_not_in)
      
      
      colnames(exon_coords)<-c("V1", "V2")
      exon_coords<- unique(rbind(exon_coords, exon_coords_keep))
      # change the second junction column to the appropriate column name to match coords_not_in
      if(nrow(exon_coords)>0){colnames(exon_coords)<-c("V1", "coords_eliminate")
      exon_coords <- anti_join(exon_coords, coords_not_in)
      colnames(exon_coords)<-c("V1", "V2")
      exon_coords<- unique(rbind(exon_coords, exon_coords_keep))
      if (nrow(exon_coords)>0){colnames(exon_coords)<-c("V1", "V2")
      Num.exons <- nrow(exon_coords)
      exon_coords<- paste0(with(exon_coords, paste(V1, V2, sep="-")), collapse = ";")
      }
      else {exon_coords<- unique(rbind(exon_coords, exon_coords_keep))
      Num.exons <-nrow(exon_coords)}
      }
      else{exon_coords<- unique(rbind(exon_coords, exon_coords_keep))
      Num.exons <-nrow(exon_coords)}
    }
    else{
      Num.junctions <- length(df_select$junc_coords)
      exon_coords <- two_sample["Exon.coords"]}
    }
    else{exon_coords <- two_sample["Exon.coords"]}
  }
  else{exon_coords <- two_sample["Exon.coords"]}
  
  #Recovering each LSV line in its original format but with updated junction and exon numbers based on the 
  #removal of low confidence junctions
  mean_dPSI <- paste0(df_select$mean_dPSI, collapse=";")
  prob_changing <- paste0(df_select$prob_changing, collapse=";")
  prob_notchanging <- paste0(df_select$prob_notchanging, collapse=";")
  sample_onePSI <- paste0(df_select$sample_onePSI, collapse=";")
  sample_twoPSI <- paste0(df_select$sample_twoPSI, collapse=";")
  denovo <- paste0(df_select$denovo, collapse=";")
  junc_coords <- paste0(df_select$junc_coords, collapse=";")
  df_select2<- data.frame(mean_dPSI, prob_changing, prob_notchanging, sample_onePSI, sample_twoPSI, denovo, junc_coords, exon_coords, Num.junctions, Num.exons)
  df_select2$Gene.Name <-two_sample["Gene.Name"]
  df_select2$Gene.ID <-two_sample["Gene.ID"]
  df_select2$LSV.ID <-two_sample["LSV.ID"]
  df_select2$Chromosome <-two_sample["Chromosome"]
  df_select2$strand <-two_sample["strand"]
  df_select2$IR_coords <- IR_coords
  rownames(df_select2)<- 1:nrow(df_select2)
  return(df_select2)
}

sample1_vs_sample2_split <- apply(sample1_vs_sample2, 1, function(X)split_values(X))

sample1_vs_sample2_split <-bind_rows(sample1_vs_sample2_split)

sample1_vs_sample2_split<- sample1_vs_sample2_split[c("Gene.Name", "Gene.ID", "LSV.ID", "Chromosome", "strand", "sample_onePSI", "sample_twoPSI", "mean_dPSI", "prob_changing", "prob_notchanging", "Num.junctions", "Num.exons", "denovo", "junc_coords", "exon_coords", "IR_coords")]

#replace blanks with NA and then remove LSV lines that no longer have at least 1 junction
sample1_vs_sample2_split<- apply(sample1_vs_sample2_split, 2, function(x) gsub("^$|^ $", NA, x))
sample1_vs_sample2_split<- data.frame(sample1_vs_sample2_split)
sample1_vs_sample2_split <- subset(sample1_vs_sample2_split, !is.na(sample1_vs_sample2_split$junc_coords))

## Combine the delta signs with the full table for each junction from the Majiq and from the SJ
split_majiq<- function(majiq_table){
  junction_coords <- strsplit(as.character(majiq_table["junc_coords"]), ";")
  sample_onePSI <- strsplit(as.character(majiq_table["sample_onePSI"]), ";")
  sample_twoPSI <- strsplit(as.character(majiq_table["sample_twoPSI"]), ";")
  mean_dPSI <- strsplit(as.character(majiq_table["mean_dPSI"]), ";")
  prob_changing <- strsplit(as.character(majiq_table["prob_changing"]), ";")
  prob_notchanging <- strsplit(as.character(majiq_table["prob_notchanging"]), ";")
  denovo <- strsplit(as.character(majiq_table["denovo"]), ";")
  coord_df <-data.frame(junction_coords, sample_onePSI, sample_twoPSI, mean_dPSI, prob_changing, prob_notchanging, denovo)
  colnames(coord_df) <-c("Junctions.coords", "sample_onePSI", "sample_twoPSI", "mean_dPSI", "prob_changing", "prob_notchanging", "denovo")
  coord_df <- cbind(coord_df, majiq_table["Gene.Name"], majiq_table["Gene.ID"], majiq_table["LSV.ID"], majiq_table["Chromosome"], majiq_table["strand"], majiq_table["Num.junctions"], majiq_table["Num.exons"], majiq_table["exon_coords"], majiq_table["IR_coords"])
  colnames(coord_df) <- c("Junctions.coords", "sample_onePSI", "sample_twoPSI", "mean_dPSI", "prob_changing", "prob_notchanging", "denovo", "Gene.Name", "Gene.ID", "LSV.ID", "Chromosome", "strand", "Num.junctions", "Num.exons", "exon_coords", "IR_coords")
  return(coord_df)
}

Majiq_expanded <- apply(sample1_vs_sample2_split, 1, function(x)split_majiq(x))
Majiq_expanded <-bind_rows(Majiq_expanded)

Majiq_expanded<- Majiq_expanded[c("Gene.Name", "Gene.ID", "LSV.ID", "Chromosome", "strand", "sample_onePSI", "sample_twoPSI", "mean_dPSI", "prob_changing", "prob_notchanging", "Num.junctions", "Num.exons", "denovo", "Junctions.coords", "exon_coords", "IR_coords")]

#PARAMETER
#assign a delta sign to the junctions to filter for low confidence changes from the SJ files
get_PSI <- function(expanded_table){
  #print(expanded_table)
  if (abs(as.numeric(expanded_table["mean_dPSI"]))<=0.15 & (as.numeric(expanded_table["prob_changing"]))<=0.6){
    delta_sign <- "-"
  }
  
  else if (abs(as.numeric(expanded_table["mean_dPSI"]))>0.15 & abs(as.numeric(expanded_table["mean_dPSI"]))<0.20 & (as.numeric(expanded_table["prob_changing"]))>0.6){
    delta_sign <- "~"
  }
  else if (abs(as.numeric(expanded_table["mean_dPSI"]))>=0.15 & abs(as.numeric(expanded_table["mean_dPSI"]))<=0.20 & (as.numeric(expanded_table["prob_changing"]))<=0.6){
    delta_sign <- "-"
  }
  else if (abs(as.numeric(expanded_table["mean_dPSI"]))>0.20 & (as.numeric(expanded_table["prob_changing"]))<0.6){
    delta_sign <- "~"
  }
  else if (abs(as.numeric(expanded_table["mean_dPSI"]))>0.20 & (as.numeric(expanded_table["prob_changing"]))>=0.6){
    delta_sign <- "+"
  }
  
  PSIs<- data.frame(expanded_table["LSV.ID"], expanded_table["Junctions.coords"], delta_sign)
  colnames(PSIs) <- c("LSV.ID", "Junctions.coords", "delta_sign")
  return(PSIs)
}
PSIs <- apply(Majiq_expanded, 1, function(x)get_PSI(x))
PSIs <-bind_rows(PSIs)

Majiq_expanded<- merge(Majiq_expanded, PSIs)
##PARAMETER
#remove entries with a delta "-" -> removing the entries where the differences in dPSI between replicates is big

merge_SJ_Majiq <-aggregate(cbind(denovo, Junctions.coords, delta_sign, sample_onePSI, sample_twoPSI, mean_dPSI, prob_changing, prob_notchanging) ~LSV.ID + Gene.Name + Gene.ID + Chromosome + strand, data=Majiq_expanded, paste, collapse=";")
merge_SJ_Majiq <- subset(merge_SJ_Majiq, merge_SJ_Majiq$delta_sign!="-")


get_cols<- Majiq_expanded[!duplicated(Majiq_expanded$LSV.ID),]
get_cols<- select(get_cols, LSV.ID, Num.junctions, Num.exons, exon_coords, IR_coords)

merge_SJ_Majiq<- left_join(merge_SJ_Majiq, get_cols, by=c("LSV.ID"))
##PARAMETER
#Filter out the junctions where there are fewer than 2 "+"s -> this gives entries where the LSV mostly 
#has junctions that have (~) with no pair of junctions with (+)
remove_tilde <- function(delta_table){
  deltas <- unlist(strsplit(as.character(delta_table["delta_sign"]), ";"))
  count_plus <- length(which(deltas =="+"))
  
  if (count_plus >=2){
    filtered <- data.frame(delta_table["LSV.ID"], delta_table["Gene.Name"], delta_table["Gene.ID"], 
                           delta_table["Chromosome"], delta_table["strand"], delta_table["denovo"],
                           delta_table["Junctions.coords"], delta_table["sample_onePSI"], delta_table["sample_twoPSI"], 
                           delta_table["mean_dPSI"], delta_table["prob_changing"], delta_table["prob_notchanging"], 
                           delta_table["Num.junctions"], delta_table["Num.exons"], delta_table["exon_coords"], delta_table["IR_coords"], delta_table["delta_sign"])
    colnames(filtered) <- c("LSV.ID", "Gene.Name", "Gene.ID", "Chromosome", "strand", "denovo", "Junctions.coords",
                            "sample_onePSI", "sample_twoPSI", "mean_dPSI", "prob_changing", "prob_notchanging", 
                            "Num.junctions", "Num.exons", "exon_coords", "IR_coords",  "delta_sign")
  }
  else {filtered <- data.frame()}
  return(filtered)
}


Majiq_in_SJ_subset <- apply(merge_SJ_Majiq, 1, function(x)remove_tilde(x))
Majiq_in_SJ_subset <- bind_rows(Majiq_in_SJ_subset)

#Remove the exons that no longer involve the junctions that we have decided to keep
exon_filter<- function(Majiq_SJ_table){
  Num.junctions <- strsplit(as.character(Majiq_SJ_table["Junctions.coords"]), ";")
  junctions_table<- data.frame(Num.junctions)
  colnames(junctions_table) <- "Num.junctions"
  Num.junctions <- length(junctions_table$Num.junctions)
  
  
  junc_coords <- strsplit(as.character(Majiq_SJ_table["Junctions.coords"]), ";|-")
  coords_keep <- unique(as.data.frame(junc_coords))
  colnames(coords_keep)<- "Junctions.coords"
  
  #split the exon coordinates
  exon_coords <- unlist(strsplit(as.character(Majiq_SJ_table["exon_coords"]), ";"))
  #need to change the 'na' values to NA to get numeric values for the coordinates
  exon_coords <-read.table(text = exon_coords, sep = "-", na.strings = c("na"), colClasses = "numeric")
  exon_coords_keep<-apply(coords_keep, 1, function(x)subset(exon_coords, between(x, exon_coords[1], exon_coords[2], NAbounds = FALSE)))
  exon_coords_keep<-unique(bind_rows(exon_coords_keep))
  
  Num.exons <- nrow(exon_coords_keep)
  exon_coords<- paste0(with(exon_coords_keep, paste(V1, V2, sep="-")), collapse = ";")
  
  df_select <- data.frame(Majiq_SJ_table["LSV.ID"], Majiq_SJ_table["Gene.Name"], Majiq_SJ_table["Gene.ID"], 
                          Majiq_SJ_table["Chromosome"], Majiq_SJ_table["strand"], Majiq_SJ_table["denovo"], 
                          Majiq_SJ_table["sample_onePSI"], Majiq_SJ_table["sample_twoPSI"], Majiq_SJ_table["mean_dPSI"], 
                          Majiq_SJ_table["prob_changing"], Majiq_SJ_table["prob_notchanging"],
                          Num.junctions, Num.exons, Majiq_SJ_table["Junctions.coords"], exon_coords, Majiq_SJ_table["IR_coords"], Majiq_SJ_table["delta_sign"])
  colnames(df_select) <- c("LSV.ID", "Gene.Name", "Gene.ID", 
                           "Chromosome", "strand", "denovo", 
                           paste0(sample_name1, "_PSI"), paste0(sample_name2, "_PSI"), "mean_dPSI", 
                           "prob_changing", "prob_notchanging", "Num.junctions", "Num.exons", 
                           "Junctions.coords","Exons.coords", "IR.coords", "delta_sign")
  return(df_select)
  
}
Majiq_in_SJ_exons <- apply(Majiq_in_SJ_subset, 1, function(x)exon_filter(x))
Majiq_in_SJ_exons <-bind_rows(Majiq_in_SJ_exons)
rownames(Majiq_in_SJ_exons) <- 1:nrow(Majiq_in_SJ_exons)

##CHANGE
write.csv(Majiq_in_SJ_exons, file = paste0(sample_name, '_Majiq_deltaPSI_filtered.csv'))

