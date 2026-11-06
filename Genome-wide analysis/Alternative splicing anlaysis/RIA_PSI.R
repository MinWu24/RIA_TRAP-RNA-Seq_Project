library("dplyr")
library("data.table")
library("varhandle")
library("stringr")
library("stringi")
library("readr")

one_tissue <- read.table(file = 'op.psi.tsv', sep = '\t', na.strings = c("", " ", "NA"), header = TRUE, stringsAsFactors = FALSE)

colnames(one_tissue)<- c("Gene.Name", "Gene.ID", "LSV.ID", "meanPSI", "stdevPSI", "LSV.type", "Num_junctions", "Num_exons", "denovo", "Chromosome", "strand", "Junc.coords", "Exon.coords", "IR.coords", "UCSC.LSV.Link")
sample_name<- "op"
#sample_name<- "pa"
#sample_name<- "gacA"
#sample_name<- "pak"

#function taking the PSI Majiq file and removing junctions from each LSV with 0.05<=meanPSI<=0.95, re-counting
# Num.junctions and Num.exons; also removing IR events (they are filtered separately)
split_values <- function(tissue){
  #print(tissue)
  #take all the relevant columns with ; split values and make a matrix data frame where each row has the
  #corresponding values from different columns. Not doing enything with the exon junctions column yet because
  #those values don't directly correspond to the other ; split values.
  meanPSI <- as.numeric(unlist(strsplit(as.character(tissue["meanPSI"]), ";")))
  denovo <- as.numeric(unlist(strsplit(as.character(tissue["denovo"]), ";")))
  junc_coords <- strsplit(as.character(tissue["Junc.coords"]), ";")
  #save the Num.junctions and the Num.exons
  Num_exons <- as.numeric(tissue["Num_exons"])
  Num_junctions <- as.numeric(tissue["Num_junctions"])
  IR_coords <- tissue["IR.coords"]
  #print(IR_coords)
  
  #have to have modified IR coords that are the last exon coords to add to list of kept junction coords
  #if the meanPSI is high enough because the IR coords don't match the exon coords
  if(!is.na(IR_coords)){
    IR_coords_add<- unlist(strsplit(as.character(tissue["IR.coords"]), ";"))
    IR_coords_add <-read.table(text = IR_coords_add, sep = "-", na.strings = c("na"), colClasses = "numeric")
    IR_coords_add$V1<- IR_coords_add$V1 - 1
    IR_coords_add$V2<- IR_coords_add$V2 + 1
    IR_coords_add<- with(IR_coords_add, paste(V1, V2, sep="-"))
    IR_coords_add<- as.data.frame(as.numeric(unlist(strsplit(as.character(IR_coords_add), "-"))))
    colnames(IR_coords_add) <- "coords_keep"
  }
  
  #data frame with all split values
  newdf <- data.frame(meanPSI, denovo, junc_coords, Num_exons, Num_junctions, IR_coords)
  colnames(newdf) <- c("meanPSI", "denovo", "junc_coords", "Num_exons", "Num_junctions", "IR_coords")
  
  ##PARAMETER
  #select for all the entries where tissues have a PSI of >=0.05 and <=0.95
  df_select <- subset(newdf, (newdf$meanPSI>=0.05 & newdf$meanPSI<=0.95))
  
  #save the new number of junctions
  Num_junctions <- length(df_select$junc_coords)
  
  ##Need the junction coordinates of the eliminated junctions to remove the irrelevant exon junctions.
  ##Find the eliminated coordinates -> find which coordinates out of these are not present in the junction coords
  ##that were kept(also want to keep those exons) -> find the exons with that junction and eliminate
  #recall the entries that were eliminated by the cutoff
  df_eliminate <- subset(newdf, (meanPSI<0.05 | meanPSI>0.95))
  
  #Remove the exons from the LSV that correspond to removed junctions
  
  #get the junction coordinates of the removed junctions
  if (nrow(df_eliminate)>0){
    coords_eliminate <- as.data.frame(as.numeric(unlist(strsplit(as.character(df_eliminate$junc_coords), "-"))))
    colnames(coords_eliminate)<- "coords_eliminate"
    #print(coords_eliminate)
    
    coords_keep <- as.data.frame(as.numeric(unlist(strsplit(as.character(df_select$junc_coords), "-"))))
    colnames(coords_keep) <- "coords_keep"
    
    #In case the kept coords correspond to the IR coords
    if(!is.na(IR_coords)){
      
      IR_coords_split <- unique(as.data.frame(as.numeric(unlist(strsplit(as.character(df_select$IR_coords), "-")))))
      colnames(IR_coords_split) <- "IR_coords_split"
      #print(IR_coords_split)
      
      #Check if the coords_keep contain the original IR_coords; if they do, add the modified IR_coords to the list
      # of coords_keep so they will match to the right exon_coord
      is_IR_coords<- subset(coords_keep, (coords_keep %in% IR_coords_split$IR_coords_split))
      #print(is_IR_coords)
      
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
      colnames(coords_keep) <- "coords_keep"}
    
    if(nrow(coords_keep)>0){colnames(coords_keep) <- "coords_keep"
    
    #get the junction coordinates to be eliminated that are not kept
    coords_not_in <- subset(coords_eliminate, !(coords_eliminate %in% coords_keep$coords_keep))
    coords_not_in<-as.data.frame(unlist(apply(coords_not_in, 2, function(x)as.numeric(x))))
    
    #print(coords_keep)
    
    if (nrow(coords_not_in)>0){
      colnames(coords_not_in) <- "coords_eliminate"
      #split the exon coordinates
      exon_coords <- unlist(strsplit(as.character(tissue["Exon.coords"]), ";"))
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
      #print(coords_not_in)
      exon_coords <- anti_join(exon_coords, coords_not_in)
      
      
      colnames(exon_coords)<-c("V1", "V2")
      exon_coords<- unique(rbind(exon_coords, exon_coords_keep))
      #print(exon_coords)
      # change the second junction column to the appropriate column name to match coords_not_in
      if(nrow(exon_coords)>0){colnames(exon_coords)<-c("V1", "coords_eliminate")
      exon_coords <- anti_join(exon_coords, coords_not_in)
      colnames(exon_coords)<-c("V1", "V2")
      exon_coords<- unique(rbind(exon_coords, exon_coords_keep))
      if (nrow(exon_coords)>0){colnames(exon_coords)<-c("V1", "V2")
      Num_exons <- nrow(exon_coords)
      exon_coords<- paste0(with(exon_coords, paste(V1, V2, sep="-")), collapse = ";")
      }
      else {exon_coords<- unique(rbind(exon_coords, exon_coords_keep))
      Num_exons <-nrow(exon_coords)}
      }
      else{exon_coords<- unique(rbind(exon_coords, exon_coords_keep))
      Num_exons <-nrow(exon_coords)}
    }
    else{
      Num_junctions <- length(df_select$junc_coords)
      exon_coords <- tissue["Exon.coords"]}
    }
    else{exon_coords <- tissue["Exon.coords"]}
  }
  else{exon_coords <- tissue["Exon.coords"]}
  
  #Recovering each LSV line in its original format but with updated junction and exon numbers based on the 
  #removal of junctions and exons
  meanPSI <- paste0(df_select$meanPSI, collapse=";")
  denovo <- paste0(df_select$denovo, collapse=";")
  junc_coords <- paste0(df_select$junc_coords, collapse=";")
  
  df_select2<- data.frame(meanPSI, denovo, junc_coords, exon_coords, Num_junctions, Num_exons)
  df_select2$Gene.Name <-tissue["Gene.Name"]
  df_select2$Gene.ID <-tissue["Gene.ID"]
  df_select2$LSV.ID <-tissue["LSV.ID"]
  df_select2$Chromosome <-tissue["Chromosome"]
  df_select2$strand <-tissue["strand"]
  df_select2$IR_coords <- IR_coords
  rownames(df_select2)<- 1:nrow(df_select2)
  return(df_select2)
}

tissue_split <- apply(one_tissue, 1, function(X)split_values(X))

tissue_split <-bind_rows(tissue_split)

tissue_split<- tissue_split[c("Gene.Name", "Gene.ID", "LSV.ID", "Chromosome", "strand", "meanPSI", "Num_junctions", "Num_exons", "denovo", "junc_coords", "exon_coords", "IR_coords")]
colnames(tissue_split)<- c("Gene.Name", "Gene.ID", "LSV.ID", "Chromosome", "strand", "meanPSI", "Num_junctions", "Num_exons", "denovo", "junc_coords", "exon_coords", "IR.coords")

#replace blanks with NA and then remove LSV lines that no longer have at least 1 junction
tissue_split<- apply(tissue_split, 2, function(x) gsub("^$|^ $", NA, x))
tissue_split<- data.frame(tissue_split)
tissue_split <- subset(tissue_split, !is.na(tissue_split$junc_coords))

##
##CHANGE
#Choose the SJ.out.tab files to do the replicate-centric analysis on


rep1 <- read.table(file = 'OP1_SJ.out.tab', sep = '\t', header = FALSE)
colnames(rep1) <- c("Chromosome", "junction_start", "junction_end", "strand", "intron_motif", "Annotated", "Num_unique_reads", "Num_multimapping_reads", "Max_spliced_alignment_overhang")

rep3 <- read.table(file = 'OP3_SJ.out.tab', sep = '\t', header = FALSE, stringsAsFactors = FALSE)
colnames(rep3) <- c("Chromosome", "junction_start", "junction_end", "strand", "intron_motif", "Annotated", "Num_unique_reads", "Num_multimapping_reads", "Max_spliced_alignment_overhang")

rep4 <- read.table(file = 'OP4_SJ.out.tab', sep = '\t', header = FALSE)
colnames(rep4) <- c("Chromosome", "junction_start", "junction_end", "strand", "intron_motif", "Annotated", "Num_unique_reads", "Num_multimapping_reads", "Max_spliced_alignment_overhang")

rep5 <- read.table(file = 'OP5_SJ.out.tab', sep = '\t', header = FALSE, stringsAsFactors = FALSE)
colnames(rep5) <- c("Chromosome", "junction_start", "junction_end", "strand", "intron_motif", "Annotated", "Num_unique_reads", "Num_multimapping_reads", "Max_spliced_alignment_overhang")

rep6 <- read.table(file = 'OP6_SJ.out.tab', sep = '\t', header = FALSE, stringsAsFactors = FALSE)
colnames(rep6) <- c("Chromosome", "junction_start", "junction_end", "strand", "intron_motif", "Annotated", "Num_unique_reads", "Num_multimapping_reads", "Max_spliced_alignment_overhang")




numpy_majiq_file1<- read.table(file = 'OP1_numpy.tsv', sep = '\t', na.strings = c("", " ", "NA"), header = TRUE, stringsAsFactors = FALSE)
colnames(numpy_majiq_file1) <- c("LSV.ID", "junction_start", "junction_end", "Num_unique_reads")

numpy_majiq_file3<- read.table(file = 'OP3_numpy.tsv', sep = '\t', na.strings = c("", " ", "NA"), header = TRUE, stringsAsFactors = FALSE)
colnames(numpy_majiq_file3) <- c("LSV.ID", "junction_start", "junction_end", "Num_unique_reads")

numpy_majiq_file4<- read.table(file = 'OP4_numpy.tsv', sep = '\t', na.strings = c("", " ", "NA"), header = TRUE, stringsAsFactors = FALSE)
colnames(numpy_majiq_file4) <- c("LSV.ID", "junction_start", "junction_end", "Num_unique_reads")

numpy_majiq_file5<- read.table(file = 'OP5_numpy.tsv', sep = '\t', na.strings = c("", " ", "NA"), header = TRUE, stringsAsFactors = FALSE)
colnames(numpy_majiq_file5) <- c("LSV.ID", "junction_start", "junction_end", "Num_unique_reads")

numpy_majiq_file6<- read.table(file = 'OP6_numpy.tsv', sep = '\t', na.strings = c("", " ", "NA"), header = TRUE, stringsAsFactors = FALSE)
colnames(numpy_majiq_file6) <- c("LSV.ID", "junction_start", "junction_end", "Num_unique_reads")


numpy_majiq_file1$LSV.ID<- str_replace_all(numpy_majiq_file1$LSV.ID, "^b", "")
numpy_majiq_file3$LSV.ID<- str_replace_all(numpy_majiq_file3$LSV.ID, "^b", "")
numpy_majiq_file4$LSV.ID<- str_replace_all(numpy_majiq_file4$LSV.ID, "^b", "")
numpy_majiq_file5$LSV.ID<- str_replace_all(numpy_majiq_file5$LSV.ID, "^b", "")
numpy_majiq_file6$LSV.ID<- str_replace_all(numpy_majiq_file6$LSV.ID, "^b", "")

#
# do a nt shift to match the counting scheme of Majiq
rep1$junction_start <- rep1$junction_start - 1
rep1$junction_end <- rep1$junction_end + 1

rep3$junction_start <- rep3$junction_start - 1
rep3$junction_end <- rep3$junction_end + 1

rep4$junction_start <- rep4$junction_start - 1
rep4$junction_end <- rep4$junction_end + 1

rep5$junction_start <- rep5$junction_start - 1
rep5$junction_end <- rep5$junction_end + 1

rep6$junction_start <- rep6$junction_start - 1
rep6$junction_end <- rep6$junction_end + 1

# #remove rows with 0 reads in the numpy file
numpy_majiq_file1<- subset(numpy_majiq_file1, !(numpy_majiq_file1$Num_unique_reads==0))
numpy_majiq_file3<- subset(numpy_majiq_file3, !(numpy_majiq_file3$Num_unique_reads==0))
numpy_majiq_file4<- subset(numpy_majiq_file4, !(numpy_majiq_file4$Num_unique_reads==0))
numpy_majiq_file5<- subset(numpy_majiq_file5, !(numpy_majiq_file5$Num_unique_reads==0))
numpy_majiq_file6<- subset(numpy_majiq_file6, !(numpy_majiq_file6$Num_unique_reads==0))



##subsets with "Chromosome", "junction_start", "junction_end", "strand", "Annotated", "Num_unique_reads"
subset_table <-function (table_to_subset){
  newtable <- table_to_subset[c("Chromosome", "junction_start", "junction_end", "strand", "Num_unique_reads")]
  return (newtable)
}
##CHANGE

replicate_list <- list(rep1, rep3, rep4, rep5, rep6)
replicate_list <-lapply(replicate_list, subset_table)



# get the junctions and the counts associated with each for both replicates
merged_tables <-Reduce(function(x,y) merge(x, y, by=c("Chromosome", "junction_start", "junction_end", "strand"), all=TRUE), replicate_list)
colnames(merged_tables)[5] <-"Num_unique_reads_T1.1"
colnames(merged_tables)[6] <-"Num_unique_reads_T1.3"
colnames(merged_tables)[7] <-"Num_unique_reads_T1.4"
colnames(merged_tables)[8] <-"Num_unique_reads_T1.5"
colnames(merged_tables)[9] <-"Num_unique_reads_T1.6"

#remove duplicates
merged_tables<- merged_tables[!duplicated(merged_tables[c("Chromosome", "junction_start", "junction_end","Num_unique_reads_T1.1", "Num_unique_reads_T1.3", "Num_unique_reads_T1.4", "Num_unique_reads_T1.5", "Num_unique_reads_T1.6")]),]

#convert NA to 0 for num reads
merged_tables[c("Num_unique_reads_T1.1", "Num_unique_reads_T1.3", "Num_unique_reads_T1.4", "Num_unique_reads_T1.5", "Num_unique_reads_T1.6")][is.na(merged_tables[c("Num_unique_reads_T1.1", "Num_unique_reads_T1.3", "Num_unique_reads_T1.4", "Num_unique_reads_T1.5", "Num_unique_reads_T1.6")])] <- 0

##PARAMETER
#add a 1 to the number of reads for every junction so there is no 0 as a denominator
merged_tables$Num_unique_reads_T1.1 <- merged_tables$Num_unique_reads_T1.1 + 1
merged_tables$Num_unique_reads_T1.3 <- merged_tables$Num_unique_reads_T1.3 + 1
merged_tables$Num_unique_reads_T1.4 <- merged_tables$Num_unique_reads_T1.4 + 1
merged_tables$Num_unique_reads_T1.5 <- merged_tables$Num_unique_reads_T1.5 + 1
merged_tables$Num_unique_reads_T1.6 <- merged_tables$Num_unique_reads_T1.6 + 1

#Get the reads for each of the junction coordinates and concantenate the junction pair coordinates together with a hyphen
#to match the junction coordinates from the Majiq file
#E.g. "Chromosome", "junction_start", "junction_end", "strand"  "Num_unique_reads_T1.1" "Num_unique_reads_T1.1" "Junctions.coords"
#      I	            2712	            3203	          2	          1	                          3	                    2712-3203
merged_tables<- cbind(merged_tables, Junctions.coords = apply(merged_tables[2:3], 1, function(x) paste(na.omit(x), collapse = "-")))
merged_tables<- merged_tables[!duplicated(merged_tables[c("Chromosome", "strand", "Junctions.coords")]),]

# ############# Do the same thing for the IR list
replicate_list_IR<- list(numpy_majiq_file1, numpy_majiq_file3, numpy_majiq_file4, numpy_majiq_file5, numpy_majiq_file6)
#same for numpy files
merged_tables_IR <-Reduce(function(x,y) merge(x, y, by=c("LSV.ID", "junction_start", "junction_end"), all=TRUE), replicate_list_IR)
colnames(merged_tables_IR)[4] <-"Num_unique_reads_T1.1"
colnames(merged_tables_IR)[5] <-"Num_unique_reads_T1.3"
colnames(merged_tables_IR)[6] <-"Num_unique_reads_T1.4"
colnames(merged_tables_IR)[7] <-"Num_unique_reads_T1.5"
colnames(merged_tables_IR)[8] <-"Num_unique_reads_T1.6"
#remove duplicates
merged_tables_IR<- merged_tables_IR[!duplicated(merged_tables_IR[c("LSV.ID", "junction_start", "junction_end","Num_unique_reads_T1.1", "Num_unique_reads_T1.3", "Num_unique_reads_T1.4", "Num_unique_reads_T1.5", "Num_unique_reads_T1.6")]),]

#convert NA to 0 for num reads

merged_tables_IR[c("Num_unique_reads_T1.1", "Num_unique_reads_T1.3", "Num_unique_reads_T1.4", "Num_unique_reads_T1.5", "Num_unique_reads_T1.6")][is.na(merged_tables_IR[c("Num_unique_reads_T1.1", "Num_unique_reads_T1.3", "Num_unique_reads_T1.4", "Num_unique_reads_T1.5", "Num_unique_reads_T1.6")])] <- 0

##PARAMETER
#add a 1 to the number of reads for every junction so there is no 0 as a denominator
merged_tables_IR$Num_unique_reads_T1.1 <- merged_tables_IR$Num_unique_reads_T1.1 + 1
merged_tables_IR$Num_unique_reads_T1.3 <- merged_tables_IR$Num_unique_reads_T1.3 + 1
merged_tables_IR$Num_unique_reads_T1.4 <- merged_tables_IR$Num_unique_reads_T1.4 + 1
merged_tables_IR$Num_unique_reads_T1.5 <- merged_tables_IR$Num_unique_reads_T1.5 + 1
merged_tables_IR$Num_unique_reads_T1.6 <- merged_tables_IR$Num_unique_reads_T1.6 + 1

#Get the reads for each of the junction coordinates and concantenate the junction pair coordinates together with a hyphen
#to match the junction coordinates from the Majiq file
#E.g. "Chromosome", "junction_start", "junction_end", "strand"  "Num_unique_reads_T1.1" "Num_unique_reads_T1.1" "Junctions.coords"
#      I	            2712	            3203	          2	          1	                          3	                    2712-3203
merged_tables_IR<- cbind(merged_tables_IR, Junctions.coords = apply(merged_tables_IR[2:3], 1, function(x) paste(na.omit(x), collapse = "-")))
merged_tables_IR<- merged_tables_IR[!duplicated(merged_tables_IR[c("LSV.ID", "Junctions.coords")]),]


#Split the Majiq LSVs into separate junction coordinates to be able to match the SJ.out files
#E.g. "Junctions.coords" "meanPSI"     "Gene.Name"        "Gene.ID"          "LSV.ID"           "Chromosome" "strand"
#       8887285-8887796	0.93220341205596	C09G9.5	         gene:WBGene00007495	gene:WBGene00007495:t:8887796-8887904	IV	+
expand_Majiq<- function(majiq_table){
  junction_coords <- strsplit(as.character(majiq_table["junc_coords"]), ";")
  meanPSI <- strsplit(as.character(majiq_table["meanPSI"]), ";")
  coord_df <-data.frame(junction_coords, meanPSI)
  colnames(coord_df) <-c("Junctions.coords", "meanPSI")
  coord_df <- cbind(coord_df, majiq_table["Gene.Name"], majiq_table["Gene.ID"], majiq_table["LSV.ID"], majiq_table["Chromosome"], majiq_table["strand"])
  colnames(coord_df) <- c("Junctions.coords", "meanPSI", "Gene.Name", "Gene.ID", "LSV.ID", "Chromosome", "strand")
  
  if (!is.na(majiq_table["IR.coords"])){
    coord_df<- coord_df[-nrow(coord_df),]
    return(coord_df)
  }
  else{return(coord_df)}
}

#expand the Majiq LSVs to just be junctions with the appropriate PSI values
expanded_tissue <- apply(tissue_split, 1, function(x)expand_Majiq(x))
expanded_tissue <-bind_rows(expanded_tissue)

#merge the expanded Majiq junctions with the PSI values to the SJ junction and read counts ->gives a table of
#the Majiq junctions with the associated read counts from the SJ.out files
Majiq_in_SJ <- left_join(expanded_tissue, merged_tables, by=c("Junctions.coords", "Chromosome"))

#get rid of junctions that do not have SJ.out support (this includes IR events)
Majiq_in_SJ <- subset(Majiq_in_SJ, !is.na(junction_start))
Majiq_in_SJ<- select(Majiq_in_SJ, Junctions.coords, meanPSI, Gene.Name, Gene.ID, LSV.ID, Chromosome, strand.x, junction_start, junction_end, Num_unique_reads_T1.1, Num_unique_reads_T1.3, Num_unique_reads_T1.4, Num_unique_reads_T1.5, Num_unique_reads_T1.6)
colnames(Majiq_in_SJ)[colnames(Majiq_in_SJ) =="strand.x"] <- "strand"

###############
expand_Majiq<- function(majiq_table){
  junction_coords <- strsplit(as.character(majiq_table["junc_coords"]), ";")
  meanPSI <- strsplit(as.character(majiq_table["meanPSI"]), ";")
  coord_df <-data.frame(junction_coords, meanPSI)
  colnames(coord_df) <-c("Junctions.coords", "meanPSI")
  coord_df <- cbind(coord_df, majiq_table["Gene.Name"], majiq_table["Gene.ID"], majiq_table["LSV.ID"], majiq_table["Chromosome"], majiq_table["strand"])
  colnames(coord_df) <- c("Junctions.coords", "meanPSI", "Gene.Name", "Gene.ID", "LSV.ID", "Chromosome", "strand")
  
  if (is.na(majiq_table["IR.coords"])){
    coord_df<- NULL
    return(coord_df)
  }
  else{
    coord_df<- coord_df[nrow(coord_df),]
    return(coord_df)
  }
}

#expand the Majiq LSVs to just be junctions with the appropriate PSI values
expanded_tissue_IR <- apply(tissue_split, 1, function(x)expand_Majiq(x))
expanded_tissue_IR <-bind_rows(expanded_tissue_IR)

#merge the expanded Majiq junctions with the PSI values to the SJ junction and read counts ->gives a table of
#the Majiq junctions with the associated read counts from the SJ.out files
Majiq_in_SJ_IR <- left_join(expanded_tissue_IR, merged_tables_IR, by=c("Junctions.coords", "LSV.ID"))

#get rid of junctions that do not have SJ.out support (this includes IR events)
Majiq_in_SJ_IR <- subset(Majiq_in_SJ_IR, !is.na(junction_start))

###

Majiq_in_SJ<- bind_rows(Majiq_in_SJ, Majiq_in_SJ_IR)
Majiq_in_SJ<- Majiq_in_SJ[order(Majiq_in_SJ$LSV.ID),]





#split the table by LSV -> this is in order to calculate the read sums, PSI values for the Majiq junctions using
#the SJ.out reads
LSVlist <-unique(Majiq_in_SJ$LSV.ID)
mergingtable <-Majiq_in_SJ %>% filter(LSV.ID %in% LSVlist)
mergingtable$LSV.ID <-factor(mergingtable$LSV.ID)
LSV_split_tables<- split(mergingtable, f=mergingtable$LSV.ID)
#gives the sum of all tissue reads for all junctions in an LSV -> use this as the denominator to calculate PSI
#E.g. Junction coordinates  Gene  Same LSV.ID                           Chr Strand    Rep1    Rep2
#     7062058-7062189	Y57E12AL.1	gene:WBGene00021956:t:7062189-7062327	V	+	          60	    136
#     7062076-7062189	Y57E12AL.1	gene:WBGene00021956:t:7062189-7062327	V	+	          35	    116
# summed the total reads for both junctions per replicate                             95      252
find_sums <- function(split_tables){
  sums_table <-aggregate(cbind(Num_unique_reads_T1.1, Num_unique_reads_T1.3, Num_unique_reads_T1.4, Num_unique_reads_T1.5, Num_unique_reads_T1.6) ~LSV.ID + Gene.Name + Gene.ID + Chromosome, data=split_tables, sum)
  colnames(sums_table) <- c("LSV.ID", "Gene.Name", "Gene.ID", "Chromosome", "T1_sum_1", "T1_sum_3", "T1_sum_4", "T1_sum_5", "T1_sum_6")
  return (sums_table)
}
LSV_junctions_summed <- lapply(LSV_split_tables, function(x) find_sums(x))
LSV_junctions_summed <-bind_rows(LSV_junctions_summed)

#add the sums to the read counts table
Majiq_in_SJ_sums <- merge(Majiq_in_SJ, LSV_junctions_summed, by=c("LSV.ID", "Gene.Name", "Gene.ID", "Chromosome"))
#E.g. Junction coordinates  Gene  Same LSV.ID                           Chr Strand    Rep1    Rep2  Sum for Rep1  Sum for Rep2
#     7062058-7062189	Y57E12AL.1	gene:WBGene00021956:t:7062189-7062327	V	+	          60	    136       95          252
#     7062076-7062189	Y57E12AL.1	gene:WBGene00021956:t:7062189-7062327	V	+	          35	    116       95          252

##PARAMETER
#at least two replicate summed counts for the LSV must be >=20
Majiq_in_SJ_sums <- Majiq_in_SJ_sums[rowSums(cbind(Majiq_in_SJ_sums$T1_sum_1 >=20, Majiq_in_SJ_sums$T1_sum_3 >=20, Majiq_in_SJ_sums$T1_sum_4 >=20, Majiq_in_SJ_sums$T1_sum_5 >=20, Majiq_in_SJ_sums$T1_sum_6 >=20)) > 1,]

#########
#Calculate PSI values from the SJ.out files by dividing Rep read count/sum for that Rep
#E.g. Junction coordinates  Gene  Same LSV.ID                           Chr Strand    Rep1    Rep2  Sum for Rep1  Sum for Rep2  PSI_Rep1      PSI_Rep2
#     7062058-7062189	Y57E12AL.1	gene:WBGene00021956:t:7062189-7062327	V	+	          60	    136       95          252         60/95=0.63  136/252=0.54
#     7062076-7062189	Y57E12AL.1	gene:WBGene00021956:t:7062189-7062327	V	+	          35	    116       95          252         0.37        0.46
#assign a delta sign to the junctions to filter for low confidence changes from the SJ files
get_PSI <- function(to_divide){
  T1_PSI_1 <- as.numeric(to_divide["Num_unique_reads_T1.1"])/ as.numeric(to_divide["T1_sum_1"])
  T1_PSI_3 <- as.numeric(to_divide["Num_unique_reads_T1.3"])/ as.numeric(to_divide["T1_sum_3"])
  T1_PSI_4 <- as.numeric(to_divide["Num_unique_reads_T1.4"])/ as.numeric(to_divide["T1_sum_4"])
  T1_PSI_5 <- as.numeric(to_divide["Num_unique_reads_T1.5"])/ as.numeric(to_divide["T1_sum_5"])
  T1_PSI_6 <- as.numeric(to_divide["Num_unique_reads_T1.6"])/ as.numeric(to_divide["T1_sum_6"])
  PSIs <- data.frame(to_divide["LSV.ID"], to_divide["Gene.Name"], to_divide["Gene.ID"], to_divide["Chromosome"], to_divide["strand"], to_divide["Junctions.coords"], to_divide["meanPSI"], 
                     to_divide["Num_unique_reads_T1.1"], to_divide["Num_unique_reads_T1.3"], 
                     to_divide["Num_unique_reads_T1.4"], to_divide["Num_unique_reads_T1.5"], 
                     to_divide["Num_unique_reads_T1.6"], 
                     to_divide["T1_sum_1"], to_divide["T1_sum_3"], 
                     to_divide["T1_sum_4"], to_divide["T1_sum_5"],
                     to_divide["T1_sum_6"],
                     T1_PSI_1, T1_PSI_3, T1_PSI_4, T1_PSI_5, T1_PSI_6)
  colnames(PSIs) <- c("LSV.ID", "Gene.Name", "Gene.ID", "Chromosome", "strand", "Junctions.coords", "meanPSI", 
                      "Num_unique_reads_T1.1", "Num_unique_reads_T1.3", 
                      "Num_unique_reads_T1.4", "Num_unique_reads_T1.5",
                      "Num_unique_reads_T1.6",
                      "T1_sum_1", "T1_sum_3", 
                      "T1_sum_4", "T1_sum_5",
                      "T1_sum_6",
                      "T1_PSI_1", "T1_PSI_3", "T1_PSI_4", "T1_PSI_5", "T1_PSI_6")
  
  #PARAMETER
  df_Tsums<- data.frame(as.numeric(to_divide["T1_sum_1"]), as.numeric(to_divide["T1_sum_3"]), 
                        as.numeric(to_divide["T1_sum_4"]), as.numeric(to_divide["T1_sum_5"]),
                        as.numeric(to_divide["T1_sum_6"]))
  df_Tsums$rows_above <- rowSums(df_Tsums >=20)
  
  PSI_values<- c()

  if (df_Tsums$rows_above == 2 | df_Tsums$rows_above == 3 | df_Tsums$rows_above == 4 | df_Tsums$rows_above == 5){

    if (as.numeric(to_divide["T1_sum_1"])>=20){
      PSI_values<- c(PSI_values, T1_PSI_1)
    }
    if (as.numeric(to_divide["T1_sum_3"])>=20){
      PSI_values<- c(PSI_values, T1_PSI_3)
    }
    if (as.numeric(to_divide["T1_sum_4"])>=20){
      PSI_values<- c(PSI_values, T1_PSI_4)
    }
    if (as.numeric(to_divide["T1_sum_5"])>=20){
      PSI_values<- c(PSI_values, T1_PSI_5)
    }
    if (as.numeric(to_divide["T1_sum_6"])>=20){
      PSI_values<- c(PSI_values, T1_PSI_6)
    }
    
    
    PSI_stdev<- sd(PSI_values)
    
    if (PSI_stdev<=0.10) {
      delta_sign <- "+"
    }
    else if (PSI_stdev<=0.15) {
      delta_sign <- "~"
    }
    else {delta_sign <- "-"}
    PSIs<- data.frame(PSIs, delta_sign)
    colnames(PSIs) <- c("LSV.ID", "Gene.Name", "Gene.ID", "Chromosome", "strand", "Junctions.coords", "meanPSI", 
                        "Num_unique_reads_T1.1", "Num_unique_reads_T1.3", 
                        "Num_unique_reads_T1.4", "Num_unique_reads_T1.5",
                        "Num_unique_reads_T1.6", 
                        "T1_sum_1", "T1_sum_3", 
                        "T1_sum_4", "T1_sum_5", 
                        "T1_sum_6", 
                        "T1_PSI_1", "T1_PSI_3", "T1_PSI_4", "T1_PSI_5", "T1_PSI_6", "delta_sign")
    return(PSIs)
  }
  
  if (df_Tsums$rows_above > 5){
    PSI_list<- c(T1_PSI_1, T1_PSI_3, T1_PSI_4, T1_PSI_5, T1_PSI_6)
    
    
    PSI_stdev<- sd(PSI_list)
    
    
    
    if (PSI_stdev<=0.10) {
      delta_sign <- "+"
    }
    else if (PSI_stdev<=0.15) {
      delta_sign <- "~"
    }
    else {delta_sign <- "-"}
    PSIs<- data.frame(PSIs, delta_sign)
    colnames(PSIs) <- c("LSV.ID", "Gene.Name", "Gene.ID", "Chromosome", "strand", "Junctions.coords", "meanPSI", 
                        "Num_unique_reads_T1.1", "Num_unique_reads_T1.3", 
                        "Num_unique_reads_T1.4", "Num_unique_reads_T1.5",
                        "Num_unique_reads_T1.6", 
                        "T1_sum_1", "T1_sum_3", 
                        "T1_sum_4", "T1_sum_5", 
                        "T1_sum_6",
                        "T1_PSI_1", "T1_PSI_3", "T1_PSI_4", "T1_PSI_5", "T1_PSI_6", "delta_sign")
    
    return(PSIs)
  }
}
PSIs <- apply(Majiq_in_SJ_sums, 1, function(x)get_PSI(x))
PSIs <-bind_rows(PSIs)

##PARAMETER

#remove entries with a delta "-" -> removing the entries where the differences in PSI between replicates is big
#PSIs <- subset(PSIs, PSIs$delta_sign!="-")
PSIs <- select(PSIs, "LSV.ID", "Gene.Name", "Gene.ID", "Chromosome", "Junctions.coords", "meanPSI", 
               "Num_unique_reads_T1.1", "Num_unique_reads_T1.3", 
               "Num_unique_reads_T1.4", "Num_unique_reads_T1.5",
               "Num_unique_reads_T1.6", 
               "T1_sum_1", "T1_sum_3", 
               "T1_sum_4", "T1_sum_5", 
               "T1_sum_6", 
               "T1_PSI_1", "T1_PSI_3", "T1_PSI_4", "T1_PSI_5", "T1_PSI_6", "delta_sign")
rownames(PSIs)<- 1:nrow(PSIs)

## Combine the delta signs with the full table for each junction from the Majiq and from the SJ, make everything ; separated
split_majiq<- function(majiq_table){
  junction_coords <- strsplit(as.character(majiq_table["junc_coords"]), ";")
  meanPSI <- strsplit(as.character(majiq_table["meanPSI"]), ";")
  denovo <- strsplit(as.character(majiq_table["denovo"]), ";")
  coord_df <-data.frame(junction_coords, meanPSI, denovo)
  colnames(coord_df) <-c("Junctions.coords", "meanPSI","denovo")
  coord_df <- cbind(coord_df, majiq_table["Gene.Name"], majiq_table["Gene.ID"], majiq_table["LSV.ID"], majiq_table["Chromosome"], majiq_table["strand"])
  colnames(coord_df) <- c("Junctions.coords", "meanPSI", "denovo", "Gene.Name", "Gene.ID", "LSV.ID", "Chromosome", "strand")
  return(coord_df)
}

Majiq_expanded <- apply(tissue_split, 1, function(x)split_majiq(x))
Majiq_expanded <-bind_rows(Majiq_expanded)

merge_SJ_Majiq<- left_join(PSIs, Majiq_expanded, by=c("Junctions.coords", "Gene.Name", "Gene.ID", "LSV.ID", "Chromosome", "meanPSI"))
merge_SJ_Majiq<- select(merge_SJ_Majiq, LSV.ID, Gene.Name, Gene.ID, Chromosome, strand, denovo, meanPSI, Junctions.coords, 
                        Num_unique_reads_T1.1, Num_unique_reads_T1.3, 
                        Num_unique_reads_T1.4, Num_unique_reads_T1.5, 
                        Num_unique_reads_T1.6,
                        T1_sum_1, T1_sum_3, T1_sum_4, T1_sum_5, T1_sum_6,  
                        T1_PSI_1, T1_PSI_3, T1_PSI_4, T1_PSI_5, T1_PSI_6, delta_sign)
merge_SJ_Majiq <-aggregate(cbind(denovo, meanPSI, Junctions.coords, 
                                 Num_unique_reads_T1.1, Num_unique_reads_T1.3, 
                                 Num_unique_reads_T1.4, Num_unique_reads_T1.5, 
                                 Num_unique_reads_T1.6,
                                 T1_sum_1, T1_sum_3, T1_sum_4, T1_sum_5, T1_sum_6,  
                                 T1_PSI_1, T1_PSI_3, T1_PSI_4, T1_PSI_5, T1_PSI_6, delta_sign) ~LSV.ID + Gene.Name + Gene.ID + Chromosome + strand, data=merge_SJ_Majiq, paste, collapse=";")
merge_SJ_Majiq<- select(merge_SJ_Majiq, LSV.ID, Gene.Name, Gene.ID, Chromosome, strand, denovo, meanPSI, Junctions.coords, 
                        Num_unique_reads_T1.1, Num_unique_reads_T1.3, 
                        Num_unique_reads_T1.4, Num_unique_reads_T1.5, 
                        Num_unique_reads_T1.6,
                        T1_sum_1, T1_sum_3, T1_sum_4, T1_sum_5, T1_sum_6,  
                        T1_PSI_1, T1_PSI_3, T1_PSI_4, T1_PSI_5, T1_PSI_6, delta_sign)
merge_SJ_Majiq <- left_join(merge_SJ_Majiq, tissue_split, by=c("LSV.ID", "Gene.Name", "Gene.ID", "Chromosome", "strand"))

#get rid of LSVs with just one junction
merge_SJ_Majiq <- subset(merge_SJ_Majiq, merge_SJ_Majiq$delta_sign!="-")
merge_SJ_Majiq <- subset(merge_SJ_Majiq, !(merge_SJ_Majiq$delta_sign == "+" | merge_SJ_Majiq$delta_sign == "~"))
#Optional for High Stringency
#merge_SJ_Majiq<- subset(merge_SJ_Majiq, !(merge_SJ_Majiq$delta_sign == "~;~" | merge_SJ_Majiq$delta_sign == "~;~;~"))

##PARAMETER
#Filter out the junctions where there are fewer than 2 "+"s -> this gives entries where the LSV mostly 
#has junctions that have (~) with no pair of junctions with (+)
remove_tilde <- function(delta_table){
  deltas <- unlist(strsplit(as.character(delta_table["delta_sign"]), ";"))
  count_plus <- length(which(deltas =="+"))
  
  if (count_plus >=2){
    filtered <- data.frame(delta_table["LSV.ID"], delta_table["Gene.Name"], delta_table["Gene.ID"], delta_table["Chromosome"], 
                           delta_table["strand"], delta_table["denovo.x"], delta_table["meanPSI.x"], delta_table["Junctions.coords"], 
                           delta_table["exon_coords"], delta_table["IR.coords"], delta_table["Num_junctions"], delta_table["Num_exons"],
                           delta_table["Num_unique_reads_T1.1"], delta_table["Num_unique_reads_T1.3"], 
                           delta_table["Num_unique_reads_T1.4"], delta_table["Num_unique_reads_T1.5"],
                           delta_table["Num_unique_reads_T1.6"], 
                           delta_table["T1_sum_1"], delta_table["T1_sum_3"], 
                           delta_table["T1_sum_4"], delta_table["T1_sum_5"], 
                           delta_table["T1_sum_6"], 
                           delta_table["T1_PSI_1"], delta_table["T1_PSI_3"], 
                           delta_table["T1_PSI_4"], delta_table["T1_PSI_5"], 
                           delta_table["T1_PSI_6"], 
                           delta_table["delta_sign"])
    
    colnames(filtered) <- c("LSV.ID", "Gene.Name", "Gene.ID", "Chromosome", 
                            "strand", "denovo.x", "meanPSI.x", "Junctions.coords",
                            "exon_coords", "IR.coords", "Num.junctions", "Num.exons",
                            "Num_unique_reads_T1.1", "Num_unique_reads_T1.3", 
                            "Num_unique_reads_T1.4", "Num_unique_reads_T1.5",
                            "Num_unique_reads_T1.6",  
                            "T1_sum_1", "T1_sum_3", 
                            "T1_sum_4", "T1_sum_5", 
                            "T1_sum_6", 
                            "T1_PSI_1", "T1_PSI_3", "T1_PSI_4", "T1_PSI_5", "T1_PSI_6", "delta_sign")
  }
  else {filtered <- data.frame()}
  return(filtered)
}


Majiq_in_SJ_subset <- apply(merge_SJ_Majiq, 1, function(x)remove_tilde(x))
Majiq_in_SJ_subset <- bind_rows(Majiq_in_SJ_subset)


#Remove the exons that no longer involve the junctions that we have decided to keep
exon_filter<- function(Majiq_SJ_table){
  Num_junctions <- strsplit(as.character(Majiq_SJ_table["Junctions.coords"]), ";")
  junctions_table<- data.frame(Num_junctions)
  colnames(junctions_table) <- "Num_junctions"
  Num_junctions <- length(junctions_table$Num_junctions)
  
  junc_coords <- as.data.frame(as.numeric(unlist(strsplit(as.character(Majiq_SJ_table["Junctions.coords"]), ";|-"))))
  coords_keep <- unique(as.data.frame(junc_coords))
  colnames(coords_keep)<- "Junctions.coords"
  
  
  
  #make modified IR.coords if there are any
  if (!is.na(Majiq_SJ_table["IR.coords"])){
    IR_coords_add<- unlist(strsplit(as.character(Majiq_SJ_table["IR.coords"]), ";"))
    IR_coords_add <-read.table(text = IR_coords_add, sep = "-", na.strings = c("na"), colClasses = "numeric")
    IR_coords_add$V1<- IR_coords_add$V1 - 1
    IR_coords_add$V2<- IR_coords_add$V2 + 1
    IR_coords_add<- with(IR_coords_add, paste(V1, V2, sep="-"))
    IR_coords_add<- as.data.frame(as.numeric(unlist(strsplit(as.character(IR_coords_add), "-"))))
    colnames(IR_coords_add)<- "Junctions.coords"
    
  }
  
  if(!is.na(Majiq_SJ_table["IR.coords"])){
    
    IR_coords_split <- as.data.frame(as.numeric(unlist(strsplit(as.character(Majiq_SJ_table["IR.coords"]), "-"))))
    colnames(IR_coords_split) <- "IR_coords_split"
    
    #Check if the coords_keep contain the original IR_coords; if they do, add the modified IR_coords to the list
    # of coords_keep so they will match to the right exon_coord
    is_IR_coords<- subset(coords_keep, (Junctions.coords %in% IR_coords_split$IR_coords_split))
    
    if (nrow(is_IR_coords)>0){
      coords_keep<- rbind(coords_keep, IR_coords_add)
    }
    else{
      junc_coords <- as.data.frame(as.numeric(unlist(strsplit(as.character(Majiq_SJ_table["Junctions.coords"]), ";|-"))))
      coords_keep<- unique(as.data.frame(junc_coords))
      colnames(coords_keep) <- "Junctions.coords"
    }
    
  }
  else{
    junc_coords <- as.data.frame(as.numeric(unlist(strsplit(as.character(Majiq_SJ_table["Junctions.coords"]), ";|-"))))
    coords_keep<- unique(as.data.frame(junc_coords))
    colnames(coords_keep) <- "Junctions.coords"
  }
  
  #split the exon coordinates
  exon_coords <- unlist(strsplit(as.character(Majiq_SJ_table["exon_coords"]), ";"))
  #need to change the 'na' values to NA to get numeric values for the coordinates
  exon_coords <-read.table(text = exon_coords, sep = "-", na.strings = c("na"), colClasses = "numeric")
  exon_coords_keep<-apply(coords_keep, 1, function(x)subset(exon_coords, between(x, exon_coords[1], exon_coords[2], NAbounds = FALSE)))
  exon_coords_keep<-unique(bind_rows(exon_coords_keep))
  
  Num_exons <- nrow(exon_coords_keep)
  exon_coords<- paste0(with(exon_coords_keep, paste(V1, V2, sep="-")), collapse = ";")
  
  
  df_select <- data.frame(Majiq_SJ_table["LSV.ID"], Majiq_SJ_table["Gene.Name"], Majiq_SJ_table["Gene.ID"], Majiq_SJ_table["Chromosome"], Majiq_SJ_table["strand"], Majiq_SJ_table["denovo.x"], Majiq_SJ_table["meanPSI.x"], Majiq_SJ_table["Junctions.coords"], exon_coords, Num_junctions, Num_exons, Majiq_SJ_table["IR.coords"], 
                          Majiq_SJ_table["Num_unique_reads_T1.1"], Majiq_SJ_table["Num_unique_reads_T1.3"], 
                          Majiq_SJ_table["Num_unique_reads_T1.4"], Majiq_SJ_table["Num_unique_reads_T1.5"],
                          Majiq_SJ_table["Num_unique_reads_T1.6"],  
                          Majiq_SJ_table["T1_sum_1"], Majiq_SJ_table["T1_sum_3"], 
                          Majiq_SJ_table["T1_sum_4"], Majiq_SJ_table["T1_sum_5"], 
                          Majiq_SJ_table["T1_sum_6"], 
                          Majiq_SJ_table["T1_PSI_1"], Majiq_SJ_table["T1_PSI_3"], 
                          Majiq_SJ_table["T1_PSI_4"], Majiq_SJ_table["T1_PSI_5"], 
                          Majiq_SJ_table["T1_PSI_6"], 
                          Majiq_SJ_table["delta_sign"])
  
  colnames(df_select) <- c("LSV.ID", "Gene.Name", "Gene.ID", "Chromosome", "strand", "denovo", "meanPSI", "Junctions.coords","Exons.coords","Num.junctions", "Num.exons", "IR.coords", 
                           "Num_unique_reads_T1.1", "Num_unique_reads_T1.3", 
                           "Num_unique_reads_T1.4", "Num_unique_reads_T1.5",
                           "Num_unique_reads_T1.6", 
                           "T1_sum_1", "T1_sum_3", 
                           "T1_sum_4", "T1_sum_5", 
                           "T1_sum_6", 
                           "T1_PSI_1", "T1_PSI_3", "T1_PSI_4", "T1_PSI_5", "T1_PSI_6", "delta_sign")
  return(df_select)
  
}
Majiq_in_SJ_exons <- apply(Majiq_in_SJ_subset, 1, function(x)exon_filter(x))

###Full table with Majiq and SJ.out information after all the filters
Majiq_in_SJ_exons <-bind_rows(Majiq_in_SJ_exons) 
Majiq_in_SJ_exons[Majiq_in_SJ_exons == ""] <- NA
Majiq_in_SJ_exons<- subset(Majiq_in_SJ_exons, !is.na(Majiq_in_SJ_exons["Exons.coords"]))

one_tissue_LSV<- select(one_tissue, LSV.ID, UCSC.LSV.Link)
Majiq_in_SJ_exons<- left_join(Majiq_in_SJ_exons, one_tissue_LSV)
##CHANGE
write.csv(Majiq_in_SJ_exons, file = paste0(sample_name, '_PSI_Majiq_SJ.csv'))
