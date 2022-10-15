data <- as.data.frame(readRDS("DATA/TCGA_BRCA/NEW_TCGABRCA_miR.rds"))
rownames(data) <- miRNA_AccessionToName(rownames(data),targetVersion = "v21")$TargetName
clinical <- readRDS("DATA/TCGA_BRCA/NEW_TCGABRCA_clinical.rds")

##extract labels er pr her2
subtype_raw_data <- matrix(0,nrow=dim(data)[2],6)
subtype_raw_data_df <- as.data.frame(subtype_raw_data)
subtype_raw_data_df[,1] <- clinical$breast_carcinoma_estrogen_receptor_status
subtype_raw_data_df[,2] <- clinical$breast_carcinoma_progesterone_receptor_status
subtype_raw_data_df[,3] <- clinical$lab_proc_her2_neu_immunohistochemistry_receptor_status
subtype_raw_data_df[,4] <- clinical$er_level_cell_percentage_category
subtype_raw_data_df[,5] <- clinical$progesterone_receptor_level_cell_percent_category
subtype_raw_data_df[,6] <- clinical$her2_erbb_pos_finding_cell_percent_category
colnames(subtype_raw_data_df) <- c("ER","PR","Her2","ERp","PRp","Her2p")
#enhancing labels based on percentage category
subtype_raw_data_df[subtype_raw_data_df$ERp=="<10%" &
                      subtype_raw_data_df$ER=="Positive","ER"] <- "Negative"
subtype_raw_data_df[subtype_raw_data_df$PRp=="<10%" &
                      subtype_raw_data_df$PR=="Positive","PR"] <- "Negative"
subtype_raw_data_df[subtype_raw_data_df$Her2p=="<10%" &
                      subtype_raw_data_df$Her2=="Positive","Her2"] <- "Negative"

#index of nomral samples and subtypes
samples <- colnames(data)
sampletypes <- gsub("(.*)-(.*)-(.*)-","",samples) #based on TCGA category

index_norm <- which(sampletypes=="11A" | sampletypes=="11B")
index_tumor <- setdiff(1:dim(data)[2],index_norm)

index_tumor_type <- list()
index_tumor_type[["LA"]] <- as.numeric(rownames(subset(subtype_raw_data_df[index_tumor,],
                                                       ER=="Positive" &
                                                         PR=="Positive" & PRp!="10-19%" &
                                                         Her2=="Negative" )))
index_tumor_type[["LBpos"]] <- as.numeric(rownames(subset(subtype_raw_data_df[index_tumor,],
                                                          ER=="Positive" &
                                                            Her2=="Positive" )))
index_tumor_type[["LBneg"]] <- as.numeric(rownames(subset(subtype_raw_data_df[index_tumor,],
                                                          ER=="Positive" &
                                                            (PR=="Negative" | PRp=="10-19%") &
                                                            Her2=="Negative" )))
index_tumor_type[["Her2over"]] <- as.numeric(rownames(subset(subtype_raw_data_df[index_tumor,],
                                                             ER=="Negative" &
                                                               PR=="Negative" &
                                                               Her2=="Positive" )))
index_tumor_type[["TNBC"]] <- as.numeric(rownames(subset(subtype_raw_data_df[index_tumor,],
                                                         ER=="Negative" &
                                                           PR=="Negative" &
                                                           Her2=="Negative" )))
index_tumor_type[["ERpos"]] <- as.numeric(rownames(subset(subtype_raw_data_df[index_tumor,],
                                                          ER=="Positive")))

# Define group samples
grs = rep("X",dim(data)[2])
grs[index_tumor] <- rep("T",length(index_tumor))
grs[index_norm] <- rep("N",length(index_norm))

#filtering lowly expressed MicroRNAs
cpm_thre <- 2.6
data_cpm <- cpm(data)
data_filtered <- data[which(rowSums(data_cpm>cpm_thre)>=10),]
data_cpm_filtered <- data_cpm[which(rowSums(data_cpm>cpm_thre)>=10),]

#removing samples without subtype label
data_lbl <- data[,which(grs!="X")]
grs_lbl <- grs[which(grs!="X")]
data_lbl_cpm <- data_cpm[,which(grs!="X")]
data_lbl_cpm_filtered <- data_cpm_filtered[,which(grs!="X")]

# retaining well expressed MicroRNAs
cpmCut2count <- function(d,thr) {
  return(apply(d,2,function(col) max(col[which(cpm(col)<thr)])))
}

CountCutter <- cpmCut2count(data,5.916)
retainInd <-  which(apply(data_filtered,1,function(row) sum(row>CountCutter)>(0.95*ncol(data_filtered))))
retainNames <- rownames(data_filtered)[retainInd]

data_lbl_cpm_filtered = data_lbl_cpm_filtered[retainInd, ]

