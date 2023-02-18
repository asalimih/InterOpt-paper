library(miRBaseConverter)
library(edgeR)
library(GEOquery)

#setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))

#---------TCGA (BRCA miRNA-Seq Data)---------
source("SOURCE/TCGA_DataPrep.R")

# exclude: at least one more than 39 and all more than 35 (these are not neither considered in normalization nor as RG)
get_exclude_idx <- function(dd) {
	excl <- which(rowSums(dd > 35) == ncol(dd) &
				  	rowSums(dd > 39) > 1)
	return(excl)
}
filterInputData = function(data, ctVal, logarithm=F){
	if(ctVal){
		exc = (rowSums(data<0)>0)
	}else{
		if(logarithm){
			exc = (rowSums(data<0)>0)
		} else{
			exc = (rowSums(data<1)>0)
		}
	}
	if(any(exc)){
		message(paste0(sum(exc),' rows from data were removed.'))
		return(data[-which(exc),])
	}else{
		return(data)
	}
}

#----------GSE78870 (BRCA qPCR miRNA)---------
# as the number of genes with all samples detected are high enough no imputation done here
pcr_ct <- read.table("DATA/GSE78870/GSE78870_non-normalized_data.txt", sep = "\t", header = T, stringsAsFactors = F)
pcr_mirs <- gsub("(.+)-(.+)", "\\1", pcr_ct$Detector)
pcr_data <- pcr_ct[, 7:ncol(pcr_ct)]

pcr_data[pcr_data == "Undetermined" | pcr_data == ""] <- NA
pcr_data_clean <- as.matrix(pcr_data)
class(pcr_data_clean) <- "numeric"
pcr_data_clean[is.na(pcr_data_clean)] <- 40
exc <- get_exclude_idx(pcr_data_clean)
pcr_data_clean <- pcr_data_clean[-exc, ]
pcr_mirs_clean <- pcr_mirs[-exc]
pcr_data_clean <- aggregate(pcr_data_clean, list(pcr_mirs_clean), mean) # aggregate duplicates by mean
rownames(pcr_data_clean) <- pcr_data_clean[, 1]
pcr_mirs_clean <- pcr_data_clean[, 1]
pcr_data_clean <- as.matrix(pcr_data_clean[, -1])
a <- miRNAVersionConvert(rownames(pcr_data_clean), targetVersion = "v21")$TargetName
keep_prev_name <- c(grep("&", a), which(is.na(a)))
a[keep_prev_name] <- rownames(pcr_data_clean)[keep_prev_name]
pcr_data_clean <- pcr_data_clean[!is.na(a), ]
rownames(pcr_data_clean) <- a[!is.na(a)]
grs_pcr <- rep("T", ncol(pcr_data_clean))

#---------GSE50013 (LIHC qPCR)---------
# rows with less than 5 undetected values are imputed
data_13_imputed <- as.matrix(read.table("DATA/GSE50013/GSE50013_non_normalized_4miss_imputed.txt", stringsAsFactors = F))
data_13_imputed <- data_13_imputed[-get_exclude_idx(data_13_imputed), ]
gr_13 <- factor(rep(c("T", "N"), 20))

#------GSE57661 (early-stage breast cancer qPCR)--------
data_57_imputed <- as.matrix(read.table("DATA/GSE57661/GSE57661_Non-normalized_data_4miss_imputed.txt", stringsAsFactors = F))
excll = get_exclude_idx(data_57_imputed)
if(length(excll)>0)
	data_57_imputed <- data_57_imputed[-excll, ]
gr_57 <- factor(gsub('(.+)\\..*','\\1',colnames(data_57_imputed)))

#------GSE59520 (embryonic tumors of testis)--------
data_59_imputed <- as.matrix(read.table("DATA/GSE59520/GSE59520_non-normalized_4miss_imputed.txt", stringsAsFactors = F))
excll = get_exclude_idx(data_59_imputed)
if(length(excll)>0)
	data_59_imputed <- data_59_imputed[-excll, ]
data_59_imputed[data_59_imputed > 40] <- 40
gr_59 <- factor(gsub('\\.','',gsub('[[:digit:]]+','',gsub('_SERUM','',colnames(data_59_imputed)))))


data_tcga_brca_filtered = filterInputData(data_lbl_cpm_filtered[, grs_lbl=='T' & subtype_raw_data_df$ER=='Positive'], ctVal=F)
grs_lbl_filtered = grs_lbl[grs_lbl=='T' & subtype_raw_data_df$ER=='Positive']

