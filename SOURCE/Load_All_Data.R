library(miRBaseConverter)
library(edgeR)
library(GEOquery)

#---------TCGA (BRCA miRNA-Seq Data)---------
source("SOURCE/TCGA_DataPrep.R")

#----------GSE78870 (BRCA qPCR miRNA)---------
pcr_ct = read.table('DATA/GSE78870/GSE78870_non-normalized_data.txt',sep = '\t', header = T, stringsAsFactors = F)
pcr_mirs = gsub("(.+)-(.+)", "\\1", pcr_ct$Detector)
pcr_data = pcr_ct[,7:ncol(pcr_ct)]

exc = which(rowSums(pcr_data=="Undetermined")>0 |
				rowSums(pcr_data=='')!=0)
pcr_data_clean = as.matrix(pcr_data[-exc,])
pcr_mirs_clean = pcr_mirs[-exc]
class(pcr_data_clean) = "numeric"
pcr_data_clean = aggregate(pcr_data_clean, list(pcr_mirs_clean), mean) #aggregate duplicates by mean
rownames(pcr_data_clean) = pcr_data_clean[,1]
pcr_mirs_clean = pcr_data_clean[,1]
pcr_data_clean = as.matrix(pcr_data_clean[,-1])
a <- miRNAVersionConvert(rownames(pcr_data_clean),targetVersion = "v21")$TargetName
keep_prev_name = c(grep('&',a), which(is.na(a)))
a[keep_prev_name] = rownames(pcr_data_clean)[keep_prev_name]
pcr_data_clean = pcr_data_clean[!is.na(a),]
rownames(pcr_data_clean) <- a[!is.na(a)]
grs_pcr = rep('T', ncol(pcr_data_clean))

#---------GSE50013 (LIHC qPCR)---------
readr::local_edition(1)
gset <- getGEO(filename = "DATA/GSE50013/GSE50013_series_matrix.txt",GSEMatrix = T,getGPL = T, destdir = 'DATA/GSE50013/')
gr_13 = factor(gset$`sample group:ch1`)
levels(gr_13) = c('N','T')
data_13 = read.delim('DATA/GSE50013/GSE50013_non_normalized.txt', stringsAsFactors = F)
data_13$ID_REF = sub('#','*',gsub('(.*)-.*','\\1',data_13$ID_REF))
rownames(data_13)  = data_13$ID_REF
data_13 = data_13[,-1]
a <- miRNAVersionConvert(rownames(data_13),targetVersion = "v21")$TargetName
keep_prev_name = c(grep('&',a), which(is.na(a)))#dont convert those having multiple targetnames
a[keep_prev_name] = rownames(data_13)[keep_prev_name]
data_13 = data_13[!is.na(a),]
rownames(data_13) <- a[!is.na(a)]
data_13 =  as.matrix(data_13)
data_13[data_13=="."] = "0"
class(data_13) = 'numeric'
data_13[is.na(data_13)] = 0 # just 32 miRNA don't have 0/NA! imputation is needed
data_13_imputed = data_13[rowSums(data_13==0)<5,]
avgs = rowMeans(data_13_imputed, na.rm = T)
have_zeros = which(rowSums(data_13_imputed==0)>0)
for(idx in have_zeros){
	data_13_imputed[idx, data_13_imputed[idx,]==0] = avgs[idx]
}

#---------GSE67075 (Colon qPCR)---------
data_67 = read.delim('DATA/GSE67075/GSE67075_Non-normalized_data.txt', stringsAsFactors = F)
data_67_norm = read.delim('DATA/GSE67075/GSE67075_Normalized_data.txt', stringsAsFactors = F)

data_67_imputed = data_67[rowSums(data_67=='Undetermined')<4,]
data_67_imputed[data_67_imputed=='Undetermined'] = NA
data_67_imputed[,2:ncol(data_67_imputed)] = sapply(data_67_imputed[,2:ncol(data_67_imputed)], as.numeric)
avgs = rowMeans(data_67_imputed[,2:ncol(data_67_imputed)], na.rm = T)
have_nas = which(rowSums(is.na(data_67_imputed[,2:ncol(data_67_imputed)]))>0)
for(idx in have_nas){
		data_67_imputed[idx, is.na(data_67_imputed[idx,])] = avgs[idx]
}
data_67_imputed = aggregate(data_67_imputed[2:ncol(data_67_imputed)], list(data_67_imputed$ID_REF), mean)
data_67_imputed[,1] = sub('#','*',gsub('(.*)-.*','\\1',data_67_imputed[,1]))
rownames(data_67_imputed) = data_67_imputed[,1]
data_67_imputed = data_67_imputed[,-1]

a <- miRNAVersionConvert(rownames(data_67_imputed),targetVersion = "v21")$TargetName
keep_prev_name = c(grep('&',a), which(is.na(a)))#dont convert those having multiple targetnames
a[keep_prev_name] = rownames(data_67_imputed)[keep_prev_name]
data_67_imputed = data_67_imputed[!is.na(a),]
rownames(data_67_imputed) <- a[!is.na(a)]
data_67_imputed =  as.matrix(data_67_imputed)
gr_67 = factor(gsub('SAMPLE\\.([^[:digit:]]+)[[:digit:]]+','\\1',colnames(data_67_imputed)))

