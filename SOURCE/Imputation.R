library('GEOquery')
library(miRBaseConverter)
library(HTqPCR)
library(nondetects)

curD <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dirname(curD))

impute_pcrarray = function(data_pcr, gr_pcr, max_miss=2, golden_ref_max_ct = 35){
	if(!is.matrix(data_pcr))
		stop('data_pcr must be a matrix')
	rows_to_impute = which(rowSums(data_pcr==40)<(max_miss+1))
	data_pcr_f = data_pcr[rows_to_impute, ]
	ref = apply(data_pcr,2, function(x) mean(x[x<golden_ref_max_ct]))
	data_pcr_f = rbind(data_pcr_f, ref) # add the golden reference to last row of data
	
	# convert data to a qPCRset class
	q.raw   <- new("qPCRset", exprs=data_pcr_f)
	featureNames(q.raw) <- rownames(data_pcr_f)
	sampleNames(q.raw) <- colnames(data_pcr_f)
	featureCategorytmp <-  as.data.frame(array("OK", dim=dim(data_pcr_f)), stringsAsFactors=F)
	featureCategorytmp[data_pcr_f==40] <- 'Undetermined'
	featureCategory(q.raw) <- featureCategorytmp
	exprs(q.raw) <- data_pcr_f
	pData(q.raw) <- data.frame(Group = gr_pcr, row.names = colnames(data_pcr))
	pData(q.raw)$sampleName <- colnames(data_pcr_f) 
	featureType(q.raw) <- c(rep("Target",nrow(data_pcr_f)-1),"Endogenous Control")
	
	# run nondedetect package imputation
	theimputed <- qpcrImpute(q.raw, groupVars="Group", iterMax=100) 

	data_pcr_imputed = data_pcr
	data_pcr_imputed[rows_to_impute,] = exprs(theimputed)[seq(length(rows_to_impute)),]
	return(data_pcr_imputed)
}

#-----GSE50013----
gr_13 = factor(rep(c('T','N'),20))
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

data_13[data_13=="."] = "40"
class(data_13) = 'numeric'
data_13[is.na(data_13)] = 40

data_13_2miss_imputed = impute_pcrarray(data_13, gr_13, max_miss=2)
data_13_4miss_imputed = impute_pcrarray(data_13, gr_13, max_miss=4)

write.table(data_13_2miss_imputed, 'DATA/GSE50013/GSE50013_non_normalized_2miss_imputed.txt', quote = T, sep = '\t')
write.table(data_13_4miss_imputed, 'DATA/GSE50013/GSE50013_non_normalized_4miss_imputed.txt', quote = T, sep = '\t')

#-----GSE57661----
data_57 = read.delim('DATA/GSE57661/GSE57661_Non-normalized_data.txt', stringsAsFactors = F)
data_57[data_57=='*'] = "40"
data_57 = data_57[,grep('Blank', colnames(data_57), invert = T)]
data_57 = data_57[grep('Blank', data_57$Sample, invert = T),]
data_57[,2:ncol(data_57)] = sapply(data_57[,2:ncol(data_57)], as.numeric)


a <- miRNAVersionConvert(data_57$Sample,targetVersion = "v21")$TargetName
keep_prev_name = c(grep('&',a), which(is.na(a)))
a[keep_prev_name] = data_57$Sample[keep_prev_name]

data_57 = aggregate(data_57[2:ncol(data_57)], list(a), mean)
rownames(data_57) = data_57[,1]
data_57 = as.matrix(data_57[,-1])
gr_57 = factor(gsub('(.+)\\..*','\\1',colnames(data_57)))

data_57_2miss_imputed = impute_pcrarray(data_57, gr_57, max_miss=2)
data_57_4miss_imputed = impute_pcrarray(data_57, gr_57, max_miss=4)

write.table(data_57, 'DATA/GSE57661/GSE57661_Non-normalized_data_not_imputed.txt', quote = T, sep = '\t')
write.table(data_57_2miss_imputed, 'DATA/GSE57661/GSE57661_Non-normalized_data_2miss_imputed.txt', quote = T, sep = '\t')
write.table(data_57_4miss_imputed, 'DATA/GSE57661/GSE57661_Non-normalized_data_4miss_imputed.txt', quote = T, sep = '\t')


#-----GSE59520----embryonic tumours of testis
data_59 = read.delim('DATA/GSE59520/GSE59520_non-normalized.txt', stringsAsFactors = F)
data_59 = data_59[data_59$Species=='Human', ]
data_59$ID_REF = data_59$miRNA_ID
data_59 = data_59[,-c(38:ncol(data_59))]
data_59[is.na(data_59)] = "40"
data_59[,2:ncol(data_59)] = sapply(data_59[,2:ncol(data_59)], as.numeric)

data_59$ID_REF = sub('#','*',data_59$ID_REF)

data_59 = aggregate(data_59[2:ncol(data_59)], list(data_59$ID_REF), mean)
rownames(data_59) = data_59[,1]
data_59 = as.matrix(data_59[,-1])
a <- miRNAVersionConvert(rownames(data_59),targetVersion = "v21")$TargetName
keep_prev_name = c(grep('&',a), which(is.na(a)))
a[keep_prev_name] = rownames(data_59)[keep_prev_name]
data_59 = data_59[!is.na(a),]
rownames(data_59) <- a[!is.na(a)]

gr_59 = factor(gsub('\\.','',gsub('[[:digit:]]+','',gsub('_SERUM','',colnames(data_59)))))

data_59_2miss_imputed = impute_pcrarray(data_59, gr_59, max_miss = 2)
data_59_4miss_imputed = impute_pcrarray(data_59, gr_59, max_miss = 4)

write.table(data_59, 'DATA/GSE59520/GSE59520_non-normalized_not_imputed.txt', quote = T, sep = '\t')
write.table(data_59_2miss_imputed, 'DATA/GSE59520/GSE59520_non-normalized_2miss_imputed.txt', quote = T, sep = '\t')
write.table(data_59_4miss_imputed, 'DATA/GSE59520/GSE59520_non-normalized_4miss_imputed.txt', quote = T, sep = '\t')




