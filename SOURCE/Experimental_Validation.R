#setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))

#------Functions------
library('CovTools')
source('SOURCE/weights_optim2.R')
# for geom_sd function input must be in log2 or ct mode
# for a combination of two it would be matrix with two rows

normalizer = function(d, i, refs, weighted = F){
	d = t(d)
	if(length(refs)>1){
		if(weighted)
			w = geom_sd(d[refs,])
		else
			w = rep(1/length(refs),length(refs))
		refx = colSums(w*d[refs,])
	}else{
		refx = d[refs,]
	}
	res = refx-d[i,]
	return(res)
}

dea_ttest = function(data, gr, gr1, gr2, paired, thetest=t.test){
	pvals = apply(data, 1, function(gene){
		a = thetest(gene[gr==gr1], gene[gr==gr2], paired = paired)
		return(a$p.value)
	})
	if(paired){
		lfcs = apply(data, 1, function(x) x[gr==gr1]-x[gr==gr2])
	}else{
		lfcs = apply(data, 1, function(x) mean(x[gr==gr1])-mean(x[gr==gr2]))
	}
	data.frame(pvals, lfcs)
}


#--------Ghanbari Data-------

paired = F
## Prepare data
library(xlsx)
ct = read.xlsx("DATA/RawCT.xlsx",sheetIndex = 1)
ct[which(!is.na(ct$SAMPLE))+1,1] = ct[which(!is.na(ct$SAMPLE)),1]
inc21 = which(!is.na(ct$Ct..miR21.))
ct_full = ct[inc21,]
colnames(ct_full) = c("sample","U48","hsa_miR_361_5p","hsa_miR_16_5p","hsa_miR_21_5p")
ct_full_m = aggregate(. ~ sample, data = ct_full, mean)
rownames(ct_full_m) = ct_full_m$sample
ct_full_m = ct_full_m[,-1]
rownames(ct_full) = make.names(ct_full$sample, unique = T)
ct_full = ct_full[,-1]

# filter outlier
excl = which(ct_full_m$hsa_miR_16_5p==40)
if(paired)
	excl = c(excl, excl-12) # to maintain just paired samples
ct_full_m = ct_full_m[-excl,]
gr_full_m = gsub("(.)(.*)","\\1",rownames(ct_full_m))


ct_full_m[["hsa_miR_21_5p_wall"]]     = normalizer(ct_full_m, 4, 1:3, T)
ct_full_m[["hsa_miR_21_5p_all"]]      = normalizer(ct_full_m, 4, 1:3)

data_full_m = t(ct_full_m) # used in plots

## Differential Expression Analysis with t test
dea_ttest(t(ct_full_m), gr_full_m, 'T', 'N', paired)

# 							   pvals      lfcs
# U48                    0.062592743 3.1605853
# hsa_miR_361_5p         0.221838036 1.8694204
# hsa_miR_16_5p          0.055941091 3.2644685
# hsa_miR_21_5p          0.046342233 4.2342147

# hsa_miR_21_5p_wall     0.001892650 2.4688539 <
# hsa_miR_21_5p_all      0.052102938 1.4693899 <

