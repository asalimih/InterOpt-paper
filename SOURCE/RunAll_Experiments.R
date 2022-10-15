source('SOURCE/InterOpt.R')
source('SOURCE/Load_All_Data.R')
source('SOURCE/run_expers_functions.R') #doesn't work on the server

#------Self Weights (Weights from Raw CT)----- 

# prepare data
data_13_filtered = data_13_imputed[rowSums(data_13_imputed>39)==0,]
pcr_data_clean_filtered = pcr_data_clean[rowSums(pcr_data_clean>39)==0,]
weight_methods = c('arith','geom', 'random','arith_cv','geom_cv','arith_sd','geom_sd','geom_sd_soft','geom_sd_hybrid','sd_simple')

res_selfweight = list()
source('SOURCE/InterOpt.R')
res_selfweight[['GSE78870']] = run_exper_pcr_selfweight(data_source = pcr_data_clean_filtered,
                                                        gr_source = grs_pcr,
                                                        k = 2,
                                                        dir = 'RESULTS/SelfWeight_3/GSE78870/', 
                                                        mc.cores=8,
                                                        weight_methods=weight_methods)


res_selfweight[['GSE50013']] = run_exper_pcr_selfweight(data_source = data_13_filtered,
                                                        gr_source = gr_13,
                                                        k = 2,
                                                        dir = 'RESULTS/SelfWeight_3/GSE50013/',
                                                        mc.cores=8,
                                                        weight_methods=weight_methods)

#------Sample Size Analysis-----

# prepare data
data_13_filtered = data_13_imputed[rowSums(data_13_imputed>39)==0,]
pcr_data_clean_filtered = pcr_data_clean[rowSums(pcr_data_clean>39)==0,]
weight_methods = c('arith','geom', 'random','arith_cv','geom_cv','arith_sd','geom_sd','geom_sd_soft','geom_sd_hybrid','sd_simple')

res_selfweight_ss = list()
source('SOURCE/InterOpt.R')
res_selfweight_ss[['GSE78870']] = run_exper_pcr_sample_num(data_source = pcr_data_clean_filtered,
                                                        gr_source = grs_pcr,
                                                        k = 2,
                                                        small_subset_n = c(10,20,30,40,50),
                                                        dir = 'RESULTS/SelfWeight_3_samplenu_rep20/GSE78870/', 
                                                        mc.cores=8,
														repeat_n = 20,
                                                        weight_methods=weight_methods)


res_selfweight_ss[['GSE50013']] = run_exper_pcr_sample_num(data_source = data_13_filtered,
                                                        gr_source = gr_13,
                                                        k = 2,
                                                        small_subset_n = c(10,20,30),
                                                        dir = 'RESULTS/SelfWeight_3_samplenum_rep20/GSE50013/',
                                                        mc.cores=8,
														repeat_n = 20,
                                                        weight_methods=weight_methods)




#------Generalization (Weights from external dataset)----- 

# prepare data
pcr_data_clean_filtered = pcr_data_clean[rowSums(pcr_data_clean>39)==0, ]
data_tcga_brca_filtered = filterInputData(data_lbl_cpm_filtered[, grs_lbl=='T' & subtype_raw_data_df$ER=='Positive'], ctVal=F)
grs_lbl_filtered = grs_lbl[grs_lbl=='T' & subtype_raw_data_df$ER=='Positive']

weight_methods = c('arith','geom', 'random','arith_cv','geom_cv','arith_sd','geom_sd','geom_sd_soft','geom_sd_hybrid','sd_simple')

res_gener = list()
source('SOURCE/InterOpt.R')

res_gener[['TCGA_GSE78870']] = run_exper_generalization(data_source = data_tcga_brca_filtered,
                                                        data_target = pcr_data_clean_filtered,
                                                        gr_source = grs_lbl_filtered,
                                                        gr_target = grs_pcr,
                                                        k = 2,
                                                        small_subset_n = c(10, 15, 20, 30, 40),
                                                        repeat_n = 20,
                                                        dir = 'RESULTS/Generalization_3/TCGA_GSE78870/', 
                                                        mc.cores=8,
                                                        weight_methods=weight_methods)

#------Iterative----- (Number of reference genes analysis)

# prepare data
data_13_filtered = data_13_imputed[rowSums(data_13_imputed>39)==0,]
pcr_data_clean_filtered = pcr_data_clean[rowSums(pcr_data_clean>39)==0,]

weight_methods = c('arith','geom', 'random','arith_cv','geom_cv','arith_sd','geom_sd','sd_simple')

res_selfweight_iter = list()


source('SOURCE/InterOpt.R')
res_selfweight_iter[['GSE50013']] = run_exper_pcr_selfweight_iter(data_source = data_13_filtered,
                                                                  gr_source = gr_13,
                                                                  k = 30,
                                                                  dir = 'RESULTS/SelfWeight_iter/GSE50013/',
                                                                  keep = 400,
                                                                  mc.cores=10,
                                                                  weight_methods = weight_methods)


res_selfweight_iter[['GSE78870']] = run_exper_pcr_selfweight_iter(data_source = pcr_data_clean_filtered,
                                                                  gr_source = grs_pcr,
                                                                  k = 30,
                                                                  keep = 400,
                                                                  dir = 'RESULTS/SelfWeight_iter/GSE78870/', 
                                                                  mc.cores=10,
                                                                  weight_methods = weight_methods)

res_k1 = list()
source('SOURCE/InterOpt.R')
res_k1[['GSE78870']] = run_exper_pcr_selfweight(data_source = pcr_data_clean_filtered,
                                                gr_source = grs_pcr,
                                                k = 1,
                                                dir = 'RESULTS/Self_k1/GSE78870/', 
                                                mc.cores=8,
                                                weight_methods = weight_methods)

res_k1[['GSE50013']] = run_exper_pcr_selfweight(data_source = data_13_filtered,
                                                gr_source = gr_13,
                                                k = 1,
                                                dir = 'RESULTS/Self_k1/GSE50013/',
                                                mc.cores=8,
                                                weight_methods = weight_methods)


#------Experimential-------
source('SOURCE/Experimental_Validation.R')
# 							   pvals      lfcs
# U48                    0.062592743 3.1605853
# hsa_miR_361_5p         0.221838036 1.8694204
# hsa_miR_16_5p          0.055941091 3.2644685
# hsa_miR_21_5p          0.046342233 4.2342147

# hsa_miR_21_5p_wall     0.001892650 2.4688539 <
# hsa_miR_21_5p_all      0.052102938 1.4693899 <