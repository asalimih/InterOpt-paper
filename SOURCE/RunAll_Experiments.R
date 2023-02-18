setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))

source('SOURCE/InterOpt.R')
source('SOURCE/Load_All_Data.R')
source('SOURCE/run_expers_functions.R')

#------Self Weights (Weights from Raw CT)----- 

weight_methods = c('arith','geom', 'random','arith_cv','geom_cv','arith_sd',
				   'geom_sd','geom_sd_soft','geom_sd_hybrid','sd_simple')

res_selfweight = list()
source('SOURCE/InterOpt.R')
res_selfweight[['GSE78870']] = run_exper_pcr_selfweight(data_source = pcr_data_clean,
														gr_source = grs_pcr,
														k = 2,
														dir = 'RESULTS/RawResults/SelfWeight/GSE78870/', 
														mc.cores=12,
														weight_methods=weight_methods)


res_selfweight[['GSE50013']] = run_exper_pcr_selfweight(data_source = data_13_imputed,
														gr_source = gr_13,
														k = 2,
														dir = 'RESULTS/RawResults/SelfWeight/GSE50013/',
														mc.cores=12,
														weight_methods=weight_methods)

res_selfweight[['GSE57661']] = run_exper_pcr_selfweight(data_source = data_57_imputed,
														gr_source = gr_57,
														k = 2,
														dir = 'RESULTS/RawResults/SelfWeight/GSE57661/',
														mc.cores=12,
														weight_methods=weight_methods)

res_selfweight[['GSE59520']] = run_exper_pcr_selfweight(data_source = data_59_imputed,
														gr_source = gr_59,
														k = 2,
														dir = 'RESULTS/RawResults/SelfWeight/GSE59520/',
														mc.cores=12,
														weight_methods=weight_methods)

#------Sample Size Analysis-----

# prepare data
weight_methods = c('arith','geom', 'random','arith_cv','geom_cv','arith_sd','geom_sd','geom_sd_soft','geom_sd_hybrid','sd_simple')

res_selfweight_ss = list()
source('SOURCE/InterOpt.R')
res_selfweight_ss[['GSE78870']] = run_exper_pcr_sample_num(data_source = pcr_data_clean,
														   gr_source = grs_pcr,
														   k = 2,
														   small_subset_n = c(10,20,30,40,50),
														   dir = 'RESULTS/RawResults/SelfWeight_samplenum_rep20/GSE78870/', 
														   mc.cores=12,
														   repeat_n = 20,
														   weight_methods=weight_methods)


res_selfweight_ss[['GSE50013']] = run_exper_pcr_sample_num(data_source = data_13_imputed,
														   gr_source = gr_13,
														   k = 2,
														   small_subset_n = c(10,20,30),
														   dir = 'RESULTS/RawResults/SelfWeight_samplenum_rep20/GSE50013/',
														   mc.cores=12,
														   repeat_n = 20,
														   weight_methods=weight_methods)

res_selfweight_ss[['GSE57661']] = run_exper_pcr_sample_num(data_source = data_57_imputed,
														   gr_source = gr_57,
														   k = 2,
														   small_subset_n = c(10,20,30),
														   dir = 'RESULTS/RawResults/SelfWeight_samplenum_rep20/GSE57661/',
														   mc.cores=12,
														   repeat_n = 20,
														   weight_methods=weight_methods)

res_selfweight_ss[['GSE59520']] = run_exper_pcr_sample_num(data_source = data_59_imputed,
														   gr_source = gr_59,
														   k = 2,
														   small_subset_n = c(10,20,30),
														   dir = 'RESULTS/RawResults/SelfWeight_samplenum_rep20/GSE59520/',
														   mc.cores=12,
														   repeat_n = 20,
														   weight_methods=weight_methods)

#------Generalization (Weights from external dataset)----- 

# prepare data
data_tcga_brca_filtered = filterInputData(data_lbl_cpm_filtered[, grs_lbl=='T' & subtype_raw_data_df$ER=='Positive'], ctVal=F)
grs_lbl_filtered = grs_lbl[grs_lbl=='T' & subtype_raw_data_df$ER=='Positive']

weight_methods = c('arith','geom', 'random','arith_cv','geom_cv','arith_sd','geom_sd','geom_sd_soft','geom_sd_hybrid','sd_simple')

res_gener = list()
source('SOURCE/InterOpt.R')
source('SOURCE/run_expers_functions.R')
res_gener[['TCGA_GSE78870']] = run_exper_generalization(data_source = data_tcga_brca_filtered,
														data_target = pcr_data_clean,
														gr_source = grs_lbl_filtered,
														gr_target = grs_pcr,
														k = 2,
														small_subset_n = c(20, 30),
														repeat_n = 20,
														dir = 'RESULTS/RawResults/Generalization/TCGA_GSE78870/', 
														mc.cores=13,
														weight_methods=weight_methods)


#------Iterative----- (Number of reference genes analysis)

# prepare data
weight_methods = c('geom_sd', 'geom')
res_selfweight_iter = list()


source('SOURCE/InterOpt.R')
res_selfweight_iter[['GSE50013']] = run_exper_pcr_selfweight_iter(data_source = data_13_imputed,
																  gr_source = gr_13,
																  k = 30,
																  dir = 'RESULTS/RawResults/SelfWeight_iter/GSE50013/',
																  keep = 400,
																  mc.cores=10,
																  weight_methods = weight_methods)


res_selfweight_iter[['GSE78870']] = run_exper_pcr_selfweight_iter(data_source = pcr_data_clean,
																  gr_source = grs_pcr,
																  k = 30,
																  keep = 400,
																  dir = 'RESULTS/RawResults/SelfWeight_iter/GSE78870/', 
																  mc.cores=10,
																  weight_methods = weight_methods)

res_selfweight_iter[['GSE57661']] = run_exper_pcr_selfweight_iter(data_source = data_57_imputed,
																  gr_source = gr_57,
																  k = 20,
																  keep = 400,
																  dir = 'RESULTS/RawResults/SelfWeight_iter/GSE57661/', 
																  mc.cores=12,
																  weight_methods = weight_methods)

res_selfweight_iter[['GSE59520']] = run_exper_pcr_selfweight_iter(data_source = data_59_imputed,
																  gr_source = gr_59,
																  k = 20,
																  keep = 400,
																  dir = 'RESULTS/RawResults/SelfWeight_iter/GSE59520/', 
																  mc.cores=12,
																  weight_methods = weight_methods)

#------Experimental-------
source('SOURCE/Experimental_Validation.R')
# 							   pvals      lfcs
# U48                    0.062592743 3.1605853
# hsa_miR_361_5p         0.221838036 1.8694204
# hsa_miR_16_5p          0.055941091 3.2644685
# hsa_miR_21_5p          0.046342233 4.2342147

# hsa_miR_21_5p_wall     0.001892650 2.4688539 <
# hsa_miR_21_5p_all      0.052102938 1.4693899 <