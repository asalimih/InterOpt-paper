curD <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dirname(curD))

library('xlsx')
library('matrixStats')
source('SOURCE/summerize_experiment.R')
source('SOURCE/pub_plot_functions.R')
saveFlag = F


##############
###  DATA  ###
##############

#----color set----
mycols =  c(ggsci::pal_npg("nrc")(10),
		"gray50","gold2","chocolate3","brown1","violetred3","gray20",
		"olivedrab3","firebrick4",ggsci::pal_nejm()(8)[2],'goldenrod1')
plot(1:length(mycols), 1:length(mycols), col=mycols, pch=19, cex=5, xlab="", ylab="")
text(1:length(mycols), 1:length(mycols),labels=c(1:length(mycols)))


#----SelfWeight (Weights from Raw CT) ----

data_ids = c('GSE78870', 'GSE50013', 'GSE57661', 'GSE59520')
measures = c('SD', 'CV', 'Genorm', 'NormFinder')
wmethod_names =c('arith','geom', 'geom(rand)','arith(cv)','geom(cv)','arith(sd)',
				 'geom(sd)','geom(sd_soft)','geom(sd+)','geom(sd_r)')
new_order = c(1:6,10,9) #used for reordering and removing unwanted wmethods
wmethod_names = wmethod_names[new_order]
s_scenarios = c('raw', 'norm') #Don't Change the Order of this
group_names = c('','Raw CT','Norm. CT')
Base = c(rep(group_names[1],3), rep(group_names[2:3], each=length(wmethod_names)-3))

res_selfweight = list()
for(data_id in data_ids){
	res_selfweight[[data_id]] = read_exper(file.path('RESULTS/RawResults/SelfWeight',data_id), s_scenarios)
	res_selfweight[[data_id]] = process_experiment(res_selfweight[[data_id]], s_scenarios, measures, wmethod_names,
												   new_order, group_names, Base, precision=3)
}

# save
if(saveFlag)
	save_summary_table(res_selfweight, 'RESULTS/InterOpt_Publication/summary_selfweight.xlsx', data_ids)


#----Sample Size Analysis----

data_ids = c('GSE78870', 'GSE50013', 'GSE57661', 'GSE59520')
measures = c('SD', 'CV', 'Genorm', 'NormFinder')
wmethod_names =c('arith','geom', 'geom(rand)','arith(cv)','geom(cv)','arith(sd)',
				 'geom(sd)','geom(sd_soft)','geom(sd+)','geom(sd_r)')
new_order = c(2,4:6,10,7,9) #used for reordering and removing unwanted wmethods
wmethod_names = wmethod_names[new_order]
ss_scenarios = c('raw_small10','raw_small20','raw_small30',
				'raw_small40','raw_small50','raw') #Don't Change the Order of this
group_names = c('No Optim.', ss_scenarios)
Base = c(rep(group_names[1],3), rep(group_names[2:6], each=length(wmethod_names)-3)) # this variable is not actually used in the plot functions

res_samples = list()
for(data_id in data_ids){
	res_samples[[data_id]] = read_exper(file.path('RESULTS/RawResults/SelfWeight_samplenum_rep20',data_id), ss_scenarios)
	res_samples[[data_id]] = process_experiment(res_samples[[data_id]], ss_scenarios, measures, wmethod_names,
												   new_order, group_names, Base, precision=3, keep_no_optims = T, rep_sds=T)
}

# save
if(saveFlag)
	save_summary_table(res_samples, 'RESULTS/InterOpt_Publication/summary_sample_analysis.xlsx', data_ids)

#----Generalizations (Weights from external dataset)----

data_ids = c('TCGA_GSE78870')
measures = c('SD', 'CV', 'Genorm', 'NormFinder')
wmethod_names =c('arith','geom', 'geom(rand)','arith(cv)','geom(cv)','arith(sd)',
				 'geom(sd)','geom(sd_soft)','geom(sd+)','geom(sd_r)')
new_order = c(1:6,10,9) #used for reordering and removing unwanted wmethods
wmethod_names = wmethod_names[new_order]
g_scenarios = c('raw_small20', 'rnaseq', 'raw_small30','raw', 'norm') #Don't Change the Order of this
group_names = c('','Raw n:20', 'RNASeq',  'Raw n:30','Raw CT','Norm. CT')
Base = c(rep(group_names[1],3), rep(group_names[2:length(group_names)], each=length(wmethod_names)-3))

res_gener = list()
for(data_id in data_ids){
	# load
	res_gener[[data_id]] = read_exper(file.path('RESULTS/RawResults/Generalization',data_id), g_scenarios)
	# we put target result in source result because next functions read results from just res_source
	res_gener[[data_id]]$expers$rnaseq$res_source = res_gener[[data_id]]$expers$rnaseq$res_target
	res_gener[[data_id]] = process_experiment(res_gener[[data_id]], g_scenarios, measures, wmethod_names,
											  new_order, group_names, Base, precision=3)
}

if(saveFlag)
	save_summary_table(res_gener, 'RESULTS/InterOpt_Publication/summary_gener.xlsx', data_ids)




#############
### PLOTS ###
#############


#----Aggregated Plots----
data_pair_infos = list()
data_pair_infos[[1]] = list(data_ids = c('GSE78870', 'GSE50013'),
						   data_labels = c('Breast Cancer (Tissue)', 'Liver Cancer (Plasma)'),
						   data_smp_nums = list(c(10,20,30,40,50,106), c(10,20,30,40)))

data_pair_infos[[3]] = list(data_ids = c('GSE57661', 'GSE59520'),
						   data_labels = c('Early Stage Breast Cancer (Plasma)', 'Embryonic Tumors of Testis (Plasma)'),
						   data_smp_nums = list( c(10,20,30,48), c(10,20,30,36)))

measures = c('SD', 'CV', 'Genorm', 'NormFinder', 'Stability')
selfweight_full = list()
sample_plot_flag = T
iter_plot_flag = T
saveFlag = T

for(data_pair_info in data_pair_infos){
	cat('Data Pair: ', data_pair_info[['data_ids']][1], ', ', data_pair_info[['data_ids']][2], '\n')
	for(MainMeasure in measures) {
		#MainMeasure = 'Stability'
		cat('  Generating ', MainMeasure, ' plots ...\n')
		
		# following R code fills the selfweight_full list
		# it's inputs are: MainMeasure, data_pair_info, sample_plot_flag, iter_plot_flag
		source('SOURCE/pub_plot_generate_panel.R')
		
		#---save plots---
		if(saveFlag){
			ggsave(sprintf("RESULTS/InterOpt_Publication/selfweight_full_%s_%s_%s.pdf",
						   data_pair_info[['data_ids']][1], data_pair_info[['data_ids']][2],MainMeasure),
				   selfweight_full[[MainMeasure]] ,
				   width = new_w, height = 27, units = "cm")
		}
	} # end of measure loop for aggregated plots 
}

#----- Generalization Plots----
gener_full = list()
measures = c('SD', 'CV', 'Genorm', 'NormFinder', 'Stability')
for(MainMeasure in measures) {
	cls_gener = mycols[c(14,15,20,2,4,3)]
	names(cls_gener) <- res_gener[[1]]$group_names
	wmethod_filter_g = res_gener[[1]]$wmethod_names[c(-1,-3)] #removes arith and geom(rand)
	wmethod_filter_g_heat = res_gener[[1]]$wmethod_names[c(2,7,8)] #removes arith and geom(rand)
	
	g1_bar = plot.measure.bar(res_gener[['TCGA_GSE78870']], MainMeasure, cls_gener, wmethod_filter_g, tip_num_size=3)+
		theme(legend.position="none", strip.text.y = element_text(size = 10))
	g1_box = plot.measure.box(res_gener[['TCGA_GSE78870']], MainMeasure, cls_gener, wmethod_filter_g, strip.text.size = 8)
	
	g1_heat1 = plot.pairwise_wilcox2(res_gener$TCGA_GSE78870, MainMeasure, wmethod_filter_g_heat, rm.norm.res = T)
	g1_heat2 = plot.pairwise_wilcox2(res_gener$TCGA_GSE78870, MainMeasure, rm.norm.res = T)
	
	gener_right_plots = plot_grid(g1_box, g1_heat1,
								  labels = c('(B)','(C)'),
								  rel_heights = c(5, 6),
								  label_size = 9, nrow = 2,hjust = 0)
	
	gener_full[[MainMeasure]] <- plot_grid(g1_bar, gener_right_plots,
										   labels = c('(A)',''),
										   rel_heights = c(9.5, 6, 6),
										   label_size = 9, nrow = 1,hjust = 0)
	
	if(saveFlag){
		ggsave(sprintf("RESULTS/InterOpt_Publication/gener_full_%s.pdf",MainMeasure),
			   gener_full[[MainMeasure]] ,
			   width = 25, height = 12, units = "cm",)
	}
} # end of measure loop for generalization



#------Experimential-------
rows = c('hsa_miR_21_5p_all', 'hsa_miR_21_5p_wall')
p_ghanbari = plot.box.grouped(data_full_m, gr_full_m, rows, F, 'miR-21-5p Expression', F, colors = mycols[c(14,3,4)])
p_ghanbari = p_ghanbari + theme(legend.position = "none")

if(saveFlag){
	ggsave("RESULTS/InterOpt_Publication/exper_ghanbari_21.pdf",
		   p_ghanbari,
		   width = 8, height = 6, units = "cm")
}
