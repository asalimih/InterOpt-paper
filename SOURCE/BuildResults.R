curD <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dirname(curD))

library('xlsx')
library('matrixStats')
source('SOURCE/summerize_experiment.R')
source('SOURCE/pub_plot_functions.R')
saveFlag = F

#----color set----
mycols =  c(ggsci::pal_npg("nrc")(10),
			"gray50","gold2","chocolate3","brown1","violetred3","gray20",
			"olivedrab3","firebrick4",ggsci::pal_nejm()(8)[2],'goldenrod1')
plot(1:length(mycols), 1:length(mycols), col=mycols, pch=19, cex=5, xlab="", ylab="")
text(1:length(mycols), 1:length(mycols),labels=c(1:length(mycols)))


#----SelfWeight (Weights from Raw CT) ----

data_ids = c('GSE78870', 'GSE50013')
measures = c('SD', 'CV', 'Genorm', 'NormFinder')#,'BestKeeper')
wmethod_names =c('arith','geom', 'geom(rand)','arith(cv)','geom(cv)','arith(sd)',
				 'geom(sd)','geom(sd_soft)','geom(sd+)','geom(sd_r)')
new_order = c(1:6,10,7) #used for reordering and removing unwanted wmethods
wmethod_names = wmethod_names[new_order]
s_scenarios = c('raw', 'norm') #Don't Change the Order of this
group_names = c('','Raw CT','Norm. CT')
cls_s = mycols[c(14,3,4)]
Base = c(rep(group_names[1],3), rep(group_names[2:3], each=5))

res_selfweight = list()
for(data_id in data_ids){
	res_selfweight[[data_id]] = read_exper(file.path('RESULTS/RawResults/SelfWeight_3',data_id), s_scenarios)
	res_selfweight[[data_id]] = process_experiment(res_selfweight[[data_id]], s_scenarios, measures, wmethod_names,
												   new_order, group_names, Base, precision=3)
}

# save
if(saveFlag)
	save_summary_table(res_selfweight, 'RESULTS/InterOpt_Publication/summary_selfweight.xlsx', data_ids)

	
#----Sample Size Analysis----

data_ids = c('GSE78870', 'GSE50013')
measures = c('SD', 'CV', 'Genorm', 'NormFinder')#,'BestKeeper')
wmethod_names =c('arith','geom', 'geom(rand)','arith(cv)','geom(cv)','arith(sd)',
				 'geom(sd)','geom(sd_soft)','geom(sd+)','geom(sd_r)')
new_order = c(2,4:6,10,7,9) #used for reordering and removing unwanted wmethods
wmethod_names = wmethod_names[new_order]
ss_scenarios = c('raw_small10','raw_small20','raw_small30',
				'raw_small40','raw_small50','raw') #Don't Change the Order of this
group_names = c('No Optim.', s_scenarios)
Base = c(group_names[1], rep(group_names[2:6], each=6))

res_samples = list()
for(data_id in data_ids){
	res_samples[[data_id]] = read_exper(file.path('RESULTS/RawResults/SelfWeight_3_samplenum',data_id), ss_scenarios)
	res_samples[[data_id]] = process_experiment(res_samples[[data_id]], s_scenarios, measures, wmethod_names,
												   new_order, group_names, Base, precision=3, keep_no_optims = T, rep_sds=T)
}

# save
if(saveFlag)
	save_summary_table(res_samples, 'RESULTS/InterOpt_Publication/summary_sample2.xlsx', data_ids)

#----Generalizations (Weights from external dataset)----

data_ids = c('TCGA_GSE78870')
measures = c('SD', 'CV', 'Genorm', 'NormFinder')#,'BestKeeper')
wmethod_names =c('arith','geom', 'geom(rand)','arith(cv)','geom(cv)','arith(sd)',
				 'geom(sd)','geom(sd_soft)','geom(sd+)','geom(sd_r)')
new_order = c(1:6,10,9) #used for reordering and removing unwanted wmethods
wmethod_names = wmethod_names[new_order]
g_scenarios = c('raw_small20', 'rnaseq', 'raw_small30','raw', 'norm') #Don't Change the Order of this
group_names = c('','Raw n:20', 'RNASeq',  'Raw n:30','Raw CT','Norm. CT')
cls_g = mycols[c(14,18,15,2,3,4)]
Base = c(rep(group_names[1],3), rep(group_names[2:length(group_names)], each=length(wmethod_names)-3))

res_gener = list()
for(data_id in data_ids){
	# load
	res_gener[[data_id]] = read_exper(file.path('RESULTS/RawResults/Generalization_3',data_id), g_scenarios)
	# we put target result in source result because next functions read results from just res_source
	res_gener[[data_id]]$expers$rnaseq$res_source = res_gener[[data_id]]$expers$rnaseq$res_target
	res_gener[[data_id]] = process_experiment(res_gener[[data_id]], g_scenarios, measures, wmethod_names,
											  new_order, group_names, Base, precision=3)
}

if(saveFlag)
	save_summary_table(res_gener, 'RESULTS/InterOpt_Publication/summary_gener.xlsx', data_ids)


#----Selfweight Plots----

wmethod_filter_s = res_selfweight$GSE78870$wmethod_names[c(-1,-3)] #removes arith and geom(rand)
measures = c('SD', 'CV', 'Genorm', 'NormFinder', 'Stability')
cls_selfweight = mycols[c(14,15,20,2,4,3)]
names(cls_selfweight) <- c('','Raw CT','Normalized CT')

# try
cls2 = mycols[c(14,15,20,2,4,3)]
cls2 = mycols[c(14,15,5,6,4,3)]
plot.measure.bar(res_selfweight$GSE50013,'Stability', cls2, wmethod_filter_s)
plot.measure.box(res_selfweight$GSE50013,'Stability', cls2, wmethod_filter_s)

if(saveFlag){
	the_dir = "RESULTS/InterOpt_Publication/measures"
	for(measure in measures){
		for(data_id in names(res_selfweight)){
			ggsave(sprintf("%s/selfweight_bar_%s_%s.pdf", the_dir, data_id, measure),
				   plot.measure.bar(res_selfweight[[data_id]], measure, cls_selfweight, wmethod_filter_s),
				   width = 13, height = 6, units = "cm")
			ggsave(sprintf("%s/selfweight_box_%s_%s.pdf", the_dir, data_id, measure),
				   plot.measure.box(res_selfweight[[data_id]], measure, cls_selfweight, wmethod_filter_s),
				   width = 10, height = 6, units = "cm")
			ggsave(sprintf("%s/selfweight_heat_%s_%s.pdf", the_dir, data_id, measure),
				   plot.pairwise_wilcox2(res_selfweight[[data_id]], measure, wmethod_filter_s),
				   width = 11, height = 7, units = "cm")
		}
	}
}

frg <- ggplot(data.frame()) + geom_point()+ theme(panel.background = element_rect(fill = 'white', colour = 'white'))

p1_bar = plot.measure.bar(res_selfweight[['GSE78870']], 'Stability', cls_selfweight, wmethod_filter_s)+ theme(legend.position="none")
p2_bar = plot.measure.bar(res_selfweight[['GSE50013']], 'Stability', cls_selfweight, wmethod_filter_s)
selfweight_bar <- plot_grid(frg, p1_bar+ggtitle('Breast Cancer (Tissue)'),
							frg, p2_bar+ggtitle('Liver Cancer (Plasma)'),
							label_size = 10, ncol = 4,hjust = 0.2 , rel_widths = c(0.02,0.4, 0.02, 0.58))

p1_box = plot.measure.box(res_selfweight[['GSE78870']], 'Stability', cls_selfweight, wmethod_filter_s)+ theme(legend.position="none")
p2_box = plot.measure.box(res_selfweight[['GSE50013']], 'Stability', cls_selfweight, wmethod_filter_s)
selfweight_box <- plot_grid(frg, p1_box, frg, frg, p2_box,
							label_size = 10, ncol = 5,hjust = 0.2 , rel_widths = c(0.06,0.31,0.02, 0.04, 0.48))

p1_heat = plot.pairwise_wilcox2(res_selfweight$GSE78870, 'Stability', wmethod_filter_s, rm.norm.res = T)+ guides(fill = FALSE)     
p2_heat = plot.pairwise_wilcox2(res_selfweight$GSE50013, 'Stability', wmethod_filter_s, rm.norm.res = T)
selfweight_heat <- plot_grid(frg, p1_heat, frg, p2_heat, frg,
							 label_size = 9, ncol = 5,hjust = 0.2 , rel_widths = c(0.14,0.26, 0.1, 0.38, 0.14))

selfweight_full <- plot_grid(selfweight_bar, selfweight_box, selfweight_heat, labels = c('(A)','(B)','(C)'),
							 rel_heights = c(0.39,0.35, 0.26),
							 label_size = 9, nrow = 3,hjust = 0)
if(saveFlag){
	ggsave("RESULTS/InterOpt_Publication/selfweight_full_v2.pdf",
		   selfweight_full,
		   width = 20, height = 14.2, units = "cm")
}

#----- Generalization Plots----
measures = c('SD', 'CV', 'Genorm', 'NormFinder', 'Stability')
#cls_gener =  mycols[c(14,3,7,9,4)]
cls_gener = mycols[c(14,15,20,2,4,3)]
names(cls_gener) <- res_gener[[1]]$group_names
wmethod_filter_g = res_gener[[1]]$wmethod_names[c(-1,-3)] #removes arith and geom(rand)
wmethod_filter_g_heat = res_gener[[1]]$wmethod_names[c(2,7,8)] #removes arith and geom(rand)
measures = c('SD', 'CV', 'Genorm', 'NormFinder', 'Stability')
the_dir = "RESULTS/InterOpt_Publication/measures"
if(saveFlag){
	for(measure in measures){
		for(data_id in names(res_gener)){
			ggsave(sprintf("%s/gener_bar_%s_%s.pdf", the_dir, data_id, measure),
				   plot.measure.bar(res_gener[[data_id]], measure, cls_gener, wmethod_filter_g),
				   width = 13, height = 6, units = "cm")
			ggsave(sprintf("%s/gener_box_%s_%s.pdf", the_dir, data_id, measure),
				   plot.measure.box(res_gener[[data_id]], measure, cls_gener, wmethod_filter_g),
				   width = 10, height = 6, units = "cm")
			ggsave(sprintf("%s/gener_heat_%s_%s.pdf", the_dir, data_id, measure),
				   plot.pairwise_wilcox2(res_gener[[data_id]], measure, wmethod_filter_g_heat),
				   width = 11, height = 7, units = "cm")
		}
	}
}

g1_bar = plot.measure.bar(res_gener[['TCGA_GSE78870']], 'Stability', cls_gener, wmethod_filter_g)+ theme(legend.position="none")
g1_box = plot.measure.box(res_gener[['TCGA_GSE78870']], 'Stability', cls_gener, wmethod_filter_g)
# g1_heat1 = plot.pairwise_wilcox2(res_gener$TCGA_GSE78870, 'Stability', wmethod_filter_g_heat)
# g1_heat2 = plot.pairwise_wilcox2(res_gener$TCGA_GSE78870, 'Stability')

g1_heat1 = plot.pairwise_wilcox2(res_gener$TCGA_GSE78870, 'Stability', wmethod_filter_g_heat, rm.norm.res = T)
g1_heat2 = plot.pairwise_wilcox2(res_gener$TCGA_GSE78870, 'Stability', rm.norm.res = T)

if(saveFlag){
	ggsave("RESULTS/InterOpt_Publication/gener_bar.pdf",
		   g1_bar,
		   width = 10, height = 9.5, units = "cm")
	
	ggsave("RESULTS/InterOpt_Publication/gener_box.pdf",
		   g1_box,
		   width = 13, height = 6, units = "cm")
	
	ggsave("RESULTS/InterOpt_Publication/gener_heat.pdf",
		   g1_heat1 + theme(axis.text.y = element_text(size=8),
		   				    axis.text.x = element_text(size=8)),
		   width = 12, height = 6, units = "cm")
	
	ggsave("RESULTS/InterOpt_Publication/gener_heat2.pdf",
		   g1_heat2,
		   width = 18, height = 12, units = "cm")
}
#-----Sample Analysis Plots----

wmethod_names =c('geom','arith(cv)','geom(cv)','arith(sd)','geom(sd_r)','geom(sd)','geom(sd_hybrid)')
cls = mycols[c(14,15,20,2,4,7,3)]
better_ylabel = ylab(expression(Mean~(SD~(RG[2]))))

data_smp_nums = c(10,20,30,40,50,106)
smp78870 = plot.sample.size(res_samples$GSE78870, data_smp_nums, 'SD',cls) +
	       ggtitle('Breast Cancer (Tissue)') + better_ylabel

data_smp_nums = c(10,20,30,40)
smp50013 = plot.sample.size(res_samples$GSE50013, data_smp_nums, 'SD',cls) +
		   ggtitle('Liver Cancer (Plasma)') + better_ylabel

if(saveFlag){
	ggsave("RESULTS/InterOpt_Publication/samples_78870.pdf", smp78870+theme(legend.position="none"), width=5, height=3.2) 
	ggsave("RESULTS/InterOpt_Publication/samples_50013.pdf", smp50013, width=5, height=4) 
}

smp_analysis <- plot_grid(smp78870+theme(legend.position="none"), smp50013, label_size = 10,
						  ncol = 1, hjust = 0.2, rel_heights = c(3.2,4) )
if(saveFlag)
	ggsave("RESULTS/InterOpt_Publication/samples_both.pdf", smp_analysis, width=5, height=7.2)

#--------Iterative (RG number analysis)---------
curD <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(sub("(.+)/(.+)","\\1",curD))


data_ids = c('GSE78870', 'GSE50013')
s_scenarios = c('raw', 'norm') #Don't Change the Order of this
res_iter = read_exper_full('RESULTS/RawResults/SelfWeight_iter', data_ids, s_scenarios)


g_iter = plot.Iterative3(data.frame(res_iter$GSE78870$norm$res_source$sd_wgmean$SD$stats)$Mean,
						 data.frame(res_iter$GSE78870$raw$res_source$sd_wgmean$SD$stats)$Mean,
						 data.frame(res_iter$GSE78870$raw$res_source$gmean$SD$stats)$Mean, rg=c(1:19))

g_iter2 = plot.Iterative3(data.frame(res_iter$GSE50013$norm$res_source$sd_wgmean$SD$stats)$Mean,
						 data.frame(res_iter$GSE50013$raw$res_source$sd_wgmean$SD$stats)$Mean,
						 data.frame(res_iter$GSE50013$raw$res_source$gmean$SD$stats)$Mean, rg=c(1:19))

if(saveFlag){
	ggsave("RESULTS/InterOpt_Publication/selfweight_iter_GSE78870.pdf",      g_iter,      width=4.5, height=2.5) 
	ggsave("RESULTS/InterOpt_Publication/selfweight_iter_GSE50013.pdf",      g_iter2,      width=4.5, height=2.5) 
}

#------Experimential-------
source('SOURCE/Experimental_Validation.R')

rows = c('hsa_miR_21_5p_all', 'hsa_miR_21_5p_wall')
p_ghanbari = plot.box.grouped(data_full_m, gr_full_m, rows, F, 'miR-21-5p Expression', F, colors = mycols[c(14,3,4)])
p_ghanbari = p_ghanbari + theme(legend.position = "none")

if(saveFlag){
	ggsave("RESULTS/InterOpt_Publication/exper_ghanbari_21.pdf",
		   p_ghanbari,
		   width = 8, height = 6, units = "cm")
}
#-----Self + Sample + Iter----

pre_w = 24
new_w = 24
scale.one.fig = function(pre_rats, pre_w=24, new_w=24, scale.which=1){
	new_rats = pre_rats*pre_w/new_w
	new_rats[scale.which] = 1-sum(new_rats)+new_rats[scale.which]
	return(new_rats)
}

selfweight_bar <- plot_grid(frg, p1_bar+ggtitle('Breast Cancer (Tissue)'),
							frg, p2_bar+ggtitle('Liver Cancer (Plasma)'),
							label_size = 10, ncol = 4,hjust = 0.2 ,
							rel_widths = scale.one.fig(c(0.18,0.33, 0.02, 0.47)))
selfweight_box <- plot_grid(p1_heat, p1_box, p2_heat+ guides(fill = FALSE), p2_box+ theme(legend.position="none"),
							label_size = 10, ncol = 4,hjust = 0.2 ,
							rel_widths = scale.one.fig(c(0.25,0.25, 0.25, 0.25)))
smp_analysis <- plot_grid(smp78870+theme(plot.title = element_blank()),
						  smp50013+theme(plot.title = element_blank()),
						  ncol = 2, hjust = 0.2 )
iter_analysis <- plot_grid(g_iter,
						   g_iter2,
						   ncol = 2, hjust = 0.2 )

selfweight_full <- plot_grid(selfweight_bar, selfweight_box, smp_analysis, iter_analysis,
							 labels = c('(A)','(B)','(C)','(D)'),
							 rel_heights = c(0.39,0.35, 0.6, 0.4),
							 label_size = 9, nrow = 4,hjust = 0)

if(saveFlag){
	ggsave("RESULTS/InterOpt_Publication/selfweight_full_v3.pdf",
		   selfweight_full,
		   width = new_w, height = 27, units = "cm")
}
