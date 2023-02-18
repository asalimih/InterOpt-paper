
data_ids = data_pair_info[['data_ids']]
data_labels = data_pair_info[['data_labels']]
data_smp_nums = data_pair_info[['data_smp_nums']]

#----Selfweight Plots----
wmethod_filter_s = res_selfweight[[data_ids[1]]]$wmethod_names[c(-1,-3)] #removes arith and geom(rand)
cls_selfweight = mycols[c(14,15,20,2,4,3)]
names(cls_selfweight) <- c('','Raw CT','Normalized CT')

frg <- ggplot(data.frame()) + geom_point()+ theme(panel.background = element_rect(fill = 'white', colour = 'white'))

p1_bar = plot.measure.bar(res_selfweight[[data_ids[1]]], MainMeasure, cls_selfweight, wmethod_filter_s)+ theme(legend.position="none")
p2_bar = plot.measure.bar(res_selfweight[[data_ids[2]]], MainMeasure, cls_selfweight, wmethod_filter_s)

p1_box = plot.measure.box(res_selfweight[[data_ids[1]]], MainMeasure, cls_selfweight, wmethod_filter_s)+ theme(legend.position="none")
p2_box = plot.measure.box(res_selfweight[[data_ids[2]]], MainMeasure, cls_selfweight, wmethod_filter_s)

p1_heat = plot.pairwise_wilcox2(res_selfweight[[data_ids[1]]], MainMeasure, wmethod_filter_s, rm.norm.res = T)+ guides(fill = FALSE)
p2_heat = plot.pairwise_wilcox2(res_selfweight[[data_ids[2]]], MainMeasure, wmethod_filter_s, rm.norm.res = T)

#-----Sample Analysis Plots----
if(sample_plot_flag){
	
	wmethod_names =c('geom','arith(cv)','geom(cv)','arith(sd)','geom(sd_r)','geom(sd)','geom(sd+)')
	cls = mycols[c(14,15,20,2,4,7,3)]
	meas = ifelse(MainMeasure=='Stability', 'SD', MainMeasure)
	better_ylabel = ylab(bquote(Mean~(.(meas)~(RG[2]))))
	
	smp_1 = plot.sample.size(res_samples[[data_ids[1]]], data_smp_nums[[1]], meas,cls) +
		ggtitle(data_labels[1]) + better_ylabel
	
	smp_2 = plot.sample.size(res_samples[[data_ids[2]]], data_smp_nums[[2]], meas,cls) +
		ggtitle(data_labels[2]) + better_ylabel
	
}else{
	smp_1 = frg
	smp_2 = frg
}

#--------Iterative (RG number analysis)---------
if(iter_plot_flag){
	
	s_scenarios = c('raw', 'norm') #Don't Change the Order of this
	res_iter = read_exper_full('RESULTS/RawResults/SelfWeight_iter', data_ids, s_scenarios)
	
	
	meas = ifelse(MainMeasure=='Stability', 'SD', MainMeasure)
	better_ylabel = ylab(bquote(Mean~(.(meas)~(RG['d']))))
	
	g_iter = plot.Iterative3(data.frame(res_iter[[data_ids[1]]]$norm$res_source$geom_sd[[meas]]$stats)$Mean,
							 data.frame(res_iter[[data_ids[1]]]$raw$res_source$geom_sd[[meas]]$stats)$Mean,
							 data.frame(res_iter[[data_ids[1]]]$raw$res_source$geom[[meas]]$stats)$Mean, rg=c(1:19)) + better_ylabel
	
	g_iter2= plot.Iterative3(data.frame(res_iter[[data_ids[2]]]$norm$res_source$geom_sd[[meas]]$stats)$Mean,
							 data.frame(res_iter[[data_ids[2]]]$raw$res_source$geom_sd[[meas]]$stats)$Mean,
							 data.frame(res_iter[[data_ids[2]]]$raw$res_source$geom[[meas]]$stats)$Mean, rg=c(1:19)) + better_ylabel
}else{
	g_iter = frg
	g_iter2 = frg
}

#-----Self + Sample + Iter----

pre_w = 24
new_w = 24
scale.one.fig = function(pre_rats, pre_w=24, new_w=24, scale.which=1){
	new_rats = pre_rats*pre_w/new_w
	new_rats[scale.which] = 1-sum(new_rats)+new_rats[scale.which]
	return(new_rats)
}

selfweight_bar <- plot_grid(frg, p1_bar+ggtitle(data_labels[1]),
							frg, p2_bar+ggtitle(data_labels[2]),
							label_size = 10, ncol = 4,hjust = 0.2 ,
							rel_widths = scale.one.fig(c(0.18,0.33, 0.02, 0.47)))
selfweight_box <- plot_grid(p1_heat, p1_box, p2_heat+ guides(fill = FALSE), p2_box+ theme(legend.position="none"),
							label_size = 10, ncol = 4,hjust = 0.2 ,
							rel_widths = scale.one.fig(c(0.25,0.25, 0.25, 0.25)))
smp_analysis <- plot_grid(smp_1+theme(plot.title = element_blank()),
						  smp_2+theme(plot.title = element_blank()),
						  ncol = 2, hjust = 0.2 )
iter_analysis <- plot_grid(g_iter,
						   g_iter2,
						   ncol = 2, hjust = 0.2 )

selfweight_full[[MainMeasure]] <- plot_grid(selfweight_bar, selfweight_box, smp_analysis, iter_analysis,
											labels = c('(A)','(B)','(C)','(D)'),
											rel_heights = c(0.39,0.35, 0.6, 0.4),
											label_size = 9, nrow = 4,hjust = 0)