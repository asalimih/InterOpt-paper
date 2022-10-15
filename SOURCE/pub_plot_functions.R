
library(pheatmap)
library(matrixStats)
library(ggplot2)
library(reshape)
library(scales)
library(cowplot)
library(ggsignif)

plot.measure.bar = function(res, measure, cls, wmethods=NULL, tip_num_size = 2.6){
	method_labels = sub('.*\\.(.*)','\\1',rownames(res$means))
	d.m = cbind(res$means, Base=res$Base, Method=method_labels)
	d.m$Method <- factor(d.m$Method, levels = res$wmethod_names)
	
	# filtering weight methods if provided
	if(!is.null(wmethods)){
		keep = which(method_labels %in% wmethods)
		d.m = d.m[keep,]
		method_labels = method_labels[keep]
	}
	
	measure_vec = d.m[[measure]]
	y_limits = c(min(measure_vec)-0.5*min(measure_vec),
				 max(measure_vec)+0.003*max(measure_vec) )
	txt_dodge = 0.01*(y_limits[2]-y_limits[1])

	names(cls) = wmethods
	p = ggplot(d.m, aes(x=Method, y=.data[[measure]], label=round(.data[[measure]],2), fill=Method)) + 
		
		geom_bar(stat = "identity", position = 'dodge', width = 1) +
		geom_text(size = tip_num_size, color='white', aes(y = .data[[measure]]-txt_dodge ),
				  hjust='right', vjust='middle', position = position_dodge(1)) +
		facet_grid(Base~., scales = "free_y", space = "free_y", switch = "y") + 
		coord_flip()+
		xlab('Optimized on')+ylab(paste('Mean',measure))+
		theme_bw()+
		theme(panel.grid.minor = element_blank(),
			  panel.grid.major = element_blank(),
			  panel.border = element_blank(),
			  axis.text.y = element_blank(),
			  axis.ticks = element_blank(),
			  plot.title = element_text(size=10),
			  strip.text.y = element_text(size = 8))+guides(shape=NA)+
		scale_fill_manual(values=cls)+
	#	scale_x_discrete(labels= method_labels) + 
		scale_y_continuous(limits = y_limits,oob = rescale_none )
	return(p)
}

plot.measure.box = function(res, measure, cls, wmethods=NULL) {
	d.m = melt(res$combs[[measure]])
	d.m$Base = res$Base[match(d.m$variable,colnames(res$combs[[measure]]))]
	colnames(d.m) = c('Method','value','Base')
	d.m$Method = sub('.*\\.(.*)','\\1',d.m$Method) # causes raw and norm be beside each other
	d.m$Method <- factor(d.m$Method, levels = res$wmethod_names)
	
	# filtering weight methods if provided
	if(!is.null(wmethods)){
		keep = which(d.m$Method %in% wmethods)
		d.m = d.m[keep,]
	}
	
	# compute lower and upper whiskers
	ylim1 = boxplot.stats(d.m$value)$stats[c(1, 5)]
	names(cls) = wmethods
	p = ggplot(d.m, aes(x=Method, y=value, label=round(value,2), fill=Method)) +
		#stat_boxplot(geom ='errorbar', width = 0.4)+
		geom_boxplot(outlier.shape = NA)+#, varwidth = T)+
		facet_grid(~Base, scales = "free_x", space = "free_x", switch = "x") + 
		xlab('weight method')+ylab(measure)+
		theme_bw()+
		theme(panel.grid.minor = element_blank(),
			  #panel.grid.major = element_blank(),
			  axis.text.x = element_blank(),
			  panel.border = element_rect(colour='gray'),
			  axis.ticks.x = element_blank())+guides(shape=NA)+
		scale_fill_manual(values=cls)+xlab('')+
		coord_cartesian(ylim = c(ylim1[1]/1.04, ylim1[2]*1.09))
	return(p)
}

cormat_to_df = function(mat, lbls){
	a = c()
	for(i in seq(nrow(mat))){
		d=data.frame(X = lbls[i],Y = lbls,p.value = mat[i,])
		a = rbind(a,d)
	}
	a$X <- factor(a$X, levels = lbls)
	a$Y <- factor(a$Y, levels = lbls)
	return(a)
}

better_label = function(x, rm.norm.res=FALSE){
	lvls = levels(x)
	if(rm.norm.res & (!any(grepl('rnaseq', lvls))))
		lvls = gsub('raw\\.(.*)','\\1',lvls)
	else
		lvls = gsub('raw\\.(.*)','Raw CT \\1',lvls)
	lvls = gsub('norm\\.(.*)','Norm. \\1',lvls)
	lvls = gsub('raw_small(.+)\\.(.*)','(n:\\1) Raw CT \\2',lvls)
	lvls = gsub('rnaseq\\.(.*)','RNA-Seq \\1',lvls)
	y = x
	levels(y) = lvls
	return(y)
}

plot.pairwise_wilcox2 = function(resy, measure, wmethods=NULL, rm.norm.res=FALSE){
	res = resy$combs
	if(rm.norm.res)
		sub.base = grep('Norm',resy$Base, invert = T)
	else
		sub.base = c(1:length(resy$Base))
	nn = length(sub.base)
	combs = combn(sub.base, 2)
	comp_res = c()
	gg = matrix(1,nn,nn) #matrix of the p values

	for(i in seq(ncol(combs))){
		wtest1 = wilcox.test(res[[measure]][[combs[2,i]]], res[[measure]][[combs[1,i]]], paired = TRUE, alternative = "less")
		wtest2 = wilcox.test(res[[measure]][[combs[1,i]]], res[[measure]][[combs[2,i]]], paired = TRUE, alternative = "less")
		gg[combs[1,i],combs[2,i]] = wtest1$p.value
		gg[combs[2,i],combs[1,i]] = wtest2$p.value
	}
	a = cormat_to_df(gg, colnames(res[[1]])[sub.base])

	a$p.value[a$p.value<0.01] = '<0.01'
	a$p.value[a$p.value!='<0.01'] = 'ns'
	a$p.value = factor(a$p.value)
	
	# filtering weight methods if provided
	if(!is.null(wmethods)){
		method_labels_x = sub('.*\\.(.*)','\\1',a$X)
		method_labels_y = sub('.*\\.(.*)','\\1',a$Y)
		keep = which((method_labels_x %in% wmethods) &
					 (method_labels_y %in% wmethods))
		a = a[keep, ]
	}
	
	a$X = better_label(a$X, rm.norm.res)
	a$Y = better_label(a$Y, rm.norm.res)
	
	ggplot(a, aes(X, Y, fill= p.value)) + 
		geom_tile(color='white') +
		scale_fill_manual(values = c(mycols[3],'lightgray'))+ 
		labs(x = NULL, y = NULL) +
		theme(axis.text.x = element_text(angle = 45, hjust=1))
}

plot.sample.size = function(res, data_smp_nums, measure, cls){
	
	sample_num = rep(as.character(data_smp_nums),each=length(res$wmethod_names))
	Method = rep(res$wmethod_names,length(data_smp_nums))
	sd_bar = res$rep_sds[[measure]]

	d2 = cbind(res$means, sample_num, Method, sd_bar)
	d2$Method = factor(d2$Method, levels=res$wmethod_names)
	d2$sample_num = factor(d2$sample_num, levels=as.character(data_smp_nums))
	ggplot(d2, aes(x=sample_num, y=.data[[measure]], colour=Method)) +
		geom_line(group = d2$Method, position = position_dodge(0.5)) + 
		geom_point(position = position_dodge(0.5))+
		geom_errorbar(aes(ymax = .data[[measure]] + sd_bar, ymin=.data[[measure]] - sd_bar),
					  width=0.3, size=0.5, position = position_dodge(0.5))+
		theme_bw()+
		ylab(measure) +
		xlab('Number of Samples') +
		scale_colour_manual(name = '', values = cls, labels = res$wmethod_names)+
		theme(panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(),
			  legend.position="bottom")
}

plot.Iterative <- function(meanSt,meanSt2,method_names=c('geom(sd)','geom'), log_flag=F) {
	if(log_flag){
		meanSt = log2(meanSt)
		meanSt2 = log2(meanSt2)
	}
	trend <- data.frame(mean=meanSt, mean2=meanSt2, iter=c(2:(length(meanSt)+1)))
	trend.m = melt(trend, measure.vars = c('mean','mean2'))
	levels(trend.m$variable) = method_names
	g <- ggplot(trend.m, aes(x=iter, y=value, colour=variable)) +
		geom_line() + 
		geom_point()+
		theme_bw()+
		ylab(ifelse(log_flag, 'log2(SD)', 'SD')) +
		xlab('Number of reference genes') +
		scale_x_continuous(breaks=seq(2,length(meanSt)+1,2))+
		scale_colour_manual(name = '',
							values = c(mycols[3], mycols[8]),
							labels = method_names) +
		theme(panel.grid.minor.x = element_blank(),
			  legend.position="bottom")
	g
}

plot.Iterative3 <- function(meanSt,meanSt2, meanSt3,
							method_names=c('geom(sd) Normalized', 'geom(sd) Raw CT', 'geom'),
							log_flag=F, rg = NULL){
	if(is.null(rg))
		rg = c(1:length(meanSt))
	
	if(log_flag){
		meanSt = log2(meanSt)
		meanSt2 = log2(meanSt2)
		meanSt3 = log2(meanSt3)
	}
	trend <- data.frame(mean=meanSt, mean2=meanSt2, mean3=meanSt3, iter=c(2:(length(meanSt)+1)))
	trend = trend[rg,]
	trend.m = melt(trend, measure.vars = c('mean','mean2','mean3'))
	levels(trend.m$variable) = method_names
	g <- ggplot(trend.m, aes(x=iter, y=value, colour=variable)) +
		geom_line() + 
		geom_point()+
		theme_bw()+
		ylab(expression(Mean~(SD~(RG['d'])))) +
		xlab('Number of reference genes (d)') +
		scale_x_continuous(breaks=seq(2,length(meanSt)+1,1))+
		scale_colour_manual(name = '',
							values = c(mycols[16],mycols[3], mycols[14]),
							labels = method_names) +
		theme(panel.grid.minor.x = element_blank(),
			  legend.position="bottom")
	g
}
 

plot.box.grouped <- function(d,gr,mirnames,logFlag=TRUE,un="CPM",jitter=T,
							 colors=mycols, oneCol=F, oneColNum=1,
							 axis_scales = "fixed", facet_labels=c('geom','geom(sd)')) {
	if(logFlag)
		exp = log2(d[mirnames,]+1)
	else
		exp = d[mirnames,]
	rownames(exp) = facet_labels
	mirexp = data.frame(t(exp),Group = gr,stringsAsFactors = F)
	colnames(mirexp) <- gsub("\\.","-",colnames(mirexp))
	colnames(mirexp) <- gsub("_","-",colnames(mirexp))
	mirexp.m = reshape2::melt(mirexp)
	colnames(mirexp.m) = c("Group","miRNA","expression")
	cls <- colors[1:length(mirnames)]
	names(cls) <- colnames(mirexp)[1:length(mirnames)]
	
	p = ggplot(mirexp.m, aes(x=Group, y=expression, fill=miRNA)) + 
		geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
		geom_signif(comparisons = list(c("T", "N")), map_signif_level=TRUE, tip_length=0, y_position=c(6.3)) +
		facet_wrap(~miRNA ,scales = axis_scales)+
		labs(y=ifelse(logFlag,paste0("log2(",un,")"),un))+theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
		theme(axis.text=element_text(size=11),
			  axis.title=element_text(size=12),
			  legend.title = element_text(size=12),
			  legend.text = element_text(size=11,face = "plain"))
	
	if (jitter) p = p+geom_jitter(position=position_jitter(width=.1, height=0), color='blue', alpha=0.5, size=1)
	if (oneCol) {cls = rep(colors[oneColNum],length(cls)); names(cls)=colnames(mirexp)[1:length(mirnames)]}
	p = p+scale_fill_manual(values=cls, labels=mirnames)
	p
}
