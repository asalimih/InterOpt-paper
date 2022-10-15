process_experiment = function(res, scenarios, measures, wmethod_names, new_order,
							  group_names, Base, precision=3, keep_no_optims=F, rep_sds=F){
	# reorder and attach group info
	res$scenarios = scenarios
	res$expers = lapply(res$expers, function(x) override_names(x,wmethod_names,new_order))
	res$wmethod_names = wmethod_names
	res$Base = Base
	res$group_names = group_names
	
	# summerize data
	res$combs   = aggregate_measures(res, measures, keep_no_optims=keep_no_optims)
	res$means   = as.data.frame(sapply(res$combs, colMeans))
	res$sds     = as.data.frame(sapply(res$combs, function(x) colSds(as.matrix(x))))
	res$summary = extract_summary_table(res$combs, precision)
	
	if(rep_sds)
		res$rep_sds = extract_rep_sds(res)
	
	return(res)
}

save_summary_table = function(res, filename, data_ids){
	for(data_id in data_ids){
		append = ifelse(data_id==data_ids[1], F, T)
		write.xlsx(res[[data_id]]$summary, filename, data_id, append = append)
		write.xlsx(res[[data_id]]$means, filename, paste0(data_id,'_'), append = T)
	}
}

override_names = function(res, new_names, new_order=NULL){

	if(!is.null(new_order)){
		res$res_source = res$res_source[new_order]
		#res$res_comb_weights = res$res_comb_weights[new_order]
		if(length(res$res_target)!=0)
			res$res_target = res$res_target[new_order]
		if('mean_list' %in% names(res)){
			res$mean_list = lapply(res$mean_list, function(a) a[new_order,])
		}
	}
	
	names(res$res_source) = new_names
	#names(res$res_comb_weights) = new_names
	if(length(res$res_target)!=0)
		names(res$res_target) = new_names
	
	return(res)
}

read_exper = function(dir, scenarios){

	res = list()
	res[['expers']] = list()
	for(sc in scenarios){
		res_rds_file = file.path(dir, sprintf('%s/res.rds',sc))
		if(file.exists(res_rds_file)){
			res[['expers']][[sc]] = readRDS(res_rds_file)
			print(res_rds_file)
		}
	}

	return(res)
}

read_exper_full = function(dir, data_ids, scenarios){ #just used in iterative 
	res = list()
	for(data_id in data_ids){
		res[[data_id]] = list()
		for(sc in scenarios){
			res[[data_id]][[sc]] = readRDS(file.path(dir, data_id, sprintf('%s/res.rds',sc)))
		}
		
	}
	return(res)
}

build_measure_df = function(res, measures){
	thetable = 'res_source'
	measure_dfs = lapply(measures, function(measure) as.data.frame(sapply(res[[thetable]], function(x) x[[measure]])))
	names(measure_dfs) = measures
	return(measure_dfs)
}

simplify_measure_df = function(measure_dfs){
	measures = names(measure_dfs[[1]])
	measure_dfs_s = lapply(measures, function(measure){
		a = lapply(measure_dfs, function(x) x[[measure]])# a list of dataframes. each is SD of all wmethods
		b = do.call(cbind,a)
		return(b)
	})
	names(measure_dfs_s) = measures
	return(measure_dfs_s)
}

calc_stability = function(res, measures){
	# for every measure the z-score is calculated based on all weight methods results
	means = sapply(res, function(x) mean(unlist(x)))
	sds   = sapply(res, function(x) sd(unlist(x)))
	res_z = mapply(function(a,m,s){(a-m)/s}, res, means, sds, SIMPLIFY = F)
	Stability = Reduce('+', res_z[measures])
	Stability = Stability+2 # this +2 is just to make numbers mroe positive
	return(Stability)
}

is_no_optim = function(x){
	return(x=='geom' | x=='arith' | grepl('rand',x))
}

aggregate_measures = function(res, measures, keep_no_optims=F){
	no_optims = which(is_no_optim(res$wmethod_names))
	
	res_dfs = lapply(res$expers, function(x) build_measure_df(x,measures))
	res_combs = simplify_measure_df(res_dfs)
	
	# because geom, arith and geom(rand) are repeated in other scenarios we remove them
	if(!keep_no_optims){
		col_del = as.vector(sapply(c(1:(length(res_dfs)-1)),
								   function(x) x * length(res$wmethod_names) + no_optims))
		
		res_combs = lapply(res_combs, function(x) {
			colnames(x)[no_optims] = res$wmethod_names[no_optims]
			return(x[,-col_del])
		})
	}
	# add aggregated Stability
	res_combs$Stability = calc_stability(res_combs, measures)
	return(res_combs)
}

extract_summary_table = function(combs, precision=3){
	means = as.data.frame(sapply(combs, colMeans))
	sds = as.data.frame(sapply(combs, function(x) colSds(as.matrix(x))))
	table = as.data.frame(mapply(function(a,b) paste0(a," ± ",b),
												round(means,precision),
												round(sds,precision), SIMPLIFY=T), stringsAsFactors = F)
	rownames(table) = rownames(means)
	return(table)
}

extract_rep_sds = function(res){
	rep_sds = c()
	for(s_sc in grep('raw_small',names(res$expers),value=T)){
		a = res$expers[[s_sc]]$mean_list
		small_sds = as.data.frame(apply(array(unlist(a), c(nrow(a[[1]]), ncol(a[[2]]), length(a))), c(1,2), sd))
		colnames(small_sds) = colnames(a[[1]])
		rep_sds = rbind(rep_sds,small_sds)
	}
	small_sds[,] = 0 # add the full sample raw scenario with 0 sd because of no repeat
	rep_sds =  rbind(rep_sds, small_sds)
	rownames(rep_sds) = rownames(res$means)
	return(rep_sds)
}
