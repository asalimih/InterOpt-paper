run_exper_pcr_selfweight = function(data_source, gr_source, k,
									weight_methods=c('arith','geom', 'random','arith_cv','geom_cv','arith_sd','geom_sd','sd_simple'),
									cv_th = 1000,
									dir = 'RESULTS/', remove_left_over=T, mc.cores=32){
	
	#data_norm = 2^(mean(data_source)-t(t(data_source) - colMeans(data_source)))
	
	#cvs = matrixStats::rowSds(data_norm)/rowMeans(data_norm)
	#data_source = data_source[cvs<cv_th,]
	# dir must have the /
  
  mirs_for_combinations = rownames(data_source)[rowSums(data_source>39)==0]
  
	res = list()
	
	#  Weights From Raw
	res[['raw']] = run_experiment(data_source = data_source,
								  gr_source = gr_source,
								  sub_names = mirs_for_combinations,
								  ctVal_source = T,
								  tmpFolder = paste0(dir,'raw'),
								  weights_from_raw= T,
								  k = k,
								  weight_methods=weight_methods,
								  remove_left_over=remove_left_over,
								  mc.cores=mc.cores)
	
	#  Weights From Normalized
	res[['norm']] = run_experiment(data_source = data_source,
								   gr_source = gr_source,
								   sub_names = mirs_for_combinations,
								   ctVal_source = T,
								   tmpFolder = paste0(dir,'norm'),
								   weights_from_raw= F,
								   k = k,
								   weight_methods=weight_methods,
								   remove_left_over=remove_left_over,
								   mc.cores=mc.cores)
	
	
	
	return(res)
}

run_rawsmall_rep = function(data_source, gr_source, k, weight_methods, sub_names=NULL, dir = 'RESULTS/',
							remove_left_over=T, mc.cores=32, s_n=10, repeat_n = NULL, norm_method='high_exp'){
	#  Weights From Raw with small sample size
	res_adder = function(x,y){
		g = x
		measure_cols = c((2*k+1):ncol(x))
		g[,measure_cols] = x[,measure_cols] + y[,measure_cols]
		return(g)
	}
	
	scenario_name = paste0('raw_small',s_n)
	if(is.null(repeat_n)) # number of repeats calculated automatically
		repeat_n = max(4,ceiling(2*ncol(data_source)/s_n))
	
	res_sum = NULL #sum of all repeats for each comb - at last divided by the number of reps
	mean_list = list() # mean of all combs for each wmethod and measure as a matrix
	for(i in 1:repeat_n){
		cat('(', i,'/',repeat_n,') -----------\n')
		sub_samples = sample_balanced2(1:ncol(data_source),s_n, gr_source)
		res_small = run_experiment(data_source = data_source,
								   gr_source = gr_source,
								   sub_names = sub_names,
								   ctVal_source = T,
								   tmpFolder = dir,
								   weights_from_raw= T,
								   sub_samples_for_weights = sub_samples,
								   k = k,
								   weight_methods=weight_methods,
								   remove_left_over=remove_left_over,
								   mc.cores=mc.cores,
								   norm_method=norm_method)
		if(!is.null(res_sum)){ # first run
			res_sum = mapply(res_adder, res_sum, res_small$res_source, SIMPLIFY = F)
		}else{
			res_sum = res_small$res_source
		}
		mean_list[[i]] = t(sapply(res_small$res_source, function(x) colMeans(x[,(2*k+1):ncol(x)])))
	} #end rep
	res_sum = lapply(res_sum, function(x) {
		a = x
		measure_cols = c((2*k+1):ncol(a))
		a[,measure_cols] = a[,measure_cols] / repeat_n
		return(a)
	})
	res_small$res_source = res_sum
	res_small$mean_list = mean_list 
	saveRDS(res_small, file.path(dir,'res.rds')) # to save the mean in the rds (the last repeat was there)
	
	return(res_small)
}

run_exper_pcr_sample_num = function(data_source, gr_source, k, small_subset_n = c(10,20,30,40), repeat_n = NULL,
									weight_methods=c('arith','geom', 'random','arith_cv','geom_cv','arith_sd','geom_sd','sd_simple'),
									cv_th = 1000,
									dir = 'RESULTS/', remove_left_over=T, mc.cores=32){
	# dir must have the /
	# data_norm = 2^(mean(data_source)-t(t(data_source) - colMeans(data_source)))
	# cvs = matrixStats::rowSds(data_norm)/rowMeans(data_norm)
	# data_source = data_source[cvs<cv_th,]
	
  mirs_for_combinations = rownames(data_source)[rowSums(data_source>39)==0]
  
	res = list()
	
	for(s_n in small_subset_n){
		scenario_name = paste0('raw_small',s_n)
		res[[scenario_name]] = run_rawsmall_rep(data_source = data_source,
												gr_source = gr_source,
												sub_names = mirs_for_combinations,
												k = k,
												s_n = s_n,
												repeat_n = repeat_n,
												dir = paste0(dir,scenario_name),
												mc.cores=mc.cores,
												weight_methods = weight_methods,
												remove_left_over = remove_left_over)
	}
	
	#  Weights From Raw full sample
	res[['raw']] = run_experiment(data_source = data_source,
								  gr_source = gr_source,
								  sub_names = mirs_for_combinations,
								  ctVal_source = T,
								  tmpFolder = paste0(dir,'raw'),
								  weights_from_raw= T,
								  k = k,
								  weight_methods=weight_methods,
								  remove_left_over=remove_left_over,
								  mc.cores=mc.cores)
	
	#  Weights From Normalized
	res[['norm']] = run_experiment(data_source = data_source,
								   gr_source = gr_source,
								   sub_names = mirs_for_combinations,
								   ctVal_source = T,
								   tmpFolder = paste0(dir,'norm'),
								   weights_from_raw= F,
								   k = k,
								   weight_methods=weight_methods,
								   remove_left_over=remove_left_over,
								   mc.cores=mc.cores)
	return(res)
}

run_exper_generalization = function(data_source, gr_source, data_target , gr_target, k, small_subset_n = c(10,15,20,30,40), repeat_n = NULL,
									weight_methods=c('arith','geom', 'random','arith_cv','geom_cv','arith_sd','geom_sd','sd_simple'), 
									norm_method = 'high_exp',dir = 'RESULTS/', remove_left_over=T, mc.cores=32){
	# dir must have the /
  
  mirs_for_combinations = rownames(data_target)[rowSums(data_target>39)==0]
	common = intersect(mirs_for_combinations, rownames(data_source))
	fixed_rep = !is.null(repeat_n)
	res = list()
	
	# Weights From Raw
	sub_samples = sample_balanced2(1:ncol(data_target),106, gr_target)
	res[['raw']] =  run_experiment(data_source = data_target,
								   gr_source = gr_target,
								   sub_names = common,
								   sub_samples_for_weights = sub_samples,
								   ctVal_source = T,
								   tmpFolder = paste0(dir,'raw'),
								   weights_from_raw= T,
								   k = k,
								   weight_methods=weight_methods,
								   mc.cores=mc.cores,
								   remove_left_over = remove_left_over,
								   norm_method=norm_method)
	
	
	#  Weights From Raw with small sample size
	for(s_n in small_subset_n){
		scenario_name = paste0('raw_small',s_n)
		res[[scenario_name]] = run_rawsmall_rep(data_source = data_target,
												gr_source = gr_target,
												sub_names = common,
												k = k,
												s_n = s_n,
												repeat_n = repeat_n,
												dir = paste0(dir,scenario_name),
												mc.cores=mc.cores,
												weight_methods = weight_methods,
												remove_left_over = remove_left_over,
												norm_method=norm_method)
	}
	
	#  Weights From Normalized
	res[['norm']] = run_experiment(data_source = data_target,
								   gr_source = gr_target,
								   sub_names = common,
								   ctVal_source = T,
								   tmpFolder = paste0(dir,'norm'),
								   weights_from_raw= F,
								   k = k,
								   weight_methods=weight_methods,
								   mc.cores=mc.cores,
								   remove_left_over = remove_left_over,
								   norm_method=norm_method)

	#  Weights From RNA-Seq
	res[['rnaseq']] = run_experiment(data_source = data_source,
									 gr_source = gr_source,
									 sub_names = common, # this is curcial to keep the order the same with previous
									 ctVal_source = F,
									 data_target = data_target,
									 gr_target = gr_target,
									 ctVal_target = T,
									 tmpFolder = paste0(dir,'rnaseq'),
									 k = k,
									 weight_methods=weight_methods,
									 mc.cores=mc.cores,
									 remove_left_over = remove_left_over,
									 norm_method=norm_method)

	
	return(res)
}

run_exper_pcr_selfweight_iter = function(data_source, gr_source, k, keep,
										 weight_methods=c('arith','geom', 'random','arith_cv','geom_cv','arith_sd','geom_sd','sd_simple'), 
										 dir = 'RESULTS/', remove_left_over=T, mc.cores=32){
	# dir must have the /
	res = list()
	mirs_for_combinations = rownames(data_source)[rowSums(data_source>39)==0]
	
	# Weights From Raw
	res[['raw']] = run_experiment(data_source = data_source,
								  gr_source = gr_source,
								  sub_names = mirs_for_combinations,
								  ctVal_source = T,
								  tmpFolder = paste0(dir,'raw'),
								  weights_from_raw= T,
								  k = k,
								  iter = T,
								  retain_iter = F,
								  keep = keep,
								  weight_methods=weight_methods,
								  mc.cores=mc.cores)
	
	#  Weights From Normalized
	res[['norm']] = run_experiment(data_source = data_source,
								   gr_source = gr_source,
								   sub_names = mirs_for_combinations,
								   ctVal_source = T,
								   tmpFolder = paste0(dir,'norm'),
								   weights_from_raw= F,
								   k = k,
								   iter = T,
								   retain_iter = F,
								   keep = keep,
								   weight_methods = weight_methods,
								   mc.cores=mc.cores)
	
	return(res)
}

sample_balanced = function(x, size, gr){
	gr = as.factor(gr)
	groups = levels(gr)
	n_gr = length(groups)
	if(n_gr==1)
		return(sample(x, size))
	smp_num = rep(floor(size/n_gr),n_gr)
	diff = size - sum(smp_num)
	if(diff>0)
		smp_num[1:diff] = smp_num[1:diff] + 1
	final = c()
	for(i_gr in seq(n_gr)){
		cands = x[which(gr==groups[i_gr])]
		final = c(final, sample(cands, smp_num[i_gr]))
	}
	return(final)
}

sample_balanced2 = function(x, size, gr){
  gr = as.factor(gr)
  groups = levels(gr)
  n_gr = length(groups)
  rat = size/length(x)
  n_grs = table(gr)
  
  if(n_gr==1)
    return(sample(x, size))
  
  smp_num = sapply(1:n_gr, function(i) floor(rat*n_grs[i]))
  diff = size - sum(smp_num)
  if(diff>0){
    groups_with_capacity = which(n_grs>smp_num)
    for(i in rep(groups_with_capacity,2)){
      if(diff>0){
        smp_num[i] = smp_num[i] + 1
        diff = diff -1
      }else
        break
    }
  }
  
  final = c()
  for(i_gr in seq(n_gr)){
    cands = x[which(gr==groups[i_gr])]
    final = c(final, sample(cands, smp_num[i_gr]))
  }
  return(final)
}

run_exper_simulation = function(data_source, gr_source, data_source_norm, combs_name_mat, k=2,
								weight_methods=c('arith','geom', 'random','arith_cv','geom_cv','arith_sd','geom_sd','sd_simple'),
								dir = 'RESULTS/', remove_left_over=T, mc.cores=8){
	# dir must have the /
	
	res = list()
	
	#  Weights From Raw
	res[['raw']] = run_experiment(data_source = data_source,
								  data_source_norm = data_source_norm,
								  combs_name_mat = combs_name_mat,
								  gr_source = gr_source,
								  ctVal_source = T,
								  algors = 'SDCV',
								  tmpFolder = paste0(dir,'raw'),
								  weights_from_raw= T,
								  k = k,
								  weight_methods=weight_methods,
								  remove_left_over=remove_left_over,
								  mc.cores=mc.cores)
	
	
	#  Weights From Normalized
	res[['norm']] = run_experiment(data_source = data_source,
								   data_source_norm = data_source_norm,
								   combs_name_mat = combs_name_mat,
								   gr_source = gr_source,
								   ctVal_source = T,
								   algors = 'SDCV',
								   tmpFolder = paste0(dir,'norm'),
								   weights_from_raw= F,
								   k = k,
								   weight_methods=weight_methods,
								   remove_left_over=remove_left_over,
								   mc.cores=mc.cores)
	
	return(res)
}
