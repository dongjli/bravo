utils::globalVariables(c("MIP", ".", "j"))

#' @title FDR using correlation threshold
#' @description Computes False Discovery Rate using SNP correlation as proximity measure.
#' @param model Integer vector of selected SNP indices.
#' @param truth Integer vector of true causal SNP indices.
#' @param x SNP matrix.
#' @param threshold Correlation threshold (default 0.9).

FDR_corrected <- function(model, truth, x, threshold = 0.9) {
  xmodel = as.matrix(x[,model])
  xtruth = as.matrix(x[,truth])
  R = abs(cor(xtruth, xmodel))
  falses = apply(R, 2, function(xx) all(xx < threshold))
  return(sum(falses)/max(length(model), 1))
}

#' @title FDR using window size
#' @description Computes False Discovery Rate using base-pair window as proximity measure.
#' @param model Integer vector of selected SNP indices.
#' @param truth Integer vector of true causal SNP indices.
#' @param mapmat Map matrix with chromosome and position columns.
#' @param winsize Window size in base pairs (default 1000).

FDR_WS <- function(model, truth, mapmat, winsize = 1000){
  model_chr = mapmat[model,2] ; model_bp = mapmat[model,3]
  truth_chr = mapmat[truth,2] ; truth_bp = mapmat[truth,3]
  if(length(model) == 0){return(0)}else{
    falses = logical(length(model))
    for(i in 1:length(model)){
      chr_check = model_chr[i] == truth_chr
      bp_check = abs(model_bp[i] - truth_bp) <= winsize
      if(sum(bp_check + chr_check == 2) < 1){falses[i] = T}else{falses[i] = F}
    }
    return(sum(falses)/max(length(model), 1))
  }
}

#' @title FPR using correlation threshold
#' @description Computes False Positive Rate using SNP correlation as proximity measure.
#' @param model Integer vector of selected SNP indices.
#' @param truth Integer vector of true causal SNP indices.
#' @param x SNP matrix.
#' @param threshold Correlation threshold (default 0.9).

FPR_corrected <- function(model, truth, x, threshold = 0.9){
  xmodel = as.matrix(x[,model])
  xtruth = as.matrix(x[,truth])
  R = abs(cor(xtruth, xmodel))
  falses = apply(R, 2, function(xx) all(xx < threshold))
  return(sum(falses)/(ncol(x) - length(truth)))
}

#' @title FPR using window size
#' @description Computes False Positive Rate using base-pair window as proximity measure.
#' @param model Integer vector of selected SNP indices.
#' @param truth Integer vector of true causal SNP indices.
#' @param mapmat Map matrix with chromosome and position columns.
#' @param winsize Window size in base pairs (default 1000).

FPR_WS <- function(model, truth, mapmat, winsize = 1000){
  model_chr = mapmat[model,2] ; model_bp = mapmat[model,3]
  truth_chr = mapmat[truth,2] ; truth_bp = mapmat[truth,3]
  if(length(model) == 0){return(0)}else{
    falses = logical(length(model))
    for(i in 1:length(model)){
      chr_check = model_chr[i] == truth_chr
      bp_check = abs(model_bp[i] - truth_bp) <= winsize
      if(sum(bp_check + chr_check == 2) < 1){falses[i] = T}else{falses[i] = F}
    }
    return(sum(falses)/(nrow(mapmat) - length(truth)))
  }
}

#' @title TPR using correlation threshold
#' @description Computes True Positive Rate using SNP correlation as proximity measure.
#' @param model Integer vector of selected SNP indices.
#' @param truth Integer vector of true causal SNP indices.
#' @param x SNP matrix.
#' @param threshold Correlation threshold (default 0.9).

TPR_corrected <- function(model, truth, x, threshold = 0.9){
  testfor <- as.matrix(x[,model])
  reference <- as.matrix(x[,truth])
  cor_check <- cor(testfor, reference)
  include <- numeric(ncol(cor_check))
  for(i in 1:ncol(cor_check)){
    if(sum(abs(cor_check[,i]) > threshold) > 0){include[i] = 1}else{include[i] = 0}
  }
  tpr = sum(include)/length(truth)
  return(tpr)
}

#' @title TPR using window size
#' @description Computes True Positive Rate using base-pair window as proximity measure.
#' @param model Integer vector of selected SNP indices.
#' @param truth Integer vector of true causal SNP indices.
#' @param mapmat Map matrix with chromosome and position columns.
#' @param winsize Window size in base pairs (default 1000).

TPR_WS <- function(model, truth, mapmat, winsize = 1000){
  if(length(model) == 0){return(0)}
  model_chr = mapmat[model,2] ; model_bp = mapmat[model,3]
  truth_chr = mapmat[truth,2] ; truth_bp = mapmat[truth,3]
  hit = logical(length(truth))
  for(i in 1:length(truth)){
    chr_check = truth_chr[i] == model_chr
    bp_check = abs(truth_bp[i] - model_bp) <= winsize
    if(sum(bp_check + chr_check == 2) < 1){hit[i] = F}else{hit[i] = T}
  }
  return(sum(hit)/length(truth))
}

#' @title Jaccard Index
#' @description Computes Jaccard index between selected and true causal SNPs.
#' @param vars Integer vector of selected SNP indices.
#' @param truth Integer vector of true causal SNP indices.

jcidx <- function(vars, truth){
  result = length(intersect(vars, truth))/length(union(vars, truth))
  return(result)
}

#' @title Clean SNP Matrix
#' @description Removes duplicate SNPs and filters low MAF variants.
#' @param SNPmat A sparse SNP matrix (dgCMatrix).
#' @param MAF_threshold Minor Allele Frequency cutoff (default 0.05).

clean <- function(SNPmat, MAF_threshold = 0.05){
  z = SNPmat
  set.seed(seed = 2441139)
  for(pass in 1:10) {
    v = round(runif(nrow(z)) %*% z, digits = 9)
    dupv = duplicated(v@x)
    ndup = sum(dupv)
    if(ndup == 0) break
    z = z[,!dupv]
    cat("pass = ", pass, "; ", ndup, " many duplicated. ncol(z) = ", ncol(z),"\n", sep="")
    pass = pass + 1
    gc()
  }
  idx <- apply(z, 2, mean)
  max_val <- max(z@x)
  if(max_val == 2){
    x <- z[, idx > 2*MAF_threshold]
  } else {
    x <- z[, idx > MAF_threshold]
  }
  return(x)
}

#' @title Create Parameter Grid
#' @description Builds a grid of lambda and w tuning parameters for SVEN.
#' @param x Cleaned SNP matrix.

create_param_mat <- function(x){
  n = nrow(x); p = ncol(x)
  lamseq = exp(c(seq(log(n/p^2), log(sqrt(n)), l = 10), log(5), log(10)))
  wseq = c(1/p, sqrt(n)/p, n/p) ; wseq = wseq[wseq < 0.5]
  if(length(wseq) == 1) wseq = c(wseq, 0.5)
  param_mat <- as.matrix(expand.grid(lambda = lamseq, w = wseq))
  return(param_mat)
}

#' @title Run SVEN Across Parameter Grid
#' @description Simulates a phenotype and evaluates all parameter combinations via Jaccard index.
#' @param k Simulation replicate index.
#' @param x Cleaned SNP matrix.
#' @param R2 Heritability (default 0.5).
#' @param nspike Number of causal SNPs to simulate (default 20).
#' @param betamax Maximum effect size magnitude.

run_all_params <- function(k, x, R2 = 0.5, nspike = 20, betamax){
  param_mat <- create_param_mat(x)
  set.seed(seed = 481 + k)
  p = ncol(x); n = nrow(x)
  nonzero = sample(1:p, nspike, replace = F)
  beta = runif(nspike, -betamax, betamax)
  xs = scale(as.matrix(x[,nonzero]))
  mu = xs %*% beta
  SSR = var(mu)
  sigma = as.numeric(sqrt(((1/R2) - 1)*SSR))
  y = mu + sigma*rnorm(n)
  m = nrow(param_mat)
  jcid = numeric(m)
  for(i in 1:m){
    params = param_mat[i,]
    set.seed(seed=481)
    object = bravo::sven(x, y, w = params[2], lam = params[1], Ntemp = 3, verbose = T)
    model = object$model.map
    jcid[i] = jcidx(model, nonzero)
  }
  return(jcid)
}

#' @title Tune SVEN Parameters
#' @description Runs 100 parallel simulations and returns the optimal (lambda, w) pair.
#' @param x Cleaned SNP matrix.
#' @param R2 Heritability (default 0.5).
#' @param ehits Expected number of causal SNPs.
#' @param betamax Maximum effect size magnitude (default 1).
#' @param n.cores Number of cores for parallel computation.

tune.sven <- function(x, R2 = 0.5, ehits = 20, betamax = 1, n.cores){
  nspike = ehits
  cl <- makeCluster(n.cores, type = "FORK")
  registerDoParallel(cl)
  output = foreach(j = 1:100, .errorhandling = "pass",
                .packages = c("bravo", "Matrix", "dplyr")) %dopar% {
                run_all_params(j, x, R2, nspike = nspike, betamax)
                            }
  stopCluster(cl)
  for(i in 1:length(output)){
    if(inherits(output[[i]], "error")){
      cat("iteration", i, "failed:", conditionMessage(output[[i]]), "\n")
    }
  }
  output = do.call(cbind, output)
  mean_jc_idx = rowMeans(output, na.rm = T)
  get_param = which.max(mean_jc_idx)
  param_mat <- create_param_mat(x)
  use_param = param_mat[get_param,]
  return(use_param)
}

#' @title Tune SVEN for All Hit Sizes
#' @description Tunes SVEN parameters for small, medium, and large expected hit sizes.
#' @param x Cleaned SNP matrix.
#' @param R2 Heritability (default 0.5).
#' @param betamax Maximum effect size magnitude (default 1).
#' @param n.cores Number of cores for parallel computation.

tune.sven.all <- function(x, R2 = 0.5, betamax = 1, n.cores){
  a1 = tune.sven(x, R2, ehits = 10, betamax, n.cores)
  a2 = tune.sven(x, R2, ehits = 20, betamax, n.cores)
  a3 = tune.sven(x, R2, ehits = 50, betamax, n.cores)
  a = list(small = a1, medium = a2, large = a3)
  return(a)
}

#' @title Run SVEN with Optimal Parameters
#' @description Fits SVEN on (x, y) using the supplied tuned parameters.
#' @param x Cleaned SNP matrix.
#' @param y Phenotype vector.
#' @param params Named numeric vector with lambda and w.
#' @param seed Random seed (default 2441139).

basic.sven.model <- function(x, y, params, seed = 2441139){
  set.seed(seed = seed)
  fit = bravo::sven(x, y, lam = params[1], w = params[2], verbose = T)
  return(fit)
}

#' @title Select Optimal Tuning Parameters
#' @description Full training step: cleans SNP matrix and tunes SVEN hyperparameters.
#' @param X Raw SNP matrix (dgCMatrix).
#' @param R2 Heritability (default 0.5).
#' @param betamax Maximum effect size magnitude (default 1).
#' @param n.cores Number of cores (default: detectCores() - 1).
#' @param hitsize One of "all", "small", "medium", "large" (default "all").
#' @param MAF_threshold MAF cutoff (default 0.05).

parameter_selection <- function(X, R2 = 0.5, betamax = 1, n.cores = max(1, parallel::detectCores() - 1), hitsize = "all", MAF_threshold = 0.05){
  bigX <- X
  X <- clean(X, MAF_threshold = MAF_threshold)
  if(hitsize == "all"){
    optimal_params <- tune.sven.all(x = X, R2 = R2, betamax = betamax, n.cores = n.cores)
  } else {
    ehits <- switch(hitsize,
                    "small"  = 10,
                    "medium" = 20,
                    "large"  = 50,
                    stop("hitsize must be one of: 'small', 'medium', 'large', or 'all'"))
    optimal_params <- tune.sven(x = X, R2 = R2, ehits = ehits, betamax = betamax, n.cores = n.cores)
  }
  result <- list(bigX = bigX, X = X, optimal_params = optimal_params, MAF_thres = MAF_threshold,
                 ehits = hitsize, n.cores = n.cores)
  class(result) <- "svenetics_trained"
  return(result)
}

#' @title UNITE: Post-process SVEN Model
#' @description Propagates MIPs from SVEN hits to the full SNP matrix via LD.
#' @param x Cleaned SNP matrix.
#' @param bigx Full (uncleaned) SNP matrix.
#' @param basic_sven_object Fitted SVEN object.
#' @param y Phenotype vector.
#' @param ehits Expected number of hits.
#' @param threshold MIP threshold (default 0).

unite.sven <- function(x, bigx, basic_sven_object, y, ehits, threshold = 0){
  if(is.null(dim(bigx))) stop(paste("bigx has no dimensions, class:", class(bigx), "length:", length(bigx)))
  bigx <- as(bigx, "dgCMatrix")  # force correct class
  object <- basic_sven_object
  if(length(object$model.map) < ehits/2){
    set.seed(seed=481)
    object <- bravo::sven(x, y, lam = sqrt(nrow(x)), w = nrow(x)/ncol(x), verbose = T)
  }
  hits = bravo::mip.sven(object, threshold = 0)
  big.sd = sqrt(colMSD_dgc(bigx, m = Matrix::colMeans(bigx)))
  R = crossprod(bigx, scale(as.matrix(bigx[,hits$Hits]))) / big.sd / (nrow(bigx) - 1)
  R = abs(as.matrix(R))
  df <- data.frame(SNP = character(0), MIP = numeric(0))
  df <- rbind(df, data.frame(MIP = {(R > 0.9)*R} %*% hits$MIP) %>%
                dplyr::filter(MIP > threshold) %>%
                dplyr::mutate(MIP = ifelse(MIP > 1, 1, MIP)) %>%
                dplyr::mutate(SNP = rownames(.), .before = 1))
  return(df)
}

#' @title BWAS: Bayesian GWAS for a Single Trait
#' @description Runs SVEN and UNITE for one trait.
#' @param x Cleaned SNP matrix.
#' @param y Phenotype vector.
#' @param params Named numeric vector with lambda and w.
#' @param bigx Full SNP matrix.
#' @param ehits Expected number of hits (default 20).

bwas <- function(x, y, params, bigx, ehits = 20){
  object <- basic.sven.model(x, y, params)
  pp_model <- unite.sven(x, bigx, object, y, ehits)
  return(pp_model)
}

#' @title Run GWAS for a Single Trait
#' @description Runs the full GWAS pipeline for the i-th trait column.
#' @param svenetics_trained_object A trained svenetics object from parameter_selection().
#' @param i Trait index.
#' @param hitsize One of "small", "medium", "large", or NULL.

pipeline_single_trait <- function(svenetics_trained_object, i, hitsize = NULL){
  obj = svenetics_trained_object
  
  # resolve params based on hitsize
  if(is.list(obj$optimal_params) && all(c("small","medium","large") %in% names(obj$optimal_params))){
    if(is.null(hitsize)) hitsize <- "medium"
    obj$optimal_params <- obj$optimal_params[[hitsize]]
    obj$ehits <- switch(hitsize, "small" = 10, "medium" = 20, "large" = 50)
  }
  
  traitfile = obj$traitfile
  snpfile = obj$X
  big_snpfile = obj$bigX
  select = intersect(rownames(traitfile), rownames(snpfile))
  snpfile = snpfile[select,]
  traitfile = traitfile[select,]
  omit = is.na(traitfile[,i])
  snpfile = snpfile[!omit,]
  trait_i = traitfile[!omit,i]
  x <- clean(snpfile, MAF_threshold = obj$MAF_thres)
  params <- obj$optimal_params
  print(tryCatch(nrow(x) == length(trait_i)))
  result <- bwas(x = x, y = trait_i, params = params, bigx = big_snpfile, ehits = obj$ehits)
  return(result)
}


#' @title Full Multi-Trait GWAS Pipeline
#' @description Runs GWAS across all traits in parallel and saves results as CSVs.
#' @param svenetics_trained_object A trained svenetics object from parameter_selection().
#' @param traitfile Data frame with sample IDs in column 1 and traits in remaining columns.
#' @param hitsizes Character vector of hit sizes per trait, or NULL for "medium" across all.
#' @export
#' @title Full Multi-Trait GWAS Pipeline
#' @description Runs GWAS across all traits and saves results as CSVs.
#' @param svenetics_trained_object A trained svenetics object from parameter_selection().
#' @param traitfile Data frame with sample IDs in column 1 and traits in remaining columns.
#' @param hitsizes Character vector of hit sizes per trait, or NULL for "medium" across all.
#' @param save_dir Directory path where results will be saved (default "~/SVENETICS_RESULTS").
svenetics_pipeline <- function(svenetics_trained_object, traitfile, hitsizes = NULL, save_dir = "~/SVENETICS_RESULTS"){
  obj = svenetics_trained_object
  trait_names <- colnames(traitfile)[-1]
  traits <- traitfile[,-1]
  rownames(traits) <- traitfile[,1]
  obj$traitfile = traits
  n.cores = obj$n.cores
  
  if(is.null(hitsizes)){
    hitsizes <- rep("medium", ncol(traits))
  }
  
  registerDoParallel(n.cores)
  output = foreach(j = 1:ncol(traits), .errorhandling = "pass",
          .packages = c("bravo", "Matrix", "dplyr")) %dopar% {
          pipeline_single_trait(obj, j, hitsize = hitsizes[j])
                               }
  
  # output <- vector("list", ncol(traits))
  # for(j in 1:ncol(traits)){
  #   output[[j]] <- pipeline_single_trait(obj, j, hitsize = hitsizes[j])
  # }
  
  save_dir <- path.expand(save_dir)
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  
  for(j in 1:length(output)){
    result_j <- output[[j]]
    result_j <- result_j[order(result_j$MIP, decreasing = TRUE),]
    output[[j]] <- result_j
    write.csv(result_j, file = file.path(save_dir, paste0(trait_names[j], "_GWAS_results.csv")), row.names = FALSE)
  }
  return(output)
}

#' @title Estimate SVEN Runtime
#' @description Times a single SVEN run on the given SNP matrix using a random phenotype.
#' @param X SNP matrix.

calc.runtime <- function(X){
  y = rnorm(nrow(X))
  start_time <- Sys.time()
  ff = bravo::sven(X, y)
  end_time <- Sys.time()
  return(end_time - start_time)
}




#' @title Convert numeric matrix to sparse matrix
#' @description Reads a numeric genotype file and converts it to a sparse matrix format. 
#' @param file.name Path to the numeric genotype file. Could be (and should be) gzipped.
#' @param num.genotypes Maximum number of genotypes to read. An upper bound is OK.
#' @param separator "\\t" or "," etc that separates the entries in a line.
#' @param progress Whether to show a progress bar (default TRUE).



dense2sparse <- function(file.name, num.genotypes, separator,progress = TRUE) {
  con <- file(file.name, open = "r")
  
  if(progress) pb = txtProgressBar(min = 0, max = num.genotypes, initial = 0,
                                   title = "Progress:",style = 3)
  
  header_line <- readLines(con, n = 1)
  col_names <- strsplit(header_line, split = separator)[[1]][-1]
  
  row_names <- character(num.genotypes)
  row_data <- vector(mode = "list", length = num.genotypes)
  row_ptr <- integer(num.genotypes+1)
  col_indices <- vector(mode = "list", length = num.genotypes)
  
  line_num <- 0
  
  repeat {
    line <- readLines(con, n = 1)
    if ({length(line) == 0} || (line_num >= num.genotypes)) {
      cat("\nNumber of genotypes read: ",line_num)
      break
    }
    line_num <- line_num + 1
    # cat("Processing line:  ",line_num,"\n")
    fields <- strsplit(line, split = separator)[[1]]
    row_names[line_num] <- fields[1]
    
    vals <- as.numeric(fields[-1])
    nz_idx <- which(vals != 0)
    col_indices[[line_num]] <- nz_idx - 1  # 0-based for dgRMatrix
    row_data[[line_num]] <- vals[nz_idx]
    row_ptr[line_num + 1] <- row_ptr[line_num] + length(nz_idx)
    if(progress) setTxtProgressBar(pb,line_num)
    
  }
  if(progress) close(pb)
  
  close(con)
  
  # Trim excess
  if(line_num < num.genotypes) {
    row_ptr = row_ptr[1:{line_num+1}]
    row_names = row_names[1:line_num]
  }
  
  gc()
  x = sparseMatrix(j=unlist(col_indices, use.names = FALSE),
                   p = row_ptr,
                   x = unlist(row_data, use.names = FALSE),
                   dims = c(line_num,length(col_names)),
                   dimnames = list(row_names,col_names),
                   index1 = FALSE,
                   repr = "C")
  return(x)
}



