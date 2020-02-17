library(bsvs)
library(compiler)
library(doParallel)
#library(MASS)

source("gen_y.R")
source("gen_xmat.R")
source("SIS.R")
source("fowardregression.R")

set.seed(27)


n = 500; p = 100000; p.true = c(9, 12, 15); m = c(3,4,5); delta = sqrt(0.1)
w = exp(seq(log(1/p), log(0.1), length.out = 20))
w = round(w, 3)
logw = log(w/(1-w))
lam = c(p/n, n/log(p), n/p)
lam = round(lam, 3)
#w2 = n*log(n)/p
trueidx <- list(1:p.true[1], 1:p.true[2], 1:p.true[3])
nonzero.value <- 2
TrueBeta = list(numeric(p), numeric(p), numeric(p))
#TrueBeta[trueidx] <- nonzero.value
for(i in 1:length(m)){
  TrueBeta[[i]][trueidx[[i]]] <- nonzero.value
}
Rsquared = c(0.3, 0.5, 0.7, 0.9)
acc_bis = 0 
tpr_bis = 0
fdr_bis = 0

registerDoParallel(8)

res <- foreach(i=1:100, .packages = c("bsvs", "compiler"), .combine = cbind)%dopar%{
  for(j in 1:length(p.true)){     #j is different model size
    for(k in 1:length(Rsquared)){       #k is different R2
      xmat.train = dtgen.group(n, p, p.true[[j]], m[j], delta)
      sigma <- calsigma(Rsquared[k], trueidx[[j]], TrueBeta[[j]], type = "group", ngroup=m[j], delta=delta)
      y.train <- gen_y(xmat.train, sigma, TrueBeta[[j]])
      
      for(l in 1:length(lam)){        #l is different lam
        bis_res = bis(xmat.train,y.train,lam=lam[l], w = 0.1, pp = F)
        model.size = length(bis_res$model.pp)
        for(ii in 1:length(logw)){      ##ii is different w
          model.pp = bis_res$postprobs + (0:model.size)*(logw[ii] - logw[length(logw)]) #model.pp for current w
          current.model.size = ifelse(any(diff(model.pp)<0), which(diff(model.pp)<0)[1], model.size)    #First index less than 0
          tpr_bis = sum(bis_res$model.pp[1:current.model.size] %in% trueidx[[j]]) #number of true idx
          acc_bis = as.numeric(tpr_bis==p.true[j])
          res <- c(Rsquared[k], m[j], lam[l], w[ii], acc_bis, current.model.size, tpr_bis)
          cat(res, '\n', file = paste0("iter",i), append = TRUE)
          #cat(length(acc_bis),'\n',file = paste0("res",i), append = TRUE)
        }
      }
    }
  }
  1
}
