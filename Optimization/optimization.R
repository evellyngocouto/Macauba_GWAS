rm(list=ls())

# Packages ----------------------------------------------------------------
library(tidyverse)
library(asreml)
library(ASRgenomics)
library(ggpubr)
library(TrainSel)
library(TSDFGS)
library(parallel)
library(doParallel)
library(foreach)
library(FactoMineR)
library(gghighlight)
library(tidytext)
library(ggtree)
library(ggpattern)
library(corehunter)

# Phenotypic data ---------------------------------------------------------
data = read.table(file = 'Data/pheno/raw_data.txt', h = T)
str(data)
data$Plant[which(data$Plant %in% 1:9)] = paste0('EC0', 1:9)
data$Plant[-which(data$Plant %in% paste0('EC0', 1:9))] = paste0('EC', data$Plant[-which(data$Plant %in% paste0('EC0', 1:9))])
dat = data[which(!is.na(data$OC) & !is.na(data$S_T) & !is.na(data$S_Polpa)),]
dat$Local[which(dat$Local == 4)] = 3

# Genotypic data ----------------------------------------------------------
## Data calling -------------
genodende = readRDS(file = 'Data/geno/genodende.RDS')
genodenovo = readRDS(file = 'Data/geno/genodenovo.RDS')
genotrans = readRDS(file = 'Data/geno/genotrans.RDS')

validgen = unique(dat$Plant)[which(unique(dat$Plant) %in% rownames(genodende) & 
                                     unique(dat$Plant) %in% rownames(genodenovo) & 
                                     unique(dat$Plant) %in% rownames(genotrans))]
genodende =  genodende[rownames(genodende) %in% validgen, ]

G.dende = G.matrix(M = genodende)$G

dat = dat[which(dat$Plant %in% validgen), ]

## PCA --------------
pca.dende = PCA(X = G.dende, scale.unit = F)

pc.df = data.frame(
  Plant = names(pca.dende$ind$coord[,1]),
  pc1 = pca.dende$ind$coord[,1],
  pc2 = pca.dende$ind$coord[,2], row.names = NULL
)  

rm(genodenovo, genotrans, validgen, data)

genodende.cent = scale(genodende, center = TRUE, scale = FALSE)
genodende[1:5, 1:5]
genodende.cent[1:5, 1:5]
PCmat = genodende.cent %*% svd(genodende.cent)$v
dim(genodende); dim(PCmat)

# Optimization ------------------------------------------------------------
TSC<-TrainSelControl()
TSC$mc.cores = 3
TSC$tolconv = 1e-5

## 50 gen: Un-targeted optimization ---------
### D-optimality --------------
dataDopt<-list(FeatureMat=PCmat)
DOPT<-function(soln, Data){
  Fmat<-Data[["FeatureMat"]]
  return(determinant(crossprod(Fmat[soln,]), logarithm=TRUE)$modulus)
}

cl = makeCluster(detectCores() - 1)
clusterEvalQ(cl, {library(TrainSel)})
clusterExport(cl, varlist = c('DOPT', 'dataDopt', 'PCmat', 'TSC'))
registerDoParallel(cl)
set.seed(500)
dopt = foreach(i = 1:50) %dopar% {
  
  TrainSel(Data=dataDopt,
           Candidates = list(1:nrow(PCmat)),
           setsizes = c(50),
           settypes = "UOS",
           Stat = DOPT, 
           control=TSC)
  
}
stopCluster(cl)

### CDMin -------
datacdmin = MakeTrainSelData(M = PCmat, K = G.dende)

cl = makeCluster(detectCores() - 1)
clusterEvalQ(cl, {library(TrainSel)})
clusterExport(cl, varlist = c('datacdmin', 'PCmat', 'TSC'))
registerDoParallel(cl)
set.seed(500)
cdmin = foreach(i = 1:50) %dopar% {
  
  TrainSel(Data=datacdmin,
           Candidates = list(1:nrow(PCmat)),
           setsizes = c(50),
           settypes = "UOS",
           Stat = NULL,
           control = TSC)
  
}
stopCluster(cl)

### Maximin criterion ----------
dataMaximin = list(DistMat = as.matrix(dist(genodende.cent)))
MaximinOPT = function(soln, Data){
  Dsoln = Data[['DistMat']][soln, soln]
  DsolnVec = Dsoln[lower.tri(Dsoln, diag = FALSE)]
  return(min(DsolnVec))
}

cl = makeCluster(detectCores() - 1)
clusterEvalQ(cl, {library(TrainSel)})
clusterExport(cl, varlist = c('dataMaximin','MaximinOPT', 'genodende', 'TSC'))
registerDoParallel(cl)
set.seed(500)
maxmin = foreach(i = 1:50) %dopar% {
  
  TrainSel(Data=dataMaximin,
           Candidates = list(1:nrow(genodende)),
           setsizes = c(50),
           settypes = "UOS",
           Stat = MaximinOPT, 
           control=TSC)
  
}
stopCluster(cl)

### r-score criterion ----------
cl = makeCluster(detectCores() - 1)
clusterEvalQ(cl, {library(TSDFGS)})
clusterExport(cl, varlist = c('PCmat'))
registerDoParallel(cl)
set.seed(500)
rscore = foreach(i = 1:50) %dopar% {
  
  optTrain(PCmat, cand = 1:nrow(PCmat), n.train = 50, method = 'rScore')
  
}

### PEV ridge (Akdemir 2015) criterion ----------
cl = makeCluster(detectCores() - 1)
clusterEvalQ(cl, {library(TSDFGS)})
clusterExport(cl, varlist = c('PCmat'))
registerDoParallel(cl)
set.seed(500)
pev = foreach(i = 1:50) %dopar% {
  
  optTrain(PCmat, cand = 1:nrow(PCmat), n.train = 50, method = 'PEV')
  
}

### CD score (Laloë 1996) criterion ----------
cl = makeCluster(detectCores() - 1)
clusterEvalQ(cl, {library(TSDFGS)})
clusterExport(cl, varlist = c('PCmat'))
registerDoParallel(cl)
set.seed(500)
cdscore = foreach(i = 1:50) %dopar% {
  
  optTrain(PCmat, cand = 1:nrow(PCmat), n.train = 50, method = 'CD')
  
}

#### Checking the convergence
all(do.call(c, lapply(dopt, function(x){
  x$maxvec[length(x$maxvec)] - x$maxvec[length(x$maxvec)-1] 
})) == 0)
all(do.call(c, lapply(cdmin, function(x){
  x$maxvec[length(x$maxvec)] - x$maxvec[length(x$maxvec)-1] 
})) == 0)
all(do.call(c, lapply(maxmin, function(x){
  x$maxvec[length(x$maxvec)] - x$maxvec[length(x$maxvec)-1] 
})) == 0)
all(do.call(c, lapply(rscore, function(x){
  x$TOPscore[length(x$TOPscore)] - x$TOPscore[length(x$TOPscore)-1] 
})) == 0)
all(do.call(c, lapply(pev, function(x){
  x$TOPscore[length(x$TOPscore)] - x$TOPscore[length(x$TOPscore)-1] 
})) == 0)
all(do.call(c, lapply(cdscore, function(x){
  x$TOPscore[length(x$TOPscore)] - x$TOPscore[length(x$TOPscore)-1] 
})) == 0)

#### Extracting the optimized training set from each criterion
trt.dopt = data.frame(table(do.call(cbind, lapply(dopt, function(x){
  rownames(genodende)[x$BestSol_int]
})))) |> arrange(desc(Freq)) |>
  mutate(crit = 'D-opt')
trt.cdmin = data.frame(table(do.call(cbind, lapply(cdmin, function(x){
  rownames(genodende)[x$BestSol_int]
})))) |> arrange(desc(Freq)) |>
  mutate(crit = 'CDmin')
trt.maxmin = data.frame(table(do.call(cbind, lapply(maxmin, function(x){
  rownames(genodende)[x$BestSol_int]
}))))  |> arrange(desc(Freq)) |>
  mutate(crit = 'MaxiMin')
trt.rscore = data.frame(table(do.call(cbind, lapply(rscore, function(x){
  rownames(genodende)[x$OPTtrain]
}))))  |> arrange(desc(Freq)) |>
  mutate(crit = 'r-score')
trt.cdscore = data.frame(table(do.call(cbind, lapply(cdscore, function(x){
  rownames(genodende)[x$OPTtrain]
})))) |> arrange(desc(Freq)) |>
  mutate(crit = 'CDmean')
trt.pev = data.frame(table(do.call(cbind, lapply(pev, function(x){
  rownames(genodende)[x$OPTtrain]
})))) |> arrange(desc(Freq)) |>
  mutate(crit = 'PEV')

pc.df.ut = rbind(
  left_join(pc.df, trt.dopt, by = c('Plant' = 'Var1')) |> 
    mutate(Freq = case_when(is.na(Freq) ~ 0, .default = as.numeric(Freq)),
           crit =  unique(na.exclude(.data$crit))),
  left_join(pc.df, trt.cdscore, by = c('Plant' = 'Var1')) |> 
    mutate(Freq = case_when(is.na(Freq) ~ 0, .default = as.numeric(Freq)),
           crit =  unique(na.exclude(.data$crit))),
  left_join(pc.df, trt.cdmin, by = c('Plant' = 'Var1')) |> 
    mutate(Freq = case_when(is.na(Freq) ~ 0, .default = as.numeric(Freq)),
           crit =  unique(na.exclude(.data$crit))),
  left_join(pc.df, trt.maxmin, by = c('Plant' = 'Var1')) |> 
    mutate(Freq = case_when(is.na(Freq) ~ 0, .default = as.numeric(Freq)),
           crit =  unique(na.exclude(.data$crit))),
  left_join(pc.df, trt.rscore, by = c('Plant' = 'Var1')) |> 
    mutate(Freq = case_when(is.na(Freq) ~ 0, .default = as.numeric(Freq)),
           crit =  unique(na.exclude(.data$crit))),
  left_join(pc.df, trt.pev, by = c('Plant' = 'Var1')) |> 
    mutate(Freq = case_when(is.na(Freq) ~ 0, .default = as.numeric(Freq)),
           crit =  unique(na.exclude(.data$crit)))
)  |> 
  mutate(smooth = Freq/max(Freq))

pc.df.ut = left_join(pc.df.ut, left_join(pc.df.ut, pc.df.ut |> group_by(crit) |> 
                                           slice_max(order_by = Freq, n = 50) |> 
                                           mutate(cand = 1) |> select(Plant, cand)) |> 
                       mutate(cand = ifelse(is.na(cand), 0, 1)) |> 
                       reframe(sel = sum(cand), .by = Plant))

pc.df.ut = pc.df.ut |> full_join(do.call(rbind, lapply(split(pc.df.ut, f = pc.df.ut$crit), function(x){
  data.frame(Plant = x[order(x$Freq, x$sel, decreasing = TRUE), 'Plant'][1:50], 
             trt = 1)
})) |> rownames_to_column('crit') |>
  mutate_at('crit', str_replace, '[.].*',''), by = c('Plant', 'crit')) |> 
  mutate(trt = ifelse(is.na(trt), 0, 1))

rm(trt.cdmin, trt.cdscore, trt.dopt, trt.maxmin, trt.pev, trt.rscore, 
   dopt, maxmin, pc.df, pev, rscore, cdmin, cdscore)

### Plots
pc.df.ut |> 
  ggplot(aes(reorder_within(Plant, -Freq, crit), Freq)) + 
  geom_bar(stat = 'identity', fill = 'steelblue') + theme_bw() + 
  theme(axis.text.x = element_blank()) + 
  scale_x_reordered() +
  facet_wrap(.~crit, scales = 'free_x') + 
  labs(x = 'Candidates', y = 'Runs')

ggplot(data = pc.df.ut, aes(x = pc1, y = pc2)) +
  geom_point(aes(colour = smooth), size = 2) + 
  xlim(-1, .9) + ylim(-1, 0.9) + 
  labs(x = paste0('PC1 (', round(pca.dende$eig[1,2],2), '%)'),
       y = paste0('PC2 (', round(pca.dende$eig[2,2],2), '%)'), 
       colour = 'Frequency') + 
  scale_colour_gradient2(low = '#2b83ba', mid = '#ffffbf', high = '#d7191c',
                         midpoint = .5) + 
  #scale_colour_viridis_c(option = 'inferno', direction = 1) + 
  theme_bw() + theme(legend.position = 'bottom') + 
  facet_wrap(.~crit) + 
  gghighlight(trt == 1,
              keep_scales = TRUE, unhighlighted_params = list(colour = NULL, alpha = .5), 
              calculate_per_facet = TRUE)

ggplot(data = pc.df.ut, aes(y = reorder(Plant, smooth), x = crit, fill = smooth)) + 
  geom_tile() + 
  scale_fill_gradient2(low = '#2b83ba', mid = '#ffffbf', high = '#d7191c',
                       midpoint = .5) + 
  theme_bw() + theme(axis.text.y = element_text(size = 6))+ 
  labs(y = 'Candidates', x = 'Criterion', fill = 'Freq')

aux = pc.df.ut |> reframe(sel = mean(sel), .by = 'Plant') |> 
  mutate(pat = ifelse(sel >= 3, 'Half', 'No'))
aux |> 
  ggplot(aes(x = reorder(Plant, -sel), y = sel)) + 
  geom_bar_pattern(stat = 'identity', aes(pattern_density = pat), width = 1,
                   pattern = 'wave', fill = 'darkred', colour = 'darkred') +
  scale_pattern_density_manual('',values = c(Half=.4, No=0), 
                               labels = c(
                                 No = paste0(round(length(which(aux$sel >=1)) / nrow(aux)*100,2), '% selected at least once'),
                                 Half = paste0(round(length(which(aux$sel >=3)) / nrow(aux)*100,2), '% selected half of the time') 
                               )) + 
  scale_pattern_fill_manual(values = c('black')) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 6), 
                     legend.position = 'bottom') +
  scale_y_continuous(breaks = seq(6)) + 
  labs(x = 'Candidates', y = 'No. optimization criteria')


### Core collection ---------------
coredata = coreHunterData(genotypes(genodende, format = "biparental"),)
obj<- objective("EN", "MR") 
core <- sampleCore(coredata, obj, size = 50)
selcore = core$sel

core$EN$MR
# diversity value = 0.3115836
# the average of the Modified Rogers distance between each 
# selected accession and the most similar other selected accession
evaluateCore(core, coredata, objective("CV")) # Allele coverage: 99%
evaluateCore(core, coredata, objective("HE")) 
#' Expected proportion of heterozygous loci in offspring produced from random crossings 
#' within the selected core: 31%


### Cross-validation with the optimized training set ------------------
G.dende.inv = G.inverse(G = G.tuneup(G = G.dende, bend = TRUE)$Gb, sparseform = TRUE)$Ginv

dat = dat |> select(Plant, Ano, Local, OC, S_T, S_Polpa) |> 
  mutate(Plant = as.factor(Plant), Ano = as.factor(Ano), Local = as.factor(Local))


### Oil content -------------
mod = asreml(fixed = OC ~ Ano + Plant, 
             data = dat, 
             na.action = na.method(x = 'exclude', y = 'exclude'))
input = predict(mod, classify = 'Plant')$pvals[,1:2] |> 
  rename(ybar = predicted.value)

rm(mod)

##### Cross-validation (random sampling)
fold.size = 50 
nrept = 100 

cv = list()
k = 1
i = 1
repeat
{
  set.seed(1 + k * 3)
  
  cvdata = cbind(input, trt = sample(c(rep(1, fold.size), rep(0, nlevels(input$Plant) - fold.size)))) |> 
    mutate(yNA = ybar) |> mutate(yNA = ifelse(trt != 1, NA, yNA))
  mod1 = tryCatch({
    asreml(fixed = yNA ~ 1, 
           random = ~ vm(Plant, G.dende.inv), 
           data = cvdata)
  }, error = function(e){cat("There was a singularity issue, sampling again...", fill = TRUE)})
  if(class(mod1) != 'asreml'){
    rm(mod1, cvdata)
    k = k + 1
    next
  }else{
    if(mod1$converge){
      res1 = as.data.frame(coef(mod1)$random + c(coef(mod1)$fixed)) |> 
        rownames_to_column('Plant') |>
        separate(col = 'Plant', into = c('aaa','Plant'), sep = '_') |>
        select(-aaa) |> right_join(cvdata, by = 'Plant') |>
        rename(yhat = effect) |>
        select(Plant, ybar, yhat, trt) |>
        filter(trt != 1)
      rm(mod1)
    }else{
      cat("The model did not converge, sampling again...", fill = TRUE)
      k = k+1
      rm(mod1)
      next
    }
  }
  
  cv[[i]] = res1
  cv[[i]]$seed = 1 + k * 3
  
  rm(res1, cvdata)
  
  k = k + 1
  message('Succeded! Repetition ',i)
  i = i + 1
  if(length(cv) == nrept) break 
}

corr = lapply(cv, function(x) cor(x$ybar, x$yhat, use = 'complete.obs'))
mspe = lapply(cv, function(x) mean((x$ybar - x$yhat)^2, na.rm = T))

rm(cv, k, i, fold.size, nrept)

##### Cross-validation (optimized training set)
cv.opt = lapply(split(pc.df.ut, f = pc.df.ut$crit), function(x){
  cvdata = left_join(input, x[, c("Plant", 'trt')], by = 'Plant')
  cvdata$yNA = cvdata$ybar; cvdata$Plant = as.factor(cvdata$Plant)
  cvdata[which(cvdata$trt != 1),'yNA'] = NA
  
  asreml.options(ai.sing = TRUE)
  mod1 = asreml(fixed = yNA ~ 1, 
                random = ~ vm(Plant, G.dende.inv), 
                data = cvdata)
  asreml.options(ai.sing = FALSE)
  
  as.data.frame(coef(mod1)$random + c(coef(mod1)$fixed)) |> 
    rownames_to_column('Plant') |>
    separate(col = 'Plant', into = c('aaa','Plant'), sep = '_') |>
    select(-aaa) |> right_join(cvdata, by = 'Plant') |>
    rename(yhat = effect) |>
    select(Plant, ybar, yhat, trt) |>
    filter(trt != 1) |> 
    reframe(corr = cor(ybar, yhat), 
            mspe = mean((yhat-ybar)^2))
})

res.un.oc = rbind(do.call(rbind, cv.opt) |> rownames_to_column('crit'),
                  data.frame(crit = 'Random', 
                             corr = do.call(rbind, corr),
                             mspe = do.call(rbind, mspe)))|> 
  pivot_longer(corr:mspe) |> mutate(val = 'val')

##### Cross-validation (core collection)
cvdata = input; cvdata$yNA = cvdata$ybar
cvdata[which(!cvdata$Plant %in% selcore),'yNA'] = NA
cvdata$Plant = as.factor(cvdata$Plant)

mod1 = asreml(fixed = yNA ~ 1, 
              random = ~ vm(Plant, G.dende.inv), 
              data = cvdata)

res.un.oc = rbind(res.un.oc, 
                  as.data.frame(coef(mod1)$random + c(coef(mod1)$fixed)) |> 
                    rownames_to_column('Plant') |>
                    separate(col = 'Plant', into = c('aaa','Plant'), sep = '_') |>
                    select(-aaa) |> right_join(cvdata, by = 'Plant') |>
                    rename(yhat = effect) |>
                    select(Plant, ybar, yhat) |>
                    filter(!Plant %in% selcore) |> 
                    reframe(corr = cor(ybar, yhat), 
                            mspe = mean((yhat-ybar)^2)) |> 
                    mutate(crit = 'Core') |> pivot_longer(corr:mspe) |>
                    mutate(val = 'val'))

facet.label = c('Correlation', 'MSPE')
names(facet.label) = unique(res.un.oc$name)

res.un.oc |> 
  ggplot(aes(y = value, x = val)) + 
  geom_violin(aes(fill = name), color = 'black', show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.1), alpha = .6, color = 'darkgrey') + 
  geom_point(data = subset(res.un.oc, subset = crit != 'Random'), 
             aes(y = value, x = val, color = crit), size = 3, 
             pch = 17) +
  facet_wrap(.~name, scales = 'free_y',labeller = labeller(.cols = facet.label)) + 
  scale_fill_manual(values = c('#edf8fb','#fef0d9')) + 
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3', '#e7298a',
                                '#a6761d','#e6ab02','#e41a1c')) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.position = 'bottom') +
  labs(x = '', y = '', color = 'Criterion')

rm(cv.opt, corr, mspe)

### Fruit dry mass -------------
mod = asreml(fixed = S_T ~ Ano + Plant, 
             data = dat, 
             na.action = na.method(x = 'exclude', y = 'exclude'))
input = predict(mod, classify = 'Plant')$pvals[,1:2] |> 
  rename(ybar = predicted.value)

rm(mod)

##### Cross-validation (random sampling)
fold.size = 50 
nrept = 100 

cv = list()
k = 1
i = 1
repeat
{
  set.seed(1 + k * 3)
  
  cvdata = cbind(input, trt = sample(c(rep(1, fold.size), rep(0, nlevels(input$Plant) - fold.size)))) |> 
    mutate(yNA = ybar) |> mutate(yNA = ifelse(trt != 1, NA, yNA))
  mod1 = tryCatch({
    asreml(fixed = yNA ~ 1, 
           random = ~ vm(Plant, G.dende.inv), 
           data = cvdata)
  }, error = function(e){cat("There was a singularity issue, sampling again...", fill = TRUE)})
  if(class(mod1) != 'asreml'){
    rm(mod1, cvdata)
    k = k + 1
    next
  }else{
    if(mod1$converge){
      res1 = as.data.frame(coef(mod1)$random + c(coef(mod1)$fixed)) |> 
        rownames_to_column('Plant') |>
        separate(col = 'Plant', into = c('aaa','Plant'), sep = '_') |>
        select(-aaa) |> right_join(cvdata, by = 'Plant') |>
        rename(yhat = effect) |>
        select(Plant, ybar, yhat, trt) |>
        filter(trt != 1)
      rm(mod1)
    }else{
      cat("The model did not converge, sampling again...", fill = TRUE)
      k = k+1
      rm(mod1)
      next
    }
  }
  
  cv[[i]] = res1
  cv[[i]]$seed = 1 + k * 3
  
  rm(res1, cvdata)
  
  k = k + 1
  message('Succeded! Repetition ',i)
  i = i + 1
  if(length(cv) == nrept) break 
}

corr = lapply(cv, function(x) cor(x$ybar, x$yhat, use = 'complete.obs'))
mspe = lapply(cv, function(x) mean((x$ybar - x$yhat)^2, na.rm = T))

rm(cv, k, i, fold.size, nrept)

##### Cross-validation (optimized training set)

cv.opt = lapply(split(pc.df.ut, f = pc.df.ut$crit), function(x){
  cvdata = left_join(input, x[, c("Plant", 'trt')], by = 'Plant')
  cvdata$yNA = cvdata$ybar; cvdata$Plant = as.factor(cvdata$Plant)
  cvdata[which(cvdata$trt != 1),'yNA'] = NA
  
  asreml.options(ai.sing = TRUE)
  mod1 = asreml(fixed = yNA ~ 1, 
                random = ~ vm(Plant, G.dende.inv), 
                data = cvdata)
  asreml.options(ai.sing = FALSE)
  
  as.data.frame(coef(mod1)$random + c(coef(mod1)$fixed)) |> 
    rownames_to_column('Plant') |>
    separate(col = 'Plant', into = c('aaa','Plant'), sep = '_') |>
    select(-aaa) |> right_join(cvdata, by = 'Plant') |>
    rename(yhat = effect) |>
    select(Plant, ybar, yhat, trt) |>
    filter(trt != 1) |> 
    reframe(corr = cor(ybar, yhat), 
            mspe = mean((yhat-ybar)^2))
})

res.un.fdm = rbind(do.call(rbind, cv.opt) |> rownames_to_column('crit'),
                   data.frame(crit = 'Random', 
                              corr = do.call(rbind, corr),
                              mspe = do.call(rbind, mspe)))|> 
  pivot_longer(corr:mspe) |> mutate(val = 'val')

##### Cross-validation (core collection)
cvdata = input; cvdata$yNA = cvdata$ybar
cvdata[which(!cvdata$Plant %in% selcore),'yNA'] = NA
cvdata$Plant = as.factor(cvdata$Plant)

asreml.options(ai.sing = TRUE)
mod1 = asreml(fixed = yNA ~ 1, 
              random = ~ vm(Plant, G.dende.inv), 
              data = cvdata)
asreml.options(ai.sing = FALSE)
res.un.fdm = rbind(res.un.fdm, 
                   as.data.frame(coef(mod1)$random + c(coef(mod1)$fixed)) |> 
                     rownames_to_column('Plant') |>
                     separate(col = 'Plant', into = c('aaa','Plant'), sep = '_') |>
                     select(-aaa) |> right_join(cvdata, by = 'Plant') |>
                     rename(yhat = effect) |>
                     select(Plant, ybar, yhat) |>
                     filter(!Plant %in% selcore) |> 
                     reframe(corr = cor(ybar, yhat), 
                             mspe = mean((yhat-ybar)^2)) |> 
                     mutate(crit = 'Core') |> pivot_longer(corr:mspe) |>
                     mutate(val = 'val'))


facet.label = c('Correlation', 'MSPE')
names(facet.label) = unique(res.un.oc$name)

res.un.fdm |> 
  ggplot(aes(y = value, x = val)) + 
  geom_violin(aes(fill = name), color = 'black', show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.1), alpha = .6, color = 'darkgrey') + 
  geom_point(data = subset(res.un.fdm, subset = crit != 'Random'), 
             aes(y = value, x = val, color = crit), size = 3, 
             pch = 17) +
  facet_wrap(.~name, scales = 'free_y',labeller = labeller(.cols = facet.label)) + 
  scale_fill_manual(values = c('#edf8fb','#fef0d9')) + 
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3', '#e7298a',
                                '#a6761d','#e6ab02','#e41a1c')) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.position = 'bottom') +
  labs(x = '', y = '', color = 'Criterion')

rm(cv.opt, corr, mspe)

### Pulp Dry mass -------------
mod = asreml(fixed = S_Polpa ~ Ano + Plant, 
             data = dat, 
             na.action = na.method(x = 'exclude', y = 'exclude'))
input = predict(mod, classify = 'Plant')$pvals[,1:2] |> 
  rename(ybar = predicted.value)

rm(mod)

##### Cross-validation (random sampling)
fold.size = 50 
nrept = 100 

cv = list()
k = 1
i = 1
repeat
{
  set.seed(1 + k * 3)
  
  cvdata = cbind(input, trt = sample(c(rep(1, fold.size), rep(0, nlevels(input$Plant) - fold.size)))) |> 
    mutate(yNA = ybar) |> mutate(yNA = ifelse(trt != 1, NA, yNA))
  mod1 = tryCatch({
    asreml(fixed = yNA ~ 1, 
           random = ~ vm(Plant, G.dende.inv), 
           data = cvdata)
  }, error = function(e){cat("There was a singularity issue, sampling again...", fill = TRUE)})
  if(class(mod1) != 'asreml'){
    rm(mod1, cvdata)
    k = k + 1
    next
  }else{
    if(mod1$converge){
      res1 = as.data.frame(coef(mod1)$random + c(coef(mod1)$fixed)) |> 
        rownames_to_column('Plant') |>
        separate(col = 'Plant', into = c('aaa','Plant'), sep = '_') |>
        select(-aaa) |> right_join(cvdata, by = 'Plant') |>
        rename(yhat = effect) |>
        select(Plant, ybar, yhat, trt) |>
        filter(trt != 1)
      rm(mod1)
    }else{
      cat("The model did not converge, sampling again...", fill = TRUE)
      k = k+1
      rm(mod1)
      next
    }
  }
  
  cv[[i]] = res1
  cv[[i]]$seed = 1 + k * 3
  
  rm(res1, cvdata)
  
  k = k + 1
  message('Succeded! Repetition ',i)
  i = i + 1
  if(length(cv) == nrept) break 
}

corr = lapply(cv, function(x) cor(x$ybar, x$yhat, use = 'complete.obs'))
mspe = lapply(cv, function(x) mean((x$ybar - x$yhat)^2, na.rm = T))

rm(cv, k, i, fold.size, nrept)

##### Cross-validation (optimized training set)

cv.opt = lapply(split(pc.df.ut, f = pc.df.ut$crit), function(x){
  cvdata = left_join(input, x[, c("Plant", 'trt')], by = 'Plant')
  cvdata$yNA = cvdata$ybar; cvdata$Plant = as.factor(cvdata$Plant)
  cvdata[which(cvdata$trt != 1),'yNA'] = NA
  
  asreml.options(ai.sing = TRUE)
  mod1 = asreml(fixed = yNA ~ 1, 
                random = ~ vm(Plant, G.dende.inv), 
                data = cvdata)
  asreml.options(ai.sing = FALSE)
  
  as.data.frame(coef(mod1)$random + c(coef(mod1)$fixed)) |> 
    rownames_to_column('Plant') |>
    separate(col = 'Plant', into = c('aaa','Plant'), sep = '_') |>
    select(-aaa) |> right_join(cvdata, by = 'Plant') |>
    rename(yhat = effect) |>
    select(Plant, ybar, yhat, trt) |>
    filter(trt != 1) |> 
    reframe(corr = cor(ybar, yhat), 
            mspe = mean((yhat-ybar)^2))
})

res.un.pdm = rbind(do.call(rbind, cv.opt) |> rownames_to_column('crit'),
                   data.frame(crit = 'Random', 
                              corr = do.call(rbind, corr),
                              mspe = do.call(rbind, mspe)))|> 
  pivot_longer(corr:mspe) |> mutate(val = 'val')

##### Cross-validation (core collection)
cvdata = input; cvdata$yNA = cvdata$ybar
cvdata[which(!cvdata$Plant %in% selcore),'yNA'] = NA
cvdata$Plant = as.factor(cvdata$Plant)

asreml.options(ai.sing = TRUE)
mod1 = asreml(fixed = yNA ~ 1, 
              random = ~ vm(Plant, G.dende.inv), 
              data = cvdata)
asreml.options(ai.sing = FALSE)

res.un.pdm = rbind(res.un.pdm, 
                   as.data.frame(coef(mod1)$random + c(coef(mod1)$fixed)) |> 
                     rownames_to_column('Plant') |>
                     separate(col = 'Plant', into = c('aaa','Plant'), sep = '_') |>
                     select(-aaa) |> right_join(cvdata, by = 'Plant') |>
                     rename(yhat = effect) |>
                     select(Plant, ybar, yhat) |>
                     filter(!Plant %in% selcore) |> 
                     reframe(corr = cor(ybar, yhat), 
                             mspe = mean((yhat-ybar)^2)) |> 
                     mutate(crit = 'Core') |> pivot_longer(corr:mspe) |>
                     mutate(val = 'val'))


facet.label = c('Correlation', 'MSPE')
names(facet.label) = unique(res.un.pdm$name)

res.un.pdm |> 
  ggplot(aes(y = value, x = val)) + 
  geom_violin(aes(fill = name), color = 'black', show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.1), alpha = .6, color = 'darkgrey') + 
  geom_point(data = subset(res.un.pdm, subset = crit != 'Random'), 
             aes(y = value, x = val, color = crit), size = 3, 
             pch = 17) +
  facet_wrap(.~name, scales = 'free_y',labeller = labeller(.cols = facet.label)) + 
  scale_fill_manual(values = c('#edf8fb','#fef0d9')) + 
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3', '#e7298a',
                                '#a6761d','#e6ab02','#e41a1c')) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.position = 'bottom') +
  labs(x = '', y = '', color = 'Criterion')

rm(cv.opt, corr, mspe)

### Wrap-up (CV) ----------
res.un = rbind(res.un.oc |> mutate(trait = 'OC'), 
               res.un.fdm |> mutate(trait = 'FDM'), 
               res.un.pdm |> mutate(trait = 'PDM')) 

res.un |> 
  ggplot(aes(y = value, x = val)) + 
  geom_violin(aes(fill = name), color = 'black', show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.1), alpha = .6, color = 'darkgrey') + 
  geom_point(data = subset(res.un, subset = crit != 'Random'), 
             aes(y = value, x = val, color = crit), size = 3, 
             pch = 17) +
  facet_wrap(name~trait, scales = 'free', ncol = 3) + 
  scale_fill_manual(values = c('#edf8fb','#fef0d9')) + 
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3', '#e7298a',
                                '#a6761d','#e6ab02','#e41a1c')) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.position = 'bottom') +
  labs(x = '', y = '', color = 'Criterion')


## Targeted optimization (1 & 2 -> 3) ---------
dat = dat |> select(Plant, Ano, Local, OC, S_T, S_Polpa) |> 
  mutate(Plant = as.factor(Plant), Ano = as.factor(Ano), Local = as.factor(Local))

target = droplevels(tapply(dat$Plant, dat$Local, unique)[[3]])

#' Plants from Locations 1 and 2 will compose the training set. The goal is to 
#' to select the best composition for predicting location 3. 

### PEV -------------
cl = makeCluster(detectCores() - 1)
clusterEvalQ(cl, {library(TSDFGS)})
clusterExport(cl, varlist = c('PCmat', 'target'))
registerDoParallel(cl)
set.seed(500)
pev = foreach(i = 1:50) %dopar% {
  
  optTrain(PCmat, cand = which(!rownames(PCmat) %in% target), 
           test = which(rownames(PCmat) %in% target),
           n.train = 50, method = 'PEV')
  
}

stopCluster(cl)

### CDmean -------------
cl = makeCluster(detectCores() - 1)
clusterEvalQ(cl, {library(TSDFGS)})
clusterExport(cl, varlist = c('PCmat', 'target'))
registerDoParallel(cl)
set.seed(500)
cdscore = foreach(i = 1:50) %dopar% {
  
  optTrain(PCmat, cand = which(!rownames(PCmat) %in% target), 
           test = which(rownames(PCmat) %in% target),
           n.train = 50, method = 'CD')
  
}

stopCluster(cl)

### r-score -------------
cl = makeCluster(detectCores() - 1)
clusterEvalQ(cl, {library(TSDFGS)})
clusterExport(cl, varlist = c('PCmat', 'target'))
registerDoParallel(cl)
set.seed(500)
rscore = foreach(i = 1:50) %dopar% {
  
  optTrain(PCmat, cand = which(!rownames(PCmat) %in% target), 
           test = which(rownames(PCmat) %in% target),
           n.train = 50, method = 'CD')
  
}

stopCluster(cl)

### Minimax -------------
dataMiniMax<-list(DistMat=as.matrix(dist(genodende.cent)))
MiniMaxOPT<-function(soln, Data){
  Dsoln<-Data[["DistMat"]][soln, which(rownames(Data[["DistMat"]]) %in% target)]
  DsolnVec<-max(c(unlist(Dsoln)))
  return(-(DsolnVec))
}

minimax<-TrainSel(Data=dataMiniMax,
                  Candidates = list(which(!rownames(dataMiniMax[["DistMat"]]) %in% target)),
                  setsizes = c(50),
                  settypes = "UOS",
                  Stat = MiniMaxOPT,
                  control=TSC)

### CDMin -------
datacdmin = MakeTrainSelData(M = PCmat, K = G.dende)

cl = makeCluster(detectCores() - 1)
clusterEvalQ(cl, {library(TrainSel)})
clusterExport(cl, varlist = c('datacdmin', 'TSC', 'target'))
registerDoParallel(cl)
set.seed(500)
cdmin = foreach(i = 1:50) %dopar% {
  
  TrainSel(Data=datacdmin,
           Candidates = list(which(!datacdmin$labels$names %in% target)),
           setsizes = c(50),
           settypes = "UOS",
           Stat = NULL,
           Target = which(datacdmin$labels$names %in% target),
           control = TSC)
  
}
stopCluster(cl)

### Multiple design criterion -------------
dataMultOpt<-list(DistMat=as.matrix(dist(genodende.cent)))

MultOPT<-function(soln, Data){
  D<-Data[["DistMat"]]
  Dsoln<-D[soln, which(rownames(Data[["DistMat"]]) %in% target)]
  DsolnVec1<- -mean(c(unlist(Dsoln)))
  Dsoln2<-D[soln,soln]
  DsolnVec2<-mean(c(unlist(Dsoln2)))
  return(c(DsolnVec2,DsolnVec1))
}

cl = makeCluster(detectCores() - 1)
clusterEvalQ(cl, {library(TrainSel)})
clusterExport(cl, varlist = c('dataMultOpt','MultOPT', 'TSC', 'target'))
registerDoParallel(cl)
set.seed(500)

mult = foreach(i = 1:50) %dopar% {
  TSOUTMultOPt<-TrainSel(Data=dataMultOpt,
                         Candidates = list(which(!rownames(dataMultOpt[["DistMat"]]) %in% target)),
                         setsizes = c(50),
                         settypes = "UOS",
                         Stat = MultOPT, 
                         nStat=2, 
                         control=TSC)
}
stopCluster(cl)

### Wrap-up ----------

#### Checking the convergence
all(do.call(c, lapply(cdmin, function(x){
  x$maxvec[length(x$maxvec)] - x$maxvec[length(x$maxvec)-1] 
})) == 0)
all(do.call(c, lapply(rscore, function(x){
  x$TOPscore[length(x$TOPscore)] - x$TOPscore[length(x$TOPscore)-1] 
})) == 0)
all(do.call(c, lapply(pev, function(x){
  x$TOPscore[length(x$TOPscore)] - x$TOPscore[length(x$TOPscore)-1] 
})) == 0)
all(do.call(c, lapply(cdscore, function(x){
  x$TOPscore[length(x$TOPscore)] - x$TOPscore[length(x$TOPscore)-1] 
})) == 0)

minimax$maxvec[length(minimax$maxvec)] - minimax$maxvec[length(minimax$maxvec)-1] == 0

#### Extracting the optimized training set from each criterion
trt.cdmin = data.frame(table(do.call(cbind, lapply(cdmin, function(x){
  rownames(genodende)[x$BestSol_int]
})))) |> arrange(desc(Freq)) |> mutate(crit = 'CDmin')
trt.minimax = data.frame(Var1 = rownames(genodende)[minimax$BestSol_int], 
                         Freq = 50, crit = 'MiniMax')
trt.rscore = data.frame(table(do.call(cbind, lapply(rscore, function(x){
  rownames(genodende)[x$OPTtrain]
}))))  |> arrange(desc(Freq)) |> mutate(crit = 'r-score')
trt.cdscore = data.frame(table(do.call(cbind, lapply(cdscore, function(x){
  rownames(genodende)[x$OPTtrain]
})))) |> arrange(desc(Freq)) |> mutate(crit = 'CDmean')
trt.pev = data.frame(table(do.call(cbind, lapply(pev, function(x){
  rownames(genodende)[x$OPTtrain]
})))) |> arrange(desc(Freq)) |> mutate(crit = 'PEV')
trt.mult = data.frame(table(do.call(cbind, lapply(mult, function(x){
  rownames(genodende)[x$BestSol_int[, sample(which((t(x$BestVal)[,1] >= 64) & 
                                                     (t(x$BestVal)[,2] >= -65.5)),1)]]
})))) |> arrange(desc(Freq)) |> mutate(crit = 'Mult')


pc.df.ut = rbind(
  left_join(pc.df |> filter(!Plant %in% target), trt.cdscore, by = c('Plant' = 'Var1')) |> 
    mutate(Freq = case_when(is.na(Freq) ~ 0, .default = as.numeric(Freq)),
           crit =  unique(na.exclude(.data$crit))),
  left_join(pc.df |> filter(!Plant %in% target), trt.cdmin, by = c('Plant' = 'Var1')) |> 
    mutate(Freq = case_when(is.na(Freq) ~ 0, .default = as.numeric(Freq)),
           crit =  unique(na.exclude(.data$crit))),
  left_join(pc.df |> filter(!Plant %in% target), trt.minimax, by = c('Plant' = 'Var1')) |> 
    mutate(Freq = case_when(is.na(Freq) ~ 0, .default = as.numeric(Freq)),
           crit =  unique(na.exclude(.data$crit))),
  left_join(pc.df |> filter(!Plant %in% target), trt.rscore, by = c('Plant' = 'Var1')) |> 
    mutate(Freq = case_when(is.na(Freq) ~ 0, .default = as.numeric(Freq)),
           crit =  unique(na.exclude(.data$crit))),
  left_join(pc.df |> filter(!Plant %in% target), trt.mult, by = c('Plant' = 'Var1')) |> 
    mutate(Freq = case_when(is.na(Freq) ~ 0, .default = as.numeric(Freq)),
           crit =  unique(na.exclude(.data$crit))),
  left_join(pc.df |> filter(!Plant %in% target), trt.pev, by = c('Plant' = 'Var1')) |> 
    mutate(Freq = case_when(is.na(Freq) ~ 0, .default = as.numeric(Freq)),
           crit =  unique(na.exclude(.data$crit)))
) |> 
  mutate(smooth = Freq/max(Freq))

pc.df.ut = left_join(pc.df.ut, left_join(pc.df.ut, pc.df.ut |> group_by(crit) |> 
                                           slice_max(order_by = Freq, n = 50) |> 
                                           mutate(cand = 1) |> select(Plant, cand)) |> 
                       mutate(cand = ifelse(is.na(cand), 0, 1)) |> 
                       reframe(sel = sum(cand), .by = Plant))

pc.df.ut = pc.df.ut |> full_join(do.call(rbind, lapply(split(pc.df.ut, f = pc.df.ut$crit), function(x){
  data.frame(Plant = x[order(x$Freq, x$sel, decreasing = TRUE), 'Plant'][1:50], 
             trt = 1)
})) |> rownames_to_column('crit') |>
  mutate_at('crit', str_replace, '[.].*',''), by = c('Plant', 'crit')) |> 
  mutate(trt = ifelse(is.na(trt), 0, 1))

rm(trt.cdmin, trt.cdscore, trt.mult, trt.minimax, trt.pev, trt.rscore, 
   mult, minimax, pev, rscore, cdmin, cdscore)

### Plots

pc.df.ut |> 
  ggplot(aes(reorder_within(Plant, -Freq, crit), Freq)) + 
  geom_bar(stat = 'identity', fill = 'steelblue') + theme_bw() + 
  theme(axis.text.x = element_blank()) + 
  scale_x_reordered() +
  facet_wrap(.~crit, scales = 'free_x') + 
  labs(x = 'Candidates', y = 'Runs')

ggplot() + 
  geom_point(data = pc.df.ut |> filter(trt == 0), 
             aes(x = pc1, y = pc2, colour = smooth, shape = 'Training'), 
             alpha = .5) +
  geom_point(data = pc.df.ut |> filter(trt == 1), 
             aes(x = pc1, y = pc2, colour = smooth, shape = 'Training')) + 
  geom_point(data = pc.df |> filter(Plant %in% target), alpha = 0.7,
             aes(x = pc1, y = pc2, shape = 'Target'), colour = 'darkgreen') +
  xlim(-1, .9) + ylim(-1, 0.9) + 
  facet_wrap(.~crit)+ 
  scale_colour_gradient2(low = '#2b83ba', mid = '#ffffbf', high = '#d7191c',
                         midpoint = .5)   + 
  theme_bw() + theme(legend.position = 'bottom') + 
  labs(x = paste0('PC1 (', round(pca.dende$eig[1,2],2), '%)'),
       y = paste0('PC2 (', round(pca.dende$eig[2,2],2), '%)'), 
       colour = 'Frequency', shape = '') 

ggplot(data = pc.df.ut, aes(y = reorder(Plant, smooth), x = crit, fill = smooth)) + 
  geom_tile() + 
  scale_fill_gradient2(low = '#2b83ba', mid = '#ffffbf', high = '#d7191c',
                       midpoint = .5) + 
  theme_bw() + theme(axis.text.y = element_text(size = 6))+ 
  labs(y = 'Candidates', x = 'Criterion', fill = 'Freq')


aux = pc.df.ut |> reframe(sel = mean(sel), .by = 'Plant') |> 
  mutate(pat = ifelse(sel >= 3, 'Half', 'No'))
aux |> 
  ggplot(aes(x = reorder(Plant, -sel), y = sel)) + 
  geom_bar_pattern(stat = 'identity', aes(pattern_density = pat), width = 1,
                   pattern = 'wave', fill = 'darkred', colour = 'darkred') +
  scale_pattern_density_manual('',values = c(Half=.4, No=0), 
                               labels = c(
                                 No = paste0(round(length(which(aux$sel >=1)) / nrow(aux)*100,2), '% selected at least once'),
                                 Half = paste0(round(length(which(aux$sel >=3)) / nrow(aux)*100,2), '% selected half of the time') 
                               )) + 
  scale_pattern_fill_manual(values = c('black')) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 6), 
                     legend.position = 'bottom') +
  scale_y_continuous(breaks = seq(6)) + 
  labs(x = 'Candidates', y = 'No. optimization criteria')

rm(aux)

### Core collection ---------------
coredata = coreHunterData(genotypes(genodende[!rownames(genodende) %in% target,],
                                    format = "biparental"),)
obj<- objective("EN", "MR") 
core <- sampleCore(coredata, obj, size = 50)
selcore = core$sel

core$EN$MR
# diversity value = 0.3115836
# the average of the Modified Rogers distance between each 
# selected accession and the most similar other selected accession
evaluateCore(core, coredata, objective("CV")) # Allele coverage: 99%
evaluateCore(core, coredata, objective("HE")) 
#' Expected proportion of heterozygous loci in offspring produced from random crossings 
#' within the selected core: 31%

### Cross-validation with the optimized training set ------------------
G.dende.inv = G.inverse(G = G.tuneup(G = G.dende, bend = TRUE)$Gb, sparseform = TRUE)$Ginv

dat = dat |> select(Plant, Ano, Local, OC, S_T, S_Polpa) |> 
  mutate(Plant = as.factor(Plant), Ano = as.factor(Ano), Local = as.factor(Local))

### Oil content -------------
mod = asreml(fixed = OC ~ Ano + Plant, 
             data = dat, 
             na.action = na.method(x = 'exclude', y = 'exclude'))
input = predict(mod, classify = 'Plant')$pvals[,1:2] |> 
  rename(ybar = predicted.value)

train = droplevels(input$Plant[which(!input$Plant %in% target)])

rm(mod)

##### Cross-validation (random sampling)
fold.size = 50 
nrept = 100 

cv = list()
k = 1
i = 1
repeat
{
  set.seed(1 + k * 3)
  
  cvdata = input |> mutate(trt = ifelse(Plant %in% sample(train, fold.size), 1, 0), 
                           yNA = ifelse(trt != 1, NA, ybar))
  mod1 = tryCatch({
    asreml(fixed = yNA ~ 1, 
           random = ~ vm(Plant, G.dende.inv), 
           data = cvdata)
  }, error = function(e){cat("There was a singularity issue, sampling again...", fill = TRUE)})
  if(class(mod1) != 'asreml'){
    rm(mod1, cvdata)
    k = k + 1
    next
  }else{
    if(mod1$converge){
      res1 = as.data.frame(coef(mod1)$random + c(coef(mod1)$fixed)) |> 
        rownames_to_column('Plant') |>
        separate(col = 'Plant', into = c('aaa','Plant'), sep = '_') |>
        select(-aaa) |> right_join(cvdata, by = 'Plant') |>
        rename(yhat = effect) |>
        select(Plant, ybar, yhat, trt) |>
        filter(Plant %in% target)
      rm(mod1)
    }else{
      cat("The model did not converge, sampling again...", fill = TRUE)
      k = k+1
      rm(mod1)
      next
    }
  }
  
  cv[[i]] = res1
  cv[[i]]$seed = 1 + k * 3
  
  rm(res1, cvdata)
  
  k = k + 1
  message('Succeded! Repetition ',i)
  i = i + 1
  if(length(cv) == nrept) break 
}

corr = lapply(cv, function(x) cor(x$ybar, x$yhat, use = 'complete.obs'))
mspe = lapply(cv, function(x) mean((x$ybar - x$yhat)^2, na.rm = T))

rm(cv, k, i, fold.size, nrept)

##### Cross-validation (optimized training set)
cv.opt = lapply(split(pc.df.ut, f = pc.df.ut$crit), function(x){
  cvdata = left_join(input, x[, c("Plant", 'trt')], by = 'Plant')
  cvdata$trt[which(is.na(cvdata$trt))] = 0
  cvdata$yNA = cvdata$ybar; cvdata$Plant = as.factor(cvdata$Plant)
  cvdata[which(cvdata$trt != 1),'yNA'] = NA
  
  asreml.options(ai.sing = TRUE)
  mod1 = asreml(fixed = yNA ~ 1, 
                random = ~ vm(Plant, G.dende.inv), 
                data = cvdata)
  asreml.options(ai.sing = FALSE)
  
  as.data.frame(coef(mod1)$random + c(coef(mod1)$fixed)) |> 
    rownames_to_column('Plant') |>
    separate(col = 'Plant', into = c('aaa','Plant'), sep = '_') |>
    select(-aaa) |> right_join(cvdata, by = 'Plant') |>
    rename(yhat = effect) |>
    select(Plant, ybar, yhat, trt) |>
    filter(Plant %in% target) |> 
    reframe(corr = cor(ybar, yhat), 
            mspe = mean((yhat-ybar)^2))
})

res.tg.oc = rbind(do.call(rbind, cv.opt) |> rownames_to_column('crit'),
                  data.frame(crit = 'Random', 
                             corr = do.call(rbind, corr),
                             mspe = do.call(rbind, mspe)))|> 
  pivot_longer(corr:mspe) |> mutate(val = 'val')

##### Cross-validation (core collection)
cvdata = input; cvdata$yNA = cvdata$ybar
cvdata[which(!cvdata$Plant %in% selcore),'yNA'] = NA
cvdata$Plant = as.factor(cvdata$Plant)

mod1 = asreml(fixed = yNA ~ 1, 
              random = ~ vm(Plant, G.dende.inv), 
              data = cvdata)

res.tg.oc = rbind(res.tg.oc, 
                  as.data.frame(coef(mod1)$random + c(coef(mod1)$fixed)) |> 
                    rownames_to_column('Plant') |>
                    separate(col = 'Plant', into = c('aaa','Plant'), sep = '_') |>
                    select(-aaa) |> right_join(cvdata, by = 'Plant') |>
                    rename(yhat = effect) |>
                    select(Plant, ybar, yhat) |>
                    filter(Plant %in% target) |> 
                    reframe(corr = cor(ybar, yhat), 
                            mspe = mean((yhat-ybar)^2)) |> 
                    mutate(crit = 'Core') |> pivot_longer(corr:mspe) |>
                    mutate(val = 'val'))


facet.label = c('Correlation', 'MSPE')
names(facet.label) = unique(res.tg.oc$name)

res.tg.oc |> 
  ggplot(aes(y = value, x = val)) + 
  geom_violin(aes(fill = name), color = 'black', show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.1), alpha = .6, color = 'darkgrey') + 
  geom_point(data = subset(res.tg.oc, subset = crit != 'Random'), 
             aes(y = value, x = val, color = crit), size = 3, 
             pch = 17) +
  facet_wrap(.~name, scales = 'free_y',labeller = labeller(.cols = facet.label)) + 
  scale_fill_manual(values = c('#edf8fb','#fef0d9')) + 
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3', '#e7298a',
                                '#a6761d','#e6ab02','#e41a1c')) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.position = 'bottom') +
  labs(x = '', y = '', color = 'Criterion')

rm(cv.opt, corr, mspe)

### Fruit dry mass -------------
mod = asreml(fixed = S_T ~ Ano + Plant, 
             data = dat, 
             na.action = na.method(x = 'exclude', y = 'exclude'))
input = predict(mod, classify = 'Plant')$pvals[,1:2] |> 
  rename(ybar = predicted.value)

train = droplevels(input$Plant[which(!input$Plant %in% target)])

rm(mod)

##### Cross-validation (random sampling)
fold.size = 50 
nrept = 100 

cv = list()
k = 1
i = 1
repeat
{
  set.seed(1 + k * 3)
  
  cvdata = input |> mutate(trt = ifelse(Plant %in% sample(train, fold.size), 1, 0), 
                           yNA = ifelse(trt != 1, NA, ybar))
  mod1 = tryCatch({
    asreml(fixed = yNA ~ 1, 
           random = ~ vm(Plant, G.dende.inv), 
           data = cvdata)
  }, error = function(e){cat("There was a singularity issue, sampling again...", fill = TRUE)})
  if(class(mod1) != 'asreml'){
    rm(mod1, cvdata)
    k = k + 1
    next
  }else{
    if(mod1$converge){
      res1 = as.data.frame(coef(mod1)$random + c(coef(mod1)$fixed)) |> 
        rownames_to_column('Plant') |>
        separate(col = 'Plant', into = c('aaa','Plant'), sep = '_') |>
        select(-aaa) |> right_join(cvdata, by = 'Plant') |>
        rename(yhat = effect) |>
        select(Plant, ybar, yhat, trt) |>
        filter(Plant %in% target)
      rm(mod1)
    }else{
      cat("The model did not converge, sampling again...", fill = TRUE)
      k = k+1
      rm(mod1)
      next
    }
  }
  
  cv[[i]] = res1
  cv[[i]]$seed = 1 + k * 3
  
  rm(res1, cvdata)
  
  k = k + 1
  message('Succeded! Repetition ',i)
  i = i + 1
  if(length(cv) == nrept) break 
}

corr = lapply(cv, function(x) cor(x$ybar, x$yhat, use = 'complete.obs'))
mspe = lapply(cv, function(x) mean((x$ybar - x$yhat)^2, na.rm = T))

rm(cv, k, i, fold.size, nrept)

##### Cross-validation (optimized training set)

cv.opt = lapply(split(pc.df.ut, f = pc.df.ut$crit), function(x){
  cvdata = left_join(input, x[, c("Plant", 'trt')], by = 'Plant')
  cvdata$trt[which(is.na(cvdata$trt))] = 0
  cvdata$yNA = cvdata$ybar; cvdata$Plant = as.factor(cvdata$Plant)
  cvdata[which(cvdata$trt != 1),'yNA'] = NA
  
  asreml.options(ai.sing = TRUE)
  mod1 = asreml(fixed = yNA ~ 1, 
                random = ~ vm(Plant, G.dende.inv), 
                data = cvdata)
  asreml.options(ai.sing = FALSE)
  
  as.data.frame(coef(mod1)$random + c(coef(mod1)$fixed)) |> 
    rownames_to_column('Plant') |>
    separate(col = 'Plant', into = c('aaa','Plant'), sep = '_') |>
    select(-aaa) |> right_join(cvdata, by = 'Plant') |>
    rename(yhat = effect) |>
    select(Plant, ybar, yhat, trt) |>
    filter(Plant %in% target) |> 
    reframe(corr = cor(ybar, yhat), 
            mspe = mean((yhat-ybar)^2))
})

res.tg.fdm = rbind(do.call(rbind, cv.opt) |> rownames_to_column('crit'),
                   data.frame(crit = 'Random', 
                              corr = do.call(rbind, corr),
                              mspe = do.call(rbind, mspe)))|> 
  pivot_longer(corr:mspe) |> mutate(val = 'val')

##### Cross-validation (core collection)
cvdata = input; cvdata$yNA = cvdata$ybar
cvdata[which(!cvdata$Plant %in% selcore),'yNA'] = NA
cvdata$Plant = as.factor(cvdata$Plant)

mod1 = asreml(fixed = yNA ~ 1, 
              random = ~ vm(Plant, G.dende.inv), 
              data = cvdata)

res.tg.fdm = rbind(res.tg.fdm, 
                   as.data.frame(coef(mod1)$random + c(coef(mod1)$fixed)) |> 
                     rownames_to_column('Plant') |>
                     separate(col = 'Plant', into = c('aaa','Plant'), sep = '_') |>
                     select(-aaa) |> right_join(cvdata, by = 'Plant') |>
                     rename(yhat = effect) |>
                     select(Plant, ybar, yhat) |>
                     filter(Plant %in% target) |> 
                     reframe(corr = cor(ybar, yhat), 
                             mspe = mean((yhat-ybar)^2)) |> 
                     mutate(crit = 'Core') |> pivot_longer(corr:mspe) |>
                     mutate(val = 'val'))

facet.label = c('Correlation', 'MSPE')
names(facet.label) = unique(res.tg.oc$name)

res.tg.fdm |> 
  ggplot(aes(y = value, x = val)) + 
  geom_violin(aes(fill = name), color = 'black', show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.1), alpha = .6, color = 'darkgrey') + 
  geom_point(data = subset(res.tg.fdm, subset = crit != 'Random'), 
             aes(y = value, x = val, color = crit), size = 3, 
             pch = 17) +
  facet_wrap(.~name, scales = 'free_y',labeller = labeller(.cols = facet.label)) + 
  scale_fill_manual(values = c('#edf8fb','#fef0d9')) + 
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3', '#e7298a',
                                '#a6761d','#e6ab02','#e41a1c')) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.position = 'bottom') +
  labs(x = '', y = '', color = 'Criterion')

rm(cv.opt, corr, mspe)

### Pulp Dry mass -------------
mod = asreml(fixed = S_Polpa ~ Ano + Plant, 
             data = dat, 
             na.action = na.method(x = 'exclude', y = 'exclude'))
input = predict(mod, classify = 'Plant')$pvals[,1:2] |> 
  rename(ybar = predicted.value)

train = droplevels(input$Plant[which(!input$Plant %in% target)])

rm(mod)

##### Cross-validation (random sampling)
fold.size = 50 
nrept = 100 

cv = list()
k = 1
i = 1
repeat
{
  set.seed(1 + k * 3)
  
  cvdata = input |> mutate(trt = ifelse(Plant %in% sample(train, fold.size), 1, 0), 
                           yNA = ifelse(trt != 1, NA, ybar))
  mod1 = tryCatch({
    asreml(fixed = yNA ~ 1, 
           random = ~ vm(Plant, G.dende.inv), 
           data = cvdata)
  }, error = function(e){cat("There was a singularity issue, sampling again...", fill = TRUE)})
  if(class(mod1) != 'asreml'){
    rm(mod1, cvdata)
    k = k + 1
    next
  }else{
    if(mod1$converge){
      res1 = as.data.frame(coef(mod1)$random + c(coef(mod1)$fixed)) |> 
        rownames_to_column('Plant') |>
        separate(col = 'Plant', into = c('aaa','Plant'), sep = '_') |>
        select(-aaa) |> right_join(cvdata, by = 'Plant') |>
        rename(yhat = effect) |>
        select(Plant, ybar, yhat, trt) |>
        filter(Plant %in% target)
      rm(mod1)
    }else{
      cat("The model did not converge, sampling again...", fill = TRUE)
      k = k+1
      rm(mod1)
      next
    }
  }
  
  cv[[i]] = res1
  cv[[i]]$seed = 1 + k * 3
  
  rm(res1, cvdata)
  
  k = k + 1
  message('Succeded! Repetition ',i)
  i = i + 1
  if(length(cv) == nrept) break 
}

corr = lapply(cv, function(x) cor(x$ybar, x$yhat, use = 'complete.obs'))
mspe = lapply(cv, function(x) mean((x$ybar - x$yhat)^2, na.rm = T))

rm(cv, k, i, fold.size, nrept)

##### Cross-validation (optimized training set)

cv.opt = lapply(split(pc.df.ut, f = pc.df.ut$crit), function(x){
  cvdata = left_join(input, x[, c("Plant", 'trt')], by = 'Plant')
  cvdata$trt[which(is.na(cvdata$trt))] = 0
  cvdata$yNA = cvdata$ybar; cvdata$Plant = as.factor(cvdata$Plant)
  cvdata[which(cvdata$trt != 1),'yNA'] = NA
  
  asreml.options(ai.sing = TRUE)
  mod1 = asreml(fixed = yNA ~ 1, 
                random = ~ vm(Plant, G.dende.inv), 
                data = cvdata)
  asreml.options(ai.sing = FALSE)
  
  as.data.frame(coef(mod1)$random + c(coef(mod1)$fixed)) |> 
    rownames_to_column('Plant') |>
    separate(col = 'Plant', into = c('aaa','Plant'), sep = '_') |>
    select(-aaa) |> right_join(cvdata, by = 'Plant') |>
    rename(yhat = effect) |>
    select(Plant, ybar, yhat, trt) |>
    filter(Plant %in% target) |> 
    reframe(corr = cor(ybar, yhat), 
            mspe = mean((yhat-ybar)^2))
})

res.tg.pdm = rbind(do.call(rbind, cv.opt) |> rownames_to_column('crit'),
                   data.frame(crit = 'Random', 
                              corr = do.call(rbind, corr),
                              mspe = do.call(rbind, mspe)))|> 
  pivot_longer(corr:mspe) |> mutate(val = 'val')

##### Cross-validation (core collection)
cvdata = input; cvdata$yNA = cvdata$ybar
cvdata[which(!cvdata$Plant %in% selcore),'yNA'] = NA
cvdata$Plant = as.factor(cvdata$Plant)

mod1 = asreml(fixed = yNA ~ 1, 
              random = ~ vm(Plant, G.dende.inv), 
              data = cvdata)

res.tg.pdm = rbind(res.tg.pdm, 
                   as.data.frame(coef(mod1)$random + c(coef(mod1)$fixed)) |> 
                     rownames_to_column('Plant') |>
                     separate(col = 'Plant', into = c('aaa','Plant'), sep = '_') |>
                     select(-aaa) |> right_join(cvdata, by = 'Plant') |>
                     rename(yhat = effect) |>
                     select(Plant, ybar, yhat) |>
                     filter(Plant %in% target) |> 
                     reframe(corr = cor(ybar, yhat), 
                             mspe = mean((yhat-ybar)^2)) |> 
                     mutate(crit = 'Core') |> pivot_longer(corr:mspe) |>
                     mutate(val = 'val'))

facet.label = c('Correlation', 'MSPE')
names(facet.label) = unique(res.tg.pdm$name)

res.tg.pdm |> 
  ggplot(aes(y = value, x = val)) + 
  geom_violin(aes(fill = name), color = 'black', show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.1), alpha = .6, color = 'darkgrey') + 
  geom_point(data = subset(res.tg.pdm, subset = crit != 'Random'), 
             aes(y = value, x = val, color = crit), size = 3, 
             pch = 17) +
  facet_wrap(.~name, scales = 'free_y',labeller = labeller(.cols = facet.label)) + 
  scale_fill_manual(values = c('#edf8fb','#fef0d9')) + 
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3', '#e7298a',
                                '#a6761d','#e6ab02','#e41a1c')) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.position = 'bottom') +
  labs(x = '', y = '', color = 'Criterion')

rm(cv.opt, corr, mspe)

### Wrap-up (CV) -----------------------------------------------------------------
res.tg = rbind(res.tg.oc |> mutate(trait = 'OC'), 
               res.tg.fdm |> mutate(trait = 'FDM'), 
               res.tg.pdm |> mutate(trait = 'PDM')) 

res.tg |> 
  ggplot(aes(y = value, x = val)) + 
  geom_violin(aes(fill = name), color = 'black', show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.1), alpha = .6, color = 'darkgrey') + 
  geom_point(data = subset(res.tg, subset = crit != 'Random'), 
             aes(y = value, x = val, color = crit), size = 3, 
             pch = 17) +
  facet_wrap(name~trait, scales = 'free', ncol = 3) + 
  scale_fill_manual(values = c('#edf8fb','#fef0d9')) + 
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3', '#e7298a',
                                '#a6761d','#e6ab02','#e41a1c')) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.position = 'bottom') +
  labs(x = '', y = '', color = 'Criterion')

tg1.2_3 = list(overview = pc.df.ut, cv = res.tg)
