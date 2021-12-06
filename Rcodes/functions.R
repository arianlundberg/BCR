library(ggcorrplot)
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(webshot)
library(htmlwidgets)
library(metafor)
library(metap)
library(purrr)
library(biomaRt)
library(ggplot2)
library(forestplot)
library(ggExtra)
library(survival)
library(tibble)
library(doParallel)
library(org.Hs.eg.db)
library(survminer)
library(ComplexHeatmap)
library(survplot)
library(plyr)
library(tableone)
library(ggthemes)
library(ggpubr)
library(dplyr)

load(file = "~/RData/META-BRCA.RData")
load(file="~/RData/META-CRC.RData")
load(file="~/RData/META-NSCLC.RData")
load(file="~/RData/META-SKCM.RData")
load(file="~/RData/sig.genes.RData")

options(stringsAsFactors = F)
registerDoParallel(cores = detectCores(all.tests = FALSE, logical = TRUE)-2)
my.corr.test <- function(x,y) {
  cbind(paste('r =',round(cor.test(x,y,method='pearson')$estimate,digits = 3)),
        ifelse(cor.test(x,y,method='pearson')$p.value < 0.001,'p. < 0.001',paste('p. =',round(cor.test(x,y,method='pearson')$p.value,digits = 3))))}

My.interaction.function <- function(ImmuneCell,Celltype.exprs,time,event, ...){

  na.omit(as.data.frame(foreach(i = 1:nrow(Celltype.exprs), .packages = 'survival', .combine = 'rbind') %dopar% {
    surv <- Surv(time, event)
    int.model <- summary(coxph(surv ~ ImmuneCell * as.numeric(Celltype.exprs[i,]), as.data.frame(Celltype.exprs)))
    c(main = int.model$coefficients[1,1],
      se = int.model$coefficients[1,3],
      z.score = int.model$coefficients[1,4],
      pvalue = int.model$coefficients[1,5],

      interaction = int.model$coefficients[nrow(int.model$coefficients),1],
      se.int = int.model$coefficients[nrow(int.model$coefficients),3],
      z.score.int = int.model$coefficients[nrow(int.model$coefficients),4],
      pvalue.int = int.model$coefficients[nrow(int.model$coefficients),5])},row.names = rownames(Celltype.exprs)))}

