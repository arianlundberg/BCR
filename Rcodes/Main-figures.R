## this code is to generate main figures of the manuscript

## load the RData all required datasets and libraries

### RData files are accessible here
# https://www.dropbox.com/sh/uggkgtelk192cwa/AADJhH6AVyooBBc2k4LGOxG4a?dl=0
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

load(file = "~/Dropbox/Public/BCM/RData/META-BRCA.RData")
load(file="~/Dropbox/Public/BCM/RData/META-CRC.RData")
load(file="~/Dropbox/Public/BCM/RData/META-NSCLC.RData")
load(file="~/Dropbox/Public/BCM/RData/META-SKCM.RData")
load(file="~/Dropbox/Public/BCM/RData/sig.genes.RData")
source(file="~/Documents/GitHub/BCM/Rcodes/functions.R")

## function to generate interaction values of all genes in the dataset
# the results are saved as a list in int.list

#### generating a list of Cytokine signalling genes interacted with B cell lineage - significant FDR

TNBC.sig.genes <- as.data.frame(foreach(i = 1:nrow(common.cytokines), .packages = 'metap', .combine = 'rbind') %dopar% {
  sumz.p <- metap::sumz(p = unlist(purrr::map(purrr::map(int.list,'pvalue.int'),i)),weights = weights)
  coef.df <- metafor::rma(unlist(purrr::map(purrr::map(int.list,'interaction'),i)),sei =unlist(purrr::map(purrr::map(int.list,'se.int'),i)),control=list(maxiter=1000),method='ML',measure = 'RR')
  c(Meta.int = round(coef.df$b,digits = 3),
    Meta.CI.l =round(coef.df$ci.lb,digits = 3),
    Meta.CI.u = round(coef.df$ci.ub,digits = 3),
    Sumz.z=round(sumz.p$z,digits = 3),
    Sumz.p=round(sumz.p$p,digits=3)
  )},row.names = common.cytokines$SYMBOL) %>% na.omit() %>% add_column(fdr = round(p.adjust(.$Sumz.p, method = 'fdr'),digits = 4),.after = 'Sumz.p')


#### Genes passing FDR < or == 0.1 criteria - Table 1

final.table.sig <- TNBC.sig.genes[TNBC.sig.genes$fdr <= 0.1,]
final.table.sig <- final.table.sig %>% add_column(SYMBOL=biomaRt::select(org.Hs.eg.db,keys = rownames(final.table.sig),columns = c('ENTREZID','SYMBOL'),keytype = 'SYMBOL')$SYMBOL,.before = 'Meta.int')
final.table.sig <- final.table.sig %>% add_column(ENTREZID=biomaRt::select(org.Hs.eg.db,keys = rownames(final.table.sig),columns = c('ENTREZID','SYMBOL'),keytype = 'SYMBOL')$ENTREZID,.before = 'Meta.int')
# dim(final.table.sig)
# 15  8
tmp.table <- final.table.sig[,c(1:3,7,8)];colnames(tmp.table)[c(3:5)] <- c('Coeff','p.','FDR')
tmp.table[c(3:5)] <- signif(tmp.table[c(3:5)],digits = 2)

tmp.table <- tmp.table[order(tmp.table$Coeff,decreasing = T),]

# Table 1

first.table.genes <- tmp.table %>% mutate(p.=replace(p.,p. <0.001,'< 0.001'),FDR=replace(FDR,FDR < 0.001,'< 0.001'))

# SYMBOL ENTREZID  Coeff    p.   FDR
# TRAF2       TRAF2     7186  0.250 0.002 0.041
# TNFRSF4   TNFRSF4     7293  0.210 0.005  0.07
# TNFRSF8   TNFRSF8      943  0.200 0.002 0.041
# IL2RG       IL2RG     3561  0.200 0.006  0.07
# LCK           LCK     3932  0.180 0.003 0.049
# BATF         BATF    10538  0.140 0.002 0.041
# TNFRSF14 TNFRSF14     8764  0.110 0.003 0.049
# CXCL13     CXCL13    10563  0.089 0.001 0.041
# IL18RAP   IL18RAP     8807  0.045 0.002 0.041
# PSMB10     PSMB10     5699  0.041 0.001 0.041
# CXCR6       CXCR6    10663  0.035 0.002 0.041
# IL23A       IL23A    51561 -0.012 0.007 0.076
# PTPN6       PTPN6     5777 -0.029 0.002 0.041
# APP           APP      351 -0.160 0.006  0.07
# NOD1         NOD1    10392 -0.170 0.006  0.07


####

#head(CTL.Bcell)
## a data.frame with subset of the data to generate figure 1A - C


# B.lineage        CTL B.bin OS  OS.year CTL.bin CTL.B.groups
# 1  0.35057272 -0.3007329     2  0       NA       1          III
# 2 -1.07839633 -0.8200763     1  0 1.849315       1            I
# 3  0.17676515 -0.4347091     2  0       NA       1          III
# 4 -0.57543003  0.1857145     1  0       NA       2           II
# 5 -0.56540467 -0.1217551     1  0 2.701370       1            I
# 6 -0.03870975 -0.0280080     2  0 4.731507       2           IV

### generating groups of B.cell Low/T.cells Low - B.cell Low/T.cells High - B.cell High/T cells low - B.cell High/T.cells High

CTL.Bcell$CTL.B.groups <- as.factor(ifelse(CTL.Bcell$B.bin %in% '1'&CTL.Bcell$CTL.bin %in% '1','I',
                                           ifelse(CTL.Bcell$B.bin %in% '1'&CTL.Bcell$CTL.bin %in% '2','II',
                                                  ifelse(CTL.Bcell$B.bin %in% '2'&CTL.Bcell$CTL.bin %in% '1','III',
                                                         ifelse(CTL.Bcell$B.bin %in% '2'&CTL.Bcell$CTL.bin %in% '2','IV','NA')))))


#### Correlation plot - Figure 1A

tmp.plot <- ggplot(CTL.Bcell, aes(x=CTL, y=B.lineage,color=CTL.B.groups)) + 
  geom_point(size = 2) + ggsci::scale_color_jco() +
  geom_smooth(method=lm, color='black',fill='lightcyan4') + 
  scale_x_continuous(name=expression(bold("T cells")),limits = c(-2.5,2.5)) + scale_y_continuous(name=expression(bold("B.lineage-cell lineage abundance score (B.cell)")),limits = c(-2.5,2.5)) +
  geom_vline(xintercept=median(CTL.Bcell$CTL), linetype="dashed", color = "hotpink4")+
  geom_hline(yintercept=median(CTL.Bcell$B.lineage), linetype="dashed", color = "hotpink4")+
  annotate("text", x = min(CTL.Bcell$CTL)+1, y = max(CTL.Bcell$B.lineage)-1, label = my.corr.test(CTL.Bcell$B.lineage, CTL.Bcell$CTL)[1],colour='black',fontface=2,size=8) +
  annotate("text", x = min(CTL.Bcell$CTL)+1, y = max(CTL.Bcell$B.lineage)-1.5, label = my.corr.test(CTL.Bcell$B.lineage, CTL.Bcell$CTL)[2],colour='black',fontface=2,size=8) +
  theme(legend.key = element_rect(fill='white'),legend.position = 'None',
        strip.text.x =element_text(size = 10),
        axis.title.x = element_text(size = 24,face='bold'),
        axis.title.y = element_text(size = 24,face='bold'),
        plot.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.text.x=element_text(size = 22,angle = 45,hjust = c(1,1),face='bold'),
        axis.text.y=element_text(size = 22,face='bold'),
        plot.margin = unit(c(0.5,0.5,0.7,0.5),'cm'),
        axis.line.x = element_line(color="gray", size = 0.5),
        axis.line.y = element_line(color="gray", size = 0.5),
        legend.key.size = unit(1.3, "cm"),legend.title = element_text(size = 15,face='bold'),legend.text = element_text(size = 15))


Figure1A <- ggExtra::ggMarginal(tmp.plot, type = 'densigram',fill='lightsteelblue2') 


###### 

# generating data for the forestplot

model.B <- coxph(Surv(CTL.Bcell$OS.year, CTL.Bcell$OS) ~ B.lineage, data= CTL.Bcell)
model.CTL <- coxph(Surv(CTL.Bcell$OS.year, CTL.Bcell$OS) ~ CTL, data= CTL.Bcell)

model <- coxph(Surv(CTL.Bcell$OS.year, CTL.Bcell$OS) ~ B.lineage+CTL, data= CTL.Bcell)

m<- c(summary(model.B)$coeff[,2],summary(model.CTL)$coeff[,2],NA,summary(model)$coeff[,2])
l<- c(signif(summary(model.B)$conf.int[,3],digits=2),signif(summary(model.CTL)$conf.int[,3],digits=2),NA,signif(summary(model)$conf.int[,3],digits = 2))
u<- c(signif(summary(model.B)$conf.int[,4],digits=),signif(summary(model.CTL)$conf.int[,4],digits=2),NA,signif(summary(model)$conf.int[,4],digits = 2))

np <- c(paste(c(length(CTL.Bcell$B.lineage),length(CTL.Bcell$B.lineage),length(CTL.Bcell$B.lineage),length(CTL.Bcell$B.lineage),length(CTL.Bcell$CTL))," (",
              signif(c(length(CTL.Bcell$B.lineage),length(CTL.Bcell$B.lineage),length(CTL.Bcell$B.lineage),length(CTL.Bcell$B.lineage),length(CTL.Bcell$CTL))/766*100,digits=1),")",sep=""))

tabletext<-cbind(c(NA,'Univariate',c(' B cells',' T.cells','Multivariate',' B cells',' T.cells')),
                 c("HR",NA,format(summary(model.B)$coeff[,2],digits=2),format(summary(model.CTL)$coeff[,2],digits=2),NA,format(summary(model)$coeff[,2],digits=2)),
                 c('95% CI',NA,paste0('[',format(l,digits = 3),' - ',format(u,digits=3),']',sep = '')),
                 c("Coefficient",NA,format(summary(model.B)$coeff[,1],digits=2),format(summary(model.CTL)$coeff[,1],digits=2),NA,format(summary(model)$coeff[,1],digits=2)),
                 c("p.",NA,format(summary(model.B)$coeff[,5],digits=2),format(summary(model.CTL)$coeff[,5],digits=2),NA,format(summary(model)$coeff[,5],digits=2)))
tabletext[5,3] <- NA
tabletext[,5] <- ifelse(tabletext[,5] < 0.001,'< 0.001',tabletext[,5])


##### Forestplot Figure 1B
forestplot::forestplot(tabletext,c(NA,NA,m),c(NA,NA,l),c(NA,NA,u),zero=1,align = c("l", "c", "c","c","c","c"),
                                   graph.pos=6, 
                                   title="Hazard Ratio",lwd.ci=2, ci.vertices=TRUE, ci.vertices.height = 0.05,lwd.zero=3,
                                   hrzl_lines = list("2" = gpar(lwd=3, columns=1:6, col = "black"),
                                                     "3" = gpar(lwd=1, columns=1, col = "black"),
                                                     "6" = gpar(lwd=1, columns=1, col = "black")),
                                   is.summary=c(TRUE,TRUE,rep(FALSE,2),TRUE,rep(FALSE,2)),boxsize = .3,
                                   colgap=unit(4,"mm"),lineheight = unit(12,"mm"),
                                   xticks=c(0.6,0.8,1.0,1.2,1.4),graphwidth=unit(7,"cm"),
                                   xlog=F,txt_gp = fpTxtGp(label = list(gpar(fontfamily = "Arial"),
                                                                        gpar(fontfamily = "")),cex = 1.75,
                                                           ticks = gpar(fontfamily = "", cex=1.5),
                                                           xlab  = gpar(fontfamily = "HersheySerif", cex = 2)),col=fpColors(box=c("blue"),lines='black',hrz_lines = "gray50"))


##### KM Figure 1C

tmp.fit <- survfit(Surv(CTL.Bcell$OS.year, CTL.Bcell$OS) ~ CTL.B.groups, data = CTL.Bcell)

Figure1C <- ggsurvplot(tmp.fit, data = CTL.Bcell, palette = "jco", pval=T,risk.table = TRUE,legend.labs=c('B.cell Low/T.cells Low','B.cell Low/T.cells High','B.cell High/T.cells Low','B.cell High/T.cells High'),
                       font.legend=20,font.x=20,font.y=20,font.tickslab=20,pval.size=10,risk.table.fontsize=8,font.main=20,xlim = c(0,10),legend='right',
                       risk.table.col = "strata",risk.table.y.text=F,legend.title='Patient groups',
                       xlab='Time (years)',ylab='Overall survival')+
  guides(colour = guide_legend(nrow = 4))




# Figure 2
# Figure 2A BCM correlation plot

#### BCM genes
tmp.BCM <- final.table.sig[final.table.sig$Meta.int > 0,]

corr.score.BCM <- cor(TNBC.exprs[,tmp.BCM$SYMBOL],method = 'pearson')
colnames(corr.score.BCM) <- paste0(tmp.BCM$SYMBOL," (",format((tmp.BCM$Meta.int),digits = 1),")")

Figure2A <- draw(Heatmap(corr.score.BCM*100,col=circlize::colorRamp2(c(-100, 0, 100), c("green", "white", "red")),cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%.1f", corr.score.BCM[i, j]*100), x, y, gp = gpar(fontsize = 10))
},column_km =2 ,row_km = 2,row_gap = unit(0.8, "mm"),heatmap_legend_param = list(direction = "horizontal"),column_title = ' ',row_title = ' ',column_gap = unit(0.8, "mm"),border = T,column_names_rot = 45,name='Pearson Correlation (%)'),heatmap_legend_side='bottom')

  
# Figure 2B BCM Forest plot

BRCA.FP.data <- merge(TNBC.info,TNBC.exprs,by='row.names',sort=F,all.x=T);rownames(BRCA.FP.data) <- BRCA.FP.data$Row.names;BRCA.FP.data <- BRCA.FP.data[,-1]
model.BRCA <- coxph(Surv(BRCA.FP.data$OS.year, BRCA.FP.data$OS) ~ B.lineage*BCM, data= BRCA.FP.data)


m.BRCA<- c(NA,NA,summary(model.BRCA)$coeff[1,2],
           NA,NA,summary(model.BRCA)$coeff[2,2],
           NA,NA,summary(model.BRCA)$coeff[3,2])

l.BRCA<- c(NA,NA,format(signif(summary(model.BRCA)$conf.int[1,3],digits = 3)),
           NA,NA,format(signif(summary(model.BRCA)$conf.int[2,3],digits = 3)),
           NA,NA,format(signif(summary(model.BRCA)$conf.int[3,3],digits = 4)))

u.BRCA<- c(NA,NA,format(signif(summary(model.BRCA)$conf.int[1,4],digits = 4)),
           NA,NA,format(signif(summary(model.BRCA)$conf.int[2,4],digits = 4)),
           NA,NA,format(round(summary(model.BRCA)$conf.int[3,4],digits = 3)))


np.BRCA <- c(paste(nrow(BRCA.FP.data)," (",'100',")",sep=""))
np.BRCA[np.BRCA=='NA (NA)'] <- NA

tabletext.BRCA<-cbind(c(c(NA,'META-BRCA','B cells','BCM','B cells * BCM')),
                      c("No. of Patients (%)",np.BRCA,NA,NA,NA), 
                      c("HR",NA,format((summary(model.BRCA)$coeff[,2]),digits=3)[1],format((summary(model.BRCA)$coeff[,2])[2],digits=3),format((summary(model.BRCA)$coeff[,2])[3],digits=4)),
                      c('95% CI',NA,paste0('[',l.BRCA[3],' - ',u.BRCA[3],']',sep = ''),
                        paste0('[',l.BRCA[6],' - ',u.BRCA[6],']',sep = ''),
                        paste0('[',l.BRCA[9],' - ',u.BRCA[9],']',sep = '')),
                      c("Coeff",NA,format((summary(model.BRCA)$coeff[,1]),digits=3)[1],format((summary(model.BRCA)$coeff[,1]),digits=3)[2],format((summary(model.BRCA)$coeff[,1])[3],digits=3)),
                      c("p.",NA,format((summary(model.BRCA)$coeff[,5]),digits=2)[1],format((summary(model.BRCA)$coeff[,5])[2],digits=2),format((summary(model.BRCA)$coeff[,5]),digits=2)[3]))



lower.FP.BRCA <- round(as.numeric(c(NA,NA,l.BRCA[c(3,6,9)])),digits=3)
upper.FP.BRCA <- round(as.numeric(c(NA,NA,u.BRCA[c(3,6,9)])),digits=3)
mean.FP.BRCA <- round(as.numeric(c(NA,NA,m.BRCA[c(3,6,9)])),digits = 3)

#generating a forestplot 
# Figure 2B
forestplot::forestplot(tabletext.BRCA,mean.FP.BRCA,lower.FP.BRCA,upper.FP.BRCA,zero=1,
                       title="Hazard Ratio",lwd.ci=2, ci.vertices=TRUE, ci.vertices.height = 0.05,
                       hrzl_lines = list("2" = gpar(lwd=3, columns=1:6, col = "black"),
                                         "3" = gpar(lwd=1, columns=1, col = "black")),lwd.zero=5,
                       clip=c(0,5),xticks=seq(from = 0.5,to = 2.5,by=0.5),graphwidth=unit(10,"cm"),
                       is.summary=c(TRUE,TRUE,rep(FALSE,3)),boxsize = .15,
                       colgap=unit(6,"mm"),lineheight = unit(16,"mm"),align = c("l", "c", "c","c","c","c"),
                       xlog=F,txt_gp = fpTxtGp(label = list(gpar(fontfamily = "Arial"),
                                                            gpar(fontfamily = "")),cex = 1.75,
                                               ticks = gpar(fontfamily = "", cex=1.5),
                                               xlab  = gpar(fontfamily = "Arial", cex = 2)),col=fpColors(box=c("blue"),lines='black',hrz_lines = "gray50"))




# Figure 2C and D BCM Kaplan Meiers


#### finding an optimal cutpoint based on HR
#### Generating a list of HR from B.lineage low group to be able to generate HR for genes in that subset. 
Gene.Low.TNBC <- foreach (i=c(5:ncol(TNBC.info))) %:%
  foreach(j=c(1:766),.combine='rbind',.errorhandling = 'pass') %dopar% {
    tmp <- as.data.frame(TNBC.info[TNBC.info$B.binary %in% 'Low',])
    gene <- summary(coxph(Surv(tmp$OS.year, tmp$OS) ~ B.lineage, 
                          data= tmp, subset = tmp[,i] < TNBC.info[,i][j]))
    c(gene=TNBC.info[,i][j],Coeff=gene$coeff[,c('exp(coef)')],p=gene$coeff[,c('Pr(>|z|)')])}
Gene.High.TNBC <- foreach (i=c(5:ncol(TNBC.info))) %:%
  foreach(j=c(1:766),.combine='rbind',.errorhandling = 'pass') %dopar% {
    tmp <- as.data.frame(TNBC.info[TNBC.info$B.binary %in% 'High',])
    gene <- summary(coxph(Surv(tmp$OS.year, tmp$OS) ~ B.lineage, 
                          data= tmp, subset = tmp[,i] < TNBC.info[,i][j]))
    c(gene=TNBC.info[,i][j],Coeff=gene$coeff[,c('exp(coef)')],p=gene$coeff[,c('Pr(>|z|)')])}

#### All genes in one data.frame
tmp.Low <- foreach (i=c(1:length(Gene.Low.TNBC)),.combine='rbind',.errorhandling = 'pass') %dopar%{
  data.frame(Coeff.Low=as.numeric(as.character(Gene.Low.TNBC[[i]][,2])),
             p.Low=as.numeric(as.character(Gene.Low.TNBC[[i]][,3])))}
tmp.High <- foreach (i=c(1:length(Gene.High.TNBC)),.combine='rbind',.errorhandling = 'pass') %dopar%{
  data.frame(gene=as.numeric(as.character(Gene.High.TNBC[[i]][,1])),
             Coeff.High=as.numeric(as.character(Gene.High.TNBC[[i]][,2])),
             p.High=as.numeric(as.character(Gene.High.TNBC[[i]][,3])))}
tmp.Gene <- cbind(tmp.High,tmp.Low)

### Convert to List
Gene.TNBC <- split(tmp.Gene, rep(1:length(Gene.Low.TNBC), length.out = nrow(tmp.Gene), each = ceiling(nrow(tmp.Gene)/length(Gene.Low.TNBC))))
names(Gene.TNBC) <- colnames(TNBC.info)[5:ncol(TNBC.info)]
tmp.cut.BCM <- data.frame(cutpoint=Gene.TNBC$BCM$gene,HR.ratio=Gene.TNBC$BCM$Coeff.High-Gene.TNBC$BCM$Coeff.Low)
tmp.cut.BCM$percent <- percent_rank(tmp.cut.BCM$HR.ratio)
cut.point.BCM <- tmp.cut.BCM[which.max(tmp.cut.BCM[tmp.cut.BCM$percent > 0.2&tmp.cut.BCM$percent < 0.8,]$HR.ratio),]$cutpoint


#### generating the plots 
# FIGURE 2C-D

par(mfrow=c(2,2),mai = c(1, 1.25, 0.75, 0.5))

foreach(i='BCM',.combine='cbind') %do% {
  tmp.Low <-   TNBC.info[TNBC.info[,i] < cut.point.BCM,]
  b.Low <- relevel(as.factor(ifelse(ntile(tmp.Low[,4],n = 2)=='1','Low','High')),ref = 'Low')
  tmp.High <-   TNBC.info[TNBC.info[,i] >cut.point.BCM,]
  b.High <- relevel(as.factor(ifelse(ntile(tmp.High[,4],n = 2)=='1','Low','High')),ref = 'Low')
  survplot(Surv(tmp.Low$OS.year, tmp.Low$OS)~b.Low,show.nrisk = T,data=tmp.Low,xlim=c(0,10),
           xlab=expression(bold('Time (years)')),ylab=expression(bold("OS")),main=paste(colnames(tmp.Low[i]),":","Low",sep = ""),
           lwd=3,col=c("royalblue3","tomato3"),legend.pos = 'bottomleft',stitle='',mark=20,cex.lab=1.5, cex.main=1.5,snames = c('B.cell < Median','B.cell > Median'))
  pval = summary(coxph(Surv(tmp.Low$OS.year, tmp.Low$OS) ~ B.lineage, data= tmp.Low))$coeff[,'Pr(>|z|)']
  pval=ifelse(formatC(format='f',pval, digits=3) < '0.001', '< 0.001',formatC(format='f',pval, digits=3))
  #plot.new()
  mtext("B cell dysfunction signature", side = 3,outer = T,cex=1.5,font = 1)
  survplot(Surv(tmp.High$OS.year, tmp.High$OS)~b.High,show.nrisk = T,data=tmp.High,xlim=c(0,10),
           xlab=expression(bold('Time (years)')),ylab=expression(bold("OS")),main=paste(colnames(tmp.High[i]),":","High",sep = ""),
           lwd=3,col=c("royalblue3","tomato3"),legend.pos = 'bottomleft',stitle='',mark=20,cex.lab=1.5, cex.main=1.5,snames = c('B.cell < Median','B.cell > Median'))
}




# Figure 3
# Figure 3 A-C-E - correlation plots

#### correlation plots
### NSCLC

LUNG.corr.score.BCM <- cor(x = t(LUNG.exprs[BCM.gene$ENTREZID,]),method = 'pearson')
rownames(LUNG.corr.score.BCM)= BCM.gene$SYMBOL;colnames(LUNG.corr.score.BCM)= BCM.gene$SYMBOL
LUNG.B.BCM <- Heatmap(LUNG.corr.score.BCM*100,col=circlize::colorRamp2(c(-100, 0, 100), c("green", "white", "red")),cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%.1f", LUNG.corr.score.BCM[i, j]*100), x, y, gp = gpar(fontsize = 7))
},heatmap_legend_param = list(direction = "horizontal"),column_title = " ",row_title = "META-NSCLC",row_gap = unit(0.8, "mm"), column_gap = unit(0.8, "mm"),border = T,column_names_rot = 45,name='Pearson Correlation (%)')

### CRC

CRC.corr.score.BCM <- cor(x = t(CRC.exprs[BCM.gene$ENTREZID,]),method = 'pearson')
rownames(CRC.corr.score.BCM)= BCM.gene$SYMBOL;colnames(CRC.corr.score.BCM)= BCM.gene$SYMBOL
CRC.B.BCM <- Heatmap(CRC.corr.score.BCM*100,col=circlize::colorRamp2(c(-100, 0, 100), c("green", "white", "red")),cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%.1f", CRC.corr.score.BCM[i, j]*100), x, y, gp = gpar(fontsize = 7))
},heatmap_legend_param = list(direction = "horizontal"),column_title = " ",row_title = "META-CRC",row_gap = unit(0.8, "mm"), column_gap = unit(0.8, "mm"),border = T,column_names_rot = 45,name='Pearson Correlation (%)')

### SKCM

MELA.corr.score.BCM <- cor(x = t(MELA.exprs[BCM.gene$ENTREZID,]),method = 'pearson')
rownames(MELA.corr.score.BCM)= BCM.gene$SYMBOL;colnames(MELA.corr.score.BCM)= BCM.gene$SYMBOL
MELA.B.BCM <- Heatmap(MELA.corr.score.BCM*100,col=circlize::colorRamp2(c(-100, 0, 100), c("green", "white", "red")),cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%.1f", MELA.corr.score.BCM[i, j]*100), x, y, gp = gpar(fontsize = 7))
},heatmap_legend_param = list(direction = "horizontal"),column_title = " ",row_title = "META-SKCM",row_gap = unit(0.8, "mm"), column_gap = unit(0.8, "mm"),border = T,column_names_rot = 45,name='Pearson Correlation (%)')


# Figure3 B-D-F KMs 

par(mfrow=c(3,4),omi=c(0, 0.9, 0, 0))


foreach(i='BCM',.combine='cbind') %do% {
  tmp.Low <-   LUNG.info[LUNG.info[,i] < cut.point.BCM.LUNG,]
  b.Low <- relevel(as.factor(ifelse(ntile(tmp.Low[,4],n = 2)=='1','Low','High')),ref = 'Low')
  tmp.High <-  LUNG.info[LUNG.info[,i] > cut.point.BCM.LUNG,]
  b.High <- relevel(as.factor(ifelse(ntile(tmp.High[,4],n = 2)=='1','Low','High')),ref = 'Low')
  
  survplot(Surv(tmp.Low$OS.year, tmp.Low$OS)~b.Low,show.nrisk = T,data=tmp.Low,stitle='',xlim=c(0,10),
           xlab=expression(bold('Time (years)')),ylab=expression(bold("OS")),main=paste(colnames(tmp.Low[i]),":","Low",sep = ""),
           lwd=3,col=c("royalblue3","tomato3"),legend.pos = 'bottomleft',mark=20,cex.lab=1.5, cex.main=1.5,snames = c('B.cell < Median','B.cell > Median'))
  pval = summary(coxph(Surv(tmp.Low$OS.year, tmp.Low$OS) ~ B.lineage, data= tmp.Low))$coeff[,'Pr(>|z|)']
  pval=ifelse(formatC(format='f',pval, digits=3) < '0.001', '< 0.001',formatC(format='f',pval, digits=3))
  mtext("B cell dysfunction (BCM) signature", side = 3,outer = T,cex=1.5,font = 1)
  
  survplot(Surv(tmp.High$OS.year, tmp.High$OS)~b.High,show.nrisk = T,data=tmp.High,stitle='',xlim=c(0,10),
           xlab=expression(bold('Time (years)')),ylab=expression(bold("OS")),main=paste(colnames(tmp.High[i]),":","High",sep = ""),
           lwd=3,col=c("royalblue3","tomato3"),legend.pos = 'bottomleft',mark=20,cex.lab=1.5, cex.main=1.5,snames = c('B.cell < Median','B.cell > Median'))
  pval = summary(coxph(Surv(tmp.High$OS.year, tmp.High$OS) ~ B.lineage, data= tmp.High))$coeff[,'Pr(>|z|)']
  pval=ifelse(formatC(format='f',pval, digits=3) < '0.001', '< 0.001',formatC(format='f',pval, digits=3))
}

plot.new()
plot.new()

foreach(i='BCM',.combine='cbind') %do% {
  tmp.Low <-   CRC.info[CRC.info[,i] < (cut.point.BCM.CRC),]
  b.Low <- relevel(as.factor(ifelse(ntile(tmp.Low[,4],n = 2)=='1','Low','High')),ref = 'Low')
  tmp.High <-  CRC.info[CRC.info[,i] > cut.point.BCM.CRC,]
  b.High <- relevel(as.factor(ifelse(ntile(tmp.High[,4],n = 2)=='1','Low','High')),ref = 'Low')
  
  survplot(Surv(tmp.Low$OS.year, tmp.Low$OS)~b.Low,show.nrisk = T,data=tmp.Low,stitle='',xlim=c(0,10),
           xlab=expression(bold('Time (years)')),ylab=expression(bold("OS")),main=paste(colnames(tmp.Low[i]),":","Low",sep = ""),
           lwd=3,col=c("royalblue3","tomato3"),legend.pos = 'bottomleft',mark=20,cex.lab=1.5, cex.main=1.5,snames = c('B.cell < Median','B.cell > Median'))
  pval = summary(coxph(Surv(tmp.Low$OS.year, tmp.Low$OS) ~ B.lineage, data= tmp.Low))$coeff[,'Pr(>|z|)']
  pval=ifelse(formatC(format='f',pval, digits=3) < '0.001', '< 0.001',formatC(format='f',pval, digits=3))
  mtext("B cell dysfunction (BCM) signature", side = 3,outer = T,cex=1.5,font = 1)
  
  survplot(Surv(tmp.High$OS.year, tmp.High$OS)~b.High,show.nrisk = T,data=tmp.High,stitle='',xlim=c(0,10),
           xlab=expression(bold('Time (years)')),ylab=expression(bold("OS")),main=paste(colnames(tmp.High[i]),":","High",sep = ""),
           lwd=3,col=c("royalblue3","tomato3"),legend.pos = 'bottomleft',mark=20,cex.lab=1.5, cex.main=1.5,snames = c('B.cell < Median','B.cell > Median'))
  pval = summary(coxph(Surv(tmp.High$OS.year, tmp.High$OS) ~ B.lineage, data= tmp.High))$coeff[,'Pr(>|z|)']
  pval=ifelse(formatC(format='f',pval, digits=3) < '0.001', '< 0.001',formatC(format='f',pval, digits=3))
}

plot.new()
plot.new()


foreach(i='BCM',.combine='cbind') %do% {
  tmp.Low <-   MELA.info[MELA.info[,i] < cut.point.BCM.MELA,]
  b.Low <- relevel(as.factor(ifelse(ntile(tmp.Low[,4],n = 2)=='1','Low','High')),ref = 'Low')
  tmp.High <-  MELA.info[MELA.info[,i] > cut.point.BCM.MELA,]
  b.High <- relevel(as.factor(ifelse(ntile(tmp.High[,4],n = 2)=='1','Low','High')),ref = 'Low')
  
  survplot(Surv(tmp.Low$OS.year, tmp.Low$OS)~b.Low,show.nrisk = T,data=tmp.Low,stitle='',xlim=c(0,10),
           xlab=expression(bold('Time (years)')),ylab=expression(bold("OS")),main=paste(colnames(tmp.Low[i]),":","Low",sep = ""),
           lwd=3,col=c("royalblue3","tomato3"),legend.pos = 'bottomleft',mark=20,cex.lab=1.5, cex.main=1.5,snames = c('B.cell < Median','B.cell > Median'))
  pval = summary(coxph(Surv(tmp.Low$OS.year, tmp.Low$OS) ~ B.lineage, data= tmp.Low))$coeff[,'Pr(>|z|)']
  pval=ifelse(formatC(format='f',pval, digits=3) < '0.001', '< 0.001',formatC(format='f',pval, digits=3))
  mtext("B cell dysfunction (BCM) signature", side = 3,outer = T,cex=1.5,font = 1)
  
  survplot(Surv(tmp.High$OS.year, tmp.High$OS)~b.High,show.nrisk = T,data=tmp.High,stitle='',xlim=c(0,10),
           xlab=expression(bold('Time (years)')),ylab=expression(bold("OS")),main=paste(colnames(tmp.High[i]),":","High",sep = ""),
           lwd=3,col=c("royalblue3","tomato3"),legend.pos = 'bottomleft',mark=20,cex.lab=1.5, cex.main=1.5,snames = c('B.cell < Median','B.cell > Median'))
  pval = summary(coxph(Surv(tmp.High$OS.year, tmp.High$OS) ~ B.lineage, data= tmp.High))$coeff[,'Pr(>|z|)']
  pval=ifelse(formatC(format='f',pval, digits=3) < '0.001', '< 0.001',formatC(format='f',pval, digits=3))
}


# Figure 3G - data-prep

LUNG1.pc <- LUNG1.info %>% dplyr::select('age at surgery:ch1','gender:ch1')
colnames(LUNG1.pc) <- c('Age','Gender');LUNG1.pc$Age <- as.numeric(LUNG1.pc$Age)
LUNG1.pc$Pstage <- NA
LUNG1.pc$Gender <- revalue(LUNG1.pc$Gender, c("M"='Male',"F"='Female'))

LUNG2.pc <- LUNG2.info %>% dplyr::select('age (years):ch1','gender:ch1','pathological stage:ch1')
colnames(LUNG2.pc) <- c('Age','Gender','Pstage');LUNG2.pc$Age <- as.numeric(LUNG2.pc$Age)
LUNG2.pc$Gender <- revalue(LUNG2.pc$Gender, c("male"='Male',"female"='Female'))

LUNG3.pc <- LUNG3.info %>% dplyr::select('age:ch1','Sex:ch1','Stage:ch1')
colnames(LUNG3.pc) <- c('Age','Gender','Pstage');LUNG3.pc$Age <- as.numeric(LUNG3.pc$Age)
LUNG3.pc$Gender <- revalue(LUNG3.pc$Gender, c("M"='Male',"F"='Female'))

LUNG4.pc <- LUNG4.info %>% dplyr::select('date of birth:ch1','date of surgery:ch1','gender:ch1','final patient stage:ch1')
LUNG4.pc <- add_column(.data =LUNG4.pc, Age=as.numeric(as.Date(LUNG4.info$`date of surgery:ch1`)-as.Date(LUNG4.info$`date of birth:ch1`))/365,.before = 'date of birth:ch1')
LUNG4.pc <- LUNG4.pc[,-c(2,3)]
colnames(LUNG4.pc) <- c('Age','Gender','Pstage')
LUNG4.pc$Gender <- revalue(LUNG4.pc$Gender, c("M"='Male',"F"='Female'))

PAN.LUNG.pc <- PAN.LUNG.info %>% dplyr::select('age_at_initial_pathologic_diagnosis','gender','ajcc_pathologic_tumor_stage')
colnames(PAN.LUNG.pc) <- c('Age','Gender','Pstage')
PAN.LUNG.pc$Gender <- revalue(PAN.LUNG.pc$Gender, c("MALE"='Male',"FEMALE"='Female'))
PAN.LUNG.pc$Pstage <- sub(".*? ", "", PAN.LUNG.pc$Pstage)
PAN.LUNG.pc$Pstage[PAN.LUNG.pc$Pstage %in% c('[Discrepancy]','Available]')] <- NA


META.LUNG.info <- rbind(LUNG1.pc,LUNG2.pc,LUNG3.pc,LUNG4.pc,PAN.LUNG.pc)
META.LUNG.info$Pstage <- factor(ifelse(META.LUNG.info$Pstage %in% c('1A','1B','I','IA','IB'),'I',
                                       ifelse(META.LUNG.info$Pstage %in% c('2A','2B','II','IIA','IIB'),'II',
                                              ifelse(META.LUNG.info$Pstage %in% c('IIIA','IIIB'),'III',
                                                     ifelse(META.LUNG.info$Pstage %in% c('IV'),'IV',NA)))),levels = c('I','II','III','IV'))

META.LUNG.info$Gender <- (as.factor(META.LUNG.info$Gender))
ST_Lung <- tableone::CreateTableOne(vars = c('Age','Gender','Pstage'),
                                    data = META.LUNG.info,includeNA = T,addOverall = T)
Lung.table <- print(ST_Lung,quote = FALSE, noSpaces = TRUE, printToggle = FALSE)


###### CRC

CRC1.pc <- CRC1.info %>% dplyr::select('age','gender','PStaging')
CRC2.pc <- CRC2.info %>% dplyr::select('age','gender','PStaging')
CRC4.pc <- CRC4.info %>% dplyr::select('age','gender:ch1','ajcc_stage:ch1')
CRC4.pc$`gender:ch1` <- revalue(CRC4.pc$`gender:ch1`, c("male"='MALE',"female"='FEMALE'));colnames(CRC4.pc) <- c('age','gender','PStaging')
CRC.PANCAN.pc <- CRC.PANCAN.info %>% dplyr::select('age','gender','PStaging')
META.CRC.info <- rbind(CRC1.pc,CRC2.pc,CRC4.pc,CRC.PANCAN.pc);colnames(META.CRC.info) <- c('Age','Gender','Pstage')
META.CRC.info$Pstage <- as.character(as.roman(as.numeric(META.CRC.info$Pstage)))

ST_CRC <- tableone::CreateTableOne(vars = c('Age','Gender','Pstage'),
                                   data = META.CRC.info,includeNA = T,addOverall = T)
CRC.table <- print(ST_CRC,quote = FALSE, noSpaces = TRUE, printToggle = FALSE)


####### SKCM cancer

patients.META2 <- colnames(MELA2.exprs[,-c(which(colnames(MELA2.exprs) %in% colnames(MELA3.exprs)))])

MELA1.pc <- data.frame(Age=rep(NA,81),Gender=rep(NA,81),Pstage=rep(NA,81),Metas=rep(NA,81),RECIST=rep(NA,81)) 

MELA2.pc <- MELA2.info[patients.META2,] %>% dplyr::select('age_start','gender','M','RECIST')
colnames(MELA2.pc) <- c('Age','Gender','Metas','RECIST');MELA2.pc$Pstage <- NA
MELA2.pc <- MELA2.pc[,c(1,2,5,3,4)];MELA2.pc$Age <- as.numeric(MELA2.pc$Age)
MELA2.pc$Gender <- revalue(MELA2.pc$Gender , c("male"='MALE',"female"='FEMALE'))
MELA2.pc$RECIST[MELA2.pc$RECIST=='X'] <- NA

MELA3.pc <- MELA3.info %>% dplyr::select('gender (Male=1, Female=0)','Mstage (IIIC=0, M1a=1, M1b=2, M1c=3)','PD1_response')
MELA3.pc$Age <- NA;MELA3.pc$Stage <- NA
MELA3.pc <- MELA3.pc[c(4,1,5,2,3)];colnames(MELA3.pc) <- c('Age','Gender','Pstage','Metas','RECIST')
MELA3.pc$Metas <- ifelse(MELA3.pc$Metas == 0,'IIIC',
                         ifelse(MELA3.pc$Metas == 1,'M1a',
                                ifelse(MELA3.pc$Metas == 2, 'M1b',
                                       ifelse(MELA3.pc$Metas == 3, 'M1c',NA))))

MELA3.pc$Gender <- ifelse(MELA3.pc$Gender == 1,'MALE',
                          ifelse(MELA3.pc$Gender == 0,'FEMALE',NA))


MELA.PANCAN.pc <- pData.PANCAN.MELA %>% dplyr::select('age_at_initial_pathologic_diagnosis','gender','ajcc_pathologic_tumor_stage')

MELA.PANCAN.pc$ajcc_pathologic_tumor_stage <- sub(".*? ", "", MELA.PANCAN.pc$ajcc_pathologic_tumor_stage)
MELA.PANCAN.pc$ajcc_pathologic_tumor_stage[MELA.PANCAN.pc$ajcc_pathologic_tumor_stage %in% c('[Discrepancy]','Available]')] <- NA
colnames(MELA.PANCAN.pc) <-  c('Age','Gender','Pstage');MELA.PANCAN.pc$Age <- as.numeric(MELA.PANCAN.pc$Age)
MELA.PANCAN.pc$Metas <- NA;MELA.PANCAN.pc$RECIST <- NA
MELA.PANCAN.pc$Pstage[MELA.PANCAN.pc$Pstage=="NOS"] <- NA
META.SKCM.info <- rbind(MELA1.pc,MELA2.pc,MELA3.pc,MELA.PANCAN.pc)
META.SKCM.info[c('Age','Gender','Pstage')]

META.SKCM.info$Pstage <- factor(ifelse(META.SKCM.info$Pstage %in% c('I','IB'),'I',
                                       ifelse(META.SKCM.info$Pstage %in% c('II','IIA','IIB','IIC'),'II',
                                              ifelse(META.SKCM.info$Pstage %in% c('III','IIIA','IIIB','IIIC'),'III',
                                                     ifelse(META.SKCM.info$Pstage %in% c('IV'),'IV',NA)))),levels = c('I','II','III','IV'))


ST_SKCM <- tableone::CreateTableOne(vars = c('Age','Gender','Pstage'),
                                    data = META.SKCM.info,includeNA = T,addOverall = T)
SKCM.table <- print(ST_SKCM,quote = FALSE, noSpaces = TRUE, printToggle = FALSE)


LUNG.FP.data <- merge(LUNG.info,META.LUNG.info,by='row.names',sort=F,all.x=T);rownames(LUNG.FP.data) <- LUNG.FP.data$Row.names;LUNG.FP.data <- LUNG.FP.data[,-1]

### NSCLC

model.NSCLC <- coxph(Surv(LUNG.FP.data$OS.year, LUNG.FP.data$OS) ~ Age+Gender+Pstage+B.lineage*BCM, data= LUNG.FP.data)
m.NSCLC<- c(NA,NA,summary(model.NSCLC)$coeff[1,2],
            NA,NA,summary(model.NSCLC)$coeff[2,2],
            NA,NA,summary(model.NSCLC)$coeff[3:5,2],
            NA,summary(model.NSCLC)$coeff[6:8,2])

l.NSCLC<- c(NA,NA,format(signif(summary(model.NSCLC)$conf.int[1,3],digits = 3)),
            NA,NA,format(signif(summary(model.NSCLC)$conf.int[2,3],digits = 3)),
            NA,NA,format(signif(summary(model.NSCLC)$conf.int[3:5,3],digits = 3)),
            NA,format(signif(summary(model.NSCLC)$conf.int[6:8,3],digits = 3)))

u.NSCLC<- c(NA,NA,format(signif(summary(model.NSCLC)$conf.int[1,4],digits = 3)),
            NA,NA,format(signif(summary(model.NSCLC)$conf.int[2,4],digits = 3)),
            NA,NA,format(signif(summary(model.NSCLC)$conf.int[3:5,4],digits = 3)),
            NA,format(signif(summary(model.NSCLC)$conf.int[6:8,4],digits = 3)))


np.NSCLC <- c(paste(c(length(LUNG.FP.data$Age),NA,table(LUNG.FP.data$Gender),NA,table(LUNG.FP.data$Pstage,useNA='always'),length(LUNG.FP.data$B.lineage),length(LUNG.FP.data$BCM),length(LUNG.FP.data$BCM))," (",
                    signif(c(length(LUNG.FP.data$Age),NA,table(LUNG.FP.data$Gender),NA,table(LUNG.FP.data$Pstage,useNA='always'),length(LUNG.FP.data$B.lineage),length(LUNG.FP.data$BCM),length(LUNG.FP.data$BCM))/nrow(LUNG.FP.data)*100,digits=1),")",sep=""))
np.NSCLC[np.NSCLC=='NA (NA)'] <- NA

tabletext.NSCLC<-cbind(c(c(NA,'META-NSCLC','Age','Gender','    Female','    Male','Tumour Stage','    I','    II','    III','    IV','    Missing cases','B cells','BCM','B cells * BCM')),
                       c("No. of Patients (%)",NA,np.NSCLC), 
                       c("HR",NA,format((summary(model.NSCLC)$coeff[,2]),digits=2)[1],NA,'Ref',format((summary(model.NSCLC)$coeff[,2])[2],digits=3),NA,'Ref',format((summary(model.NSCLC)$coeff[,2])[c(3:5)],digits=3),NA,format((summary(model.NSCLC)$coeff[,2])[c(6:8)],digits=3)),
                       c('95% CI',NA,paste0('[',l.NSCLC[3],' - ',u.NSCLC[3],']',sep = ''),
                         NA,'-',paste0('[',l.NSCLC[6],' - ',u.NSCLC[6],']',sep = ''),
                         NA,'-',paste0('[',l.NSCLC[c(9:11)],' - ',u.NSCLC[c(9:11)],']',sep = ''),NA,
                         paste0('[',l.NSCLC[c(13:15)],' - ',u.NSCLC[c(13:15)],']',sep = '')),
                       c("Coeff",NA,format((summary(model.NSCLC)$coeff[,1]),digits=2)[1],NA,'-',format((summary(model.NSCLC)$coeff[,1])[2],digits=3),NA,'-',format((summary(model.NSCLC)$coeff[,1])[c(3:5)],digits=3),NA,format((summary(model.NSCLC)$coeff[,1])[c(6:8)],digits=2)),
                       c("p.",NA,format((summary(model.NSCLC)$coeff[,5]),digits=2)[1],NA,'-',format((summary(model.NSCLC)$coeff[,5])[2],digits=3),NA,'-',format((summary(model.NSCLC)$coeff[,5])[c(3:5)],digits=2),NA,format((summary(model.NSCLC)$coeff[,5])[c(6:8)],digits=1)))
tabletext.NSCLC[as.numeric(tabletext.NSCLC[,6]) < 0.001,6] <- '< 0.001'

########## CRC

CRC.FP.data <- merge(CRC.info,META.CRC.info,by='row.names',sort=F,all.x=T);rownames(CRC.FP.data) <- CRC.FP.data$Row.names;CRC.FP.data <- CRC.FP.data[,-1]

model.CRC <- coxph(Surv(CRC.FP.data$OS.year, CRC.FP.data$OS) ~ Age+Gender+Pstage+B.lineage*BCM, data= CRC.FP.data)
m.CRC<- c(NA,NA,summary(model.CRC)$coeff[1,2],
          NA,NA,summary(model.CRC)$coeff[2,2],
          NA,NA,summary(model.CRC)$coeff[3:6,2],
          NA,summary(model.CRC)$coeff[7:9,2])

l.CRC<- c(NA,NA,format(signif(summary(model.CRC)$conf.int[1,3],digits = 3)),
          NA,NA,format(signif(summary(model.CRC)$conf.int[2,3],digits = 3)),
          NA,NA,format(signif(summary(model.CRC)$conf.int[3:6,3],digits = 3)),
          NA,format(signif(summary(model.CRC)$conf.int[7:9,3],digits = 3)))

u.CRC<- c(NA,NA,format(signif(summary(model.CRC)$conf.int[1,4],digits = 3)),
          NA,NA,format(signif(summary(model.CRC)$conf.int[2,4],digits = 3)),
          NA,NA,format(round(summary(model.CRC)$conf.int[3:6,4],digits = 3)),
          NA,format(round(summary(model.CRC)$conf.int[7:9,4],digits = 3)))


np.CRC <- c(paste(c(length(CRC.FP.data$Age),NA,table(CRC.FP.data$Gender),NA,table(CRC.FP.data$Pstage,useNA='always'),length(CRC.FP.data$B.lineage),length(CRC.FP.data$BCM),length(CRC.FP.data$BCM))," (",
                  signif(c(length(CRC.FP.data$Age),NA,table(CRC.FP.data$Gender),NA,table(CRC.FP.data$Pstage,useNA='always'),length(CRC.FP.data$B.lineage),length(CRC.FP.data$BCM),length(CRC.FP.data$BCM))/nrow(CRC.FP.data)*100,digits=1),")",sep=""))
np.CRC[np.CRC=='NA (NA)'] <- NA

tabletext.CRC<-cbind(c(c(NA,'META-CRC','Age','Gender','    Female','    Male','Tumour Stage','    I','    II','    III','    IV','    V','    Missing cases','B cells','BCM','B cells * BCM')),
                     c("No. of Patients (%)",NA,np.CRC), 
                     c("HR",NA,format((summary(model.CRC)$coeff[,2]),digits=3)[1],NA,'Ref',format((summary(model.CRC)$coeff[,2])[2],digits=3),NA,'Ref',format((summary(model.CRC)$coeff[,2])[c(3:6)],digits=3),NA,format((summary(model.CRC)$coeff[,2])[c(7:9)],digits=3)),
                     c('95% CI',NA,paste0('[',l.CRC[3],' - ',u.CRC[3],']',sep = ''),
                       NA,'-',paste0('[',l.CRC[6],' - ',u.CRC[6],']',sep = ''),
                       NA,'-',paste0('[',l.CRC[c(9:12)],' - ',u.CRC[c(9:12)],']',sep = ''),NA,
                       paste0('[',l.CRC[c(14:16)],' - ',u.CRC[c(14:16)],']',sep = '')),
                     c("Coeff",NA,format((summary(model.CRC)$coeff[,1]),digits=3)[1],NA,'-',format((summary(model.CRC)$coeff[,1])[2],digits=3),NA,'-',format((summary(model.CRC)$coeff[,1])[c(3:6)],digits=3),NA,format((summary(model.CRC)$coeff[,1])[c(7:9)],digits=2)),
                     c("p.",NA,format((summary(model.CRC)$coeff[,5]),digits=3)[1],NA,'-',format((summary(model.CRC)$coeff[,5])[2],digits=3),NA,'-',format((summary(model.CRC)$coeff[,5])[c(3:6)],digits=2),NA,format((summary(model.CRC)$coeff[,5])[c(7:9)],digits=3)))
tabletext.CRC[as.numeric(tabletext.CRC[,6]) < 0.001,6] <- '< 0.001'


########## SKCM

SKCM.FP.data <- merge(MELA.info,META.SKCM.info,by='row.names',sort=F,all.x=T);rownames(SKCM.FP.data) <- SKCM.FP.data$Row.names;SKCM.FP.data <- SKCM.FP.data[,-1]
model.SKCM <- coxph(Surv(SKCM.FP.data$OS.year, SKCM.FP.data$OS) ~ B.lineage*BCM, data= SKCM.FP.data)

m.SKCM<- c(NA,summary(model.SKCM)$coeff[1,2],
           summary(model.SKCM)$coeff[2:3,2])

l.SKCM<- c(format(signif(summary(model.SKCM)$conf.int[1,3],digits = 3)),
           format(signif(summary(model.SKCM)$conf.int[2:3,3],digits = 3)))

u.SKCM<- c('1.000',
           format(round(summary(model.SKCM)$conf.int[2:3,4],digits = 3)))

np.SKCM <- c(paste(c(length(SKCM.FP.data$B.lineage),length(SKCM.FP.data$BCM),length(SKCM.FP.data$BCM))," (",
                   signif(c(length(SKCM.FP.data$B.lineage),length(SKCM.FP.data$BCM),length(SKCM.FP.data$BCM))/nrow(SKCM.FP.data)*100,digits=1),")",sep=""))
np.SKCM[np.SKCM=='NA (NA)'] <- NA

tabletext.SKCM<-cbind(c(c(NA,'META-SKCM','B cells','BCM','B cells * BCM')),
                      c("No. of Patients (%)",NA,np.SKCM), 
                      c("HR",NA,format((summary(model.SKCM)$coeff[,2])[1],digits=3),format((summary(model.SKCM)$coeff[,2])[c(2:3)],digits=3)),
                      c('95% CI',NA,paste0('[',l.SKCM[1],' - ',u.SKCM[1],']',sep = ''),
                        paste0('[',l.SKCM[c(2:3)],' - ',u.SKCM[c(2:3)],']',sep = '')),
                      c("Coeff",NA,format((summary(model.SKCM)$coeff[,1])[1],digits=3),format((summary(model.SKCM)$coeff[,1])[c(2:3)],digits=2)),
                      c("p.",NA,signif((summary(model.SKCM)$coeff[,5])[1],digits=2),signif((summary(model.SKCM)$coeff[,5])[c(2:3)],digits=2)))
tabletext.SKCM[as.numeric(tabletext.SKCM[,6]) < 0.001,6] <- '< 0.001'



tabletext.FP <- rbind(tabletext.NSCLC,tabletext.CRC[-1,],tabletext.SKCM[-1,])
lower.FP <- c(l.NSCLC,l.CRC[-c(1)],NA,l.SKCM)
upper.FP <- c(u.NSCLC,u.CRC[-c(1)],NA,u.SKCM)
mean.FP <- c(m.NSCLC,m.CRC[-c(1)],m.SKCM)

tabletext.FP.SS <- tabletext.FP[c(1:2,13:15,16,28:29,30:34),]
lower.FP.SS <- lower.FP[c(1:2,13:15,16,28:29,30:34)]
upper.FP.SS <- upper.FP[c(1:2,13:15,16,28:29,30:34)]
mean.FP.SS <- mean.FP[c(1:2,13:15,16,28:29,30:34)]

tabletext.FP.SS[c(4:10,13),6] <- format(round(as.numeric(tabletext.FP.SS[c(4:10,13),6]),digits=3))
tabletext.FP.SS[,6][tabletext.FP.SS[,6]=="   NA"] <- NA
tabletext.FP.SS[-1,5] <- format(signif(as.numeric(tabletext.FP.SS[-1,5]),digits=1))
tabletext.FP.SS[,5][tabletext.FP.SS[,5]=="    NA"] <- NA
tabletext.FP.SS[,2] <- c("No. of Patients (%)","1247 (100)",rep(NA,3),"1247 (100)",rep(NA,3),"325 (100)",rep(NA,3))


### Figure 3C - forest plot

forestplot::forestplot(tabletext.FP.SS,mean.FP.SS,as.numeric(lower.FP.SS),as.numeric(upper.FP.SS),zero=1,
                       title="Hazard Ratio",lwd.ci=3, ci.vertices=TRUE, ci.vertices.height = 0.05,
                       hrzl_lines = list("2" = gpar(lwd=3, columns=1:7, col = "black"),
                                         "3" = gpar(lwd=1, columns=1, col = "black"),
                                         "7" = gpar(lwd=1, columns=1, col = "black"),
                                         "11" = gpar(lwd=1, columns=1, col = "black"),
                                         "14" = gpar(lwd=1, columns=1:6, col = "black")),lwd.zero=5,
                       
                       clip=c(0,5),xticks=seq(from = 0.5,to = 1.5,by=0.2),graphwidth=unit(24,"cm"),
                       is.summary=c(rep(TRUE,2),rep(FALSE,3),TRUE,rep(FALSE,3),TRUE,rep(FALSE,6)),boxsize = .3,
                       colgap=unit(25,"mm"),lineheight = 'auto',align = c("l", "c", "c","c","c","c"),
                       xlog=F,txt_gp = fpTxtGp(label = list(gpar(fontfamily = "Arial"),
                                                            gpar(fontfamily = "")),cex = 3,
                                               ticks = gpar(fontfamily = "", cex=3),
                                               xlab  = gpar(fontfamily = "HersheySerif", cex = 4)),col=fpColors(box=c("blue"),lines='black',hrz_lines = "gray50"))
grid.text('Note: Models were adjusted for clinical variables (Age, Gender, Pathological stage) whenever available.',
          x = unit(0.34, 'npc'),gp=gpar(col="black", fontsize=33),
          y = unit(1, 'lines'))




# Figure 4


MELA2.info$BCM <- base::rowMeans(t(MELA2.exprs[BCM.gene$ENTREZID,]))
MELA2.info$BCM.bin <- relevel(as.factor(ifelse(ntile(MELA2.info$BCM,n = 2)=='1','Low','High')),ref = 'Low')

MELA2.info$CTLA4.2group <- as.factor(ifelse(MELA2.info$response=='nonresponse','PD',
                                            ifelse(MELA2.info$response%in%c('response','long-survival'),'non-PD',NA)))


# waterfall plot and KM - anti-CTLA4

tmp.water.CTLA4 <- MELA2.info[c('response','CTLA4.2group','BCM','BCM.bin','OS','OS.year')]
tmp.water.CTLA4$BCM.rescaled  <- scale(tmp.water.CTLA4$BCM, center = T,scale=T)
tmp.water.CTLA4 <- tmp.water.CTLA4[order(tmp.water.CTLA4$BCM.rescaled),]
tmp.water.CTLA4$id <- as.factor(1:nrow(tmp.water.CTLA4))
tmp.bp <- data.frame(table(tmp.water.CTLA4$response));tmp.bp$total <- 'All patients'
tmp.bp.2group <- data.frame(table(tmp.water.CTLA4$CTLA4.2group));tmp.bp.2group$total <- 'All patients'

##  Testing association of BCM with worse ICB response
tmp.water.CTLA4$id2 <- as.factor(c(rep('btm',20),rep('top',20)))

tmp.bplot.CTLA4.2group <- ggplot(tmp.bp.2group, aes(fill=Var1, y=Freq, x=total,label=Freq)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(labels= c('Responders','Progressors'),
                    values = c("springgreen4","firebrick3"))+
  geom_text(size = 6, position = position_stack(vjust = 0.5))+
  theme_few(base_size = 12)+
  labs(fill='',
       x ="", y = "All patients (n= 40)")+
  theme(axis.line.x  = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.text.x=element_blank(),axis.text.y=element_blank(),
        plot.title = element_text(),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"))


CTLA4.wp.2 <- ggplot(tmp.water.CTLA4, aes(x=id,y=na.omit(BCM.rescaled), fill = CTLA4.2group)) + 
  geom_bar(stat='identity')+
  labs(x = "Patients", y = "BCM score",fill='Anti-CTLA4 response')+
  geom_segment(aes(x= 10,xend=10,y=0.7,yend=0.7))+
  geom_segment(aes(x= 30,xend=30,y=0.7,yend=0.7))+
  scale_fill_manual(labels= c('Responders','Progressors'),
                    values = c("springgreen4","firebrick3"))+
  geom_vline(xintercept = 20.5, colour="violetred4", linetype = "longdash")+theme_few(base_size = 12)+
  ggtitle('')+
  theme(axis.line.x  = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x=element_blank(),
        plot.title = element_text(),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"))

CTLA4.boxplot.2 <- ggboxplot(tmp.water.CTLA4,x='CTLA4.2group',y='BCM.rescaled',
                             bxp.errorbar = TRUE,bxp.errorbar.width=0.15,
                             fill='CTLA4.2group',error.plot = 'errorbar',
                             size = 0.2,
                             palette = c("springgreen4","firebrick3"),
                             add = "jitter", 
                             ggtheme = theme_few())+
  stat_compare_means(method = 't.test',tip.length = 0.01,label.x = 2.1,label='p.format',size=5,label.y=2.2)+
  scale_x_discrete(labels=c('Responders','Progressors')) +
  labs(x='Anti-CTLA4 response',y='BCM score',fill='') +
  theme(axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 15, face = "bold",hjust = c(0.5,0.5)),
        axis.text.x =element_text(size=10,face = 'bold'),
        axis.text.y =element_text(size=10,face = 'bold',hjust = c(1,1)),
        strip.background = element_rect(fill='whitesmoke'),strip.text = element_text(face = 'bold',size = 10),
        legend.key = element_rect(fill='white'),legend.position = 'none',
        legend.box.background = element_rect(fill='black'),legend.key.size =unit(0.4, "cm"))


# PERFORMANCE TEST WITH WILCOXON and AUC
# Compare BCM with IPS, IFGN, CD8, PDL1, CRMA, Total count of non-synonymous mutations in a tumor

IFGN.sig <-  biomaRt::select(org.Hs.eg.db,keys =c('IFNG','STAT1','IDO1','CXCL10','CXCL9','HLA-DRA'),
                             columns = c('ENTREZID','SYMBOL'),keytype = 'SYMBOL')
CD8.sig <- biomaRt::select(org.Hs.eg.db,keys =c('CD8A','CD8B'),
                           columns = c('ENTREZID','SYMBOL'),keytype = 'SYMBOL')

PDL1.sig <- biomaRt::select(org.Hs.eg.db,keys =c('CD274'),
                            columns = c('ENTREZID','SYMBOL'),keytype = 'SYMBOL')

CRMA.sig <- biomaRt::select(org.Hs.eg.db,keys =c('MAGEA2','MAGEA2B','MAGEA3','MAGEA6','MAGEA12'),
                            columns = c('ENTREZID','SYMBOL'),keytype = 'SYMBOL')

IPS.sig <- biomaRt::select(org.Hs.eg.db,keys =c('GZMA','PRF1'),
                           columns = c('ENTREZID','SYMBOL'),keytype = 'SYMBOL')

MELA2.info$IFGN <- base::rowMeans(t(MELA2.exprs[which(rownames(MELA2.exprs) %in% IFGN.sig$ENTREZID),]))
MELA2.info$CD8 <- base::rowMeans(t(MELA2.exprs[which(rownames(MELA2.exprs) %in% CD8.sig$ENTREZID),]))
MELA2.info$PDL1 <- c(t(MELA2.exprs[which(rownames(MELA2.exprs) %in% PDL1.sig$ENTREZID),]))
MELA2.info$CRMA <- base::rowMeans(t(MELA2.exprs[which(rownames(MELA2.exprs) %in% CRMA.sig$ENTREZID),]))
MELA2.info$IPS <- base::rowMeans(t(MELA2.exprs[which(rownames(MELA2.exprs) %in% IPS.sig$ENTREZID),]))
MELA2.info$MUT <- MELA2.info$nonsynonymous
MELA2.info$T.cells <- as.numeric(MCPcounter::MCPcounter.estimate(MELA2.exprs,featuresType = 'ENTREZ_ID')['T cells',])


CTLA4.table <- MELA2.info[c('response','CTLA4.2group','BCM.bin','OS','OS.year','BCM','T.cells','IFGN','CD8','PDL1','CRMA','IPS','MUT')]
CTLA4.table$BCM.rescaled  <- scale(CTLA4.table$BCM, center = T,scale=T)
CTLA4.table <- CTLA4.table[order(CTLA4.table$BCM.rescaled),]
CTLA4.table$id <- as.factor(1:nrow(CTLA4.table))


#### run through all signatures column 6-12
MELA2.CTLA4.performance <- data.frame(foreach(i=c(colnames(CTLA4.table)[7:13]),.combine='rbind') %do% {
  w.df <- wilcox.test(CTLA4.table[,i]~CTLA4.table$CTLA4.2group,alternative = "two.sided")$p.value
  w.df <- ifelse(w.df < 0.001,'< 0.001',round(w.df,digits=3))
  auc.df <-  round(as.numeric(pROC::roc(CTLA4.table$CTLA4.2group,CTLA4.table[,i])$auc),digits=4)
  auc.df.l <- as.numeric(pROC::ci.auc(CTLA4.table$CTLA4.2group,CTLA4.table[,i])[1])
  auc.df.u <- as.numeric(pROC::ci.auc(CTLA4.table$CTLA4.2group,CTLA4.table[,i])[3])
  c(patients='All',Responders='18',Progressors='22',w.p=w.df,AUC=as.numeric(auc.df),AUC.l=as.numeric(auc.df.l),AUC.u=as.numeric(auc.df.u))},row.names = colnames(CTLA4.table)[7:13])


MELA2.CTLA4.onlyBCM.performance <- data.frame(t(foreach(i=c('BCM'),.combine='rbind') %do% {
  w.df <- t.test(CTLA4.table[,i]~CTLA4.table$CTLA4.2group,alternative = "two.sided")$p.value
  w.df <- ifelse(w.df < 0.001,'< 0.001',round(w.df,digits=3))
  auc.df <-  round(as.numeric(pROC::roc(CTLA4.table$CTLA4.2group,CTLA4.table[,i])$auc),digits=4)
  auc.df.l <- as.numeric(pROC::ci.auc(CTLA4.table$CTLA4.2group,CTLA4.table[,i])[1])
  auc.df.u <- as.numeric(pROC::ci.auc(CTLA4.table$CTLA4.2group,CTLA4.table[,i])[3])
  c(patients='All',Responders='18',Progressors='22',w.p=w.df,AUC=as.numeric(auc.df),AUC.l=as.numeric(auc.df.l),AUC.u=as.numeric(auc.df.u))}),row.names = c('BCM'))

MELA2.CTLA4.BCM.performance <- rbind(MELA2.CTLA4.onlyBCM.performance,MELA2.CTLA4.performance)
MELA2.CTLA4.BCM.performance <- MELA2.CTLA4.BCM.performance %>% add_column(Signatures=colnames(CTLA4.table)[6:13],.before='patients')
MELA2.CTLA4.BCM.performance$AUC <- as.numeric(MELA2.CTLA4.BCM.performance$AUC)
MELA2.CTLA4.BCM.performance$AUC.l <- as.numeric(MELA2.CTLA4.BCM.performance$AUC.l)
MELA2.CTLA4.BCM.performance$AUC.u <- as.numeric(MELA2.CTLA4.BCM.performance$AUC.u)

# Signatures patients Responders Progressors   w.p    AUC     AUC.l     AUC.u
# BCM         BCM      All         18          22 0.024 0.7323 0.5737328 0.8909137
# T.cells T.cells      All         18          22 0.013 0.7298 0.5699435 0.8896525
# IFGN       IFGN      All         18          22 0.032 0.6995 0.5314384 0.8675515
# CD8         CD8      All         18          22 0.045 0.6869 0.5203285 0.8534089
# PDL1       PDL1      All         18          22  0.18 0.6263 0.4481426 0.8043826
# CRMA       CRMA      All         18          22 0.034 0.6970 0.5257480 0.8681914
# IPS         IPS      All         18          22 0.022 0.7121 0.5495980 0.8746444
# MUT         MUT      All         18          22 0.064 0.6730 0.5008654 0.8450942

CTLA.AUC.bar <- ggbarplot(MELA2.CTLA4.BCM.performance,x='Signatures',y='AUC',fill='Signatures',
                          color='white',palette='jco',lab.size = 4, font.tickslab = c(10,'bold'),font.x = c(12, "bold"),
                          font.y = c(12, "bold"),ylim=c(0.25,0.9),label=format(MELA2.CTLA4.BCM.performance$AUC,digits=3),
                          x.text.angle = 90,xlab = 'Gene signatures',ggtheme = theme_few())+
  geom_errorbar(aes(ymin=MELA2.CTLA4.BCM.performance$AUC.l, ymax=MELA2.CTLA4.BCM.performance$AUC.u), width=.15,alpha=0.25,color='black',
                position=position_dodge()) +
  labs(fill='')+
  geom_hline(yintercept=0.5, linetype="dashed", 
             color = "black", size=0.5)+theme(legend.position = "none")

##################### Anti-PD1 analyses 
MELA3.info$BCM <- base::rowMeans(t(MELA3.exprs[BCM.gene$ENTREZID,]))
MELA3.info$BCM.bin <- relevel(as.factor(ifelse(ntile(MELA3.info$BCM,n = 2)=='1','Low','High')),ref = 'Low')


# Waterfall plot - PD1 test

tmp.water.PD1 <- MELA3.info[c('PD1_response','PD1.2group','priorCTLA4','BCM','BCM.bin','OS','OS.year')]
tmp.water.PD1$BCM.rescaled  <- scale(tmp.water.PD1$BCM, center = T,scale=T)
tmp.water.PD1 <- tmp.water.PD1[order(tmp.water.PD1$BCM.rescaled),]
tmp.water.PD1$id <- as.factor(1:nrow(tmp.water.PD1))

tmp.water.2group <- tmp.water.PD1[!is.na(tmp.water.PD1$PD1.2group),]
tmp.water.2group$id <- as.factor(1:nrow(tmp.water.2group))
tmp.bp <- data.frame(table(tmp.water.PD1$PD1_response));tmp.bp$total <- 'All patients'
tmp.bp.2group <- data.frame(table(tmp.water.2group$PD1.2group));tmp.bp.2group$total <- 'All patients'
tmp.water.2group$BCM.bin <- relevel(as.factor(ifelse(ntile(tmp.water.2group$BCM.rescaled,n = 2)=='1','Low','High')),ref = 'Low')


tmp.bplot.PD1.2group <- ggplot(tmp.bp.2group, aes(fill=Var1, y=Freq, x=total,label=Freq)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(labels= c('Responders','Progressors'),
                    values = c("springgreen4","firebrick3"))+
  geom_text(size = 6, position = position_stack(vjust = 0.5))+
  theme_few(base_size = 12)+
  labs(fill='',
       x ="", y = "All patients (n= 103)")+
  theme(axis.line.x  = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.text.x=element_blank(),axis.text.y=element_blank(),
        plot.title = element_text(),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"))


PD1.wp.2 <- ggplot(tmp.water.2group, aes(x=id,y=(BCM.rescaled), fill = PD1.2group)) + 
  geom_bar(stat='identity')+
  labs(x = "Patients", y = "BCM score",fill='Anti-PD1 response')+
  geom_segment(aes(x= 103/4,xend=103/4,y=2,yend=2))+
  geom_segment(aes(x= 103/1.25,xend=103/1.25,y=2,yend=2))+
  scale_fill_manual(labels= c('Responders','Progressors'),
                    values = c("springgreen4","firebrick3"))+
  geom_vline(xintercept = 52.5, colour="violetred4", linetype = "longdash")+theme_few(base_size = 12)+
  ggtitle('')+
  theme(axis.line.x  = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x=element_blank(),
        plot.title = element_text(),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"))


PD1.boxplot.2 <- ggboxplot(tmp.water.2group,x='PD1.2group',y='BCM.rescaled',
                           bxp.errorbar = TRUE,bxp.errorbar.width=0.15,
                           fill='PD1.2group',error.plot = 'errorbar',
                           size = 0.2,
                           palette = c("springgreen4","firebrick3"),
                           add = "jitter", 
                           ggtheme = theme_few())+
  stat_compare_means(method = 't.test',tip.length = 0.01,label.x = 2.1,label='p.format',size=5,label.y = 2.5)+
  scale_x_discrete(labels=c('Responders','Progressors')) +
  labs(x='Anti-PD1 response',y='BCM score',fill='') +
  theme(axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 15, face = "bold",hjust = c(0.5,0.5)),
        axis.text.x =element_text(size=10,face = 'bold'),
        axis.text.y =element_text(size=10,face = 'bold',hjust = c(1,1)),
        strip.background = element_rect(fill='whitesmoke'),strip.text = element_text(face = 'bold',size = 10),
        legend.key = element_rect(fill='white'),legend.position = 'none',
        legend.box.background = element_rect(fill='black'),legend.key.size =unit(0.4, "cm"))


MELA3.info$IFGN <- base::rowMeans(t(MELA3.exprs[which(rownames(MELA3.exprs) %in% IFGN.sig$ENTREZID),]))
MELA3.info$CD8 <- base::rowMeans(t(MELA3.exprs[which(rownames(MELA3.exprs) %in% CD8.sig$ENTREZID),]))
MELA3.info$PDL1 <- c(t(MELA3.exprs[which(rownames(MELA3.exprs) %in% PDL1.sig$ENTREZID),]))
MELA3.info$CRMA <- base::rowMeans(t(MELA3.exprs[which(rownames(MELA3.exprs) %in% CRMA.sig$ENTREZID),]))
MELA3.info$IPS <- base::rowMeans(t(MELA3.exprs[which(rownames(MELA3.exprs) %in% IPS.sig$ENTREZID),]))
MELA3.info$MUT <- MELA3.info$nonsyn_muts
MELA3.info$T.cells <- as.numeric(MCPcounter::MCPcounter.estimate(MELA3.exprs,featuresType = 'ENTREZ_ID')['T cells',])

tmp <- MELA3.info[c('PD1_response','PD1.2group','priorCTLA4','BCM.bin','OS','OS.year','BCM','T.cells','IFGN','CD8','PDL1','CRMA','IPS','MUT')]
tmp$BCM.rescaled  <- scale(tmp$BCM, center = T,scale=T)
tmp <- tmp[order(tmp$BCM.rescaled),]
tmp$id <- as.factor(1:nrow(tmp))

PD1.table <- tmp[!is.na(tmp$PD1.2group),]
PD1.table$id <- as.factor(1:nrow(PD1.table))

MELA3.PD1.performance.all <- data.frame(foreach(i=c(colnames(PD1.table)[8:14]),.combine='rbind') %do% {
  w.df <- wilcox.test(PD1.table[,i]~PD1.table$PD1.2group,alternative = "two.sided")$p.value
  w.df <- ifelse(w.df < 0.001,'< 0.001',round(w.df,digits=3))
  auc.df <-  round(as.numeric(pROC::roc(PD1.table$PD1.2group,PD1.table[,i])$auc),digits=4)
  auc.df.l <- as.numeric(pROC::ci.auc(PD1.table$PD1.2group,PD1.table[,i])[1])
  auc.df.u <- as.numeric(pROC::ci.auc(PD1.table$PD1.2group,PD1.table[,i])[3])
  
  c(patients='All',Responders='47',Progressors='56',w.p=w.df,auc=auc.df,auc.l=auc.df.l,auc.u=auc.df.u)},row.names = colnames(PD1.table)[8:14])

MELA3.PD1.performance.ipi.naive <- data.frame(foreach(i=c(colnames(PD1.table)[8:14]),.combine='rbind') %do% {
  tmp.ipi.naive <- PD1.table[PD1.table$priorCTLA4==0,]
  w.df <- wilcox.test(tmp.ipi.naive[,i]~tmp.ipi.naive$PD1.2group,alternative = "two.sided")$p.value
  w.df <- ifelse(w.df < 0.001,'< 0.001',round(w.df,digits=4))
  auc.df <-  round(as.numeric(pROC::roc(tmp.ipi.naive$PD1.2group,tmp.ipi.naive[,i])$auc),digits=4)
  c(patients='ipi-naive',Responders='31',Progressors='33',w.p=w.df,auc=auc.df)},row.names = colnames(PD1.table)[8:14])

MELA3.PD1.performance.ipi.treated <- data.frame(foreach(i=c(colnames(PD1.table)[8:14]),.combine='rbind') %do% {
  tmp.ipi.treated <- PD1.table[PD1.table$priorCTLA4==1,]
  w.df <- wilcox.test(tmp.ipi.treated[,i]~tmp.ipi.treated$PD1.2group,alternative = "two.sided")$p.value
  w.df <- ifelse(w.df < 0.001,'< 0.001',round(w.df,digits=3))
  auc.df <-  round(as.numeric(pROC::roc(tmp.ipi.treated$PD1.2group,tmp.ipi.treated[,i])$auc),digits=4)
  c(patients='ipi-treated',Responders='16',Progressors='23',w.p=w.df,auc=auc.df)},row.names = colnames(PD1.table)[8:14])


MELA3.PD1.onlyBCM.performance.all <- data.frame(t(foreach(i='BCM',.combine='rbind') %do% {
  w.df <- t.test(PD1.table[,i]~PD1.table$PD1.2group,alternative = "two.sided")$p.value
  w.df <- ifelse(w.df < 0.001,'< 0.001',round(w.df,digits=3))
  auc.df <-  round(as.numeric(pROC::roc(PD1.table$PD1.2group,PD1.table[,i])$auc),digits=4)
  c(patients='All',Responders='47',Progressors='56',w.p=w.df,auc=auc.df,auc.l=auc.df.l,auc.u=auc.df.u)}),row.names ='BCM')

MELA3.PD1.onlyBCM.performance.ipi.naive <-  data.frame(t(foreach(i='BCM',.combine='rbind') %do% {
  tmp.ipi.naive <- PD1.table[PD1.table$priorCTLA4==0,]
  w.df <- t.test(tmp.ipi.naive[,i]~tmp.ipi.naive$PD1.2group,alternative = "two.sided")$p.value
  w.df <- ifelse(w.df < 0.001,'< 0.001',round(w.df,digits=3))
  auc.df <-  round(as.numeric(pROC::roc(tmp.ipi.naive$PD1.2group,tmp.ipi.naive[,i])$auc),digits=4)
  c(patients='ipi-naive',Responders='31',Progressors='33',w.p=w.df,auc=auc.df)}),row.names ='BCM')

MELA3.PD1.onlyBCM.performance.ipi.treated <- data.frame(t(foreach(i='BCM',.combine='rbind') %do% {
  tmp.ipi.treated <- PD1.table[PD1.table$priorCTLA4==1,]
  w.df <- t.test(tmp.ipi.treated[,i]~tmp.ipi.treated$PD1.2group,alternative = "two.sided")$p.value
  w.df <- ifelse(w.df < 0.001,'< 0.001',round(w.df,digits=3))
  auc.df <-  round(as.numeric(pROC::roc(tmp.ipi.treated$PD1.2group,tmp.ipi.treated[,i])$auc),digits=4)
  c(patients='ipi-treated',Responders='16',Progressors='23',w.p=w.df,auc=auc.df)}),row.names ='BCM')


MELA3.PD1.BCM.performance.all <- rbind(MELA3.PD1.onlyBCM.performance.all,MELA3.PD1.performance.all)
MELA3.PD1.BCM.performance.ipi.naive <- rbind(MELA3.PD1.onlyBCM.performance.ipi.naive,MELA3.PD1.performance.ipi.naive)
MELA3.PD1.BCM.performance.ipi.treated <- rbind(MELA3.PD1.onlyBCM.performance.ipi.treated,MELA3.PD1.performance.ipi.treated)

PD1.BCM <- cbind(MELA3.PD1.BCM.performance.all[,-c(6,7)],MELA3.PD1.BCM.performance.ipi.naive,MELA3.PD1.BCM.performance.ipi.treated)
PD1.BCM <- PD1.BCM %>% add_column(Signatures=colnames(PD1.table)[7:14],.before='patients')
bar.comparisons <- list(c('BCM','T.cells'),c('BCM','IFGN'),c('BCM','CD8'),c('BCM','PDL1'),c('BCM','CRMA'),c('BCM','IPS'),c('BCM','MUT'))

PD1.BCM$auc <- signif(as.numeric(PD1.BCM$auc),digits = 3)
PD1.AUC.bar <- ggbarplot(PD1.BCM,x='Signatures',y='auc',fill='Signatures',
                         color='white',palette='jco',lab.size = 4, font.tickslab = c(10,'bold'),font.x = c(12, "bold"),
                         font.y = c(12, "bold"),ylim=c(0.25,0.8),ylab='AUC',
                         x.text.angle = 90,xlab = 'Gene signatures',label = T,ggtheme = theme_few(),lab.nb.digits = 3)+
  geom_errorbar(aes(ymin=as.numeric(MELA3.PD1.BCM.performance.all$auc.l), ymax=as.numeric(MELA3.PD1.BCM.performance.all$auc.u)), width=.15,alpha=0.25,color='black',
                position=position_dodge())+
  labs(fill='')+
  geom_hline(yintercept=0.5, linetype="dashed", 
             color = "black", size=0.5)+theme(legend.position = "none")





# Main figure 4 

# Figure 4 A-F


ggpubr::ggarrange(tmp.bplot.PD1.2group,PD1.wp.2,PD1.boxplot.2,PD1.AUC.bar,
                  tmp.bplot.CTLA4.2group,CTLA4.wp.2,CTLA4.boxplot.2,CTLA.AUC.bar,common.legend = T,
                  ncol = 4, nrow = 2,  align = "h",widths=c(0.3,1.5,1.5,1.5),heights = c(1,1,1,1))



# Figure 4 KMs G-H

par(mfrow=c(1,3),omi=c(0, 0, 0, 0),tcl=-0.5,mai=c(1,1,0.5,0.5))

survplot(Surv(tmp.water.2group$OS.year, tmp.water.2group$OS)~BCM.bin,show.nrisk = T,data=tmp.water.2group,stitle='',subset = tmp.water.2group$OS.year <= 2.6,
         xlab=expression(bold('Time (years)')),ylab=expression(bold("OS")),main='BCM: Response to Anti-PD1 therapy',
         lwd=3,col=c("royalblue3","tomato3"),legend.pos = 'bottomleft',mark=20,cex.lab=1.5, cex.main=1.5,snames = c('BCM < Median','BCM > Median'))
pval = summary(coxph(Surv(tmp.Low$OS.year, tmp.Low$OS) ~ B.lineage, data= tmp.Low))$coeff[,'Pr(>|z|)']
pval=ifelse(formatC(format='f',pval, digits=3) < '0.001', '< 0.001',formatC(format='f',pval, digits=3))


survplot(Surv(tmp.water.CTLA4$OS.year, tmp.water.CTLA4$OS)~BCM.bin,show.nrisk = T,data=tmp.water.CTLA4,stitle='',
         xlab=expression(bold('Time (years)')),ylab=expression(bold("OS")),main='BCM: Response to Anti-CTLA4 therapy',
         lwd=3,col=c("royalblue3","tomato3"),legend.pos = 'bottomleft',mark=20,cex.lab=1.5, cex.main=1.5,snames = c('BCM < Median','BCM > Median'))
pval = summary(coxph(Surv(tmp.Low$OS.year, tmp.Low$OS) ~ B.lineage, data= tmp.Low))$coeff[,'Pr(>|z|)']
pval=ifelse(formatC(format='f',pval, digits=3) < '0.001', '< 0.001',formatC(format='f',pval, digits=3))


















