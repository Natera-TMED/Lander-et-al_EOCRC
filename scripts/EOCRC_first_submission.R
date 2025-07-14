source('/Users/cfungtammasan/gitlab/tmp_code/scripts/functions_ML.R')
library(data.table)
require(vautils)
library(ggplot2)
library(ggExtra)
library(patchwork)
library(plotly)
library(tidyverse)
library(gtsummary)
library(officer)
library(flextable)
library(Polychrome)
library(ggsignif)

###############
## load data ##
###############

## load all patient age data for making distribution
dt.master=fread('deid.dt.master.csv')
dt.master=dt.master[Pass_legal_check==TRUE]

## load selected cohort patient-level data
dt.pat.ori=fread('deid.dt.pat.csv')
dt.pat.ori[AgeGroup=='Early',AgeGroup:='EOCRC']
dt.pat.ori[AgeGroup=='Average',AgeGroup:='AOCRC']
dt.pat.ori$AgeGroup=factor(dt.pat.ori$AgeGroup,levels=c('EOCRC','AOCRC'))
dt.pat.ori=dt.pat.ori[Pass_legal_check==TRUE]

## load selected cohort genomics data
maf.pat=fread('deid.maf.pat.csv')
maf.pat=maf.pat[Pass_legal_check==TRUE]

###############################
## making tables and figures ##
###############################

## Table 1
dt.pat=dt.pat.ori

dt.pat[!(SEX %in% c('MALE','FEMALE')),SEX:='Unknown']
dt.pat.table.source=dt.pat[,.(Age,SEX,Cancer.Type,AgeGroup,Stage,MSI_TMB,MSI.score,TMB)]

reset_gtsummary_theme()
table <- dt.pat.table.source |>
  tbl_summary(
    by=AgeGroup,
    type = Age ~ "continuous2",
    label = list(Cancer.Type ~ "Cancer type",SEX ~ "Sex")
  ) 
table
# table(dt.pat[,.(AgeGroup,MSI_TMB)])
# prop.table(table(dt.pat[,.(AgeGroup,MSI_TMB)]),margin=1)
table_flextable=as_flex_table(table)

doc <- read_docx() %>%
  body_add_flextable(table_flextable) 

print(doc,target='gtsummary_table.docx')


## Figure 1
## 1a
dt.master[Age<=50,AgeGroup:="EOCRC"]
dt.master[Age>50 & Age<60,AgeGroup:="Mid"]
dt.master[Age>=60 ,AgeGroup:="AOCRC"]
dt.master$AgeGroup=factor(dt.master$AgeGroup,levels=c('EOCRC','Mid','AOCRC'))

ggplot(dt.master,aes(Age,fill=AgeGroup))+geom_histogram(bins=200,size=0.1,color='gray')+
  theme_bw()+
  labs(y='Number of participant',fill= 'Age Group')+
  scale_fill_manual(values=c(AOCRC='#56B4E9',Mid='gray',EOCRC='#E69F00'))

dt.master[,.(N=.N,percentage=round(100*length(Age)/dim(dt.master)[1],1),median=median(Age)),by=AgeGroup]


## 1b
dt.pat=dt.pat.ori
table(dt.pat[,.(Stage,Cancer.Stage)])
table(dt.pat$Stage)

subset_dt.pat=dt.pat[!is.na(Stage),.(Stage,AgeGroup)]
subset_dt.pat_table=table(subset_dt.pat)
subset_dt.pat=subset_dt.pat[,.N,by=.(Stage,AgeGroup)]
subset_dt.pat[,sum_N:=sum(N),by=.(AgeGroup)]
subset_dt.pat[,percentage:=100*N/sum_N]
subset_dt.pat
ggplot(subset_dt.pat,aes(x=AgeGroup,y=percentage,fill=Stage))+
  geom_bar(position="fill", stat="identity")+
  scale_y_continuous(labels=scales::percent_format(scale=100),expand=expansion(add=c(0,0.1)),breaks=c(0,.25,.50,.75,1))+
  geom_text(aes(label=paste0(N," (",round(percentage,1),"%)")),position=position_fill(vjust=0.5))+
  labs(y="Percentage",x="Age Group", fill= 'Stage')+
  scale_fill_brewer(palette="Set2")+
  theme_bw()+
  geom_signif(comparisons = list(c('EOCRC','AOCRC')),
                annotations=paste("p =", format(fisher.test(subset_dt.pat_table,workspace=2e8)$p.value,scientific=T,digit=2),'***'),
                y_position=-0.85,
                tip_length=0.001,
                map_signif_level=T)+
  NULL

## 1c
dt.pat=dt.pat.ori
table(dt.pat$Cancer.Type)
table(dt.pat$AgeGroup)

subset_dt.pat=dt.pat[!is.na(Cancer.Type),.(Cancer.Type,AgeGroup,Stage.Group)]
subset_dt.pat=subset_dt.pat[,.N,by=.(Cancer.Type,AgeGroup,Stage.Group)]
subset_dt.pat[,sum_N:=sum(N),by=.(AgeGroup,Stage.Group)]
subset_dt.pat[,percentage:=100*N/sum_N]
subset_dt.pat[,.(N),by=.(Stage.Group,Cancer.Type,AgeGroup)]
p_values <- subset_dt.pat[, {
  contingency_table <- dcast(.SD, Cancer.Type ~ AgeGroup, value.var = "N", fill = 0)
  print(contingency_table)
  test_result <- fisher.test(contingency_table[, -1],workspace=2e8) 
  list(p_value = test_result$p.value)
}, by = Stage.Group]
subset_dt.pat <- merge(subset_dt.pat, p_values, by = "Stage.Group", all.x = TRUE)
subset_dt.pat[, significance_label := fifelse(p_value < 0.001, "***",
                                              fifelse(p_value < 0.01, "**",
                                                      fifelse(p_value < 0.05, "*", "ns")))]
subset_dt.pat
facet_label=unique(subset_dt.pat[,.(p_value,significance_label),by=Stage.Group])
facet_label[,Sig:=paste("p =", format(p_value,scientific=T,digit=2),significance_label)]
facet_label[,start:='EOCRC']
facet_label[,end:='AOCRC']
facet_label[,position:=1.05]

ggplot(subset_dt.pat,aes(x=AgeGroup,y=percentage,fill=Cancer.Type))+
  geom_bar(position="fill", stat="identity")+
  scale_y_continuous(labels=scales::percent_format(scale=100),expand=expansion(add=c(0,0.1)),breaks=c(0,.25,.50,.75,1))+
  scale_fill_manual(values=c(Colon='forestgreen',Rectal='chartreuse3'))+
  geom_text(aes(label=paste0(N,"\n(",round(percentage,1),"%)")),position=position_fill(vjust=0.5))+
  theme_bw()+
  labs(y="Percentage",x="Cancer Type", fill= 'Age Group')+
  geom_signif(data=facet_label,
              # comparisons = list(c('Colon','Rectal')),
              # aes(annotations),
              aes(xmin = start, xmax = end, annotations = Sig, y_position = position),
              tip_length=0.0005,
              inherit.aes=FALSE,
              manual = TRUE)+
  facet_grid(cols=vars(Stage.Group))+
  NULL


## 1d, S1a
dt.pat=dt.pat.ori[Stage %in% c('I','II','III')]
# dt.pat=dt.pat.ori[Stage %in% c('IV')]

scatter=ggplot(dt.pat[!(is.na(TMB)) & !(is.na(MSI.score)) ],aes(x=TMB,y=MSI.score,color=AgeGroup))+
  geom_point(alpha=0.5,size=1)+
  scale_y_log10()+
  scale_x_log10(breaks = c(0.1,1,10,100,1000),labels= c(0.1,1,10,100,1000))+
  scale_color_manual(values=c(AOCRC='#56B4E9',EOCRC='#E69F00'))+
  labs(y="MSI Score",x="TMB", color= 'Age Group')+
  theme_bw()+
  theme(legend.position=c(1.1,1.07))+
  geom_vline(xintercept=20,linetype="dotted",color="red",linewidth=1) +
  geom_hline(yintercept=6,linetype="dotted",color="red",linewidth=1) +
  NULL
scatter
ggExtra::ggMarginal(scatter,type="histogram",fill="lightblue",margins="both")

## 1e, S1b
dt.pat=dt.pat.ori[Stage %in% c('I','II','III')]
# dt.pat=dt.pat.ori[Stage %in% c('IV')]
subset_dt.pat=dt.pat[!is.na(MSI_TMB) ,.(MSI_TMB,AgeGroup)]
subset_dt.pat_table=table(subset_dt.pat)
subset_dt.pat=subset_dt.pat[,.N,by=.(MSI_TMB,AgeGroup)]
subset_dt.pat[,sum_N:=sum(N),by=.(AgeGroup)]
subset_dt.pat[,percentage:=100*N/sum_N]
fill_order=subset_dt.pat[order(percentage),unique(MSI_TMB)]
subset_dt.pat$MSI_TMB=factor(subset_dt.pat$MSI_TMB,levels=fill_order)

ggplot(subset_dt.pat,aes(x=AgeGroup,y=percentage,fill=MSI_TMB))+
  geom_bar(position="fill", stat="identity")+
  scale_y_continuous(labels=scales::percent_format(scale=100),expand=expansion(add=c(0,0.1)),breaks=c(0,.25,.50,.75,1))+
  # scale_fill_manual(values=c(AOCRC='#00BFC4',EOCRC='#F8766D'))+
  geom_text(aes(label=paste0(N," (",round(percentage,1),"%)")),position=position_fill(vjust=0.5))+
  labs(y="Percentage",x="Age Group", fill= 'MSI and TMB status')+
  scale_fill_manual(values=c(`MSS/TMB-Low`="#FFCD00",`MSI/TMB-High`= "steelblue4" , `MSS/TMB-High`= "#C4D600"))+
  # guides(fill=guide_legend(reverse=T))+
  geom_signif(comparisons = list(c('EOCRC','AOCRC')),
              # annotations=paste("p =", format(chisq.test(subset_dt.pat_table)$p.value,scientific=T,digit=2),'***'),
              annotations=paste("p =", format(fisher.test(subset_dt.pat_table)$p.value,scientific=T,workspace=2e8,digit=2),'**'),
              # y_position=-3.05,
              y_position=-3.45,
              tip_length=0.0005,
              map_signif_level=T)+
  theme_bw()+
  NULL


## S2
dt.pat=dt.pat.ori
subset_dt.pat=dt.pat[!is.na(MSI_TMB) & !is.na(Cancer.Type),.(MSI_TMB,AgeGroup,Cancer.Type,Stage.Group)]
subset_dt.pat=subset_dt.pat[,.N,by=.(MSI_TMB,AgeGroup,Cancer.Type,Stage.Group)]
subset_dt.pat[,sum_N:=sum(N),by=.(AgeGroup,Cancer.Type,Stage.Group)]
subset_dt.pat[,percentage:=100*N/sum_N]
fill_order=subset_dt.pat[order(percentage),unique(MSI_TMB)]
subset_dt.pat$MSI_TMB=factor(subset_dt.pat$MSI_TMB,levels=fill_order)
facet_label <- subset_dt.pat[, {
  contingency_table <- dcast(.SD, MSI_TMB ~ AgeGroup, value.var = "N", fill = 0)
  print(contingency_table)
  test_result <- fisher.test(contingency_table[, -1],workspace=2e8) 
  list(p_value = test_result$p.value)
}, by = .(Stage.Group,Cancer.Type)]
facet_label[, significance_label := fifelse(p_value < 0.001, "***",
                                              fifelse(p_value < 0.01, "**",
                                                      fifelse(p_value < 0.05, "*", "ns")))]

facet_label[,Sig:=paste("p =", format(p_value,scientific=T,digit=2),significance_label)]
facet_label[,start:='EOCRC']
facet_label[,end:='AOCRC']
facet_label[,position:=1.05]

ggplot(subset_dt.pat,aes(x=AgeGroup,y=percentage,fill=MSI_TMB))+
  geom_bar(position="fill", stat="identity")+
  scale_y_continuous(labels=scales::percent_format(scale=100),expand=expansion(add=c(0,0.1)),breaks=c(0,.25,.50,.75,1))+
  geom_text_repel(aes(label=ifelse(percentage>0,paste0(N," (",round(percentage,1),"%)"),NA)),position=position_fill(vjust=0.5))+
  labs(y="Percentage",x="Age Group", fill= 'MSI and TMB status')+
  scale_fill_manual(values=c(`MSS/TMB-Low`="#FFCD00",`MSI/TMB-High`= "steelblue4" , `MSS/TMB-High`= "#C4D600"))+
  theme_bw()+
  facet_grid(rows=vars(Cancer.Type),cols=vars(Stage.Group))+
  geom_signif(data=facet_label,
              # comparisons = list(c('Colon','Rectal')),
              # aes(annotations),
              aes(xmin = start, xmax = end, annotations = Sig, y_position = position),
              tip_length=0.0002,
              inherit.aes=FALSE,
              manual = TRUE)+
  NULL

## 2, S3
# dt.pat=dt.pat.ori[Stage %in% c('I','II','III')]
dt.pat=dt.pat.ori[Stage %in% c('IV')]
dim(dt.pat)
subset_maf.pat=maf.pat
subset_maf.pat=subset_maf.pat[Bundling.ID %in% dt.pat$Bundling.ID]
list_A <- dt.pat[vcf.intersect==TRUE & AgeGroup=="AOCRC" , Tumor_Sample_Barcode]
list_B <- dt.pat[vcf.intersect==TRUE & AgeGroup=="EOCRC" , Tumor_Sample_Barcode]
subset_maf.pat[Bundling.ID %in% list_A,AgeGroup:='AOCRC']
subset_maf.pat[Bundling.ID %in% list_B,AgeGroup:='EOCRC']
subset_maf.pat$AgeGroup=factor(subset_maf.pat$AgeGroup,levels=c('EOCRC','AOCRC'))
subset_maf.pat=subset_maf.pat[!is.na(Hugo_Symbol) & FLAGS==FALSE & !is.na(AgeGroup) & IMPACT %in% c('HIGH','MODERATE') ,.(Bundling.ID,AgeGroup,Hugo_Symbol)]
subset_maf.pat=unique(subset_maf.pat)
subset_maf.pat[,N:=.N,by=.(Hugo_Symbol)]
subset_maf.pat[,participant_N:=length(unique(Bundling.ID))]
subset_maf.pat[,percentage:=100*N/participant_N]
subset_maf.pat=subset_maf.pat[order(-percentage)]
subset_maf.pat[,N_by_group:=.N,by=.(Hugo_Symbol,AgeGroup)]
subset_maf.pat[,participant_N_by_group:=length(unique(Bundling.ID)),by=AgeGroup]
subset_maf.pat[,percentage_by_group:=100*N_by_group/participant_N_by_group]
subset_maf.pat$Bundling.ID=NULL
subset_maf.pat=unique(subset_maf.pat)
fill_order=subset_maf.pat[order(-percentage),unique(Hugo_Symbol )]
subset_maf.pat$Hugo_Symbol =factor(subset_maf.pat$Hugo_Symbol ,levels=fill_order)
dim(subset_maf.pat)
ggplot(subset_maf.pat[1:50,],aes(x=Hugo_Symbol,y=percentage_by_group,fill=AgeGroup))+
  geom_bar(position='dodge', stat="identity")+
  labs(x='HUGO Symbol of gene',y='Percentage')+
  scale_fill_manual(values=c(AOCRC='#56B4E9',EOCRC='#E69F00'))+
  geom_text(data=subset_maf.pat[1:50,][AgeGroup=='AOCRC',],aes(label=round(percentage_by_group,0)),vjust=-0.5,hjust=0)+
  geom_text(data=subset_maf.pat[1:50,][AgeGroup=='EOCRC',],aes(label=round(percentage_by_group,0)),vjust=-0.5,hjust=1)+
  expand_limits(y=max(subset_maf.pat$percentage_by_group*1.05))+
  NULL
  

## 3, S4

dt.pat=dt.pat.ori[Stage %in% c('I','II','III')]
# dt.pat=dt.pat.ori[Stage %in% c('IV')]
list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-Low"  & AgeGroup=="AOCRC" , Tumor_Sample_Barcode]
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-Low"  & AgeGroup=="EOCRC" , Tumor_Sample_Barcode]
name_B <- "AOCRC_MSS/TMB-Lo"
name_A <- "EOCRC_MSS/TMB-Lo"
compare_genomic_groups_genes(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,min_gene_Freq=1,Plot_n_genes=15,Drivers=T)
# maf_forest(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B)
list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSI/TMB-High"  & AgeGroup=="AOCRC" , Tumor_Sample_Barcode]
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSI/TMB-High"  & AgeGroup=="EOCRC" , Tumor_Sample_Barcode]
name_B <- "AOCRC_MSI/TMB-Hi"
name_A <- "EOCRC_MSI/TMB-Hi"
compare_genomic_groups_genes(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,min_gene_Freq=1,Plot_n_genes=15,Drivers=T)
list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-High"  & AgeGroup=="AOCRC" , Tumor_Sample_Barcode]
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-High"  & AgeGroup=="EOCRC" , Tumor_Sample_Barcode]
name_B <- "AOCRC_MSS/TMB-Hi"
name_A <- "EOCRC_MSS/TMB-Hi"
compare_genomic_groups_genes(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,min_gene_Freq=1,Plot_n_genes=15,Drivers=T)


## 4, S5

dt.pat=dt.pat.ori[Stage %in% c('I','II','III')]
# dt.pat=dt.pat.ori[Stage %in% c('IV')]
list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-Low"  & AgeGroup=="AOCRC" , Tumor_Sample_Barcode]
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-Low"  & AgeGroup=="EOCRC" , Tumor_Sample_Barcode]
name_B <- "AOCRC_MSS/TMB-Lo"
name_A <- "EOCRC_MSS/TMB-Lo"
compare_genomic_groups_variants(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,Plot_n_variants=15,Drivers=T)
# maf_forest(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B)
list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSI/TMB-High"  & AgeGroup=="AOCRC" , Tumor_Sample_Barcode]
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSI/TMB-High"  & AgeGroup=="EOCRC" , Tumor_Sample_Barcode]
name_B <- "AOCRC_MSI/TMB-Hi"
name_A <- "EOCRC_MSI/TMB-Hi"
compare_genomic_groups_variants(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,Plot_n_variants=15,Drivers=T)
list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-High"  & AgeGroup=="AOCRC" , Tumor_Sample_Barcode]
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-High"  & AgeGroup=="EOCRC" , Tumor_Sample_Barcode]
name_B <- "AOCRC_MSS/TMB-Hi"
name_A <- "EOCRC_MSS/TMB-Hi"
compare_genomic_groups_variants(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,Plot_n_variants=15,Drivers=T)



## S6a
dt.pat=dt.pat.ori
subset_dt.pat=dt.pat[!is.na(SEX) & SEX %in% c('MALE','FEMALE') ,.(SEX,AgeGroup,Stage.Group)]
subset_dt.pat=subset_dt.pat[,.N,by=.(SEX,AgeGroup,Stage.Group)]
subset_dt.pat[,sum_N:=sum(N),by=.(SEX,Stage.Group)]
subset_dt.pat[,percentage:=100*N/sum_N]
subset_dt.pat
facet_label <- subset_dt.pat[, {
  contingency_table <- dcast(.SD, SEX ~ AgeGroup, value.var = "N", fill = 0)
  print(contingency_table)
  test_result <- fisher.test(contingency_table[, -1],workspace=2e8) 
  list(p_value = test_result$p.value)
}, by = .(Stage.Group)]
facet_label[, significance_label := fifelse(p_value < 0.001, "***",
                                            fifelse(p_value < 0.01, "**",
                                                    fifelse(p_value < 0.05, "*", "ns")))]

facet_label[,Sig:=paste("p =", format(p_value,scientific=T,digit=2),significance_label)]
facet_label[,start:='FEMALE']
facet_label[,end:='MALE']
facet_label[,position:=1.05]
ggplot(subset_dt.pat,aes(x=SEX,y=percentage,fill=AgeGroup))+
  geom_bar(position="fill", stat="identity")+
  scale_y_continuous(labels=scales::percent_format(scale=100),expand=expansion(add=c(0,0.1)),breaks=c(0,.25,.50,.75,1))+
  scale_fill_manual(values=c(AOCRC='#56B4E9',EOCRC='#E69F00'))+
  geom_text(aes(label=paste0(N,"\n(",round(percentage,1),"%)")),position=position_fill(vjust=0.5))+
  theme_bw()+
  facet_grid(cols=vars(Stage.Group))+
  labs(y="Percentage",x="Sex", fill= 'Age Group')+
  geom_signif(data=facet_label,
              aes(xmin = start, xmax = end, annotations = Sig, y_position = position),
              tip_length=0.0002,
              inherit.aes=FALSE,
              manual = TRUE)+
  NULL

## S6b
dt.pat=dt.pat.ori

subset_dt.pat=dt.pat[!is.na(SEX) & SEX %in% c('MALE','FEMALE') & !is.na(AgeGroup) & !is.na(MSI_TMB) ,.(SEX,AgeGroup,MSI_TMB,Stage.Group)]
subset_dt.pat=subset_dt.pat[,.N,by=.(SEX,AgeGroup,MSI_TMB,Stage.Group)]
subset_dt.pat[,sum_N:=sum(N),by=.(AgeGroup,MSI_TMB,Stage.Group)]
subset_dt.pat[,percentage:=100*N/sum_N]
subset_dt.pat$MSI_TMB=factor(subset_dt.pat$MSI_TMB,levels=c('MSS/TMB-Low','MSI/TMB-High','MSS/TMB-High'))
subset_dt.pat
facet_label <- subset_dt.pat[, {
  contingency_table <- dcast(.SD, SEX ~ AgeGroup, value.var = "N", fill = 0)
  print(contingency_table)
  test_result <- fisher.test(contingency_table[, -1],workspace=2e8) 
  list(p_value = test_result$p.value)
}, by = .(Stage.Group,MSI_TMB)]
facet_label[, significance_label := fifelse(p_value < 0.001, "***",
                                            fifelse(p_value < 0.01, "**",
                                                    fifelse(p_value < 0.05, "*", "ns")))]

facet_label[,Sig:=paste("p =", format(p_value,scientific=T,digit=2),significance_label)]
facet_label[,start:='EOCRC']
facet_label[,end:='AOCRC']
facet_label[,position:=1.05]
ggplot(subset_dt.pat,aes(x=AgeGroup,y=percentage,fill=SEX))+
  geom_bar(position="fill", stat="identity")+
  scale_y_continuous(labels=scales::percent_format(scale=100),expand=expansion(add=c(0,0.1)),breaks=c(0,.25,.50,.75,1))+
  geom_text(aes(label=paste0(N,"\n(",round(percentage,1),"%)")),position=position_fill(vjust=0.5))+
  labs(y="Percentage",x="Age Group", fill= 'Sex')+
  scale_fill_manual(values=c(MALE='#009E73',FEMALE='#F0E442'))+
  # scale_fill_manual(values=c('TRUE'='red','FALSE'='green'))+
  geom_signif(data=facet_label,
              aes(xmin = start, xmax = end, annotations = Sig, y_position = position),
              tip_length=0.0002,
              inherit.aes=FALSE,
              manual = TRUE)+
  theme_bw()+
  facet_grid(cols=vars(MSI_TMB),rows=vars(Stage.Group))


## S7-S8
dt.pat=dt.pat.ori[Stage %in% c('I','II','III')]
# dt.pat=dt.pat.ori[Stage %in% c('IV')]
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-Low"  & SEX=="MALE"   & AgeGroup=="AOCRC", Tumor_Sample_Barcode]
list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-Low"  & SEX=="FEMALE" & AgeGroup=="AOCRC", Tumor_Sample_Barcode]
name_A <- "Male_AOCRC_MSS/TMB-Lo"
name_B <- "Female_AOCRC_MSS/TMB-Lo"
compare_genomic_groups_genes(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,min_gene_Freq=1,Plot_n_genes=15,Drivers=T)
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSI/TMB-High" & SEX=="MALE"   & AgeGroup=="AOCRC", Tumor_Sample_Barcode]
list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSI/TMB-High" & SEX=="FEMALE" & AgeGroup=="AOCRC", Tumor_Sample_Barcode]
name_A <- "Male_AOCRC_MSI/TMB-Hi"
name_B <- "Female_AOCRC_MSI/TMB-Hi"
compare_genomic_groups_genes(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,min_gene_Freq=1,Plot_n_genes=15,Drivers=T)
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-High"  & SEX=="MALE"  & AgeGroup=="AOCRC", Tumor_Sample_Barcode]
list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-High"  & SEX=="FEMALE" & AgeGroup=="AOCRC", Tumor_Sample_Barcode]
name_A <- "Male_AOCRC_MSS/TMB-Hi"
name_B <- "Female_AOCRC_MSS/TMB-Hi"
compare_genomic_groups_genes(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,min_gene_Freq=1,Plot_n_genes=15,Drivers=T)

list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-Low"  & SEX=="MALE"   & AgeGroup=="EOCRC", Tumor_Sample_Barcode]
list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-Low"  & SEX=="FEMALE" & AgeGroup=="EOCRC", Tumor_Sample_Barcode]
name_A <- "Male_EOCRC_MSS/TMB-Lo"
name_B <- "Female_EOCRC_MSS/TMB-Lo"
compare_genomic_groups_genes(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,min_gene_Freq=1,Plot_n_genes=15,Drivers=T)
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSI/TMB-High" & SEX=="MALE"   & AgeGroup=="EOCRC", Tumor_Sample_Barcode]
list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSI/TMB-High" & SEX=="FEMALE" & AgeGroup=="EOCRC", Tumor_Sample_Barcode]
name_A <- "Male_EOCRC_MSI/TMB-Hi"
name_B <- "Female_EOCRC_MSI/TMB-Hi"
compare_genomic_groups_genes(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,min_gene_Freq=1,Plot_n_genes=15,Drivers=T)
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-High"  & SEX=="MALE"  & AgeGroup=="EOCRC", Tumor_Sample_Barcode]
list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-High"  & SEX=="FEMALE" & AgeGroup=="EOCRC", Tumor_Sample_Barcode]
name_A <- "Male_EOCRC_MSS/TMB-Hi"
name_B <- "Female_EOCRC_MSS/TMB-Hi"
compare_genomic_groups_genes(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,min_gene_Freq=1,Plot_n_genes=15,Drivers=T)


## S9-S10
dt.pat=dt.pat.ori[Stage %in% c('I','II','III')]
# dt.pat=dt.pat.ori[Stage %in% c('IV')]
list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-Low" & SEX=="MALE"  & AgeGroup=="AOCRC" , Tumor_Sample_Barcode]
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-Low" & SEX=="MALE"  & AgeGroup=="EOCRC" , Tumor_Sample_Barcode]
name_B <- "Male_AOCRC_MSS/TMB-Lo"
name_A <- "Male_EOCRC_MSS/TMB-Lo"
compare_genomic_groups_genes(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,min_gene_Freq=1,Plot_n_genes=15,Drivers=T)
list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSI/TMB-High" & SEX=="MALE"  & AgeGroup=="AOCRC" , Tumor_Sample_Barcode]
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSI/TMB-High" & SEX=="MALE"  & AgeGroup=="EOCRC" , Tumor_Sample_Barcode]
name_B <- "Male_AOCRC_MSI/TMB-Hi"
name_A <- "Male_EOCRC_MSI/TMB-Hi"
compare_genomic_groups_genes(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,min_gene_Freq=1,Plot_n_genes=15,Drivers=T)
list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-High" & SEX=="MALE"  & AgeGroup=="AOCRC" , Tumor_Sample_Barcode]
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-High" & SEX=="MALE"  & AgeGroup=="EOCRC" , Tumor_Sample_Barcode]
name_B <- "Male_AOCRC_MSS/TMB-Hi"
name_A <- "Male_EOCRC_MSS/TMB-Hi"
compare_genomic_groups_genes(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,min_gene_Freq=1,Plot_n_genes=15,Drivers=T)

list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-Low"  & SEX=="FEMALE"& AgeGroup=="AOCRC" , Tumor_Sample_Barcode]
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-Low"  & SEX=="FEMALE"& AgeGroup=="EOCRC" , Tumor_Sample_Barcode]
name_B <- "Female_AOCRC_MSS/TMB-Lo"
name_A <- "Female_EOCRC_MSS/TMB-Lo"
compare_genomic_groups_genes(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,min_gene_Freq=1,Plot_n_genes=15,Drivers=T)
list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSI/TMB-High" & SEX=="FEMALE" & AgeGroup=="AOCRC" , Tumor_Sample_Barcode]
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSI/TMB-High" & SEX=="FEMALE" & AgeGroup=="EOCRC" , Tumor_Sample_Barcode]
name_B <- "Female_AOCRC_MSI/TMB-Hi"
name_A <- "Female_EOCRC_MSI/TMB-Hi"
compare_genomic_groups_genes(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,min_gene_Freq=1,Plot_n_genes=15,Drivers=T)
list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-High" & SEX=="FEMALE" & AgeGroup=="AOCRC" , Tumor_Sample_Barcode]
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-High" & SEX=="FEMALE" & AgeGroup=="EOCRC" , Tumor_Sample_Barcode]
name_B <- "Female_AOCRC_MSS/TMB-Hi"
name_A <- "Female_EOCRC_MSS/TMB-Hi"
compare_genomic_groups_genes(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,min_gene_Freq=1,Plot_n_genes=15,Drivers=T)


color_signatures <- c(`SIG-clock-like`="#FFCD00", 
                      `SIG-MMR deficiency`= "steelblue4" , 
                      `SIG-Thiopurine chemotherapy`= "#E87722", 
                      `SIG-POLE exo domain mutation`= "#C4D600", 
                      `SIG-Damage by ROS`= "forestgreen", 
                      `SIG-Unknown`= "gray",`SIG-Other`= "black", 
                      `SIG-APOBEC activity`= "maroon", 
                      `SIG-Colibactin exposure`="pink",
                      `SIG-Defective POLD1 proofreading`="#E31A1C",
                      `SIG-MMR deficiency + temozolomide`="deeppink",
                      `SIG-MMR deficiency + POLD1 mutation`="orchid1",     
                      `SIG-Aristolochic acid exposure`="khaki2",          
                      `SIG-Aflatoxin exposure`="darkturquoise",                 
                      `SIG-Unknown chemotherapy`="#FB9A99",             
                      `SIG-HR deficiency`="palegreen2",            
                      `SIG-BER deficiency`="dodgerblue2",                   
                      `SIG-Platinum chemotherapy`="gold1",               
                      `SIG-Azathioprine exposure`="#CAB2D6",               
                      `SIG-Haloalkanes exposure`="darkorange4",                
                      `SIG-Possible_sequencing_artefacts`="#6A3D9A",       
                      `SIG-UV light exposure`="white",                   
                      `SIG-AID activity`="#911eb4",                        
                      `SIG-Polymerase eta somatic hypermutation`="#E31A1C")



## S12
# run_maftools_signatures(maf.pat, clinicalData=NULL, name=NULL) # generated mut_mat.Rdat, mut_mat.csv, cos_matrix.Rdat
load('uncorrected_cos_matrix.Rdat')

dt.pat=dt.pat.ori
patient_label=dt.pat[,.(Tumor_Sample_Barcode=Bundling.ID,AgeGroup,MSI_TMB,Stage,TMB)]

merged_cos.matrix=patient_label[cos.matrix,on='Tumor_Sample_Barcode']
merged_cos.matrix=merged_cos.matrix[Stage %in% c('I','II','III')]
# merged_cos.matrix=merged_cos.matrix[Stage %in% c('IV')]

subset_cos.matrix=merged_cos.matrix[!is.na(MSI_TMB) & !is.na(SBS_group) & SBS_group!='NA',.(SBS_group,AgeGroup,MSI_TMB)]
subset_cos.matrix=subset_cos.matrix[,.N,by=.(SBS_group,AgeGroup,MSI_TMB)]
subset_cos.matrix[,sum_N:=sum(N),by=.(AgeGroup,MSI_TMB)]
subset_cos.matrix[,percentage:=100*N/sum_N]
subset_cos.matrix$MSI_TMB=factor(subset_cos.matrix$MSI_TMB,levels=c('MSS/TMB-Low','MSI/TMB-High','MSS/TMB-High'))
fill_order=subset_cos.matrix[order(percentage),unique(SBS_group)]
subset_cos.matrix$SBS_group=factor(subset_cos.matrix$SBS_group,levels=fill_order)
facet_label <- subset_cos.matrix[, {
  contingency_table <- dcast(.SD, SBS_group ~ AgeGroup, value.var = "N", fill = 0)
  print(contingency_table)
  test_result <- chisq.test(contingency_table[, -1]) 
  list(p_value = test_result$p.value)
}, by = .(MSI_TMB)]
facet_label[, significance_label := fifelse(p_value < 0.001, "***",
                                            fifelse(p_value < 0.01, "**",
                                                    fifelse(p_value < 0.05, "*", "ns")))]


facet_label[,Sig:=paste("p =", format(p_value,scientific=T,digit=2),significance_label)]
facet_label[,start:='EOCRC']
facet_label[,end:='AOCRC']
facet_label[,position:=1.05]

ggplot(subset_cos.matrix,aes(x=AgeGroup,y=percentage,fill=SBS_group))+
  geom_bar(position="fill", stat="identity")+
  scale_y_continuous(labels=scales::percent_format(scale=100),expand=expansion(add=c(0,0.1)),breaks=c(0,.25,.50,.75,1))+
  geom_text(aes(label=ifelse(percentage>1,paste0(N," (",round(percentage,1),"%)"),NA)),position=position_fill(vjust=0.5))+
  theme_bw()+
  scale_fill_manual(values=color_signatures)+
  labs(y="Percentage",x="Age Group", fill= 'SBS Group')+
  guides(fill=guide_legend(reverse=T))+
  facet_grid(cols=vars(MSI_TMB))+
  geom_signif(data=facet_label,
              aes(xmin = start, xmax = end, annotations = Sig, y_position = position),
              tip_length=0.0002,
              inherit.aes=FALSE,
              manual = TRUE)+
  NULL


## Table S1
post_correction_profile=fread('post_correction_profile.csv')
cos.matrix=run_maftools_signatures2(post_correction_profile, clinicalData=NULL, name=NULL)

cos.matrix[,count_ind:=.N,by=SBS_group]
name_signature=unique(cos.matrix[,.(SBS_group,COSMIC.SBS.max,count_ind)])
name_signature[,COSMIC.SBS.max_group_number:=as.integer(gsub('[^0-9]','',COSMIC.SBS.max))]
table=name_signature[order(-count_ind,COSMIC.SBS.max_group_number,COSMIC.SBS.max)][,toString(COSMIC.SBS.max),by=SBS_group]
table_flextable=flextable(table)

doc <- read_docx() %>%
  body_add_flextable(table_flextable) 

print(doc,target='tableS2.docx')


cos.matrix$new_max=apply(cos.matrix[,3:80],1,function(row){
  sorted_index=order(row,decreasing=T)
  colnames(cos.matrix)[3:80][sorted_index[1]]
})
cos.matrix$new_second_max=apply(cos.matrix[,3:80],1,function(row){
  sorted_index=order(row,decreasing=T)
  colnames(cos.matrix)[3:80][sorted_index[2]]
})
mapping_SBS=unique(cos.matrix[,.(COSMIC.SBS.max,SBS_group)])
names(mapping_SBS)[1]='new_second_max'
names(mapping_SBS)[2]='SBS_group2'
mapping_SBS=mapping_SBS[!is.na(`SBS_group2`)]
cos.matrix2=mapping_SBS[cos.matrix,on='new_second_max']
write.csv(cos.matrix2,'corrected_cos_matrix2.csv',row.names=F,quote=F)

## 5,S11


# dt.pat=dt.pat.ori[Stage %in% c('I','II','III')]
dt.pat=dt.pat.ori[Stage %in% c('IV')]

patient_label=dt.pat[,.(Tumor_Sample_Barcode=Bundling.ID,AgeGroup,MSI_TMB,Stage)]
merged_cos.matrix=patient_label[cos.matrix,on='Tumor_Sample_Barcode']


subset_cos.matrix=merged_cos.matrix[!is.na(MSI_TMB) & !is.na(SBS_group) & SBS_group!='NA',.(SBS_group,AgeGroup,MSI_TMB)]
subset_cos.matrix=subset_cos.matrix[,.N,by=.(SBS_group,AgeGroup,MSI_TMB)]
subset_cos.matrix[,sum_N:=sum(N),by=.(AgeGroup,MSI_TMB)]
subset_cos.matrix[,percentage:=100*N/sum_N]
subset_cos.matrix$MSI_TMB=factor(subset_cos.matrix$MSI_TMB,levels=c('MSS/TMB-Low','MSI/TMB-High','MSS/TMB-High'))
fill_order=subset_cos.matrix[order(percentage),unique(SBS_group)]
subset_cos.matrix$SBS_group=factor(subset_cos.matrix$SBS_group,levels=fill_order)
subset_cos.matrix
facet_label <- subset_cos.matrix[, {
  contingency_table <- dcast(.SD, SBS_group ~ AgeGroup, value.var = "N", fill = 0)
  print(contingency_table)
  test_result <- chisq.test(contingency_table[, -1]) 
  list(p_value = test_result$p.value)
}, by = .(MSI_TMB)]
facet_label[, significance_label := fifelse(p_value < 0.001, "***",
                                            fifelse(p_value < 0.01, "**",
                                                    fifelse(p_value < 0.05, "*", "ns")))]


facet_label[,Sig:=paste("p =", format(p_value,scientific=T,digit=2),significance_label)]
facet_label[,start:='EOCRC']
facet_label[,end:='AOCRC']
facet_label[,position:=1.05]
ggplot(subset_cos.matrix,aes(x=AgeGroup,y=percentage,fill=SBS_group))+
  geom_bar(position="fill", stat="identity")+
  scale_y_continuous(labels=scales::percent_format(scale=100),expand=expansion(add=c(0,0.1)),breaks=c(0,.25,.50,.75,1))+
  geom_text(aes(label=ifelse(percentage>2,paste0(N," (",round(percentage,1),"%)"),NA)),position=position_fill(vjust=0.5))+
  theme_bw()+
  scale_fill_manual(values=color_signatures)+
  labs(y="Percentage",x="Age Group", fill= 'SBS Group')+
  guides(fill=guide_legend(reverse=T))+
  facet_grid(cols=vars(MSI_TMB))+
  geom_signif(data=facet_label,
              aes(xmin = start, xmax = end, annotations = Sig, y_position = position),
              tip_length=0.0002,
              inherit.aes=FALSE,
              manual = TRUE)+
  NULL


## S13

dt.pat=dt.pat.ori
patient_label=dt.pat[,.(Tumor_Sample_Barcode=Bundling.ID,AgeGroup,MSI_TMB,Stage)]
merged_cos.matrix=patient_label[cos.matrix2,on='Tumor_Sample_Barcode']

# merged_cos.matrix=merged_cos.matrix[Stage %in% c('I','II','III')]
merged_cos.matrix=merged_cos.matrix[Stage %in% c('IV')]

subset_cos.matrix=merged_cos.matrix[!is.na(MSI_TMB) & !is.na(SBS_group) & SBS_group!='NA' & !is.na(SBS_group2) & SBS_group2!='NA',.(SBS_group,SBS_group2,AgeGroup,MSI_TMB)]
subset_cos.matrix=subset_cos.matrix[,.N,by=.(SBS_group,SBS_group2,AgeGroup,MSI_TMB)]
subset_cos.matrix[,sum_N:=sum(N),by=.(AgeGroup,MSI_TMB)]
subset_cos.matrix[,percentage:=100*N/sum_N]
subset_cos.matrix$MSI_TMB=factor(subset_cos.matrix$MSI_TMB,levels=c('MSS/TMB-Low','MSI/TMB-High','MSS/TMB-High'))
valid_combination=subset_cos.matrix[,.(total_pct=sum(percentage)),by=.(SBS_group,SBS_group2)]
valid_combination=valid_combination[total_pct>=1]
subset_cos.matrix=subset_cos.matrix[valid_combination,on=.(SBS_group,SBS_group2)]
fill_order=unique(subset_cos.matrix[order(-percentage),SBS_group])
fill_order=rev(fill_order)
subset_cos.matrix$SBS_group=factor(subset_cos.matrix$SBS_group,levels=fill_order)
subset_cos.matrix$SBS_group2=factor(subset_cos.matrix$SBS_group2,levels=fill_order)
subset_cos.matrix$SBS_group=droplevels(subset_cos.matrix$SBS_group)
subset_cos.matrix$SBS_group2=droplevels(subset_cos.matrix$SBS_group2)
subset_cos.matrix
subset_cos.matrix=na.omit(subset_cos.matrix)
write.csv(subset_cos.matrix,'forsam.csv',row.names=F)
ggplot(subset_cos.matrix,aes(x=SBS_group2,y=SBS_group,fill=percentage))+
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", percentage)), size = 3) +
  scale_fill_gradient(low = "white", high = "blue", name = "Frequency (%)") +
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  theme_bw()+
  labs(x="Second highest probability SBS group", y= 'Highest probability SBS group')+
  facet_grid(cols=vars(MSI_TMB),rows=vars(AgeGroup))+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 10),
    legend.text = element_text(size = 5)
  )
  



## S14
merged_cos.matrix=patient_label[cos.matrix2,on='Tumor_Sample_Barcode']
merged_cos.matrix=merged_cos.matrix[Stage %in% c('I','II','III')]
subset_cos.matrix=merged_cos.matrix[!is.na(MSI_TMB) & !is.na(COSMIC.SBS.max ) & COSMIC.SBS.max !='NA',.(COSMIC.SBS.max ,AgeGroup,MSI_TMB)]
subset_cos.matrix=subset_cos.matrix[,.N,by=.(COSMIC.SBS.max ,AgeGroup,MSI_TMB)]
subset_cos.matrix[,sum_N:=sum(N),by=.(AgeGroup,MSI_TMB)]
subset_cos.matrix[,percentage:=100*N/sum_N]
subset_cos.matrix$MSI_TMB=factor(subset_cos.matrix$MSI_TMB,levels=c('MSS/TMB-Low','MSI/TMB-High','MSS/TMB-High'))
fill_order=subset_cos.matrix[order(percentage),unique(COSMIC.SBS.max )]
subset_cos.matrix$COSMIC.SBS.max =factor(subset_cos.matrix$COSMIC.SBS.max ,levels=fill_order)
subset_cos.matrix
facet_label <- subset_cos.matrix[, {
  contingency_table <- dcast(.SD, COSMIC.SBS.max ~ AgeGroup, value.var = "N", fill = 0)
  print(contingency_table)
  test_result <- chisq.test(contingency_table[, -1]) 
  list(p_value = test_result$p.value)
}, by = .(MSI_TMB)]
facet_label[, significance_label := fifelse(p_value < 0.001, "***",
                                            fifelse(p_value < 0.01, "**",
                                                    fifelse(p_value < 0.05, "*", "ns")))]


facet_label[,Sig:=paste("p =", format(p_value,scientific=T,digit=2),significance_label)]
facet_label[,start:='EOCRC']
facet_label[,end:='AOCRC']
facet_label[,position:=1.05]
ggplot(subset_cos.matrix,aes(x=AgeGroup,y=percentage,fill=COSMIC.SBS.max ))+
  geom_bar(position="fill", stat="identity")+
  scale_y_continuous(labels=scales::percent_format(scale=100),expand=expansion(add=c(0,0.1)),breaks=c(0,.25,.50,.75,1))+
  geom_text(aes(label=ifelse(percentage>2,paste0(COSMIC.SBS.max,' ',N," (",round(percentage,1),"%)"),NA)),position=position_fill(vjust=0.5))+
  theme_bw()+
  # scale_fill_manual(values=P100)+
  scale_fill_manual(values=P54_color)+
  labs(y="Percentage",x="Age Group", fill= 'SBS')+
  guides(fill=guide_legend(reverse=T))+
  facet_grid(cols=vars(MSI_TMB))+
  geom_signif(data=facet_label,
              aes(xmin = start, xmax = end, annotations = Sig, y_position = position),
              tip_length=0.0002,
              inherit.aes=FALSE,
              manual = TRUE)+
  NULL



## S15-S17

cos.matrix2_subset=cos.matrix2[,.(Tumor_Sample_Barcode,SBS_group,COSMIC.SBS.max)]
dt.pat=dt.pat.ori[Stage %in% c('I','II','III')]
# dt.pat=dt.pat.ori[Stage %in% c('IV')]
# dt.pat=dt.pat.ori[Stage %in% c('I','II','III','IV')]

dt.pat=cos.matrix2_subset[dt.pat,on='Tumor_Sample_Barcode']
list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-Low"  & SBS_group=="SIG-clock-like" & AgeGroup=="AOCRC" , Tumor_Sample_Barcode]
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-Low"  & SBS_group=="SIG-clock-like" & AgeGroup=="EOCRC" , Tumor_Sample_Barcode]
name_B <- "AOCRC_MSS/TMB-Lo"
name_A <- "EOCRC_MSS/TMB-Lo"
compare_genomic_groups_genes(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,min_gene_Freq=1,Plot_n_genes=15,Drivers=T)
# maf_forest(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B)
list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSI/TMB-High"  & SBS_group=="SIG-MMR deficiency" & AgeGroup=="AOCRC" , Tumor_Sample_Barcode]
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSI/TMB-High"  & SBS_group=="SIG-MMR deficiency" & AgeGroup=="EOCRC" , Tumor_Sample_Barcode]
name_B <- "AOCRC_MSI/TMB-Hi"
name_A <- "EOCRC_MSI/TMB-Hi"
compare_genomic_groups_genes(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,min_gene_Freq=1,Plot_n_genes=15,Drivers=T)
list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-High"  & SBS_group=="SIG-POLE exo domain mutation" & AgeGroup=="AOCRC" , Tumor_Sample_Barcode]
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-High"  & SBS_group=="SIG-POLE exo domain mutation" & AgeGroup=="EOCRC" , Tumor_Sample_Barcode]
name_B <- "AOCRC_MSS/TMB-Hi"
name_A <- "EOCRC_MSS/TMB-Hi"
compare_genomic_groups_genes(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,min_gene_Freq=1,Plot_n_genes=15,Drivers=T)



## S18
dt.pat=dt.pat.ori[Stage %in% c('I','II','III','IV')]
subset_dt.pat=dt.pat[hasMRD_new==T & !is.na(MSI_TMB),.(posMRD_new=as.character(posMRD_new),AgeGroup,Stage,MSI_TMB)]
length(subset_dt.pat[,posMRD_new])
table(subset_dt.pat[,posMRD_new])
prop.table(table(subset_dt.pat[,posMRD_new]))
subset_dt.pat=subset_dt.pat[,.N,by=.(posMRD_new,AgeGroup,Stage,MSI_TMB)]
subset_dt.pat[,sum_N:=sum(N),by=.(AgeGroup,Stage,MSI_TMB)]
subset_dt.pat[,percentage:=100*N/sum_N]
subset_dt.pat$MSI_TMB=factor(subset_dt.pat$MSI_TMB,levels=c('MSS/TMB-Low','MSI/TMB-High','MSS/TMB-High'))
subset_dt.pat[,posMRD_new:=ifelse(posMRD_new=='TRUE','Positive','Negative')]
subset_dt.pat
facet_label <- subset_dt.pat[, {
  contingency_table <- dcast(.SD, posMRD_new ~ AgeGroup, value.var = "N", fill = 0)
  print(contingency_table)
  test_result <- fisher.test(contingency_table[, -1],workspace=2e8) 
  list(p_value = test_result$p.value)
}, by = .(Stage,MSI_TMB)]
facet_label[, significance_label := fifelse(p_value < 0.001, "***",
                                            fifelse(p_value < 0.01, "**",
                                                    fifelse(p_value < 0.05, "*", "ns")))]

facet_label[,Sig:=paste("p =", format(p_value,scientific=T,digit=2),significance_label)]
facet_label[,start:='EOCRC']
facet_label[,end:='AOCRC']
facet_label[,position:=1.05]
ggplot(subset_dt.pat,aes(x=AgeGroup,y=percentage,fill=posMRD_new))+
  geom_bar(position="fill", stat="identity")+
  scale_y_continuous(labels=scales::percent_format(scale=100),expand=expansion(add=c(0,0.15)),breaks=c(0,.25,.50,.75,1))+
  scale_fill_manual(values=c('Positive'='coral','Negative'='green'))+

    geom_text(aes(label=paste0(N," (",round(percentage,1),"%)")),position=position_fill(vjust=0.5))+
  theme_bw()+
  labs(y="Percentage",x="Age Group", fill= 'MRD positivity')+
  facet_grid(cols=vars(Stage),rows=vars(MSI_TMB))+
  geom_signif(data=facet_label,
              aes(xmin = start, xmax = end, annotations = Sig, y_position = position),
              tip_length=0.0002,
              inherit.aes=FALSE,
              manual = TRUE)+
  NULL


  
## S19
dt.pat=dt.pat.ori[Stage %in% c('I','II','III')]
dt.pat=dt.pat.ori[Stage %in% c('IV')]

#"MSS/TMB-Low" 
list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-Low"  & AgeGroup=="EOCRC" & posMRD_new==FALSE, Tumor_Sample_Barcode]
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-Low"  & AgeGroup=="EOCRC" & posMRD_new==TRUE, Tumor_Sample_Barcode]
name_B <- "EOCRC_MSS/TMB-Lo_negMRD"
name_A <- "EOCRC_MSS/TMB-Lo_posMRD"
compare_genomic_groups_genes(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,min_gene_Freq=1,Plot_n_genes=15,Drivers=T)

list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-Low"  & AgeGroup=="AOCRC" & posMRD_new==FALSE, Tumor_Sample_Barcode]
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-Low"  & AgeGroup=="AOCRC" & posMRD_new==TRUE, Tumor_Sample_Barcode]
name_B <- "AOCRC_MSS/TMB-Lo_negMRD"
name_A <- "AOCRC_MSS/TMB-Lo_posMRD"
compare_genomic_groups_genes(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,min_gene_Freq=1,Plot_n_genes=15,Drivers=T)

"MSI/TMB-High"  
list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSI/TMB-High"  & AgeGroup=="EOCRC" & posMRD_new==FALSE, Tumor_Sample_Barcode]
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSI/TMB-High"  & AgeGroup=="EOCRC" & posMRD_new==TRUE, Tumor_Sample_Barcode]
name_B <- "EOCRC_MSI/TMB-Hi_negMRD"
name_A <- "EOCRC_MSI/TMB-Hi_posMRD"
compare_genomic_groups_genes(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,min_gene_Freq=1,Plot_n_genes=15,Drivers=T)

list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSI/TMB-High"  & AgeGroup=="AOCRC" & posMRD_new==FALSE, Tumor_Sample_Barcode]
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSI/TMB-High"  & AgeGroup=="AOCRC" & posMRD_new==TRUE, Tumor_Sample_Barcode]
name_B <- "AOCRC_MSI/TMB-Hi_negMRD"
name_A <- "AOCRC_MSI/TMB-Hi_posMRD"
compare_genomic_groups_genes(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,min_gene_Freq=1,Plot_n_genes=15,Drivers=T)

"MSS/TMB-High"  
list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-High"  & AgeGroup=="EOCRC" & posMRD_new==FALSE, Tumor_Sample_Barcode]
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-High"  & AgeGroup=="EOCRC" & posMRD_new==TRUE, Tumor_Sample_Barcode]
name_B <- "EOCRC_MSS/TMB-Hi_negMRD"
name_A <- "EOCRC_MSS/TMB-Hi_posMRD"
compare_genomic_groups_genes(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,min_gene_Freq=1,Plot_n_genes=15,Drivers=T)

list_B <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-High"  & AgeGroup=="AOCRC" & posMRD_new==FALSE, Tumor_Sample_Barcode]
list_A <- dt.pat[vcf.intersect==TRUE & MSI_TMB=="MSS/TMB-High"  & AgeGroup=="AOCRC" & posMRD_new==TRUE, Tumor_Sample_Barcode]
name_B <- "AOCRC_MSS/TMB-Hi_negMRD"
name_A <- "AOCRC_MSS/TMB-Hi_posMRD"
compare_genomic_groups_genes(maf.pat, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.1,min_gene_Freq=1,Plot_n_genes=15,Drivers=T)

