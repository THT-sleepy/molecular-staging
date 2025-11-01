##loading packages
library(dplyr)
library(survminer)
library(survival)

##loading files
load("RData/clin_mut.Rdata")

##先看看详细分期本身能不能分开
folder_path <- "~/molecular_staging/output"
filename_os <- paste0("all","_OS stage.pdf")
filename_os <- file.path(folder_path,filename_os)
filename_rfs <- paste0("all","_RFS stage.pdf")
filename_rfs <- file.path(folder_path,filename_rfs)
#os
pdf(filename_os,width = 8,height = 8,onefile = F)
km_os<-survfit(Surv(OS,OS_E)~Stage1,data = clin_mut)
p <- ggsurvplot(km_os, 
                conf.int=F, #是否显示生存率的95%CI
                risk.table=F,#显示风险表
                palette="jco", #柳叶刀配色
                title=paste0(" Result IASLC 9th Edition TNM Pathological Stage"), #大标题
                #risk.table.height = 0.15,#风险表的高度比例
                tables.height = 0.17,
                pval=T,  
                pval.method=T,
                surv.median.line="v",
                ylab="OS(%)",
                xlab="Months"
)
print(p)
dev.off()
#rfs
pdf(filename_rfs,width = 8,height = 8,onefile = F)
km_rfs<-survfit(Surv(RFS,RFS_E)~Stage1,data = clin_mut)
p <- ggsurvplot(km_rfs, 
                conf.int=F, #是否显示生存率的95%CI
                risk.table=F,#显示风险表
                palette="jco", #柳叶刀配色
                title=paste0(" Result IASLC 9th Edition TNM Pathological Stage"), #大标题
                #risk.table.height = 0.15,#风险表的高度比例
                tables.height = 0.17,
                pval=T,  
                pval.method=T,
                surv.median.line="v",
                ylab="RFS(%)",
                xlab="Months"
)
print(p)
dev.off()

##IA升IB
folder_path <- "~/molecular_staging/output"
filename_os <- paste0("IA-IB","_OS stage.pdf")
filename_os <- file.path(folder_path,filename_os)
filename_rfs <- paste0("IA-IB","_RFS stage.pdf")
filename_rfs <- file.path(folder_path,filename_rfs)
clin_mut1 <- clin_mut %>%
  filter(Stage1 %in% c("IA","IB")) %>%
  mutate(Category =ifelse(TP53_category=="TP53 Wt" & Stage1=="IA","IA TP53 Wt",
                          ifelse(TP53_category=="TP53 Mut" & Stage1=="IA","IA TP53 Mut","IB")))
pdf(filename_os,width = 8,height = 6,onefile = F)
km_os<-survfit(Surv(OS,OS_E)~Category,data = clin_mut1)
p <- ggsurvplot(km_os, 
                conf.int=F, #是否显示生存率的95%CI
                risk.table=T,#显示风险表
                palette="jco", #柳叶刀配色
                title=paste0(" Result IASLC 9th Edition TNM Pathological Stage"), #大标题
                risk.table.height = 0.25,#风险表的高度比例
                #tables.height = 0.17,
                pval=F,  
                pval.method=T,
                surv.median.line="v",
                ylab="OS(%)",
                xlab="Months"
)
print(p)
dev.off()

pdf(filename_rfs,width = 8,height = 6,onefile = F)
km_rfs<-survfit(Surv(RFS,RFS_E)~Category,data = clin_mut1)
p <- ggsurvplot(km_rfs, 
                conf.int=F, #是否显示生存率的95%CI
                risk.table=T,#显示风险表
                palette="jco", #柳叶刀配色
                title=paste0(" Result IASLC 9th Edition TNM Pathological Stage"), #大标题
                risk.table.height = 0.25,#风险表的高度比例
                #tables.height = 0.17,
                pval=F,  
                pval.method=T,
                surv.median.line="v",
                ylab="OS(%)",
                xlab="Months"
)
print(p)
dev.off()
















##后面的升期没有意义，因为我们的数据本来就分不开，所以就不分了
if(0){
  clin_mut$Stage3 <- ifelse(clin_mut$Stage1=="IA","IA",
                            ifelse(clin_mut$Stage1=="IB","IB",
                                   ifelse(clin_mut$Stage1 %in% c("IIA","IIB"),"II",
                                          ifelse(clin_mut$Stage1 %in% c("IIIA","IIIB","IVA","IVB"),"III/IV",NA))))
  km_os<-survfit(Surv(OS,OS_E)~Stage3,data = clin_mut)
  ggsurvplot(km_os, 
             conf.int=F, #是否显示生存率的95%CI
             risk.table=F,#显示风险表
             palette="jco", #柳叶刀配色
             title=paste0(" Result IASLC 9th Edition TNM Pathological Stage"), #大标题
             #risk.table.height = 0.15,#风险表的高度比例
             tables.height = 0.17,
             pval=T,  
             pval.method=T,
             surv.median.line="v",
             ylab="OS(%)",
             xlab="Months"
  )
  #IB升II
  clin_mut$Stage3 <- ifelse(clin_mut$Stage1=="IA","IA",
                            ifelse(clin_mut$Stage1=="IB","IB",
                                   ifelse(clin_mut$Stage1 %in% c("IIA","IIB"),"II",
                                          ifelse(clin_mut$Stage1 %in% c("IIIA","IIIB","IVA","IVB"),"III/IV",NA))))
  clin_mut1 <- clin_mut %>%
    filter(Stage3 %in% c("IB","II")) %>%
    mutate(Category =ifelse(TP53_category=="TP53 Wt" & Stage3=="IB","IB TP53 Wt",
                            ifelse(TP53_category=="TP53 Mut" & Stage3=="IB","IB TP53 Mut","II")))
  km_os<-survfit(Surv(OS,OS_E)~Category,data = clin_mut1)
  p <- ggsurvplot(km_os, 
                  conf.int=F, #是否显示生存率的95%CI
                  risk.table=T,#显示风险表
                  palette="jco", #柳叶刀配色
                  title=paste0(" Result IASLC 9th Edition TNM Pathological Stage"), #大标题
                  #risk.table.height = 0.15,#风险表的高度比例
                  #tables.height = 0.17,
                  pval=F,  
                  pval.method=T,
                  surv.median.line="v",
                  ylab="OS(%)",
                  xlab="Months"
  )
  print(p)
  #II升III
  clin_mut$Stage3 <- ifelse(clin_mut$Stage1=="IA","IA",
                            ifelse(clin_mut$Stage1=="IB","IB",
                                   ifelse(clin_mut$Stage1 %in% c("IIA","IIB"),"II",
                                          ifelse(clin_mut$Stage1 %in% c("IIIA","IIIB","IVA","IVB"),"III/IV",NA))))
  clin_mut1 <- clin_mut %>%
    filter(Stage3 %in% c("II","III/IV")) %>%
    mutate(Category =ifelse(TP53_category=="TP53 Wt" & Stage3=="II","II TP53 Wt",
                            ifelse(TP53_category=="TP53 Mut" & Stage3=="II","II TP53 Mut","III/IV")))
  km_os<-survfit(Surv(OS,OS_E)~Category,data = clin_mut1)
  p <- ggsurvplot(km_os, 
                  conf.int=F, #是否显示生存率的95%CI
                  risk.table=T,#显示风险表
                  palette="jco", #柳叶刀配色
                  title=paste0(" Result IASLC 9th Edition TNM Pathological Stage"), #大标题
                  #risk.table.height = 0.15,#风险表的高度比例
                  #tables.height = 0.17,
                  pval=F,  
                  pval.method=T,
                  surv.median.line="v",
                  ylab="OS(%)",
                  xlab="Months"
  )
  print(p)
  ##TK
  clin_mut$Stage3 <- ifelse(clin_mut$Stage1=="IA","IA",
                            ifelse(clin_mut$Stage1=="IB","IB",
                                   ifelse(clin_mut$Stage1 %in% c("IIA","IIB"),"II",
                                          ifelse(clin_mut$Stage1 %in% c("IIIA","IIIB","IVA","IVB"),"III/IV",NA))))
  km_os<-survfit(Surv(OS,OS_E)~Stage3,data = clin_mut)
  ggsurvplot(km_os, 
             conf.int=F, #是否显示生存率的95%CI
             risk.table=F,#显示风险表
             palette="jco", #柳叶刀配色
             title=paste0(" Result IASLC 9th Edition TNM Pathological Stage"), #大标题
             #risk.table.height = 0.15,#风险表的高度比例
             tables.height = 0.17,
             pval=T,  
             pval.method=T,
             surv.median.line="v",
             ylab="OS(%)",
             xlab="Months"
  )
  #IA升IB
  clin_mut1 <- clin_mut %>%
    filter(Stage3 %in% c("IA","IB")) %>%
    mutate(TK_category=ifelse(TP53_category=="TP53 Mut" | KRAS_category=="KRAS Canonical","TK Mut","TK Wt"))%>%
    mutate(Category =ifelse(TK_category=="TK Wt" & Stage3=="IA","IA TK Wt",
                            ifelse(TK_category=="TK Mut" & Stage3=="IA","IA TK Mut","IB")))
  km_os<-survfit(Surv(OS,OS_E)~Category,data = clin_mut1)
  p <- ggsurvplot(km_os, 
                  conf.int=F, #是否显示生存率的95%CI
                  risk.table=T,#显示风险表
                  palette="jco", #柳叶刀配色
                  title=paste0(" Result IASLC 9th Edition TNM Pathological Stage"), #大标题
                  #risk.table.height = 0.15,#风险表的高度比例
                  #tables.height = 0.17,
                  pval=F,  
                  pval.method=T,
                  surv.median.line="v",
                  ylab="OS(%)",
                  xlab="Months"
  )
  print(p)
  #IB升II
  clin_mut1 <- clin_mut %>%
    filter(Stage3 %in% c("IB","II")) %>%
    mutate(TK_category=ifelse(TP53_category=="TP53 Mut" | KRAS_category=="KRAS Canonical","TK Mut","TK Wt"))%>%
    mutate(Category =ifelse(TK_category=="TK Wt" & Stage3=="IB","IB TK Wt",
                            ifelse(TK_category=="TK Mut" & Stage3=="IB","IB TK Mut","II")))
  km_os<-survfit(Surv(OS,OS_E)~Category,data = clin_mut1)
  p <- ggsurvplot(km_os, 
                  conf.int=F, #是否显示生存率的95%CI
                  risk.table=T,#显示风险表
                  palette="jco", #柳叶刀配色
                  title=paste0(" Result IASLC 9th Edition TNM Pathological Stage"), #大标题
                  #risk.table.height = 0.15,#风险表的高度比例
                  #tables.height = 0.17,
                  pval=F,  
                  pval.method=T,
                  surv.median.line="v",
                  ylab="OS(%)",
                  xlab="Months"
  )
  print(p)
  #II升III/IV
  clin_mut1 <- clin_mut %>%
    filter(Stage3 %in% c("II","III/IV")) %>%
    mutate(TK_category=ifelse(TP53_category=="TP53 Mut" | KRAS_category=="KRAS Canonical","TK Mut","TK Wt"))%>%
    mutate(Category =ifelse(TK_category=="TK Wt" & Stage3=="II","II TK Wt",
                            ifelse(TK_category=="TK Mut" & Stage3=="II","II TK Mut","III/IV")))
  km_os<-survfit(Surv(OS,OS_E)~Category,data = clin_mut1)
  p <- ggsurvplot(km_os, 
                  conf.int=F, #是否显示生存率的95%CI
                  risk.table=T,#显示风险表
                  palette="jco", #柳叶刀配色
                  title=paste0(" Result IASLC 9th Edition TNM Pathological Stage"), #大标题
                  #risk.table.height = 0.15,#风险表的高度比例
                  #tables.height = 0.17,
                  pval=F,  
                  pval.method=T,
                  surv.median.line="v",
                  ylab="OS(%)",
                  xlab="Months"
  )
  print(p)
  
  
  
  
  
   
  
  
#IB升IIA
clin_mut1 <- clin_mut %>%
  filter(Stage1 %in% c("IB","IIA")) %>%
  mutate(Category =ifelse(TP53_category=="TP53 Wt" & Stage1=="IB","IB TP53 Wt",
                          ifelse(TP53_category=="TP53 Mut" & Stage1=="IB","IB TP53 Mut","IIA")))
km_os<-survfit(Surv(OS,OS_E)~Category,data = clin_mut1)
p <- ggsurvplot(km_os, 
                conf.int=F, #是否显示生存率的95%CI
                risk.table=T,#显示风险表
                palette="jco", #柳叶刀配色
                title=paste0(" Result IASLC 9th Edition TNM Pathological Stage"), #大标题
                #risk.table.height = 0.15,#风险表的高度比例
                tables.height = 0.17,
                pval=F,  
                pval.method=T,
                surv.median.line="v",
                ylab="OS(%)",
                xlab="Months"
)
print(p)
clin_mut1 <- clin_mut %>%
  filter(Stage1 %in% c("IA","IB")) %>%
  mutate(Category =ifelse(TP53_category=="TP53 Wt" & Stage1=="IA","IA TP53 Wt",
                          ifelse(TP53_category=="TP53 Mut" & Stage1=="IA","IA TP53 Mut",
                                 ifelse(TP53_category=="TP53 Wt" & Stage1=="IB","IB TP53 Wt",
                                        ifelse(TP53_category=="TP53 Mut" & Stage1=="IB","IB TP53 Mut","others")))))
km_os<-survfit(Surv(OS,OS_E)~Category,data = clin_mut1)
p <- ggsurvplot(km_os, 
                conf.int=F, #是否显示生存率的95%CI
                risk.table=F,#显示风险表
                palette="jco", #柳叶刀配色
                title=paste0(" Result IASLC 9th Edition TNM Pathological Stage"), #大标题
                #risk.table.height = 0.15,#风险表的高度比例
                tables.height = 0.17,
                pval=F,  
                pval.method=T,
                surv.median.line="v",
                ylab="OS(%)",
                xlab="Months"
)
print(p)
#IIA升IIB
clin_mut1 <- clin_mut %>%
  filter(Stage1 %in% c("IIA","IIB")) %>%
  mutate(Category =ifelse(TP53_category=="TP53 Wt" & Stage1=="IIA","IIA TP53 Wt",
                          ifelse(TP53_category=="TP53 Mut" & Stage1=="IIA","IIA TP53 Mut","IIB")))
km_os<-survfit(Surv(OS,OS_E)~Category,data = clin_mut1)
p <- ggsurvplot(km_os, 
                conf.int=F, #是否显示生存率的95%CI
                risk.table=T,#显示风险表
                palette="jco", #柳叶刀配色
                title=paste0(" Result IASLC 9th Edition TNM Pathological Stage"), #大标题
                #risk.table.height = 0.15,#风险表的高度比例
                tables.height = 0.17,
                pval=F,  
                pval.method=T,
                surv.median.line="v",
                ylab="OS(%)",
                xlab="Months"
)
print(p)

##IIB升IIIA
clin_mut1 <- clin_mut %>%
  filter(Stage1 %in% c("IIB","IIIA")) %>%
  mutate(Category =ifelse(TP53_category=="TP53 Wt" & Stage1=="IIB","IIB TP53 Wt",
                          ifelse(TP53_category=="TP53 Mut" & Stage1=="IIB","IIB TP53 Mut","IIIA")))
km_os<-survfit(Surv(OS,OS_E)~Category,data = clin_mut1)
p <- ggsurvplot(km_os, 
                conf.int=F, #是否显示生存率的95%CI
                risk.table=T,#显示风险表
                palette="jco", #柳叶刀配色
                title=paste0(" Result IASLC 9th Edition TNM Pathological Stage"), #大标题
                #risk.table.height = 0.15,#风险表的高度比例
                tables.height = 0.17,
                pval=F,  
                pval.method=T,
                surv.median.line="v",
                ylab="OS(%)",
                xlab="Months"
)
print(p)
}

