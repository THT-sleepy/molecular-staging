##loading packages
library(dplyr)
library(survminer)
library(survival)
library(MASS) # 用于contr.sdif
library(tibble)

##loading files
load("RData/clin_mut.Rdata")

##TP53
gene <- "TP53"
#I
if(1){
#I-KM
clin_mut1 <- clin_mut %>%
  filter(Stage2=="I") %>%
  mutate(Category=ifelse(TP53_category=="TP53 Mut","TP53 MUT","TP53 WT")
         )
folder_path <- "~/molecular_staging/output"
filename_os <- paste0(gene,"_OS stage I.pdf")
filename_os <- file.path(folder_path,filename_os)
filename_rfs <- paste0(gene,"_RFS stage I.pdf")
filename_rfs <- file.path(folder_path,filename_rfs)
#os
pdf(filename_os,width = 8,height = 8,onefile = F)
km_os<-survfit(Surv(OS,OS_E)~Category,data = clin_mut1)
summary_median <- summary(km_os)


# 从 summary_median 对象中提取我们需要的信息
# summary_median$table 是一个包含详细统计量的矩阵
median_data <- as.data.frame(summary_median$table) %>%
  rownames_to_column("Group") %>% # 将行名（分组名）转为列
  dplyr::select(Group, median, `0.95LCL`, `0.95UCL`) %>% # 选择需要的列
  rename(
    `Median Survival (days)` = median,
    `Median 95% CI Lower` = `0.95LCL`,
    `Median 95% CI Upper` = `0.95UCL`
  )

# 定义5年对应的天数
five_year_days <- 5 * 12

# 使用 summary 函数在指定时间点（5年）提取生存率和95%CI
# 这是最推荐的方法，因为它会进行线性插值
summary_5yr <- summary(km_os, times = five_year_days)

# 从 summary 结果中构建最终的整洁表格
# summary_5yr$surv, lower, upper 是按分组顺序排列的向量
result_table <- tibble(
  Group = c("TP53 Mut","TP53 Wt"),
  `5-year Survival (%)` = round(summary_5yr$surv * 100, 1),
  `95% CI` = paste0(
    round(summary_5yr$lower * 100, 1),
    " - ",
    round(summary_5yr$upper * 100, 1)
  )
)
info_table <- cbind(median_data,result_table[,c(2,3)])
info_table[is.na(info_table)] <- "."
info_table <- info_table %>%
  mutate(`Median in Months`= paste0(`Median Survival (days)`,"(",`Median 95% CI Lower`,"-",`Median 95% CI Upper`,")")) %>%
  mutate(`5-Year Estimate`= paste0(`5-year Survival (%)`,"(",`95% CI`,")"))%>%
  dplyr::select(Group,`Median in Months`,`5-Year Estimate`)

p <- ggsurvplot(km_os, 
                conf.int=F, #是否显示生存率的95%CI
                risk.table=T,#显示风险表
                palette="jco", #柳叶刀配色
                title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage I"), #大标题
                #risk.table.height = 0.15,#风险表的高度比例
                tables.height = 0.17,
                pval=T,  
                pval.method=T,
                surv.median.line="v",
                cumevents=T,
                cumcensor = T,
                ylab="OS(%)",
                xlab="Months"
)
print(p)
dev.off()

#rfs
pdf(filename_rfs,width = 8,height = 8,onefile = F)
km_rfs<-survfit(Surv(RFS,RFS_E)~Category,data = clin_mut1)
p <- ggsurvplot(km_rfs, 
                conf.int=F, #是否显示生存率的95%CI
                risk.table=T,#显示风险表
                palette="jco", #柳叶刀配色
                title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage I"), #大标题
                #risk.table.height = 0.15,#风险表的高度比例
                tables.height = 0.17,
                pval=T,  
                pval.method=T,
                surv.median.line="v",
                cumevents=T,
                cumcensor = T,
                ylab="RFS(%)",
                xlab="Months"
)
print(p)
dev.off()
#I-Cox
#多因素
clin_mut1$Differentiation <- factor(clin_mut1$Differentiation, levels = c("Well", "Mod", "Poor"), ordered = TRUE)
contrasts(clin_mut1$Differentiation) <- contr.sdif(3)
#clin_mut1$Age_E <- factor(clin_mut1$Age_E,levels=c("<40","40-49","50-59","60-69",">=70"))


clin_mut1$Category <- factor(clin_mut1$Category,levels = c("TP53 Wt","TP53 Mut"))
clin_mut1$Smoking_E <- factor(clin_mut1$Smoking_E,levels=c("Never","Ever"))

multicox <- coxph(Surv(OS,OS_E)~Category+Sex+Age_E+Smoking_E+Differentiation, data = clin_mut1)
multisum <- summary(multicox)




#提取所有基因的多因素COX回归分析结果至multiresult对象中
factor <-  c("TP53 Mut(vs Wt)","Male(vs Female)","Age>=60 (vs <60)","Ever Smoking(vs Never)","Differentiation Mod(vs Well)","Differentiation Poor(vs Mod)")
HR <- multisum$coefficients[,2]
L95CI <- multisum$conf.int[,3]
H95CI <- multisum$conf.int[,4]
pvalue <- multisum$coefficients[,5]
proportion <- clin_mut1 %>%
  summarise(
    Mut_n = sum(Category == "TP53 Mut", na.rm = TRUE),
    Male_n = sum(Sex == "Male", na.rm = TRUE),
    Age_n = sum(Age_E  == ">=60",na.rm = TRUE),
    Smoking_n =sum(Smoking_E== "Ever",na.rm = TRUE),
    Differentiation_mod_n = sum(Differentiation=="Mod",na.rm = TRUE),
    Differentiation_poor_n = sum(Differentiation=="Poor",na.rm = TRUE),
    Total_N = n(),
    Proportion_mut = paste0(Mut_n, "/", Total_N, " (", round(Mut_n / Total_N * 100, 1), "%)"),
    Proportion_male = paste0(Male_n, "/", Total_N, " (", round(Male_n / Total_N * 100, 1), "%)"),
    Proportion_age = paste0(Age_n, "/", Total_N, " (", round(Age_n / Total_N * 100, 1), "%)"),
    Proportion_smoking = paste0(Smoking_n, "/", Total_N, " (", round(Smoking_n / Total_N * 100, 1), "%)"),
    Proportion_mod = paste0(Differentiation_mod_n, "/", Total_N, " (", round(Differentiation_mod_n / Total_N * 100, 1), "%)"),
    Proportion_poor = paste0(Differentiation_poor_n, "/", Total_N, " (", round(Differentiation_poor_n / Total_N * 100, 1), "%)")
  ) 
multiresult <- data.frame(factor=factor,
                          proportion=as.vector(t(proportion[1,c(8:13)])),
                          HR=HR,
                          L95CI=L95CI,
                          H95CI=H95CI,
                          pvalue=pvalue)
multiresult$HR <- round(multiresult$HR,digits = 2)
multiresult$L95CI <- round(multiresult$L95CI,digits = 2)
multiresult$H95CI <- round(multiresult$H95CI,digits = 2)
multiresult$pvalue  <-  ifelse(
  multiresult$pvalue < 0.0001,        # 条件：p值是否小于0.0001
  "<0.0001",              # 如果是，则返回此字符串
  sprintf("%.4f",multiresult$pvalu) # 如果否，则格式化为保留四位小数的字符串
)
 
multiresult$`HR(95% CI)` <- paste0(multiresult$HR,"(",multiresult$L95CI,"-",multiresult$H95CI,")")
multiresult <- multiresult[,c(1,2,7,6)]
print(multiresult)
}
#II
if(1){
#II-KM
clin_mut1 <- clin_mut %>%
  filter(Stage2=="II") %>%
  mutate(Category=TP53_category)
folder_path <- "~/molecular_staging/output"
filename_os <- paste0(gene,"_OS stage II.pdf")
filename_os <- file.path(folder_path,filename_os)
filename_rfs <- paste0(gene,"_RFS stage II.pdf")
filename_rfs <- file.path(folder_path,filename_rfs)
#os
pdf(filename_os,width = 8,height = 8,onefile = F)
km_os<-survfit(Surv(OS,OS_E)~Category,data = clin_mut1)
summary_median <- summary(km_os)

# 从 summary_median 对象中提取我们需要的信息
# summary_median$table 是一个包含详细统计量的矩阵
median_data <- as.data.frame(summary_median$table) %>%
  rownames_to_column("Group") %>% # 将行名（分组名）转为列
  dplyr::select(Group, median, `0.95LCL`, `0.95UCL`) %>% # 选择需要的列
  rename(
    `Median Survival (days)` = median,
    `Median 95% CI Lower` = `0.95LCL`,
    `Median 95% CI Upper` = `0.95UCL`
  )

# 定义5年对应的天数
five_year_days <- 5 * 12

# 使用 summary 函数在指定时间点（5年）提取生存率和95%CI
# 这是最推荐的方法，因为它会进行线性插值
summary_5yr <- summary(km_os, times = five_year_days)

# 从 summary 结果中构建最终的整洁表格
# summary_5yr$surv, lower, upper 是按分组顺序排列的向量
result_table <- tibble(
  Group = c("TP53 Mut","TP53 Wt"),
  `Median Survival (days)` = median_surv,
  `5-year Survival (%)` = round(summary_5yr$surv * 100, 1),
  `95% CI` = paste0(
    round(summary_5yr$lower * 100, 1),
    " - ",
    round(summary_5yr$upper * 100, 1)
  )
)
info_table <- cbind(median_data,result_table[,c(2,3)])
info_table[is.na(info_table)] <- "."
info_table <- info_table %>%
  mutate(`Median in Months`= paste0(`Median Survival (days)`,"(",`Median 95% CI Lower`,"-",`Median 95% CI Upper`,")")) %>%
  mutate(`5-Year Estimate`= paste0(`5-year Survival (%)`,"(",`95% CI`,")"))%>%
  dplyr::select(Group,`Median in Months`,`5-Year Estimate`)
p <- ggsurvplot(km_os, 
                conf.int=F, #是否显示生存率的95%CI
                risk.table=T,#显示风险表
                palette="jco", #柳叶刀配色
                title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage II"), #大标题
                #risk.table.height = 0.15,#风险表的高度比例
                tables.height = 0.17,
                pval=T,  
                pval.method=T,
                surv.median.line="v",
                cumevents=T,
                cumcensor = T,
                ylab="OS(%)",
                xlab="Months"
)
print(p)
dev.off()

#rfs
pdf(filename_rfs,width = 8,height = 8,onefile = F)
km_rfs<-survfit(Surv(RFS,RFS_E)~Category,data = clin_mut1)
p <- ggsurvplot(km_rfs, 
                conf.int=F, #是否显示生存率的95%CI
                risk.table=T,#显示风险表
                palette="jco", #柳叶刀配色
                title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage II"), #大标题
                #risk.table.height = 0.15,#风险表的高度比例
                tables.height = 0.17,
                pval=T,  
                pval.method=T,
                surv.median.line="v",
                cumevents=T,
                cumcensor = T,
                ylab="RFS(%)",
                xlab="Months"
)
print(p)
dev.off()
#II-Cox
#多因素
clin_mut1$Differentiation <- factor(clin_mut1$Differentiation, levels = c("Well", "Mod", "Poor"), ordered = TRUE)
contrasts(clin_mut1$Differentiation) <- contr.sdif(3)
#clin_mut1$Age_E <- factor(clin_mut1$Age_E,levels=c("<40","40-49","50-59","60-69",">=70"))


clin_mut1$Category <- factor(clin_mut1$Category,levels = c("TP53 Wt","TP53 Mut"))
clin_mut1$Smoking_E <- factor(clin_mut1$Smoking_E,levels=c("Never","Ever"))

multicox <- coxph(Surv(OS,OS_E)~Category+Sex+Age_E+Smoking_E+Differentiation, data = clin_mut1)
multisum <- summary(multicox)




#提取所有基因的多因素COX回归分析结果至multiresult对象中
factor <-  c("TP53 Mut(vs Wt)","Male(vs Female)","Age>=60 (vs <60)","Ever Smoking(vs Never)","Differentiation Mod(vs Well)","Differentiation Poor(vs Mod)")
HR <- multisum$coefficients[,2]
L95CI <- multisum$conf.int[,3]
H95CI <- multisum$conf.int[,4]
pvalue <- multisum$coefficients[,5]
proportion <- clin_mut1 %>%
  summarise(
    Mut_n = sum(Category == "TP53 Mut", na.rm = TRUE),
    Male_n = sum(Sex == "Male", na.rm = TRUE),
    Age_n = sum(Age_E  == ">=60",na.rm = TRUE),
    Smoking_n =sum(Smoking_E== "Ever",na.rm = TRUE),
    Differentiation_mod_n = sum(Differentiation=="Mod",na.rm = TRUE),
    Differentiation_poor_n = sum(Differentiation=="Poor",na.rm = TRUE),
    Total_N = n(),
    Proportion_mut = paste0(Mut_n, "/", Total_N, " (", round(Mut_n / Total_N * 100, 1), "%)"),
    Proportion_male = paste0(Male_n, "/", Total_N, " (", round(Male_n / Total_N * 100, 1), "%)"),
    Proportion_age = paste0(Age_n, "/", Total_N, " (", round(Age_n / Total_N * 100, 1), "%)"),
    Proportion_smoking = paste0(Smoking_n, "/", Total_N, " (", round(Smoking_n / Total_N * 100, 1), "%)"),
    Proportion_mod = paste0(Differentiation_mod_n, "/", Total_N, " (", round(Differentiation_mod_n / Total_N * 100, 1), "%)"),
    Proportion_poor = paste0(Differentiation_poor_n, "/", Total_N, " (", round(Differentiation_poor_n / Total_N * 100, 1), "%)")
  ) 
multiresult <- data.frame(factor=factor,
                          proportion=as.vector(t(proportion[1,c(8:13)])),
                          HR=HR,
                          L95CI=L95CI,
                          H95CI=H95CI,
                          pvalue=pvalue)
multiresult$HR <- round(multiresult$HR,digits = 2)
multiresult$L95CI <- round(multiresult$L95CI,digits = 2)
multiresult$H95CI <- round(multiresult$H95CI,digits = 2)
multiresult$pvalue  <-  ifelse(
  multiresult$pvalue < 0.0001,        # 条件：p值是否小于0.0001
  "<0.0001",              # 如果是，则返回此字符串
  sprintf("%.4f",multiresult$pvalu) # 如果否，则格式化为保留四位小数的字符串
)

multiresult$`HR(95% CI)` <- paste0(multiresult$HR,"(",multiresult$L95CI,"-",multiresult$H95CI,")")
multiresult <- multiresult[,c(1,2,7,6)]
print(multiresult)
}
#III&IV
if(1){
  #III-KM
  clin_mut1 <- clin_mut %>%
    filter(Stage2 %in% c("III","IV")) %>%
    mutate(Category=TP53_category)
  folder_path <- "~/molecular_staging/output"
  filename_os <- paste0(gene,"_OS stage III&IV.pdf")
  filename_os <- file.path(folder_path,filename_os)
  filename_rfs <- paste0(gene,"_RFS stage III&IV.pdf")
  filename_rfs <- file.path(folder_path,filename_rfs)
  #os
  pdf(filename_os,width = 8,height = 8,onefile = F)
  km_os<-survfit(Surv(OS,OS_E)~Category,data = clin_mut1)
  summary_median <- summary(km_os)
  
  # 从 summary_median 对象中提取我们需要的信息
  # summary_median$table 是一个包含详细统计量的矩阵
  median_data <- as.data.frame(summary_median$table) %>%
    rownames_to_column("Group") %>% # 将行名（分组名）转为列
    dplyr::select(Group, median, `0.95LCL`, `0.95UCL`) %>% # 选择需要的列
    rename(
      `Median Survival (days)` = median,
      `Median 95% CI Lower` = `0.95LCL`,
      `Median 95% CI Upper` = `0.95UCL`
    )
  
  # 定义5年对应的天数
  five_year_days <- 5 * 12
  
  # 使用 summary 函数在指定时间点（5年）提取生存率和95%CI
  # 这是最推荐的方法，因为它会进行线性插值
  summary_5yr <- summary(km_os, times = five_year_days)
  
  # 从 summary 结果中构建最终的整洁表格
  # summary_5yr$surv, lower, upper 是按分组顺序排列的向量
  result_table <- tibble(
    Group = c("TP53 Mut","TP53 Wt"),
    `Median Survival (days)` = median_surv,
    `5-year Survival (%)` = round(summary_5yr$surv * 100, 1),
    `95% CI` = paste0(
      round(summary_5yr$lower * 100, 1),
      " - ",
      round(summary_5yr$upper * 100, 1)
    )
  )
  info_table <- cbind(median_data,result_table[,c(2,3)])
  info_table[is.na(info_table)] <- "."
  info_table <- info_table %>%
    mutate(`Median in Months`= paste0(`Median Survival (days)`,"(",`Median 95% CI Lower`,"-",`Median 95% CI Upper`,")")) %>%
    mutate(`5-Year Estimate`= paste0(`5-year Survival (%)`,"(",`95% CI`,")"))%>%
    dplyr::select(Group,`Median in Months`,`5-Year Estimate`)
  
    
    
  
  p <- ggsurvplot(km_os, 
                  conf.int=F, #是否显示生存率的95%CI
                  risk.table=T,#显示风险表
                  palette="jco", #柳叶刀配色
                  title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage III&IV"), #大标题
                  #risk.table.height = 0.15,#风险表的高度比例
                  tables.height = 0.17,
                  pval=T,  
                  pval.method=T,
                  surv.median.line="v",
                  cumevents=T,
                  cumcensor = T,
                  ylab="OS(%)",
                  xlab="Months"
  )
  print(p)
  dev.off()
  
  #rfs
  pdf(filename_rfs,width = 8,height = 8,onefile = F)
  km_rfs<-survfit(Surv(RFS,RFS_E)~Category,data = clin_mut1)
  p <- ggsurvplot(km_rfs, 
                  conf.int=F, #是否显示生存率的95%CI
                  risk.table=T,#显示风险表
                  palette="jco", #柳叶刀配色
                  title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage III&IV"), #大标题
                  #risk.table.height = 0.15,#风险表的高度比例
                  tables.height = 0.17,
                  pval=T,  
                  pval.method=T,
                  surv.median.line = "v",
                  cumevents=T,
                  cumcensor = T,
                  ylab="RFS(%)",
                  xlab="Months"
  )
  print(p)
  dev.off()
  #III&IV-Cox
  #多因素
  clin_mut1$Differentiation <- factor(clin_mut1$Differentiation, levels = c("Well", "Mod", "Poor"), ordered = TRUE)
  contrasts(clin_mut1$Differentiation) <- contr.sdif(3)
  #clin_mut1$Age_E <- factor(clin_mut1$Age_E,levels=c("<40","40-49","50-59","60-69",">=70"))
  
  
  clin_mut1$Category <- factor(clin_mut1$Category,levels = c("TP53 Wt","TP53 Mut"))
  clin_mut1$Smoking_E <- factor(clin_mut1$Smoking_E,levels=c("Never","Ever"))
  
  multicox <- coxph(Surv(OS,OS_E)~Category+Sex+Age_E+Smoking_E+Differentiation, data = clin_mut1)
  multisum <- summary(multicox)
  
  
  
  
  #提取所有基因的多因素COX回归分析结果至multiresult对象中
  factor <-  c("TP53 Mut(vs Wt)","Male(vs Female)","Age>=60 (vs <60)","Ever Smoking(vs Never)","Differentiation Mod(vs Well)","Differentiation Poor(vs Mod)")
  HR <- multisum$coefficients[,2]
  L95CI <- multisum$conf.int[,3]
  H95CI <- multisum$conf.int[,4]
  pvalue <- multisum$coefficients[,5]
  proportion <- clin_mut1 %>%
    summarise(
      Mut_n = sum(Category == "TP53 Mut", na.rm = TRUE),
      Male_n = sum(Sex == "Male", na.rm = TRUE),
      Age_n = sum(Age_E  == ">=60",na.rm = TRUE),
      Smoking_n =sum(Smoking_E== "Ever",na.rm = TRUE),
      Differentiation_mod_n = sum(Differentiation=="Mod",na.rm = TRUE),
      Differentiation_poor_n = sum(Differentiation=="Poor",na.rm = TRUE),
      Total_N = n(),
      Proportion_mut = paste0(Mut_n, "/", Total_N, " (", round(Mut_n / Total_N * 100, 1), "%)"),
      Proportion_male = paste0(Male_n, "/", Total_N, " (", round(Male_n / Total_N * 100, 1), "%)"),
      Proportion_age = paste0(Age_n, "/", Total_N, " (", round(Age_n / Total_N * 100, 1), "%)"),
      Proportion_smoking = paste0(Smoking_n, "/", Total_N, " (", round(Smoking_n / Total_N * 100, 1), "%)"),
      Proportion_mod = paste0(Differentiation_mod_n, "/", Total_N, " (", round(Differentiation_mod_n / Total_N * 100, 1), "%)"),
      Proportion_poor = paste0(Differentiation_poor_n, "/", Total_N, " (", round(Differentiation_poor_n / Total_N * 100, 1), "%)")
    ) 
  multiresult <- data.frame(factor=factor,
                            proportion=as.vector(t(proportion[1,c(8:13)])),
                            HR=HR,
                            L95CI=L95CI,
                            H95CI=H95CI,
                            pvalue=pvalue)
  multiresult$HR <- round(multiresult$HR,digits = 2)
  multiresult$L95CI <- round(multiresult$L95CI,digits = 2)
  multiresult$H95CI <- round(multiresult$H95CI,digits = 2)
  multiresult$pvalue  <-  ifelse(
    multiresult$pvalue < 0.0001,        # 条件：p值是否小于0.0001
    "<0.0001",              # 如果是，则返回此字符串
    sprintf("%.4f",multiresult$pvalu) # 如果否，则格式化为保留四位小数的字符串
  )
  
  multiresult$`HR(95% CI)` <- paste0(multiresult$HR,"(",multiresult$L95CI,"-",multiresult$H95CI,")")
  multiresult <- multiresult[,c(1,2,7,6)]
  print(multiresult)
}

##EGFR
gene <- "EGFR"
#I
if(1){
  #I-KM
  clin_mut1 <- clin_mut %>%
    filter(Stage2=="I") %>%
    mutate(Category=EGFR_category) %>%
    filter(Category %in% c("EGFR Canonical","EGFR WT"))
  folder_path <- "~/molecular_staging/output"
  filename_os <- paste0(gene,"_OS stage I.pdf")
  filename_os <- file.path(folder_path,filename_os)
  filename_rfs <- paste0(gene,"_RFS stage I.pdf")
  filename_rfs <- file.path(folder_path,filename_rfs)
  #os
  pdf(filename_os,width = 8,height = 8,onefile = F)
  km_os<-survfit(Surv(OS,OS_E)~Category,data = clin_mut1)
  summary_median <- summary(km_os)
  
  # 从 summary_median 对象中提取我们需要的信息
  # summary_median$table 是一个包含详细统计量的矩阵
  median_data <- as.data.frame(summary_median$table) %>%
    rownames_to_column("Group") %>% # 将行名（分组名）转为列
    dplyr::select(Group, median, `0.95LCL`, `0.95UCL`) %>% # 选择需要的列
    rename(
      `Median Survival (days)` = median,
      `Median 95% CI Lower` = `0.95LCL`,
      `Median 95% CI Upper` = `0.95UCL`
    )
  
  # 定义5年对应的天数
  five_year_days <- 5 * 12
  
  # 使用 summary 函数在指定时间点（5年）提取生存率和95%CI
  # 这是最推荐的方法，因为它会进行线性插值
  summary_5yr <- summary(km_os, times = five_year_days)
  
  # 从 summary 结果中构建最终的整洁表格
  # summary_5yr$surv, lower, upper 是按分组顺序排列的向量
  result_table <- tibble(
    Group = c("EGFR Mut","EGFR Wt"),
    `Median Survival (days)` = median_surv,
    `5-year Survival (%)` = round(summary_5yr$surv * 100, 1),
    `95% CI` = paste0(
      round(summary_5yr$lower * 100, 1),
      " - ",
      round(summary_5yr$upper * 100, 1)
    )
  )
  info_table <- cbind(median_data,result_table[,c(2,3)])
  info_table[is.na(info_table)] <- "."
  info_table <- info_table %>%
    mutate(`Median in Months`= paste0(`Median Survival (days)`,"(",`Median 95% CI Lower`,"-",`Median 95% CI Upper`,")")) %>%
    mutate(`5-Year Estimate`= paste0(`5-year Survival (%)`,"(",`95% CI`,")"))%>%
    dplyr::select(Group,`Median in Months`,`5-Year Estimate`)
  
  p <- ggsurvplot(km_os, 
                  conf.int=F, #是否显示生存率的95%CI
                  risk.table=T,#显示风险表
                  palette="jco", #柳叶刀配色
                  title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage I"), #大标题
                  #risk.table.height = 0.15,#风险表的高度比例
                  tables.height = 0.17,
                  pval=T,  
                  pval.method=T,
                  surv.median.line="v",
                  cumevents=T,
                  cumcensor = T,
                  ylab="OS(%)",
                  xlab="Months"
  )
  print(p)
  dev.off()
  
  #rfs
  pdf(filename_rfs,width = 8,height = 8,onefile = F)
  km_rfs<-survfit(Surv(RFS,RFS_E)~Category,data = clin_mut1)
  p <- ggsurvplot(km_rfs, 
                  conf.int=F, #是否显示生存率的95%CI
                  risk.table=T,#显示风险表
                  palette="jco", #柳叶刀配色
                  title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage I"), #大标题
                  #risk.table.height = 0.15,#风险表的高度比例
                  tables.height = 0.17,
                  pval=T,  
                  pval.method=T,
                  surv.median.line="v",
                  cumevents=T,
                  cumcensor = T,
                  ylab="RFS(%)",
                  xlab="Months"
  )
  print(p)
  dev.off()
  #I-Cox
  #多因素
  clin_mut1$Differentiation <- factor(clin_mut1$Differentiation, levels = c("Well", "Mod", "Poor"), ordered = TRUE)
  contrasts(clin_mut1$Differentiation) <- contr.sdif(3)
  #clin_mut1$Age_E <- factor(clin_mut1$Age_E,levels=c("<40","40-49","50-59","60-69",">=70"))
  
  
  clin_mut1$Category <- factor(clin_mut1$Category,levels = c("EGFR WT","EGFR Canonical"))
  clin_mut1$Smoking_E <- factor(clin_mut1$Smoking_E,levels=c("Never","Ever"))
  
  multicox <- coxph(Surv(OS,OS_E)~Category+Sex+Age_E+Smoking_E+Differentiation, data = clin_mut1)
  multisum <- summary(multicox)
  
  
  
  
  #提取所有基因的多因素COX回归分析结果至multiresult对象中
  factor <-  c("EGFR Mut(vs Wt)","Male(vs Female)","Age>=60 (vs <60)","Ever Smoking(vs Never)","Differentiation Mod(vs Well)","Differentiation Poor(vs Mod)")
  HR <- multisum$coefficients[,2]
  L95CI <- multisum$conf.int[,3]
  H95CI <- multisum$conf.int[,4]
  pvalue <- multisum$coefficients[,5]
  proportion <- clin_mut1 %>%
    summarise(
      Mut_n = sum(Category == "EGFR Canonical", na.rm = TRUE),
      Male_n = sum(Sex == "Male", na.rm = TRUE),
      Age_n = sum(Age_E  == ">=60",na.rm = TRUE),
      Smoking_n =sum(Smoking_E== "Ever",na.rm = TRUE),
      Differentiation_mod_n = sum(Differentiation=="Mod",na.rm = TRUE),
      Differentiation_poor_n = sum(Differentiation=="Poor",na.rm = TRUE),
      Total_N = n(),
      Proportion_mut = paste0(Mut_n, "/", Total_N, " (", round(Mut_n / Total_N * 100, 1), "%)"),
      Proportion_male = paste0(Male_n, "/", Total_N, " (", round(Male_n / Total_N * 100, 1), "%)"),
      Proportion_age = paste0(Age_n, "/", Total_N, " (", round(Age_n / Total_N * 100, 1), "%)"),
      Proportion_smoking = paste0(Smoking_n, "/", Total_N, " (", round(Smoking_n / Total_N * 100, 1), "%)"),
      Proportion_mod = paste0(Differentiation_mod_n, "/", Total_N, " (", round(Differentiation_mod_n / Total_N * 100, 1), "%)"),
      Proportion_poor = paste0(Differentiation_poor_n, "/", Total_N, " (", round(Differentiation_poor_n / Total_N * 100, 1), "%)")
    ) 
  multiresult <- data.frame(factor=factor,
                            proportion=as.vector(t(proportion[1,c(8:13)])),
                            HR=HR,
                            L95CI=L95CI,
                            H95CI=H95CI,
                            pvalue=pvalue)
  multiresult$HR <- round(multiresult$HR,digits = 2)
  multiresult$L95CI <- round(multiresult$L95CI,digits = 2)
  multiresult$H95CI <- round(multiresult$H95CI,digits = 2)
  multiresult$pvalue  <-  ifelse(
    multiresult$pvalue < 0.0001,        # 条件：p值是否小于0.0001
    "<0.0001",              # 如果是，则返回此字符串
    sprintf("%.4f",multiresult$pvalu) # 如果否，则格式化为保留四位小数的字符串
  )
  
  multiresult$`HR(95% CI)` <- paste0(multiresult$HR,"(",multiresult$L95CI,"-",multiresult$H95CI,")")
  multiresult <- multiresult[,c(1,2,7,6)]
  print(multiresult)
}
#II
if(1){
  #I-KM
  clin_mut1 <- clin_mut %>%
    filter(Stage2=="II") %>%
    mutate(Category=EGFR_category) %>%
    filter(Category %in% c("EGFR Canonical","EGFR WT"))
  folder_path <- "~/molecular_staging/output"
  filename_os <- paste0(gene,"_OS stage II.pdf")
  filename_os <- file.path(folder_path,filename_os)
  filename_rfs <- paste0(gene,"_RFS stage II.pdf")
  filename_rfs <- file.path(folder_path,filename_rfs)
  #os
  pdf(filename_os,width = 8,height = 8,onefile = F)
  km_os<-survfit(Surv(OS,OS_E)~Category,data = clin_mut1)
  summary_median <- summary(km_os)
  
  # 从 summary_median 对象中提取我们需要的信息
  # summary_median$table 是一个包含详细统计量的矩阵
  median_data <- as.data.frame(summary_median$table) %>%
    rownames_to_column("Group") %>% # 将行名（分组名）转为列
    dplyr::select(Group, median, `0.95LCL`, `0.95UCL`) %>% # 选择需要的列
    rename(
      `Median Survival (days)` = median,
      `Median 95% CI Lower` = `0.95LCL`,
      `Median 95% CI Upper` = `0.95UCL`
    )
  
  # 定义5年对应的天数
  five_year_days <- 5 * 12
  
  # 使用 summary 函数在指定时间点（5年）提取生存率和95%CI
  # 这是最推荐的方法，因为它会进行线性插值
  summary_5yr <- summary(km_os, times = five_year_days)
  
  # 从 summary 结果中构建最终的整洁表格
  # summary_5yr$surv, lower, upper 是按分组顺序排列的向量
  result_table <- tibble(
    Group = c("EGFR Mut","EGFR Wt"),
    `Median Survival (days)` = median_surv,
    `5-year Survival (%)` = round(summary_5yr$surv * 100, 1),
    `95% CI` = paste0(
      round(summary_5yr$lower * 100, 1),
      " - ",
      round(summary_5yr$upper * 100, 1)
    )
  )
  info_table <- cbind(median_data,result_table[,c(2,3)])
  info_table[is.na(info_table)] <- "."
  info_table <- info_table %>%
    mutate(`Median in Months`= paste0(`Median Survival (days)`,"(",`Median 95% CI Lower`,"-",`Median 95% CI Upper`,")")) %>%
    mutate(`5-Year Estimate`= paste0(`5-year Survival (%)`,"(",`95% CI`,")"))%>%
    dplyr::select(Group,`Median in Months`,`5-Year Estimate`)
  
  p <- ggsurvplot(km_os, 
                  conf.int=F, #是否显示生存率的95%CI
                  risk.table=T,#显示风险表
                  palette="jco", #柳叶刀配色
                  title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage II"), #大标题
                  #risk.table.height = 0.15,#风险表的高度比例
                  tables.height = 0.17,
                  pval=T,  
                  pval.method=T,
                  surv.median.line="v",
                  cumevents=T,
                  cumcensor = T,
                  ylab="OS(%)",
                  xlab="Months"
  )
  print(p)
  dev.off()
  
  #rfs
  pdf(filename_rfs,width = 8,height = 8,onefile = F)
  km_rfs<-survfit(Surv(RFS,RFS_E)~Category,data = clin_mut1)
  p <- ggsurvplot(km_rfs, 
                  conf.int=F, #是否显示生存率的95%CI
                  risk.table=T,#显示风险表
                  palette="jco", #柳叶刀配色
                  title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage II"), #大标题
                  #risk.table.height = 0.15,#风险表的高度比例
                  tables.height = 0.17,
                  pval=T,  
                  pval.method=T,
                  surv.median.line="v",
                  cumevents=T,
                  cumcensor = T,
                  ylab="RFS(%)",
                  xlab="Months"
  )
  print(p)
  dev.off()
  #I-Cox
  #多因素
  clin_mut1$Differentiation <- factor(clin_mut1$Differentiation, levels = c("Well", "Mod", "Poor"), ordered = TRUE)
  contrasts(clin_mut1$Differentiation) <- contr.sdif(3)
  #clin_mut1$Age_E <- factor(clin_mut1$Age_E,levels=c("<40","40-49","50-59","60-69",">=70"))
  
  
  clin_mut1$Category <- factor(clin_mut1$Category,levels = c("EGFR WT","EGFR Canonical"))
  clin_mut1$Smoking_E <- factor(clin_mut1$Smoking_E,levels=c("Never","Ever"))
  
  multicox <- coxph(Surv(OS,OS_E)~Category+Sex+Age_E+Smoking_E+Differentiation, data = clin_mut1)
  multisum <- summary(multicox)
  
  
  
  
  #提取所有基因的多因素COX回归分析结果至multiresult对象中
  factor <-  c("EGFR Mut(vs Wt)","Male(vs Female)","Age>=60 (vs <60)","Ever Smoking(vs Never)","Differentiation Mod(vs Well)","Differentiation Poor(vs Mod)")
  HR <- multisum$coefficients[,2]
  L95CI <- multisum$conf.int[,3]
  H95CI <- multisum$conf.int[,4]
  pvalue <- multisum$coefficients[,5]
  proportion <- clin_mut1 %>%
    summarise(
      Mut_n = sum(Category == "EGFR Canonical", na.rm = TRUE),
      Male_n = sum(Sex == "Male", na.rm = TRUE),
      Age_n = sum(Age_E  == ">=60",na.rm = TRUE),
      Smoking_n =sum(Smoking_E== "Ever",na.rm = TRUE),
      Differentiation_mod_n = sum(Differentiation=="Mod",na.rm = TRUE),
      Differentiation_poor_n = sum(Differentiation=="Poor",na.rm = TRUE),
      Total_N = n(),
      Proportion_mut = paste0(Mut_n, "/", Total_N, " (", round(Mut_n / Total_N * 100, 1), "%)"),
      Proportion_male = paste0(Male_n, "/", Total_N, " (", round(Male_n / Total_N * 100, 1), "%)"),
      Proportion_age = paste0(Age_n, "/", Total_N, " (", round(Age_n / Total_N * 100, 1), "%)"),
      Proportion_smoking = paste0(Smoking_n, "/", Total_N, " (", round(Smoking_n / Total_N * 100, 1), "%)"),
      Proportion_mod = paste0(Differentiation_mod_n, "/", Total_N, " (", round(Differentiation_mod_n / Total_N * 100, 1), "%)"),
      Proportion_poor = paste0(Differentiation_poor_n, "/", Total_N, " (", round(Differentiation_poor_n / Total_N * 100, 1), "%)")
    ) 
  multiresult <- data.frame(factor=factor,
                            proportion=as.vector(t(proportion[1,c(8:13)])),
                            HR=HR,
                            L95CI=L95CI,
                            H95CI=H95CI,
                            pvalue=pvalue)
  multiresult$HR <- round(multiresult$HR,digits = 2)
  multiresult$L95CI <- round(multiresult$L95CI,digits = 2)
  multiresult$H95CI <- round(multiresult$H95CI,digits = 2)
  multiresult$pvalue  <-  ifelse(
    multiresult$pvalue < 0.0001,        # 条件：p值是否小于0.0001
    "<0.0001",              # 如果是，则返回此字符串
    sprintf("%.4f",multiresult$pvalu) # 如果否，则格式化为保留四位小数的字符串
  )
  
  multiresult$`HR(95% CI)` <- paste0(multiresult$HR,"(",multiresult$L95CI,"-",multiresult$H95CI,")")
  multiresult <- multiresult[,c(1,2,7,6)]
  print(multiresult)
}
#III&IV
if(1){
  #III&IV-KM
  clin_mut1 <- clin_mut %>%
    filter(Stage2 %in% c("III","IV")) %>%
    mutate(Category=EGFR_category) %>%
    filter(Category %in% c("EGFR Canonical","EGFR WT"))
  folder_path <- "~/molecular_staging/output"
  filename_os <- paste0(gene,"_OS stage III&IV.pdf")
  filename_os <- file.path(folder_path,filename_os)
  filename_rfs <- paste0(gene,"_RFS stage III&IV.pdf")
  filename_rfs <- file.path(folder_path,filename_rfs)
  #os
  pdf(filename_os,width = 8,height = 8,onefile = F)
  km_os<-survfit(Surv(OS,OS_E)~Category,data = clin_mut1)
  summary_median <- summary(km_os)
  
  # 从 summary_median 对象中提取我们需要的信息
  # summary_median$table 是一个包含详细统计量的矩阵
  median_data <- as.data.frame(summary_median$table) %>%
    rownames_to_column("Group") %>% # 将行名（分组名）转为列
    dplyr::select(Group, median, `0.95LCL`, `0.95UCL`) %>% # 选择需要的列
    rename(
      `Median Survival (days)` = median,
      `Median 95% CI Lower` = `0.95LCL`,
      `Median 95% CI Upper` = `0.95UCL`
    )
  
  # 定义5年对应的天数
  five_year_days <- 5 * 12
  
  # 使用 summary 函数在指定时间点（5年）提取生存率和95%CI
  # 这是最推荐的方法，因为它会进行线性插值
  summary_5yr <- summary(km_os, times = five_year_days)
  
  # 从 summary 结果中构建最终的整洁表格
  # summary_5yr$surv, lower, upper 是按分组顺序排列的向量
  result_table <- tibble(
    Group = c("EGFR Mut","EGFR Wt"),
    `Median Survival (days)` = median_surv,
    `5-year Survival (%)` = round(summary_5yr$surv * 100, 1),
    `95% CI` = paste0(
      round(summary_5yr$lower * 100, 1),
      " - ",
      round(summary_5yr$upper * 100, 1)
    )
  )
  info_table <- cbind(median_data,result_table[,c(2,3)])
  info_table[is.na(info_table)] <- "."
  info_table <- info_table %>%
    mutate(`Median in Months`= paste0(`Median Survival (days)`,"(",`Median 95% CI Lower`,"-",`Median 95% CI Upper`,")")) %>%
    mutate(`5-Year Estimate`= paste0(`5-year Survival (%)`,"(",`95% CI`,")"))%>%
    dplyr::select(Group,`Median in Months`,`5-Year Estimate`)
  
  p <- ggsurvplot(km_os, 
                  conf.int=F, #是否显示生存率的95%CI
                  risk.table=T,#显示风险表
                  palette="jco", #柳叶刀配色
                  title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage III&IV"), #大标题
                  #risk.table.height = 0.15,#风险表的高度比例
                  tables.height = 0.17,
                  pval=T,  
                  pval.method=T,
                  surv.median.line="v",
                  cumevents=T,
                  cumcensor = T,
                  ylab="OS(%)",
                  xlab="Months"
  )
  print(p)
  dev.off()
  
  #rfs
  pdf(filename_rfs,width = 8,height = 8,onefile = F)
  km_rfs<-survfit(Surv(RFS,RFS_E)~Category,data = clin_mut1)
  p <- ggsurvplot(km_rfs, 
                  conf.int=F, #是否显示生存率的95%CI
                  risk.table=T,#显示风险表
                  palette="jco", #柳叶刀配色
                  title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage III&IV"), #大标题
                  #risk.table.height = 0.15,#风险表的高度比例
                  tables.height = 0.17,
                  pval=T,  
                  pval.method=T,
                  surv.median.line="v",
                  cumevents=T,
                  cumcensor = T,
                  ylab="RFS(%)",
                  xlab="Months"
  )
  print(p)
  dev.off()
  #I-Cox
  #多因素
  clin_mut1$Differentiation <- factor(clin_mut1$Differentiation, levels = c("Well", "Mod", "Poor"), ordered = TRUE)
  contrasts(clin_mut1$Differentiation) <- contr.sdif(3)
  #clin_mut1$Age_E <- factor(clin_mut1$Age_E,levels=c("<40","40-49","50-59","60-69",">=70"))
  
  
  clin_mut1$Category <- factor(clin_mut1$Category,levels = c("EGFR WT","EGFR Canonical"))
  clin_mut1$Smoking_E <- factor(clin_mut1$Smoking_E,levels=c("Never","Ever"))
  
  multicox <- coxph(Surv(OS,OS_E)~Category+Sex+Age_E+Smoking_E+Differentiation, data = clin_mut1)
  multisum <- summary(multicox)
  
  
  
  
  #提取所有基因的多因素COX回归分析结果至multiresult对象中
  factor <-  c("EGFR Mut(vs Wt)","Male(vs Female)","Age>=60 (vs <60)","Ever Smoking(vs Never)","Differentiation Mod(vs Well)","Differentiation Poor(vs Mod)")
  HR <- multisum$coefficients[,2]
  L95CI <- multisum$conf.int[,3]
  H95CI <- multisum$conf.int[,4]
  pvalue <- multisum$coefficients[,5]
  proportion <- clin_mut1 %>%
    summarise(
      Mut_n = sum(Category == "EGFR Canonical", na.rm = TRUE),
      Male_n = sum(Sex == "Male", na.rm = TRUE),
      Age_n = sum(Age_E  == ">=60",na.rm = TRUE),
      Smoking_n =sum(Smoking_E== "Ever",na.rm = TRUE),
      Differentiation_mod_n = sum(Differentiation=="Mod",na.rm = TRUE),
      Differentiation_poor_n = sum(Differentiation=="Poor",na.rm = TRUE),
      Total_N = n(),
      Proportion_mut = paste0(Mut_n, "/", Total_N, " (", round(Mut_n / Total_N * 100, 1), "%)"),
      Proportion_male = paste0(Male_n, "/", Total_N, " (", round(Male_n / Total_N * 100, 1), "%)"),
      Proportion_age = paste0(Age_n, "/", Total_N, " (", round(Age_n / Total_N * 100, 1), "%)"),
      Proportion_smoking = paste0(Smoking_n, "/", Total_N, " (", round(Smoking_n / Total_N * 100, 1), "%)"),
      Proportion_mod = paste0(Differentiation_mod_n, "/", Total_N, " (", round(Differentiation_mod_n / Total_N * 100, 1), "%)"),
      Proportion_poor = paste0(Differentiation_poor_n, "/", Total_N, " (", round(Differentiation_poor_n / Total_N * 100, 1), "%)")
    ) 
  multiresult <- data.frame(factor=factor,
                            proportion=as.vector(t(proportion[1,c(8:13)])),
                            HR=HR,
                            L95CI=L95CI,
                            H95CI=H95CI,
                            pvalue=pvalue)
  multiresult$HR <- round(multiresult$HR,digits = 2)
  multiresult$L95CI <- round(multiresult$L95CI,digits = 2)
  multiresult$H95CI <- round(multiresult$H95CI,digits = 2)
  multiresult$pvalue  <-  ifelse(
    multiresult$pvalue < 0.0001,        # 条件：p值是否小于0.0001
    "<0.0001",              # 如果是，则返回此字符串
    sprintf("%.4f",multiresult$pvalu) # 如果否，则格式化为保留四位小数的字符串
  )
  
  multiresult$`HR(95% CI)` <- paste0(multiresult$HR,"(",multiresult$L95CI,"-",multiresult$H95CI,")")
  multiresult <- multiresult[,c(1,2,7,6)]
  print(multiresult)
}

##KRAS
gene <- "KRAS"
#I
if(1){
  #I-KM
  clin_mut1 <- clin_mut %>%
    filter(Stage2=="I") %>%
    mutate(Category=KRAS_category) %>%
    filter(Category %in% c("KRAS Canonical","KRAS WT"))
  folder_path <- "~/molecular_staging/output"
  filename_os <- paste0(gene,"_OS stage I.pdf")
  filename_os <- file.path(folder_path,filename_os)
  filename_rfs <- paste0(gene,"_RFS stage I.pdf")
  filename_rfs <- file.path(folder_path,filename_rfs)
  #os
  pdf(filename_os,width = 8,height = 8,onefile = F)
  km_os<-survfit(Surv(OS,OS_E)~Category,data = clin_mut1)
  summary_median <- summary(km_os)
  
  # 从 summary_median 对象中提取我们需要的信息
  # summary_median$table 是一个包含详细统计量的矩阵
  median_data <- as.data.frame(summary_median$table) %>%
    rownames_to_column("Group") %>% # 将行名（分组名）转为列
    dplyr::select(Group, median, `0.95LCL`, `0.95UCL`) %>% # 选择需要的列
    rename(
      `Median Survival (days)` = median,
      `Median 95% CI Lower` = `0.95LCL`,
      `Median 95% CI Upper` = `0.95UCL`
    )
  
  # 定义5年对应的天数
  five_year_days <- 5 * 12
  
  # 使用 summary 函数在指定时间点（5年）提取生存率和95%CI
  # 这是最推荐的方法，因为它会进行线性插值
  summary_5yr <- summary(km_os, times = five_year_days)
  
  # 从 summary 结果中构建最终的整洁表格
  # summary_5yr$surv, lower, upper 是按分组顺序排列的向量
  result_table <- tibble(
    Group = c("KRAS Mut","KRAS Wt"),
    `Median Survival (days)` = median_surv,
    `5-year Survival (%)` = round(summary_5yr$surv * 100, 1),
    `95% CI` = paste0(
      round(summary_5yr$lower * 100, 1),
      " - ",
      round(summary_5yr$upper * 100, 1)
    )
  )
  info_table <- cbind(median_data,result_table[,c(2,3)])
  info_table[is.na(info_table)] <- "."
  info_table <- info_table %>%
    mutate(`Median in Months`= paste0(`Median Survival (days)`,"(",`Median 95% CI Lower`,"-",`Median 95% CI Upper`,")")) %>%
    mutate(`5-Year Estimate`= paste0(`5-year Survival (%)`,"(",`95% CI`,")"))%>%
    dplyr::select(Group,`Median in Months`,`5-Year Estimate`)
  
  p <- ggsurvplot(km_os, 
                  conf.int=F, #是否显示生存率的95%CI
                  risk.table=T,#显示风险表
                  palette="jco", #柳叶刀配色
                  title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage I"), #大标题
                  #risk.table.height = 0.15,#风险表的高度比例
                  tables.height = 0.17,
                  pval=T,  
                  pval.method=T,
                  surv.median.line="v",
                  cumevents=T,
                  cumcensor = T,
                  ylab="OS(%)",
                  xlab="Months"
  )
  print(p)
  dev.off()
  
  #rfs
  pdf(filename_rfs,width = 8,height = 8,onefile = F)
  km_rfs<-survfit(Surv(RFS,RFS_E)~Category,data = clin_mut1)
  p <- ggsurvplot(km_rfs, 
                  conf.int=F, #是否显示生存率的95%CI
                  risk.table=T,#显示风险表
                  palette="jco", #柳叶刀配色
                  title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage I"), #大标题
                  #risk.table.height = 0.15,#风险表的高度比例
                  tables.height = 0.17,
                  pval=T,  
                  pval.method=T,
                  surv.median.line="v",
                  cumevents=T,
                  cumcensor = T,
                  ylab="RFS(%)",
                  xlab="Months"
  )
  print(p)
  dev.off()
  #I-Cox
  #多因素
  clin_mut1$Differentiation <- factor(clin_mut1$Differentiation, levels = c("Well", "Mod", "Poor"), ordered = TRUE)
  contrasts(clin_mut1$Differentiation) <- contr.sdif(3)
  #clin_mut1$Age_E <- factor(clin_mut1$Age_E,levels=c("<40","40-49","50-59","60-69",">=70"))
  
  
  clin_mut1$Category <- factor(clin_mut1$Category,levels = c("KRAS WT","KRAS Canonical"))
  clin_mut1$Smoking_E <- factor(clin_mut1$Smoking_E,levels=c("Never","Ever"))
  
  multicox <- coxph(Surv(OS,OS_E)~Category+Sex+Age_E+Smoking_E+Differentiation, data = clin_mut1)
  multisum <- summary(multicox)
  
  
  
  
  #提取所有基因的多因素COX回归分析结果至multiresult对象中
  factor <-  c("KRAS Mut(vs Wt)","Male(vs Female)","Age>=60 (vs <60)","Ever Smoking(vs Never)","Differentiation Mod(vs Well)","Differentiation Poor(vs Mod)")
  HR <- multisum$coefficients[,2]
  L95CI <- multisum$conf.int[,3]
  H95CI <- multisum$conf.int[,4]
  pvalue <- multisum$coefficients[,5]
  proportion <- clin_mut1 %>%
    summarise(
      Mut_n = sum(Category == "KRAS Canonical", na.rm = TRUE),
      Male_n = sum(Sex == "Male", na.rm = TRUE),
      Age_n = sum(Age_E  == ">=60",na.rm = TRUE),
      Smoking_n =sum(Smoking_E== "Ever",na.rm = TRUE),
      Differentiation_mod_n = sum(Differentiation=="Mod",na.rm = TRUE),
      Differentiation_poor_n = sum(Differentiation=="Poor",na.rm = TRUE),
      Total_N = n(),
      Proportion_mut = paste0(Mut_n, "/", Total_N, " (", round(Mut_n / Total_N * 100, 1), "%)"),
      Proportion_male = paste0(Male_n, "/", Total_N, " (", round(Male_n / Total_N * 100, 1), "%)"),
      Proportion_age = paste0(Age_n, "/", Total_N, " (", round(Age_n / Total_N * 100, 1), "%)"),
      Proportion_smoking = paste0(Smoking_n, "/", Total_N, " (", round(Smoking_n / Total_N * 100, 1), "%)"),
      Proportion_mod = paste0(Differentiation_mod_n, "/", Total_N, " (", round(Differentiation_mod_n / Total_N * 100, 1), "%)"),
      Proportion_poor = paste0(Differentiation_poor_n, "/", Total_N, " (", round(Differentiation_poor_n / Total_N * 100, 1), "%)")
    ) 
  multiresult <- data.frame(factor=factor,
                            proportion=as.vector(t(proportion[1,c(8:13)])),
                            HR=HR,
                            L95CI=L95CI,
                            H95CI=H95CI,
                            pvalue=pvalue)
  multiresult$HR <- round(multiresult$HR,digits = 2)
  multiresult$L95CI <- round(multiresult$L95CI,digits = 2)
  multiresult$H95CI <- round(multiresult$H95CI,digits = 2)
  multiresult$pvalue  <-  ifelse(
    multiresult$pvalue < 0.0001,        # 条件：p值是否小于0.0001
    "<0.0001",              # 如果是，则返回此字符串
    sprintf("%.4f",multiresult$pvalu) # 如果否，则格式化为保留四位小数的字符串
  )
  
  multiresult$`HR(95% CI)` <- paste0(multiresult$HR,"(",multiresult$L95CI,"-",multiresult$H95CI,")")
  multiresult <- multiresult[,c(1,2,7,6)]
  print(multiresult)
}
#II
if(1){
  #I-KM
  clin_mut1 <- clin_mut %>%
    filter(Stage2=="II") %>%
    mutate(Category=KRAS_category) %>%
    filter(Category %in% c("KRAS Canonical","KRAS WT"))
  folder_path <- "~/molecular_staging/output"
  filename_os <- paste0(gene,"_OS stage II.pdf")
  filename_os <- file.path(folder_path,filename_os)
  filename_rfs <- paste0(gene,"_RFS stage II.pdf")
  filename_rfs <- file.path(folder_path,filename_rfs)
  #os
  pdf(filename_os,width = 8,height = 8,onefile = F)
  km_os<-survfit(Surv(OS,OS_E)~Category,data = clin_mut1)
  summary_median <- summary(km_os)
  
  # 从 summary_median 对象中提取我们需要的信息
  # summary_median$table 是一个包含详细统计量的矩阵
  median_data <- as.data.frame(summary_median$table) %>%
    rownames_to_column("Group") %>% # 将行名（分组名）转为列
    dplyr::select(Group, median, `0.95LCL`, `0.95UCL`) %>% # 选择需要的列
    rename(
      `Median Survival (days)` = median,
      `Median 95% CI Lower` = `0.95LCL`,
      `Median 95% CI Upper` = `0.95UCL`
    )
  
  # 定义5年对应的天数
  five_year_days <- 5 * 12
  
  # 使用 summary 函数在指定时间点（5年）提取生存率和95%CI
  # 这是最推荐的方法，因为它会进行线性插值
  summary_5yr <- summary(km_os, times = five_year_days)
  
  # 从 summary 结果中构建最终的整洁表格
  # summary_5yr$surv, lower, upper 是按分组顺序排列的向量
  result_table <- tibble(
    Group = c("KRAS Mut","KRAS Wt"),
    `Median Survival (days)` = median_surv,
    `5-year Survival (%)` = round(summary_5yr$surv * 100, 1),
    `95% CI` = paste0(
      round(summary_5yr$lower * 100, 1),
      " - ",
      round(summary_5yr$upper * 100, 1)
    )
  )
  info_table <- cbind(median_data,result_table[,c(2,3)])
  info_table[is.na(info_table)] <- "."
  info_table <- info_table %>%
    mutate(`Median in Months`= paste0(`Median Survival (days)`,"(",`Median 95% CI Lower`,"-",`Median 95% CI Upper`,")")) %>%
    mutate(`5-Year Estimate`= paste0(`5-year Survival (%)`,"(",`95% CI`,")"))%>%
    dplyr::select(Group,`Median in Months`,`5-Year Estimate`)
  
  p <- ggsurvplot(km_os, 
                  conf.int=F, #是否显示生存率的95%CI
                  risk.table=T,#显示风险表
                  palette="jco", #柳叶刀配色
                  title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage II"), #大标题
                  #risk.table.height = 0.15,#风险表的高度比例
                  tables.height = 0.17,
                  pval=T,  
                  pval.method=T,
                  surv.median.line="v",
                  cumevents=T,
                  cumcensor = T,
                  ylab="OS(%)",
                  xlab="Months"
  )
  print(p)
  dev.off()
  
  #rfs
  pdf(filename_rfs,width = 8,height = 8,onefile = F)
  km_rfs<-survfit(Surv(RFS,RFS_E)~Category,data = clin_mut1)
  p <- ggsurvplot(km_rfs, 
                  conf.int=F, #是否显示生存率的95%CI
                  risk.table=T,#显示风险表
                  palette="jco", #柳叶刀配色
                  title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage II"), #大标题
                  #risk.table.height = 0.15,#风险表的高度比例
                  tables.height = 0.17,
                  pval=T,  
                  pval.method=T,
                  surv.median.line="v",
                  cumevents=T,
                  cumcensor = T,
                  ylab="RFS(%)",
                  xlab="Months"
  )
  print(p)
  dev.off()
  #I-Cox
  #多因素
  clin_mut1$Differentiation <- factor(clin_mut1$Differentiation, levels = c("Well", "Mod", "Poor"), ordered = TRUE)
  contrasts(clin_mut1$Differentiation) <- contr.sdif(3)
  #clin_mut1$Age_E <- factor(clin_mut1$Age_E,levels=c("<40","40-49","50-59","60-69",">=70"))
  
  
  clin_mut1$Category <- factor(clin_mut1$Category,levels = c("KRAS WT","KRAS Canonical"))
  clin_mut1$Smoking_E <- factor(clin_mut1$Smoking_E,levels=c("Never","Ever"))
  
  multicox <- coxph(Surv(OS,OS_E)~Category+Sex+Age_E+Smoking_E+Differentiation, data = clin_mut1)
  multisum <- summary(multicox)
  
  
  
  
  #提取所有基因的多因素COX回归分析结果至multiresult对象中
  factor <-  c("KRAS Mut(vs Wt)","Male(vs Female)","Age>=60 (vs <60)","Ever Smoking(vs Never)","Differentiation Mod(vs Well)","Differentiation Poor(vs Mod)")
  HR <- multisum$coefficients[,2]
  L95CI <- multisum$conf.int[,3]
  H95CI <- multisum$conf.int[,4]
  pvalue <- multisum$coefficients[,5]
  proportion <- clin_mut1 %>%
    summarise(
      Mut_n = sum(Category == "KRAS Canonical", na.rm = TRUE),
      Male_n = sum(Sex == "Male", na.rm = TRUE),
      Age_n = sum(Age_E  == ">=60",na.rm = TRUE),
      Smoking_n =sum(Smoking_E== "Ever",na.rm = TRUE),
      Differentiation_mod_n = sum(Differentiation=="Mod",na.rm = TRUE),
      Differentiation_poor_n = sum(Differentiation=="Poor",na.rm = TRUE),
      Total_N = n(),
      Proportion_mut = paste0(Mut_n, "/", Total_N, " (", round(Mut_n / Total_N * 100, 1), "%)"),
      Proportion_male = paste0(Male_n, "/", Total_N, " (", round(Male_n / Total_N * 100, 1), "%)"),
      Proportion_age = paste0(Age_n, "/", Total_N, " (", round(Age_n / Total_N * 100, 1), "%)"),
      Proportion_smoking = paste0(Smoking_n, "/", Total_N, " (", round(Smoking_n / Total_N * 100, 1), "%)"),
      Proportion_mod = paste0(Differentiation_mod_n, "/", Total_N, " (", round(Differentiation_mod_n / Total_N * 100, 1), "%)"),
      Proportion_poor = paste0(Differentiation_poor_n, "/", Total_N, " (", round(Differentiation_poor_n / Total_N * 100, 1), "%)")
    ) 
  multiresult <- data.frame(factor=factor,
                            proportion=as.vector(t(proportion[1,c(8:13)])),
                            HR=HR,
                            L95CI=L95CI,
                            H95CI=H95CI,
                            pvalue=pvalue)
  multiresult$HR <- round(multiresult$HR,digits = 2)
  multiresult$L95CI <- round(multiresult$L95CI,digits = 2)
  multiresult$H95CI <- round(multiresult$H95CI,digits = 2)
  multiresult$pvalue  <-  ifelse(
    multiresult$pvalue < 0.0001,        # 条件：p值是否小于0.0001
    "<0.0001",              # 如果是，则返回此字符串
    sprintf("%.4f",multiresult$pvalu) # 如果否，则格式化为保留四位小数的字符串
  )
  
  multiresult$`HR(95% CI)` <- paste0(multiresult$HR,"(",multiresult$L95CI,"-",multiresult$H95CI,")")
  multiresult <- multiresult[,c(1,2,7,6)]
  print(multiresult)
}
#III&IV
if(1){
  #III&IV-KM
  clin_mut1 <- clin_mut %>%
    filter(Stage2 %in% c("III","IV")) %>%
    mutate(Category=KRAS_category) %>%
    filter(Category %in% c("KRAS Canonical","KRAS WT"))
  folder_path <- "~/molecular_staging/output"
  filename_os <- paste0(gene,"_OS stage III&IV.pdf")
  filename_os <- file.path(folder_path,filename_os)
  filename_rfs <- paste0(gene,"_RFS stage III&IV.pdf")
  filename_rfs <- file.path(folder_path,filename_rfs)
  #os
  pdf(filename_os,width = 8,height = 8,onefile = F)
  km_os<-survfit(Surv(OS,OS_E)~Category,data = clin_mut1)
  summary_median <- summary(km_os)
  
  # 从 summary_median 对象中提取我们需要的信息
  # summary_median$table 是一个包含详细统计量的矩阵
  median_data <- as.data.frame(summary_median$table) %>%
    rownames_to_column("Group") %>% # 将行名（分组名）转为列
    dplyr::select(Group, median, `0.95LCL`, `0.95UCL`) %>% # 选择需要的列
    rename(
      `Median Survival (days)` = median,
      `Median 95% CI Lower` = `0.95LCL`,
      `Median 95% CI Upper` = `0.95UCL`
    )
  
  # 定义5年对应的天数
  five_year_days <- 5 * 12
  
  # 使用 summary 函数在指定时间点（5年）提取生存率和95%CI
  # 这是最推荐的方法，因为它会进行线性插值
  summary_5yr <- summary(km_os, times = five_year_days)
  
  # 从 summary 结果中构建最终的整洁表格
  # summary_5yr$surv, lower, upper 是按分组顺序排列的向量
  result_table <- tibble(
    Group = c("KRAS Mut","KRAS Wt"),
    `Median Survival (days)` = median_surv,
    `5-year Survival (%)` = round(summary_5yr$surv * 100, 1),
    `95% CI` = paste0(
      round(summary_5yr$lower * 100, 1),
      " - ",
      round(summary_5yr$upper * 100, 1)
    )
  )
  info_table <- cbind(median_data,result_table[,c(2,3)])
  info_table[is.na(info_table)] <- "."
  info_table <- info_table %>%
    mutate(`Median in Months`= paste0(`Median Survival (days)`,"(",`Median 95% CI Lower`,"-",`Median 95% CI Upper`,")")) %>%
    mutate(`5-Year Estimate`= paste0(`5-year Survival (%)`,"(",`95% CI`,")"))%>%
    dplyr::select(Group,`Median in Months`,`5-Year Estimate`)
  
  p <- ggsurvplot(km_os, 
                  conf.int=F, #是否显示生存率的95%CI
                  risk.table=T,#显示风险表
                  palette="jco", #柳叶刀配色
                  title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage III&IV"), #大标题
                  #risk.table.height = 0.15,#风险表的高度比例
                  tables.height = 0.17,
                  pval=T,  
                  pval.method=T,
                  surv.median.line="v",
                  cumevents=T,
                  cumcensor = T,
                  ylab="OS(%)",
                  xlab="Months"
  )
  print(p)
  dev.off()
  
  #rfs
  pdf(filename_rfs,width = 8,height = 8,onefile = F)
  km_rfs<-survfit(Surv(RFS,RFS_E)~Category,data = clin_mut1)
  p <- ggsurvplot(km_rfs, 
                  conf.int=F, #是否显示生存率的95%CI
                  risk.table=T,#显示风险表
                  palette="jco", #柳叶刀配色
                  title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage III&IV"), #大标题
                  #risk.table.height = 0.15,#风险表的高度比例
                  tables.height = 0.17,
                  pval=T,  
                  pval.method=T,
                  surv.median.line="v",
                  cumevents=T,
                  cumcensor = T,
                  ylab="RFS(%)",
                  xlab="Months"
  )
  print(p)
  dev.off()
  #I-Cox
  #多因素
  clin_mut1$Differentiation <- factor(clin_mut1$Differentiation, levels = c("Well", "Mod", "Poor"), ordered = TRUE)
  contrasts(clin_mut1$Differentiation) <- contr.sdif(3)
  #clin_mut1$Age_E <- factor(clin_mut1$Age_E,levels=c("<40","40-49","50-59","60-69",">=70"))
  
  
  clin_mut1$Category <- factor(clin_mut1$Category,levels = c("KRAS WT","KRAS Canonical"))
  clin_mut1$Smoking_E <- factor(clin_mut1$Smoking_E,levels=c("Never","Ever"))
  
  multicox <- coxph(Surv(OS,OS_E)~Category+Sex+Age_E+Smoking_E+Differentiation, data = clin_mut1)
  multisum <- summary(multicox)
  
  
  
  
  #提取所有基因的多因素COX回归分析结果至multiresult对象中
  factor <-  c("KRAS Mut(vs Wt)","Male(vs Female)","Age>=60 (vs <60)","Ever Smoking(vs Never)","Differentiation Mod(vs Well)","Differentiation Poor(vs Mod)")
  HR <- multisum$coefficients[,2]
  L95CI <- multisum$conf.int[,3]
  H95CI <- multisum$conf.int[,4]
  pvalue <- multisum$coefficients[,5]
  proportion <- clin_mut1 %>%
    summarise(
      Mut_n = sum(Category == "KRAS Canonical", na.rm = TRUE),
      Male_n = sum(Sex == "Male", na.rm = TRUE),
      Age_n = sum(Age_E  == ">=60",na.rm = TRUE),
      Smoking_n =sum(Smoking_E== "Ever",na.rm = TRUE),
      Differentiation_mod_n = sum(Differentiation=="Mod",na.rm = TRUE),
      Differentiation_poor_n = sum(Differentiation=="Poor",na.rm = TRUE),
      Total_N = n(),
      Proportion_mut = paste0(Mut_n, "/", Total_N, " (", round(Mut_n / Total_N * 100, 1), "%)"),
      Proportion_male = paste0(Male_n, "/", Total_N, " (", round(Male_n / Total_N * 100, 1), "%)"),
      Proportion_age = paste0(Age_n, "/", Total_N, " (", round(Age_n / Total_N * 100, 1), "%)"),
      Proportion_smoking = paste0(Smoking_n, "/", Total_N, " (", round(Smoking_n / Total_N * 100, 1), "%)"),
      Proportion_mod = paste0(Differentiation_mod_n, "/", Total_N, " (", round(Differentiation_mod_n / Total_N * 100, 1), "%)"),
      Proportion_poor = paste0(Differentiation_poor_n, "/", Total_N, " (", round(Differentiation_poor_n / Total_N * 100, 1), "%)")
    ) 
  multiresult <- data.frame(factor=factor,
                            proportion=as.vector(t(proportion[1,c(8:13)])),
                            HR=HR,
                            L95CI=L95CI,
                            H95CI=H95CI,
                            pvalue=pvalue)
  multiresult$HR <- round(multiresult$HR,digits = 2)
  multiresult$L95CI <- round(multiresult$L95CI,digits = 2)
  multiresult$H95CI <- round(multiresult$H95CI,digits = 2)
  multiresult$pvalue  <-  ifelse(
    multiresult$pvalue < 0.0001,        # 条件：p值是否小于0.0001
    "<0.0001",              # 如果是，则返回此字符串
    sprintf("%.4f",multiresult$pvalu) # 如果否，则格式化为保留四位小数的字符串
  )
  
  multiresult$`HR(95% CI)` <- paste0(multiresult$HR,"(",multiresult$L95CI,"-",multiresult$H95CI,")")
  multiresult <- multiresult[,c(1,2,7,6)]
  print(multiresult)
}

if(0) #ALK
{
  ## ALK
  load("~/molecular_staging/RData/merge_long.Rdata")
  gene <- "ALK"
  # I
  ALK_mut <- merge_long %>%
    mutate(ALK_category=ifelse(Hugo_Symbol=="ALK" & ONCOGENIC %in% c("Oncogenic","Likely Oncogenic"),"ALK Mut","ALK Wt")) %>%
    dplyr::select(Patient_ID,ALK_category) %>%
    distinct() %>%
    group_by(Patient_ID) %>%
    arrange(ALK_category) %>%  # 按"ALK Mut"优先排序（Mut在Wt之前）
    slice(1) %>%  # 保留每组的第一条记录（即优先保留Mut）
    ungroup()  # 取消分组
  ALK_mut$Patient_ID <- as.character(ALK_mut$Patient_ID)
  if(1){
    # I-KM
    clin_mut1 <- clin_mut %>%
      left_join(ALK_mut,by=c("Patient_ID")) %>%
      filter(Stage2=="I") %>%
      mutate(Category=ALK_category)  # 替换TP53_category为ALK_category
    folder_path <- "~/molecular_staging/output"
    filename_os <- paste0(gene,"_OS stage I.pdf")
    filename_os <- file.path(folder_path,filename_os)
    filename_rfs <- paste0(gene,"_RFS stage I.pdf")
    filename_rfs <- file.path(folder_path,filename_rfs)
    # os
    pdf(filename_os,width = 8,height = 8,onefile = F)
    km_os<-survfit(Surv(OS,OS_E)~Category,data = clin_mut1)
    summary_median <- summary(km_os)
    
    
    # 从 summary_median 对象中提取我们需要的信息
    # summary_median$table 是一个包含详细统计量的矩阵
    median_data <- as.data.frame(summary_median$table) %>%
      rownames_to_column("Group") %>% # 将行名（分组名）转为列
      dplyr::select(Group, median, `0.95LCL`, `0.95UCL`) %>% # 选择需要的列
      rename(
        `Median Survival (days)` = median,
        `Median 95% CI Lower` = `0.95LCL`,
        `Median 95% CI Upper` = `0.95UCL`
      )
    
    # 定义5年对应的天数
    five_year_days <- 5 * 12
    
    # 使用 summary 函数在指定时间点（5年）提取生存率和95%CI
    # 这是最推荐的方法，因为它会进行线性插值
    summary_5yr <- summary(km_os, times = five_year_days)
    
    # 从 summary 结果中构建最终的整洁表格
    # summary_5yr$surv, lower, upper 是按分组顺序排列的向量
    result_table <- tibble(
      Group = c("ALK Mut","ALK Wt"),  # 替换TP53 Mut/TP53 Wt为ALK Mut/ALK Wt
      `5-year Survival (%)` = round(summary_5yr$surv * 100, 1),
      `95% CI` = paste0(
        round(summary_5yr$lower * 100, 1),
        " - ",
        round(summary_5yr$upper * 100, 1)
      )
    )
    info_table <- cbind(median_data,result_table[,c(2,3)])
    info_table[is.na(info_table)] <- "."
    info_table <- info_table %>%
      mutate(`Median in Months`= paste0(`Median Survival (days)`,"(",`Median 95% CI Lower`,"-",`Median 95% CI Upper`,")")) %>%
      mutate(`5-Year Estimate`= paste0(`5-year Survival (%)`,"(",`95% CI`,")"))%>%
      dplyr::select(Group,`Median in Months`,`5-Year Estimate`)
    
    p <- ggsurvplot(km_os, 
                    conf.int=F, #是否显示生存率的95%CI
                    risk.table=T,#显示风险表
                    palette="jco", #柳叶刀配色
                    title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage I"), #大标题
                    #risk.table.height = 0.15,#风险表的高度比例
                    tables.height = 0.17,
                    pval=T,  
                    pval.method=T,
                    surv.median.line="v",
                    cumevents=T,
                    cumcensor = T,
                    ylab="OS(%)",
                    xlab="Months"
    )
    print(p)
    dev.off()
    
    # rfs
    pdf(filename_rfs,width = 8,height = 8,onefile = F)
    km_rfs<-survfit(Surv(RFS,RFS_E)~Category,data = clin_mut1)
    p <- ggsurvplot(km_rfs, 
                    conf.int=F, #是否显示生存率的95%CI
                    risk.table=T,#显示风险表
                    palette="jco", #柳叶刀配色
                    title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage I"), #大标题
                    #risk.table.height = 0.15,#风险表的高度比例
                    tables.height = 0.17,
                    pval=T,  
                    pval.method=T,
                    surv.median.line="v",
                    cumevents=T,
                    cumcensor = T,
                    ylab="RFS(%)",
                    xlab="Months"
    )
    print(p)
    dev.off()
    # I-Cox
    # 多因素
    clin_mut1$Differentiation <- factor(clin_mut1$Differentiation, levels = c("Well", "Mod", "Poor"), ordered = TRUE)
    contrasts(clin_mut1$Differentiation) <- contr.sdif(3)
    # clin_mut1$Age_E <- factor(clin_mut1$Age_E,levels=c("<40","40-49","50-59","60-69",">=70"))
    
    
    clin_mut1$Category <- factor(clin_mut1$Category,levels = c("ALK Wt","ALK Mut"))  # 替换TP53 Wt/TP53 Mut为ALK Wt/ALK Mut
    clin_mut1$Smoking_E <- factor(clin_mut1$Smoking_E,levels=c("Never","Ever"))
    
    multicox <- coxph(Surv(OS,OS_E)~Category+Sex+Age_E+Smoking_E+Differentiation, data = clin_mut1)
    multisum <- summary(multicox)
    
    
    
    
    # 提取所有基因的多因素COX回归分析结果至multiresult对象中
    factor <-  c("ALK Mut(vs Wt)","Male(vs Female)","Age>=60 (vs <60)","Ever Smoking(vs Never)","Differentiation Mod(vs Well)","Differentiation Poor(vs Mod)")  # 替换TP53 Mut(vs Wt)为ALK Mut(vs Wt)
    HR <- multisum$coefficients[,2]
    L95CI <- multisum$conf.int[,3]
    H95CI <- multisum$conf.int[,4]
    pvalue <- multisum$coefficients[,5]
    proportion <- clin_mut1 %>%
      summarise(
        Mut_n = sum(Category == "ALK Mut", na.rm = TRUE),  # 替换TP53 Mut为ALK Mut
        Male_n = sum(Sex == "Male", na.rm = TRUE),
        Age_n = sum(Age_E  == ">=60",na.rm = TRUE),
        Smoking_n =sum(Smoking_E== "Ever",na.rm = TRUE),
        Differentiation_mod_n = sum(Differentiation=="Mod",na.rm = TRUE),
        Differentiation_poor_n = sum(Differentiation=="Poor",na.rm = TRUE),
        Total_N = n(),
        Proportion_mut = paste0(Mut_n, "/", Total_N, " (", round(Mut_n / Total_N * 100, 1), "%)"),
        Proportion_male = paste0(Male_n, "/", Total_N, " (", round(Male_n / Total_N * 100, 1), "%)"),
        Proportion_age = paste0(Age_n, "/", Total_N, " (", round(Age_n / Total_N * 100, 1), "%)"),
        Proportion_smoking = paste0(Smoking_n, "/", Total_N, " (", round(Smoking_n / Total_N * 100, 1), "%)"),
        Proportion_mod = paste0(Differentiation_mod_n, "/", Total_N, " (", round(Differentiation_mod_n / Total_N * 100, 1), "%)"),
        Proportion_poor = paste0(Differentiation_poor_n, "/", Total_N, " (", round(Differentiation_poor_n / Total_N * 100, 1), "%)")
      ) 
    multiresult <- data.frame(factor=factor,
                              proportion=as.vector(t(proportion[1,c(8:13)])),
                              HR=HR,
                              L95CI=L95CI,
                              H95CI=H95CI,
                              pvalue=pvalue)
    multiresult$HR <- round(multiresult$HR,digits = 2)
    multiresult$L95CI <- round(multiresult$L95CI,digits = 2)
    multiresult$H95CI <- round(multiresult$H95CI,digits = 2)
    multiresult$pvalue  <-  ifelse(
      multiresult$pvalue < 0.0001,        # 条件：p值是否小于0.0001
      "<0.0001",              # 如果是，则返回此字符串
      sprintf("%.4f",multiresult$pvalue) # 修正原代码拼写错误：multiresult$pvalu改为multiresult$pvalue
    )
    
    multiresult$`HR(95% CI)` <- paste0(multiresult$HR,"(",multiresult$L95CI,"-",multiresult$H95CI,")")
    multiresult <- multiresult[,c(1,2,7,6)]
    print(multiresult)
  }
  # II
  if(1){
    # II-KM
    clin_mut1 <- clin_mut %>%
      left_join(ALK_mut,by=c("Patient_ID")) %>%
      filter(Stage2=="II") %>%
      mutate(Category=ALK_category)  # 替换TP53_category为ALK_category
    folder_path <- "~/molecular_staging/output"
    filename_os <- paste0(gene,"_OS stage II.pdf")
    filename_os <- file.path(folder_path,filename_os)
    filename_rfs <- paste0(gene,"_RFS stage II.pdf")
    filename_rfs <- file.path(folder_path,filename_rfs)
    # os
    pdf(filename_os,width = 8,height = 8,onefile = F)
    km_os<-survfit(Surv(OS,OS_E)~Category,data = clin_mut1)
    summary_median <- summary(km_os)
    
    # 从 summary_median 对象中提取我们需要的信息
    # summary_median$table 是一个包含详细统计量的矩阵
    median_data <- as.data.frame(summary_median$table) %>%
      rownames_to_column("Group") %>% # 将行名（分组名）转为列
      dplyr::select(Group, median, `0.95LCL`, `0.95UCL`) %>% # 选择需要的列
      rename(
        `Median Survival (days)` = median,
        `Median 95% CI Lower` = `0.95LCL`,
        `Median 95% CI Upper` = `0.95UCL`
      )
    
    # 定义5年对应的天数
    five_year_days <- 5 * 12
    
    # 使用 summary 函数在指定时间点（5年）提取生存率和95%CI
    # 这是最推荐的方法，因为它会进行线性插值
    summary_5yr <- summary(km_os, times = five_year_days)
    
    # 从 summary 结果中构建最终的整洁表格
    # summary_5yr$surv, lower, upper 是按分组顺序排列的向量
    # 修正原代码错误：删除未定义的median_surv变量
    result_table <- tibble(
      Group = c("ALK Mut","ALK Wt"),  # 替换TP53 Mut/TP53 Wt为ALK Mut/ALK Wt
      `5-year Survival (%)` = round(summary_5yr$surv * 100, 1),
      `95% CI` = paste0(
        round(summary_5yr$lower * 100, 1),
        " - ",
        round(summary_5yr$upper * 100, 1)
      )
    )
    info_table <- cbind(median_data,result_table[,c(2,3)])
    info_table[is.na(info_table)] <- "."
    info_table <- info_table %>%
      mutate(`Median in Months`= paste0(`Median Survival (days)`,"(",`Median 95% CI Lower`,"-",`Median 95% CI Upper`,")")) %>%
      mutate(`5-Year Estimate`= paste0(`5-year Survival (%)`,"(",`95% CI`,")"))%>%
      dplyr::select(Group,`Median in Months`,`5-Year Estimate`)
    p <- ggsurvplot(km_os, 
                    conf.int=F, #是否显示生存率的95%CI
                    risk.table=T,#显示风险表
                    palette="jco", #柳叶刀配色
                    title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage II"), #大标题
                    #risk.table.height = 0.15,#风险表的高度比例
                    tables.height = 0.17,
                    pval=T,  
                    pval.method=T,
                    surv.median.line="v",
                    cumevents=T,
                    cumcensor = T,
                    ylab="OS(%)",
                    xlab="Months"
    )
    print(p)
    dev.off()
    
    # rfs
    pdf(filename_rfs,width = 8,height = 8,onefile = F)
    km_rfs<-survfit(Surv(RFS,RFS_E)~Category,data = clin_mut1)
    p <- ggsurvplot(km_rfs, 
                    conf.int=F, #是否显示生存率的95%CI
                    risk.table=T,#显示风险表
                    palette="jco", #柳叶刀配色
                    title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage II"), #大标题
                    #risk.table.height = 0.15,#风险表的高度比例
                    tables.height = 0.17,
                    pval=T,  
                    pval.method=T,
                    surv.median.line="v",
                    cumevents=T,
                    cumcensor = T,
                    ylab="RFS(%)",
                    xlab="Months"
    )
    print(p)
    dev.off()
    # II-Cox
    # 多因素
    clin_mut1$Differentiation <- factor(clin_mut1$Differentiation, levels = c("Well", "Mod", "Poor"), ordered = TRUE)
    contrasts(clin_mut1$Differentiation) <- contr.sdif(3)
    # clin_mut1$Age_E <- factor(clin_mut1$Age_E,levels=c("<40","40-49","50-59","60-69",">=70"))
    
    
    clin_mut1$Category <- factor(clin_mut1$Category,levels = c("ALK Wt","ALK Mut"))  # 替换TP53 Wt/TP53 Mut为ALK Wt/ALK Mut
    clin_mut1$Smoking_E <- factor(clin_mut1$Smoking_E,levels=c("Never","Ever"))
    
    multicox <- coxph(Surv(OS,OS_E)~Category+Sex+Age_E+Smoking_E+Differentiation, data = clin_mut1)
    multisum <- summary(multicox)
    
    
    
    
    # 提取所有基因的多因素COX回归分析结果至multiresult对象中
    factor <-  c("ALK Mut(vs Wt)","Male(vs Female)","Age>=60 (vs <60)","Ever Smoking(vs Never)","Differentiation Mod(vs Well)","Differentiation Poor(vs Mod)")  # 替换TP53 Mut(vs Wt)为ALK Mut(vs Wt)
    HR <- multisum$coefficients[,2]
    L95CI <- multisum$conf.int[,3]
    H95CI <- multisum$conf.int[,4]
    pvalue <- multisum$coefficients[,5]
    proportion <- clin_mut1 %>%
      summarise(
        Mut_n = sum(Category == "ALK Mut", na.rm = TRUE),  # 替换TP53 Mut为ALK Mut
        Male_n = sum(Sex == "Male", na.rm = TRUE),
        Age_n = sum(Age_E  == ">=60",na.rm = TRUE),
        Smoking_n =sum(Smoking_E== "Ever",na.rm = TRUE),
        Differentiation_mod_n = sum(Differentiation=="Mod",na.rm = TRUE),
        Differentiation_poor_n = sum(Differentiation=="Poor",na.rm = TRUE),
        Total_N = n(),
        Proportion_mut = paste0(Mut_n, "/", Total_N, " (", round(Mut_n / Total_N * 100, 1), "%)"),
        Proportion_male = paste0(Male_n, "/", Total_N, " (", round(Male_n / Total_N * 100, 1), "%)"),
        Proportion_age = paste0(Age_n, "/", Total_N, " (", round(Age_n / Total_N * 100, 1), "%)"),
        Proportion_smoking = paste0(Smoking_n, "/", Total_N, " (", round(Smoking_n / Total_N * 100, 1), "%)"),
        Proportion_mod = paste0(Differentiation_mod_n, "/", Total_N, " (", round(Differentiation_mod_n / Total_N * 100, 1), "%)"),
        Proportion_poor = paste0(Differentiation_poor_n, "/", Total_N, " (", round(Differentiation_poor_n / Total_N * 100, 1), "%)")
      ) 
    multiresult <- data.frame(factor=factor,
                              proportion=as.vector(t(proportion[1,c(8:13)])),
                              HR=HR,
                              L95CI=L95CI,
                              H95CI=H95CI,
                              pvalue=pvalue)
    multiresult$HR <- round(multiresult$HR,digits = 2)
    multiresult$L95CI <- round(multiresult$L95CI,digits = 2)
    multiresult$H95CI <- round(multiresult$H95CI,digits = 2)
    multiresult$pvalue  <-  ifelse(
      multiresult$pvalue < 0.0001,        # 条件：p值是否小于0.0001
      "<0.0001",              # 如果是，则返回此字符串
      sprintf("%.4f",multiresult$pvalue) # 修正原代码拼写错误：multiresult$pvalu改为multiresult$pvalue
    )
    
    multiresult$`HR(95% CI)` <- paste0(multiresult$HR,"(",multiresult$L95CI,"-",multiresult$H95CI,")")
    multiresult <- multiresult[,c(1,2,7,6)]
    print(multiresult)
  }
  # III&IV
  if(1){
    # III-KM
    clin_mut1 <- clin_mut %>%
      left_join(ALK_mut,by=c("Patient_ID")) %>%
      filter(Stage2 %in% c("III","IV")) %>%
      mutate(Category=ALK_category)
    folder_path <- "~/molecular_staging/output"
    filename_os <- paste0(gene,"_OS stage III&IV.pdf")
    filename_os <- file.path(folder_path,filename_os)
    filename_rfs <- paste0(gene,"_RFS stage III&IV.pdf")
    filename_rfs <- file.path(folder_path,filename_rfs)
    # os
    pdf(filename_os,width = 8,height = 8,onefile = F)
    km_os<-survfit(Surv(OS,OS_E)~Category,data = clin_mut1)
    summary_median <- summary(km_os)
    
    # 从 summary_median 对象中提取我们需要的信息
    # summary_median$table 是一个包含详细统计量的矩阵
    median_data <- as.data.frame(summary_median$table) %>%
      rownames_to_column("Group") %>% # 将行名（分组名）转为列
      dplyr::select(Group, median, `0.95LCL`, `0.95UCL`) %>% # 选择需要的列
      rename(
        `Median Survival (days)` = median,
        `Median 95% CI Lower` = `0.95LCL`,
        `Median 95% CI Upper` = `0.95UCL`
      )
    
    # 定义5年对应的天数
    five_year_days <- 5 * 12
    
    # 使用 summary 函数在指定时间点（5年）提取生存率和95%CI
    # 这是最推荐的方法，因为它会进行线性插值
    summary_5yr <- summary(km_os, times = five_year_days)
    
    # 从 summary 结果中构建最终的整洁表格
    # summary_5yr$surv, lower, upper 是按分组顺序排列的向量
    # 修正原代码错误：删除未定义的median_surv变量
    result_table <- tibble(
      Group = c("ALK Mut","ALK Wt"),  # 替换TP53 Mut/TP53 Wt为ALK Mut/ALK Wt
      `5-year Survival (%)` = round(summary_5yr$surv * 100, 1),
      `95% CI` = paste0(
        round(summary_5yr$lower * 100, 1),
        " - ",
        round(summary_5yr$upper * 100, 1)
      )
    )
    info_table <- cbind(median_data,result_table[,c(2,3)])
    info_table[is.na(info_table)] <- "."
    info_table <- info_table %>%
      mutate(`Median in Months`= paste0(`Median Survival (days)`,"(",`Median 95% CI Lower`,"-",`Median 95% CI Upper`,")")) %>%
      mutate(`5-Year Estimate`= paste0(`5-year Survival (%)`,"(",`95% CI`,")"))%>%
      dplyr::select(Group,`Median in Months`,`5-Year Estimate`)
    
    
    
    
    p <- ggsurvplot(km_os, 
                    conf.int=F, #是否显示生存率的95%CI
                    risk.table=T,#显示风险表
                    palette="jco", #柳叶刀配色
                    title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage III&IV"), #大标题
                    #risk.table.height = 0.15,#风险表的高度比例
                    tables.height = 0.17,
                    pval=T,  
                    pval.method=T,
                    surv.median.line="v",
                    cumevents=T,
                    cumcensor = T,
                    ylab="OS(%)",
                    xlab="Months"
    )
    print(p)
    dev.off()
    
    # rfs
    pdf(filename_rfs,width = 8,height = 8,onefile = F)
    km_rfs<-survfit(Surv(RFS,RFS_E)~Category,data = clin_mut1)
    p <- ggsurvplot(km_rfs, 
                    conf.int=F, #是否显示生存率的95%CI
                    risk.table=T,#显示风险表
                    palette="jco", #柳叶刀配色
                    title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage III&IV"), #大标题
                    #risk.table.height = 0.15,#风险表的高度比例
                    tables.height = 0.17,
                    pval=T,  
                    pval.method=T,
                    surv.median.line = "v",
                    cumevents=T,
                    cumcensor = T,
                    ylab="RFS(%)",
                    xlab="Months"
    )
    print(p)
    dev.off()
    # III&IV-Cox
    # 多因素
    clin_mut1$Differentiation <- factor(clin_mut1$Differentiation, levels = c("Well", "Mod", "Poor"), ordered = TRUE)
    contrasts(clin_mut1$Differentiation) <- contr.sdif(3)
    # clin_mut1$Age_E <- factor(clin_mut1$Age_E,levels=c("<40","40-49","50-59","60-69",">=70"))
    
    
    clin_mut1$Category <- factor(clin_mut1$Category,levels = c("ALK Wt","ALK Mut"))  # 替换TP53 Wt/TP53 Mut为ALK Wt/ALK Mut
    clin_mut1$Smoking_E <- factor(clin_mut1$Smoking_E,levels=c("Never","Ever"))
    
    multicox <- coxph(Surv(OS,OS_E)~Category+Sex+Age_E+Smoking_E+Differentiation, data = clin_mut1)
    multisum <- summary(multicox)
    
    
    
    
    # 提取所有基因的多因素COX回归分析结果至multiresult对象中
    factor <-  c("ALK Mut(vs Wt)","Male(vs Female)","Age>=60 (vs <60)","Ever Smoking(vs Never)","Differentiation Mod(vs Well)","Differentiation Poor(vs Mod)")  # 替换TP53 Mut(vs Wt)为ALK Mut(vs Wt)
    HR <- multisum$coefficients[,2]
    L95CI <- multisum$conf.int[,3]
    H95CI <- multisum$conf.int[,4]
    pvalue <- multisum$coefficients[,5]
    proportion <- clin_mut1 %>%
      summarise(
        Mut_n = sum(Category == "ALK Mut", na.rm = TRUE),  # 替换TP53 Mut为ALK Mut
        Male_n = sum(Sex == "Male", na.rm = TRUE),
        Age_n = sum(Age_E  == ">=60",na.rm = TRUE),
        Smoking_n =sum(Smoking_E== "Ever",na.rm = TRUE),
        Differentiation_mod_n = sum(Differentiation=="Mod",na.rm = TRUE),
        Differentiation_poor_n = sum(Differentiation=="Poor",na.rm = TRUE),
        Total_N = n(),
        Proportion_mut = paste0(Mut_n, "/", Total_N, " (", round(Mut_n / Total_N * 100, 1), "%)"),
        Proportion_male = paste0(Male_n, "/", Total_N, " (", round(Male_n / Total_N * 100, 1), "%)"),
        Proportion_age = paste0(Age_n, "/", Total_N, " (", round(Age_n / Total_N * 100, 1), "%)"),
        Proportion_smoking = paste0(Smoking_n, "/", Total_N, " (", round(Smoking_n / Total_N * 100, 1), "%)"),
        Proportion_mod = paste0(Differentiation_mod_n, "/", Total_N, " (", round(Differentiation_mod_n / Total_N * 100, 1), "%)"),
        Proportion_poor = paste0(Differentiation_poor_n, "/", Total_N, " (", round(Differentiation_poor_n / Total_N * 100, 1), "%)")
      ) 
    multiresult <- data.frame(factor=factor,
                              proportion=as.vector(t(proportion[1,c(8:13)])),
                              HR=HR,
                              L95CI=L95CI,
                              H95CI=H95CI,
                              pvalue=pvalue)
    multiresult$HR <- round(multiresult$HR,digits = 2)
    multiresult$L95CI <- round(multiresult$L95CI,digits = 2)
    multiresult$H95CI <- round(multiresult$H95CI,digits = 2)
    multiresult$pvalue  <-  ifelse(
      multiresult$pvalue < 0.0001,        # 条件：p值是否小于0.0001
      "<0.0001",              # 如果是，则返回此字符串
      sprintf("%.4f",multiresult$pvalue) # 修正原代码拼写错误：multiresult$pvalu改为multiresult$pvalue
    )
    
    multiresult$`HR(95% CI)` <- paste0(multiresult$HR,"(",multiresult$L95CI,"-",multiresult$H95CI,")")
    multiresult <- multiresult[,c(1,2,7,6)]
    print(multiresult)
  }
}
if(0){ #TP53/KRAS
  ## TP53_KRAS
  gene <- "TK"
  # I
  if(1){
    # I-KM
    clin_mut1 <- clin_mut %>%
      filter(Stage2=="I") %>% 
      # 恢复TP53/KRAS分组逻辑：TP53突变或KRAS典型突变归为Mut，否则为Wt
      mutate(Category=ifelse(TP53_category=="TP53 Mut" | KRAS_category=="KRAS Canonical","TP53_KRAS Mut","TP53_KRAS Wt"))  
    folder_path <- "~/molecular_staging/output"
    filename_os <- paste0(gene,"_OS stage I.pdf")
    filename_os <- file.path(folder_path,filename_os)
    filename_rfs <- paste0(gene,"_RFS stage I.pdf")
    filename_rfs <- file.path(folder_path,filename_rfs)
    # os
    pdf(filename_os,width = 8,height = 8,onefile = F)
    km_os<-survfit(Surv(OS,OS_E)~Category,data = clin_mut1)
    summary_median <- summary(km_os)
    
    
    # 从 summary_median 对象中提取我们需要的信息
    # summary_median$table 是一个包含详细统计量的矩阵
    median_data <- as.data.frame(summary_median$table) %>%
      rownames_to_column("Group") %>% # 将行名（分组名）转为列
      dplyr::select(Group, median, `0.95LCL`, `0.95UCL`) %>% # 选择需要的列
      rename(
        `Median Survival (days)` = median,
        `Median 95% CI Lower` = `0.95LCL`,
        `Median 95% CI Upper` = `0.95UCL`
      )
    
    # 定义5年对应的天数（原代码逻辑：5年=5*12个月，此处单位与x轴“Months”一致，无需修改）
    five_year_days <- 5 * 12
    
    # 使用 summary 函数在指定时间点（5年）提取生存率和95%CI
    summary_5yr <- summary(km_os, times = five_year_days)
    
    # 从 summary 结果中构建最终的整洁表格
    result_table <- tibble(
      Group = c("TP53_KRAS Mut","TP53_KRAS Wt"),  # 修正为TP53_KRAS
      `5-year Survival (%)` = round(summary_5yr$surv * 100, 1),
      `95% CI` = paste0(
        round(summary_5yr$lower * 100, 1),
        " - ",
        round(summary_5yr$upper * 100, 1)
      )
    )
    info_table <- cbind(median_data,result_table[,c(2,3)])
    info_table[is.na(info_table)] <- "."
    info_table <- info_table %>%
      mutate(`Median in Months`= paste0(`Median Survival (days)`,"(",`Median 95% CI Lower`,"-",`Median 95% CI Upper`,")")) %>%
      mutate(`5-Year Estimate`= paste0(`5-year Survival (%)`,"(",`95% CI`,")"))%>%
      dplyr::select(Group,`Median in Months`,`5-Year Estimate`)
    
    p <- ggsurvplot(km_os, 
                    conf.int=F, # 是否显示生存率的95%CI
                    risk.table=T,# 显示风险表
                    palette="jco", # 柳叶刀配色
                    title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage I"), 
                    tables.height = 0.17,
                    pval=T,  
                    pval.method=T,
                    surv.median.line="v",
                    cumevents=T,
                    cumcensor = T,
                    ylab="OS(%)",
                    xlab="Months"
    )
    print(p)
    dev.off()
    
    # rfs
    pdf(filename_rfs,width = 8,height = 8,onefile = F)
    km_rfs<-survfit(Surv(RFS,RFS_E)~Category,data = clin_mut1)
    p <- ggsurvplot(km_rfs, 
                    conf.int=F, 
                    risk.table=T,
                    palette="jco", 
                    title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage I"), 
                    tables.height = 0.17,
                    pval=T,  
                    pval.method=T,
                    surv.median.line="v",
                    cumevents=T,
                    cumcensor = T,
                    ylab="RFS(%)",
                    xlab="Months"
    )
    print(p)
    dev.off()
    # I-Cox
    # 多因素
    clin_mut1$Differentiation <- factor(clin_mut1$Differentiation, levels = c("Well", "Mod", "Poor"), ordered = TRUE)
    contrasts(clin_mut1$Differentiation) <- contr.sdif(3)
    
    clin_mut1$Category <- factor(clin_mut1$Category,levels = c("TP53_KRAS Wt","TP53_KRAS Mut"))  # 修正为TP53_KRAS
    clin_mut1$Smoking_E <- factor(clin_mut1$Smoking_E,levels=c("Never","Ever"))
    
    multicox <- coxph(Surv(OS,OS_E)~Category+Sex+Age_E+Smoking_E+Differentiation, data = clin_mut1)
    multisum <- summary(multicox)
    
    # 提取所有基因的多因素COX回归分析结果至multiresult对象中
    factor <-  c("TP53_KRAS Mut(vs Wt)","Male(vs Female)","Age>=60 (vs <60)","Ever Smoking(vs Never)","Differentiation Mod(vs Well)","Differentiation Poor(vs Mod)")  # 修正为TP53_KRAS
    HR <- multisum$coefficients[,2]
    L95CI <- multisum$conf.int[,3]
    H95CI <- multisum$conf.int[,4]
    pvalue <- multisum$coefficients[,5]
    proportion <- clin_mut1 %>%
      summarise(
        Mut_n = sum(Category == "TP53_KRAS Mut", na.rm = TRUE),  # 修正为TP53_KRAS
        Male_n = sum(Sex == "Male", na.rm = TRUE),
        Age_n = sum(Age_E  == ">=60",na.rm = TRUE),
        Smoking_n =sum(Smoking_E== "Ever",na.rm = TRUE),
        Differentiation_mod_n = sum(Differentiation=="Mod",na.rm = TRUE),
        Differentiation_poor_n = sum(Differentiation=="Poor",na.rm = TRUE),
        Total_N = n(),
        Proportion_mut = paste0(Mut_n, "/", Total_N, " (", round(Mut_n / Total_N * 100, 1), "%)"),
        Proportion_male = paste0(Male_n, "/", Total_N, " (", round(Male_n / Total_N * 100, 1), "%)"),
        Proportion_age = paste0(Age_n, "/", Total_N, " (", round(Age_n / Total_N * 100, 1), "%)"),
        Proportion_smoking = paste0(Smoking_n, "/", Total_N, " (", round(Smoking_n / Total_N * 100, 1), "%)"),
        Proportion_mod = paste0(Differentiation_mod_n, "/", Total_N, " (", round(Differentiation_mod_n / Total_N * 100, 1), "%)"),
        Proportion_poor = paste0(Differentiation_poor_n, "/", Total_N, " (", round(Differentiation_poor_n / Total_N * 100, 1), "%)")
      ) 
    multiresult <- data.frame(factor=factor,
                              proportion=as.vector(t(proportion[1,c(8:13)])),
                              HR=HR,
                              L95CI=L95CI,
                              H95CI=H95CI,
                              pvalue=pvalue)
    multiresult$HR <- round(multiresult$HR,digits = 2)
    multiresult$L95CI <- round(multiresult$L95CI,digits = 2)
    multiresult$H95CI <- round(multiresult$H95CI,digits = 2)
    multiresult$pvalue  <-  ifelse(
      multiresult$pvalue < 0.0001,        
      "<0.0001",              
      sprintf("%.4f",multiresult$pvalue)
    )
    
    multiresult$`HR(95% CI)` <- paste0(multiresult$HR,"(",multiresult$L95CI,"-",multiresult$H95CI,")")
    multiresult <- multiresult[,c(1,2,7,6)]
    print(multiresult)
  }
  # II
  if(1){
    # II-KM
    clin_mut1 <- clin_mut %>%
      filter(Stage2=="II") %>%
      # 恢复TP53/KRAS分组逻辑：TP53突变或KRAS典型突变归为Mut，否则为Wt
      mutate(Category=ifelse(TP53_category=="TP53 Mut" | KRAS_category=="KRAS Canonical","TP53_KRAS Mut","TP53_KRAS Wt"))  
    folder_path <- "~/molecular_staging/output"
    filename_os <- paste0(gene,"_OS stage II.pdf")
    filename_os <- file.path(folder_path,filename_os)
    filename_rfs <- paste0(gene,"_RFS stage II.pdf")
    filename_rfs <- file.path(folder_path,filename_rfs)
    # os
    pdf(filename_os,width = 8,height = 8,onefile = F)
    km_os<-survfit(Surv(OS,OS_E)~Category,data = clin_mut1)
    summary_median <- summary(km_os)
    
    # 从 summary_median 对象中提取我们需要的信息
    median_data <- as.data.frame(summary_median$table) %>%
      rownames_to_column("Group") %>% 
      dplyr::select(Group, median, `0.95LCL`, `0.95UCL`) %>% 
      rename(
        `Median Survival (days)` = median,
        `Median 95% CI Lower` = `0.95LCL`,
        `Median 95% CI Upper` = `0.95UCL`
      )
    
    # 定义5年对应的天数
    five_year_days <- 5 * 12
    
    # 使用 summary 函数在指定时间点（5年）提取生存率和95%CI
    summary_5yr <- summary(km_os, times = five_year_days)
    
    # 从 summary 结果中构建最终的整洁表格
    result_table <- tibble(
      Group = c("TP53_KRAS Mut","TP53_KRAS Wt"),  # 修正为TP53_KRAS
      `5-year Survival (%)` = round(summary_5yr$surv * 100, 1),
      `95% CI` = paste0(
        round(summary_5yr$lower * 100, 1),
        " - ",
        round(summary_5yr$upper * 100, 1)
      )
    )
    info_table <- cbind(median_data,result_table[,c(2,3)])
    info_table[is.na(info_table)] <- "."
    info_table <- info_table %>%
      mutate(`Median in Months`= paste0(`Median Survival (days)`,"(",`Median 95% CI Lower`,"-",`Median 95% CI Upper`,")")) %>%
      mutate(`5-Year Estimate`= paste0(`5-year Survival (%)`,"(",`95% CI`,")"))%>%
      dplyr::select(Group,`Median in Months`,`5-Year Estimate`)
    p <- ggsurvplot(km_os, 
                    conf.int=F, 
                    risk.table=T,
                    palette="jco", 
                    title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage II"), 
                    tables.height = 0.17,
                    pval=T,  
                    pval.method=T,
                    surv.median.line="v",
                    cumevents=T,
                    cumcensor = T,
                    ylab="OS(%)",
                    xlab="Months"
    )
    print(p)
    dev.off()
    
    # rfs
    pdf(filename_rfs,width = 8,height = 8,onefile = F)
    km_rfs<-survfit(Surv(RFS,RFS_E)~Category,data = clin_mut1)
    p <- ggsurvplot(km_rfs, 
                    conf.int=F, 
                    risk.table=T,
                    palette="jco", 
                    title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage II"), 
                    tables.height = 0.17,
                    pval=T,  
                    pval.method=T,
                    surv.median.line="v",
                    cumevents=T,
                    cumcensor = T,
                    ylab="RFS(%)",
                    xlab="Months"
    )
    print(p)
    dev.off()
    # II-Cox
    # 多因素
    clin_mut1$Differentiation <- factor(clin_mut1$Differentiation, levels = c("Well", "Mod", "Poor"), ordered = TRUE)
    contrasts(clin_mut1$Differentiation) <- contr.sdif(3)
    
    clin_mut1$Category <- factor(clin_mut1$Category,levels = c("TP53_KRAS Wt","TP53_KRAS Mut"))  # 修正为TP53_KRAS
    clin_mut1$Smoking_E <- factor(clin_mut1$Smoking_E,levels=c("Never","Ever"))
    
    multicox <- coxph(Surv(OS,OS_E)~Category+Sex+Age_E+Smoking_E+Differentiation, data = clin_mut1)
    multisum <- summary(multicox)
    
    # 提取所有基因的多因素COX回归分析结果至multiresult对象中
    factor <-  c("TP53_KRAS Mut(vs Wt)","Male(vs Female)","Age>=60 (vs <60)","Ever Smoking(vs Never)","Differentiation Mod(vs Well)","Differentiation Poor(vs Mod)")  # 修正为TP53_KRAS
    HR <- multisum$coefficients[,2]
    L95CI <- multisum$conf.int[,3]
    H95CI <- multisum$conf.int[,4]
    pvalue <- multisum$coefficients[,5]
    proportion <- clin_mut1 %>%
      summarise(
        Mut_n = sum(Category == "TP53_KRAS Mut", na.rm = TRUE),  # 修正为TP53_KRAS
        Male_n = sum(Sex == "Male", na.rm = TRUE),
        Age_n = sum(Age_E  == ">=60",na.rm = TRUE),
        Smoking_n =sum(Smoking_E== "Ever",na.rm = TRUE),
        Differentiation_mod_n = sum(Differentiation=="Mod",na.rm = TRUE),
        Differentiation_poor_n = sum(Differentiation=="Poor",na.rm = TRUE),
        Total_N = n(),
        Proportion_mut = paste0(Mut_n, "/", Total_N, " (", round(Mut_n / Total_N * 100, 1), "%)"),
        Proportion_male = paste0(Male_n, "/", Total_N, " (", round(Male_n / Total_N * 100, 1), "%)"),
        Proportion_age = paste0(Age_n, "/", Total_N, " (", round(Age_n / Total_N * 100, 1), "%)"),
        Proportion_smoking = paste0(Smoking_n, "/", Total_N, " (", round(Smoking_n / Total_N * 100, 1), "%)"),
        Proportion_mod = paste0(Differentiation_mod_n, "/", Total_N, " (", round(Differentiation_mod_n / Total_N * 100, 1), "%)"),
        Proportion_poor = paste0(Differentiation_poor_n, "/", Total_N, " (", round(Differentiation_poor_n / Total_N * 100, 1), "%)")
      ) 
    multiresult <- data.frame(factor=factor,
                              proportion=as.vector(t(proportion[1,c(8:13)])),
                              HR=HR,
                              L95CI=L95CI,
                              H95CI=H95CI,
                              pvalue=pvalue)
    multiresult$HR <- round(multiresult$HR,digits = 2)
    multiresult$L95CI <- round(multiresult$L95CI,digits = 2)
    multiresult$H95CI <- round(multiresult$H95CI,digits = 2)
    multiresult$pvalue  <-  ifelse(
      multiresult$pvalue < 0.0001,        
      "<0.0001",              
      sprintf("%.4f",multiresult$pvalue)
    )
    
    multiresult$`HR(95% CI)` <- paste0(multiresult$HR,"(",multiresult$L95CI,"-",multiresult$H95CI,")")
    multiresult <- multiresult[,c(1,2,7,6)]
    print(multiresult)
  }
  # III&IV
  if(1){
    # III-KM
    clin_mut1 <- clin_mut %>%
      filter(Stage2 %in% c("III","IV")) %>%
      # 恢复TP53/KRAS分组逻辑：TP53突变或KRAS典型突变归为Mut，否则为Wt
      mutate(Category=ifelse(TP53_category=="TP53 Mut" | KRAS_category=="KRAS Canonical","TP53_KRAS Mut","TP53_KRAS Wt"))  
    folder_path <- "~/molecular_staging/output"
    filename_os <- paste0(gene,"_OS stage III&IV.pdf")
    filename_os <- file.path(folder_path,filename_os)
    filename_rfs <- paste0(gene,"_RFS stage III&IV.pdf")
    filename_rfs <- file.path(folder_path,filename_rfs)
    # os
    pdf(filename_os,width = 8,height = 8,onefile = F)
    km_os<-survfit(Surv(OS,OS_E)~Category,data = clin_mut1)
    summary_median <- summary(km_os)
    
    # 从 summary_median 对象中提取我们需要的信息
    median_data <- as.data.frame(summary_median$table) %>%
      rownames_to_column("Group") %>% 
      dplyr::select(Group, median, `0.95LCL`, `0.95UCL`) %>% 
      rename(
        `Median Survival (days)` = median,
        `Median 95% CI Lower` = `0.95LCL`,
        `Median 95% CI Upper` = `0.95UCL`
      )
    
    # 定义5年对应的天数
    five_year_days <- 5 * 12
    
    # 使用 summary 函数在指定时间点（5年）提取生存率和95%CI
    summary_5yr <- summary(km_os, times = five_year_days)
    
    # 从 summary 结果中构建最终的整洁表格
    result_table <- tibble(
      Group = c("TP53_KRAS Mut","TP53_KRAS Wt"),  # 修正为TP53_KRAS
      `5-year Survival (%)` = round(summary_5yr$surv * 100, 1),
      `95% CI` = paste0(
        round(summary_5yr$lower * 100, 1),
        " - ",
        round(summary_5yr$upper * 100, 1)
      )
    )
    info_table <- cbind(median_data,result_table[,c(2,3)])
    info_table[is.na(info_table)] <- "."
    info_table <- info_table %>%
      mutate(`Median in Months`= paste0(`Median Survival (days)`,"(",`Median 95% CI Lower`,"-",`Median 95% CI Upper`,")")) %>%
      mutate(`5-Year Estimate`= paste0(`5-year Survival (%)`,"(",`95% CI`,")"))%>%
      dplyr::select(Group,`Median in Months`,`5-Year Estimate`)
    
    p <- ggsurvplot(km_os, 
                    conf.int=F, 
                    risk.table=T,
                    palette="jco", 
                    title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage III&IV"), 
                    tables.height = 0.17,
                    pval=T,  
                    pval.method=T,
                    surv.median.line="v",
                    cumevents=T,
                    cumcensor = T,
                    ylab="OS(%)",
                    xlab="Months"
    )
    print(p)
    dev.off()
    
    # rfs
    pdf(filename_rfs,width = 8,height = 8,onefile = F)
    km_rfs<-survfit(Surv(RFS,RFS_E)~Category,data = clin_mut1)
    p <- ggsurvplot(km_rfs, 
                    conf.int=F, 
                    risk.table=T,
                    palette="jco", 
                    title=paste0(gene," Result IASLC 9th Edition TNM Pathological Stage III&IV"), 
                    tables.height = 0.17,
                    pval=T,  
                    pval.method=T,
                    surv.median.line = "v",
                    cumevents=T,
                    cumcensor = T,
                    ylab="RFS(%)",
                    xlab="Months"
    )
    print(p)
    dev.off()
    # III&IV-Cox
    # 多因素
    clin_mut1$Differentiation <- factor(clin_mut1$Differentiation, levels = c("Well", "Mod", "Poor"), ordered = TRUE)
    contrasts(clin_mut1$Differentiation) <- contr.sdif(3)
    
    clin_mut1$Category <- factor(clin_mut1$Category,levels = c("TP53_KRAS Wt","TP53_KRAS Mut"))  # 修正为TP53_KRAS
    clin_mut1$Smoking_E <- factor(clin_mut1$Smoking_E,levels=c("Never","Ever"))
    
    multicox <- coxph(Surv(OS,OS_E)~Category+Sex+Age_E+Smoking_E+Differentiation, data = clin_mut1)
    multisum <- summary(multicox)
    
    # 提取所有基因的多因素COX回归分析结果至multiresult对象中
    factor <-  c("TP53_KRAS Mut(vs Wt)","Male(vs Female)","Age>=60 (vs <60)","Ever Smoking(vs Never)","Differentiation Mod(vs Well)","Differentiation Poor(vs Mod)")  # 修正为TP53_KRAS
    HR <- multisum$coefficients[,2]
    L95CI <- multisum$conf.int[,3]
    H95CI <- multisum$conf.int[,4]
    pvalue <- multisum$coefficients[,5]
    proportion <- clin_mut1 %>%
      summarise(
        Mut_n = sum(Category == "TP53_KRAS Mut", na.rm = TRUE),  # 修正为TP53_KRAS
        Male_n = sum(Sex == "Male", na.rm = TRUE),
        Age_n = sum(Age_E  == ">=60",na.rm = TRUE),
        Smoking_n =sum(Smoking_E== "Ever",na.rm = TRUE),
        Differentiation_mod_n = sum(Differentiation=="Mod",na.rm = TRUE),
        Differentiation_poor_n = sum(Differentiation=="Poor",na.rm = TRUE),
        Total_N = n(),
        Proportion_mut = paste0(Mut_n, "/", Total_N, " (", round(Mut_n / Total_N * 100, 1), "%)"),
        Proportion_male = paste0(Male_n, "/", Total_N, " (", round(Male_n / Total_N * 100, 1), "%)"),
        Proportion_age = paste0(Age_n, "/", Total_N, " (", round(Age_n / Total_N * 100, 1), "%)"),
        Proportion_smoking = paste0(Smoking_n, "/", Total_N, " (", round(Smoking_n / Total_N * 100, 1), "%)"),
        Proportion_mod = paste0(Differentiation_mod_n, "/", Total_N, " (", round(Differentiation_mod_n / Total_N * 100, 1), "%)"),
        Proportion_poor = paste0(Differentiation_poor_n, "/", Total_N, " (", round(Differentiation_poor_n / Total_N * 100, 1), "%)")
      ) 
    multiresult <- data.frame(factor=factor,
                              proportion=as.vector(t(proportion[1,c(8:13)])),
                              HR=HR,
                              L95CI=L95CI,
                              H95CI=H95CI,
                              pvalue=pvalue)
    multiresult$HR <- round(multiresult$HR,digits = 2)
    multiresult$L95CI <- round(multiresult$L95CI,digits = 2)
    multiresult$H95CI <- round(multiresult$H95CI,digits = 2)
    multiresult$pvalue  <-  ifelse(
      multiresult$pvalue < 0.0001,        
      "<0.0001",              
      sprintf("%.4f",multiresult$pvalue)
    )
    
    multiresult$`HR(95% CI)` <- paste0(multiresult$HR,"(",multiresult$L95CI,"-",multiresult$H95CI,")")
    multiresult <- multiresult[,c(1,2,7,6)]
    print(multiresult)
  }
  }
