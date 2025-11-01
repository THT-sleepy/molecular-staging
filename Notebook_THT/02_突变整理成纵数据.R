##loading packages
library(tidyr)
library(dplyr)
library(stringr)

##loading files
NGS_raw_table <- read.csv("~/molecular_staging/input/NGS用于R语言输入.csv") 
NGS_mutoncokb <- read.delim("~/molecular_staging/input/NGS.oncokb.txt")
NGS_mutoncokb <- NGS_mutoncokb[!(NGS_mutoncokb$HGVSp %in% c("","?")),]
NGS_fusiononcokb <- read.delim("~/molecular_staging/input/NGS_fusions.oncokb.txt")
#NGS_fusiononcokb <- NGS_fusiononcokb[-c(141),c(1,2,8)]
NGS_fusiononcokb <- NGS_fusiononcokb[,c(1,2,8)]
colnames(NGS_fusiononcokb)[3] <- "oncogenic"
load("~/molecular_staging/RData/clin_merge.Rdata")




##处理数据框

#提取出突变相应的列
NGS_mut <- NGS_raw_table[,c(1:10)]
#检查各列包含的项数一致
# 创建检查逗号数量的函数
check_comma_count <- function(data) {
  # 对每行检查第4到10列的逗号数量是否一致
  data$"各列数量一致否" <- apply(data[, 4:10], 1, function(row) {
    # 计算每列中的逗号数量
    comma_counts <- sapply(row, function(x) {
      if (is.na(x)) {
        return(NA)  # 处理NA值
      } else {
        return(lengths(regmatches(x, gregexpr(",", x))))
      }
    })
    
    # 检查所有逗号数量是否相同(忽略NA)
    if (length(unique(na.omit(comma_counts))) == 1) {
      return("相同")
    } else {
      return("不同")
    }
  })
  
  return(data)
}
NGS_mut <- check_comma_count(NGS_mut)
#拆成纵数据
NGS_mut_split <- NGS_mut[,c(1:10)] %>%
  separate_rows(4:10, sep = ",") %>%
  group_by(across(c(1:3))) %>%  # 假设前3列是ID或其他不需要拆分的列
  ungroup() %>%
  mutate(across(4:10, ~ str_replace_all(., "\\s+", ""))) %>%
  filter(across(4:10, ~ . != ""))%>%
  mutate("融合.基因A.基因B."=NA,"融合位点"=NA,
         "拷贝数变异基因"=NA,"染色体位置"=NA,
         "拷贝数趋势"=NA,"绝对拷贝数"=NA,Gene=突变基因) %>%
  select(1:3,Gene,4:17)

#弄一个用于oncokb注释的mut表格
#NGS_mut_split_foroncokb <- NGS_mut_split %>%
#  mutate(NCBI_Build="GRCh38",Hugo_Symbol=Gene,Patient_ID=住院号) %>%
#  separate(HGVS_P, into = c("prefix", "suffix"), sep = "\\.", remove = FALSE, extra = "merge") %>%
#  mutate(HGVSp = suffix) %>%
#  select(-prefix, -suffix) %>%
#  select(Patient_ID,NCBI_Build,Hugo_Symbol,HGVSp)
#write.table(NGS_mut_split_foroncokb,"~/molecular_staging/2.analysis/oncokb-annotator/NGS.input",
#           row.names = F,quote = F,sep = "\t")



#提取出融合相应的列
NGS_fusion <- NGS_raw_table[,c(1:3,11:12)]
#检查各列包含的项数一致
check_comma_count <- function(data) {
  # 对每行检查第4到10列的逗号数量是否一致
  data$"各列数量一致否" <- apply(data[, 4:5], 1, function(row) {
    # 计算每列中的逗号数量
    comma_counts <- sapply(row, function(x) {
      if (is.na(x)) {
        return(NA)  # 处理NA值
      } else {
        return(lengths(regmatches(x, gregexpr(",", x))))
      }
    })
    
    # 检查所有逗号数量是否相同(忽略NA)
    if (length(unique(na.omit(comma_counts))) == 1) {
      return("相同")
    } else {
      return("不同")
    }
  })
  
  return(data)
}
NGS_fusion <- check_comma_count(NGS_fusion)
#拆成纵数据
NGS_fusion_split <- NGS_fusion[,c(1:5)] %>%
  separate_rows(4:5, sep = ",") %>%
  group_by(across(c(1:3))) %>%  # 假设前3列是ID或其他不需要拆分的列
  ungroup() %>%
  mutate(across(4:5, ~ str_replace_all(., "\\s+", ""))) %>%
  filter(across(4:5, ~ . != "")) %>%
  mutate("突变基因"=NA,"RefSeq转录本ID"=NA,   
         "外显子序号"=NA,"突变类型"="Fusion",
         "HGVS_C"=NA,"HGVS_P"=NA,
         "突变丰度"=NA,"拷贝数变异基因"=NA,
         "染色体位置"=NA,"拷贝数趋势"=NA,       
         "绝对拷贝数"=NA) %>%
  mutate(Gene = case_when(
    grepl("ALK", `融合.基因A.基因B.`) ~ "ALK",
    grepl("RET", `融合.基因A.基因B.`) ~ "RET",
    grepl("ROS1", `融合.基因A.基因B.`) ~ "ROS1",
    grepl("NRG1", `融合.基因A.基因B.`) ~ "NRG1",
    grepl("NTRK1", `融合.基因A.基因B.`) ~ "NTRK1",
    TRUE ~ "Others"
  ))
#弄一个用于oncokb注释的fusion表格
#NGS_fusion_split_foroncokb <- NGS_fusion_split %>%
#  mutate(SAMPLE_ID=住院号,Fusion=融合.基因A.基因B.) %>%
#  select(SAMPLE_ID,Fusion)
#write.table(NGS_fusion_split_foroncokb,"~/molecular_staging/2.analysis/oncokb-annotator/NGS.fusion.input",
#            row.names = F,quote = F,sep = "\t")



#提取出拷贝数变异相应的列
NGS_cnv <- NGS_raw_table[,c(1:3,13:16)]
#检查各列包含的项数一致
check_comma_count <- function(data) {
  # 对每行检查第4到10列的逗号数量是否一致
  data$"各列数量一致否" <- apply(data[, 4:7], 1, function(row) {
    # 计算每列中的逗号数量
    comma_counts <- sapply(row, function(x) {
      if (is.na(x)) {
        return(NA)  # 处理NA值
      } else {
        return(lengths(regmatches(x, gregexpr(",", x))))
      }
    })
    
    # 检查所有逗号数量是否相同(忽略NA)
    if (length(unique(na.omit(comma_counts))) == 1) {
      return("相同")
    } else {
      return("不同")
    }
  })
  
  return(data)
}
NGS_cnv <- check_comma_count(NGS_cnv)
#拆成纵数据
NGS_cnv_split <- NGS_cnv[,c(1:7)] %>%
  separate_rows(4:7, sep = ",") %>%
  group_by(across(c(1:3))) %>%  # 假设前3列是ID或其他不需要拆分的列
  ungroup() %>%
  mutate(across(4:7, ~ str_replace_all(., "\\s+", ""))) %>%
  filter(across(4:7, ~ . != "")) %>%
  mutate("突变基因"=NA,"RefSeq转录本ID"=NA,
         "外显子序号"=NA,"突变类型"=paste0("拷贝数",拷贝数趋势),
         "HGVS_C"=NA,"HGVS_P"=NA,"突变丰度"=NA,         
         "融合.基因A.基因B."=NA,"融合位点"=NA) %>%
  mutate(Gene=拷贝数变异基因)

NGS_long <- rbind(NGS_mut_split,NGS_fusion_split,NGS_cnv_split)
NGS_long <- NGS_long[,-c(5,14,16)]
NGS_long <- NGS_long %>%
  arrange(住院号)

##添加oncokb注释
NGS_long <- NGS_long %>%
  mutate(SAMPLE_ID=住院号,Fusion=融合.基因A.基因B.,Hugo_Symbol=Gene,HGVS_P1=HGVS_P) %>%
  separate(HGVS_P, into = c("prefix", "suffix"), sep = "\\.", remove = FALSE, extra = "merge") %>%
  mutate(HGVSp = suffix) %>%
  select(-prefix, -suffix) %>%
  left_join(NGS_mutoncokb,by=c("HGVSp","Hugo_Symbol","SAMPLE_ID")) %>%
  left_join(NGS_fusiononcokb,by=c("SAMPLE_ID","Fusion"))
NGS_long <- NGS_long %>%
  mutate(oncokb_oncogenic=ifelse(!is.na(NGS_long$ONCOGENIC),NGS_long$ONCOGENIC,
                                 ifelse(!is.na(NGS_long$oncogenic),NGS_long$oncogenic,NA))) %>%
  select(1:14,oncokb_oncogenic)

##修正MET exon14的ONCOGENIC
NGS_long$oncokb_oncogenic[NGS_long$Gene=="MET" & NGS_long$住院号 %in% 
                            c(397775,494165,540426,
                              543054,543140,563803,
                              586239,589709,593391,
                              606046,610431,618765,
                              623141,627633,637452,
                              648461,648529,649627,653103)] <- "Likely Oncogenic" 

save(NGS_long,file="~/molecular_staging/RData/NGS_long.Rdata")

##loading files
WGS_maf <- read.delim("~/molecular_staging/input/new_mut_fus_combined.txt")

##突变频率柱状图
Mut_Genes <- c("AKT1","ALK","APC",
               "ATM","AXL","BRAF",
               "BRCA1","BRCA2","CCND1",
               "CDK4","CDK6","CDKN2A",
               "CTNNB1","DDR2","EGFR",
               "ERBB2","ERBB3","ERBB4",
               "ESR1","FGF19","FGF3",
               "FGF4","FGFR1","FGFR2",
               "FGFR3","FLT3","HRAS",
               "IDH1","IDH2","IGF1R",
               "JAK1","JAK2","KDR",
               "KIT","KRAS","MAP2K1",
               "MET","MTOR","MYC",
               "NF1","NOTCH1","NRAS",
               "NTRK1","NTRK2","NTRK3",
               "PDGFRA","PIK3CA",
               "PTCH1","PTEN","RAF1","RB1",
               "RET","ROS1","SMAD4","SMO",
               "STK11","TOP2A","TP53",
               "TSC1","TSC2","AR",
               "ARAF","BCL2L11","CD74","NRG1")
Fusion_Genes <- c("ALK","RET","ROS1","CD74","FGFR3","NRG1","NTRK1")
exon_skipping_gene <- "MET"

##加上外显子序号和突变类型
WGS_long <- WGS_maf %>%
  mutate(method="WGS",Exon_Number=str_extract(Exon_Number, "^[^/]+")) %>%
  filter(Patient_ID %in% clin_merge$Patient_ID) %>%
  select(Patient_ID,Hugo_Symbol,method,HGVSp,ONCOGENIC,Exon_Number,Variant_Classification)
NGS_long <- NGS_long %>%
  mutate(method="NGS",Patient_ID=住院号,Hugo_Symbol=Gene,ONCOGENIC=oncokb_oncogenic)%>%
  filter(Patient_ID %in% clin_merge$Patient_ID) %>%
  mutate(Exon_Number=外显子序号,Variant_Classification=突变类型) %>%
  separate(HGVS_P, into = c("prefix", "suffix"), sep = "\\.", remove = FALSE, extra = "merge") %>%
  mutate(HGVSp = suffix) %>%
  select(-prefix, -suffix) %>%
  select(Patient_ID,Hugo_Symbol,method,HGVSp,ONCOGENIC,Exon_Number,Variant_Classification) 
merge_long <- rbind(WGS_long,NGS_long)
save(merge_long,file="~/molecular_staging/RData/merge_long.Rdata")
