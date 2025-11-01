setwd("~/molecular_staging/")

##loading packages
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)

##loading files
load("RData/clin_merge.Rdata")

load("RData/merge_long.Rdata")

##step1 把突变整合到临床信息中
#EGFR

#整理表格
EGFR_mut_long <- merge_long %>%
  filter(Hugo_Symbol=="EGFR",ONCOGENIC %in% c("Oncogenic","Likely Oncogenic"))

EGFR_mut_long <- unique(EGFR_mut_long)
EGFR_mut_long$HGVSp <- gsub("^p\\.", "", EGFR_mut_long$HGVSp)
#Varianr Classification太杂了，重新弄一下
# 第一步：创建辅助数据框，包含每个HGVSp对应的英文分类
english_terms <- EGFR_mut_long %>%
  # 筛选英文行（假设英文只包含字母和数字，不包含中文）
  filter(grepl("^[A-Za-z0-9_]+$", Variant_Classification)) %>%
  group_by(HGVSp) %>%
  # 检查每个HGVSp对应的英文分类是否唯一
  summarise(
    eng_class = if(n_distinct(Variant_Classification) == 1) 
      first(Variant_Classification) else NA_character_,
    .groups = "drop"
  )

# 第二步：创建新列，中文分类替换为对应的英文分类
EGFR_mut_long <- EGFR_mut_long %>%
  left_join(english_terms, by = "HGVSp") %>%
  mutate(
    Variant_Classification1 = ifelse(
      # 如果是中文且能找到唯一对应的英文分类
      grepl("[\\p{Han}]", Variant_Classification, perl = TRUE) & !is.na(eng_class),
      eng_class,  # 替换为英文
      Variant_Classification  # 否则保持原样
    )
  ) %>%
  select(-eng_class)  # 移除辅助列

#框内插入突变/框内缺失突变/错义突变
#错义突变有几个不对的
#543162 Glu746_Ala750delinsGlnPro "In_Frame_Del"
#599697 Glu746_Leu747delinsValPro "Others"
#626955 Glu746_Leu747delinsAspPro "Others

EGFR_mut_long$Variant_Classification1 <- ifelse(EGFR_mut_long$Variant_Classification1 == "框内插入突变","In_Frame_Ins",
                                                ifelse(EGFR_mut_long$Variant_Classification1 == "框内缺失突变","In_Frame_Del",
                                                       ifelse(EGFR_mut_long$Variant_Classification1 == "错义突变","Missense_Mutation",
                                                              EGFR_mut_long$Variant_Classification1)))

EGFR_mut_long$Variant_Classification1[EGFR_mut_long$Patient_ID==543162] <- "In_Frame_Del"
EGFR_mut_long$Variant_Classification1[EGFR_mut_long$Patient_ID==599697] <- "Others"
EGFR_mut_long$Variant_Classification1[EGFR_mut_long$Patient_ID==626955] <- "Others"



#复杂突变 del和ins数量相同
EGFR_mut_long$Variant_Classification1[EGFR_mut_long$Variant_Classification1=="复杂突变"] <- "Others"

#插入缺失突变/框内插入缺失突变/缺失突变/非移码突变/框内移码突变
#589893 Ser768_Val769delinsIleLeu "Others"
#590154 Glu746_Leu747delinsAspPro "Others"
#563560 His773delinsAsnProTyr "In_Frame_Ins"

#其它都是 In_Frame_Del
EGFR_mut_long$Variant_Classification1[EGFR_mut_long$Variant_Classification1 %in% 
                                        c("插入缺失突变","框内插入缺失突变","缺失突变","非移码突变","框内移码突变")] <- "In_Frame_Del"
EGFR_mut_long$Variant_Classification1[EGFR_mut_long$Patient_ID %in% c(589893,590154)] <- "Others"
EGFR_mut_long$Variant_Classification1[EGFR_mut_long$Patient_ID == 563560] <- "In_Frame_Ins"

EGFR_mut_long$type <- paste0("exon",EGFR_mut_long$Exon_Number)
EGFR_mut_long$type <- paste(EGFR_mut_long$type,EGFR_mut_long$Variant_Classification1)

EGFR_mut_long$category <- ifelse((EGFR_mut_long$HGVSp == "Leu858Arg" | EGFR_mut_long$type=="exon19 In_Frame_Del"),"EGFR Canonical","EGFR Noncanonical")
EGFR_mut_long <- EGFR_mut_long[,c("Patient_ID","category")]
save(EGFR_mut_long,file = "RData/EGFR突变详细情况.Rdata")


EGFR_mut_long$Patient_ID <- as.character(EGFR_mut_long$Patient_ID)
EGFR_mut_long <- unique(EGFR_mut_long)
EGFR_mut_long <- EGFR_mut_long %>%
  group_by(Patient_ID) %>%
  arrange(Patient_ID, factor(category, levels = c("EGFR_Canonical", "EGFR_Noncanonical"))) %>%
  slice(1) %>%
  ungroup()

clin_mut <- clin_merge %>%
  left_join(EGFR_mut_long,by=c("Patient_ID"))
clin_mut$EGFR_category <- ifelse(is.na(clin_mut$category),"EGFR WT",clin_mut$category)
clin_mut$category <- NULL

#KRAS
##整理表格
KRAS_mut_long <- merge_long %>%
  filter(Hugo_Symbol=="KRAS",ONCOGENIC %in% c("Oncogenic","Likely Oncogenic"))
KRAS_mut_long <- unique(KRAS_mut_long)
KRAS_mut_long$HGVSp <- gsub("^p\\.", "", KRAS_mut_long$HGVSp)
KRAS_mut_long$Variant_Classification[KRAS_mut_long$Variant_Classification=="错义突变"] <- "Missense_Mutation"
KRAS_mut_long$Variant_Classification[KRAS_mut_long$Patient_ID %in% c(571225,649024)] <- "In_Frame_Ins"
aa_table <- data.frame(
  tri = c("Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", 
          "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", 
          "Tyr", "Val"),
  single = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", 
             "I", "L", "K", "M", "F", "P", "S", "T", "W", 
             "Y", "V")
)
replace_aa_tri_to_single <- function(input_string, aa_table) {
  # 找出所有三字母氨基酸缩写
  tri_matches <- str_extract_all(input_string, "[A-Z][a-z]{2}")[[1]]
  
  # 只保留在aa_table中存在的有效氨基酸缩写
  valid_tri <- tri_matches[tri_matches %in% aa_table$tri]
  
  # 逐个替换
  for(tri in valid_tri) {
    single <- aa_table$single[aa_table$tri == tri]
    input_string <- str_replace_all(input_string, tri, single)
  }
  
  return(input_string)
}
KRAS_mut_long$HGVSp_short <- sapply(KRAS_mut_long$HGVSp, function(x) replace_aa_tri_to_single(x, aa_table))
KRAS_mut_long <- KRAS_mut_long %>%
      mutate(category1=ifelse(str_detect(HGVSp_short,"12"),"G12X",
      ifelse(str_detect(HGVSp_short,"13"),"G13X",
      ifelse(str_detect(HGVSp_short,"61"),"Q61X","others")))) %>%
      mutate(category=ifelse(category1 %in% c("G12X","G13X","G61X"),"KRAS Canonical","KRAS Noncanonical"))
save(KRAS_mut_long,file = "RData/KRAS突变详细情况.Rdata")

KRAS_mut_long <- KRAS_mut_long[,c("Patient_ID","category")]
KRAS_mut_long$Patient_ID <- as.character(KRAS_mut_long$Patient_ID)
KRAS_mut_long <- unique(KRAS_mut_long)
KRAS_mut_long <- KRAS_mut_long %>%
  group_by(Patient_ID) %>%
  arrange(Patient_ID, factor(category, levels = c("KRAS_Canonical", "KRAS_Noncanonical"))) %>%
  slice(1) %>%
  ungroup()

clin_mut <- clin_mut %>%
  left_join(KRAS_mut_long,by=c("Patient_ID"))
clin_mut$KRAS_category <- ifelse(is.na(clin_mut$category),"KRAS WT",clin_mut$category)
clin_mut$category <- NULL

#TP53
TP53_mut_long <- merge_long %>% 
  filter(Hugo_Symbol=="TP53" & ONCOGENIC %in% c("Oncogenic","Likely Oncogenic")) %>%
  mutate(category = "TP53 Mut")
TP53_mut_long <- unique(TP53_mut_long[,c("Patient_ID","category")])
TP53_mut_long$Patient_ID <- as.character(TP53_mut_long$Patient_ID)
clin_mut <- clin_mut %>%
  left_join(TP53_mut_long,by=c("Patient_ID"))
clin_mut$TP53_category <- ifelse(is.na(clin_mut$category),"TP53 Wt",clin_mut$category)
clin_mut$category <- NULL
save(clin_mut,file = "RData/clin_mut.Rdata")

##step2 把想要的临床-突变信息表格整理出来 
# 定义临床变量和突变变量
clinical_vars <- c("Sex", "Age_E", "Smoking_E", "Stage1", "Subtype_short")
mutation_vars <- c("EGFR_category", "KRAS_category", "TP53_category")

# 步骤一：生成包含突变信息的宽表格 (与之前相同)
summary_table_mutations <- clin_mut %>%
  mutate(across(all_of(clinical_vars), as.character)) %>%
  select(all_of(c(clinical_vars, mutation_vars))) %>%
  pivot_longer(cols = all_of(clinical_vars), names_to = "Clinical_Variable", values_to = "Clinical_Level") %>%
  pivot_longer(cols = all_of(mutation_vars), names_to = "Mutation_Gene", values_to = "Mutation_Status") %>%
  group_by(Clinical_Variable, Clinical_Level, Mutation_Gene, Mutation_Status) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Mutation_Gene, Mutation_Status, Clinical_Variable) %>%
  mutate(
    Percentage = round(Count / sum(Count) * 100, 1),
    Count_Percent = paste0(Count, " (", Percentage, "%)")
  ) %>%
  ungroup() %>%
  select(-Count, -Percentage) %>%
  pivot_wider(
    names_from = c(Mutation_Gene, Mutation_Status),
    values_from = Count_Percent,
    values_fill = "0 (0.0%)"
  )

# 步骤二：计算每个临床水平的总数和百分比
summary_table_totals <-clin_mut %>%
  mutate(across(all_of(clinical_vars), as.character)) %>%
  pivot_longer(cols = all_of(clinical_vars), names_to = "Clinical_Variable", values_to = "Clinical_Level") %>%
  group_by(Clinical_Variable, Clinical_Level) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Clinical_Variable) %>%
  mutate(
    Total_Percentage = round(Count / sum(Count) * 100, 1),
    Total = paste0(Count, " (", Total_Percentage, "%)")
  ) %>%
  ungroup() %>%
  select(-Count, -Total_Percentage)

# 步骤三：合并两个表格
final_summary_table <- summary_table_mutations %>%
  left_join(summary_table_totals, by = c("Clinical_Variable", "Clinical_Level")) %>%
  arrange(Clinical_Variable, Clinical_Level)

write.csv(final_summary_table,file="output/table1_临床病理信息.csv")
