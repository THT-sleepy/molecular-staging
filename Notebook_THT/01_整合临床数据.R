##loading packages
library(dplyr)

##loading files
WGS_clin <- read.csv("~/molecular_staging/input/lc902_clin.csv")
NGS_clin <- read.csv("~/molecular_staging/input/08-21_FUSCC_LUAD_NGS_clin.csv")

##得到合并的临床信息表格
clin_wgs <- WGS_clin %>%
  mutate(method="WGS",Stage=TNM9,Smoking=ifelse(Smoking_history_brief=="Never",0,1)) %>%
  select(Patient_ID,method,Sex,Age,Smoking,Stage,RFS,RFS_E,OS,OS_E,Subtype,Subtype_short)
clin_ngs <- NGS_clin %>%
  mutate(method="Panel",Stage=TNM9_Stage,
         Sex=ifelse(gender==1,"Male","Female"),
         RFS_E=RFSS,
         OS_E=OSS,Subtype=ifelse(SUBGROUP1 %in% c("AIS","MIA"),0,
                                 ifelse(SUBGROUP1 == "APA","Acinar",
                                        ifelse(SUBGROUP1 == "APA/PPA","Acinar/Papillary",
                                               ifelse(SUBGROUP1 == "IMA","Invasive_mucinous",
                                                      ifelse(SUBGROUP1 == "LPA","Lepidic",
                                                             ifelse(SUBGROUP1 == "MPA","Micropapillary",
                                                                    ifelse(SUBGROUP1 == "PPA","Papillary",
                                                                           ifelse(SUBGROUP1 == "SPA","Solid","Unreported")))))))),
         Subtype_short=ifelse(Subtype=="Lepidic","Lepidic",
                              ifelse(Subtype %in% c("Acinar","Papillary","Invasive_mucinous","Acinar/Papillary"),"Acinar/Papillary/IMA",
                                     ifelse(Subtype %in% c("Micropapillary","Solid"),"Micropapillary/Solid",
                                            ifelse(Subtype==0,"AAH/AIS/MIA","Unreported")))))%>%
  select(Patient_ID,method,Sex,Age,Smoking,Stage,RFS,RFS_E,OS,OS_E,Subtype,Subtype_short)
clin_merge <- rbind(clin_ngs,clin_wgs)
# 步骤2: 对于Patient_ID有重复的行(即接受了两种检测)，只保留OS最大的行
# 如果OS相同，则保留第一个出现的行
clin_merge <- clin_merge %>%
  group_by(Patient_ID) %>%
  filter(OS == max(OS)) %>%
  # 如果OS相同，保留第一个出现的行
  filter(row_number() == 1) %>%
  ungroup()
clin_merge$Patient_ID <- as.character(clin_merge$Patient_ID)
clin_merge$Subtype_short <- ifelse(clin_merge$Subtype_short %in% c("0","AAH/AIS/MIA"),"MIA",clin_merge$Subtype_short)
clin_merge$Age_E <- ifelse(clin_merge$Age<60 ,"<60",">=60")
clin_merge$Smoking_E <- ifelse(clin_merge$Smoking==1,"Ever","Never")
clin_merge$Stage1 <- ifelse(clin_merge$Stage=="MIA","IA",clin_merge$Stage)
clin_merge$Stage2 <- ifelse(clin_merge$Stage1 %in% c("IA","IB"),"I",
                     ifelse(clin_merge$Stage1 %in% c("IIA","IIB"),"II",
                     ifelse(clin_merge$Stage1 %in% c("IIIA","IIIB"),"III",
                     ifelse(clin_merge$Stage1 %in% c("IVA"),"IV","?"))))
clin_merge$Differentiation <- ifelse(clin_merge$Subtype_short %in% c("MIA","Lepidic"),"Well",
                              ifelse(clin_merge$Subtype_short =="Acinar/Papillary/IMA","Mod",
                              ifelse(clin_merge$Subtype_short=="Micropapillary/Solid","Poor",NA)))

save(clin_merge,file = "~/molecular_staging/RData/clin_merge.Rdata")
##522977只有OS没有RFS



