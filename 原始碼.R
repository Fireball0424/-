#-----library and read data--------
library(dplyr)
library(ggplot2)
library(tidyverse)
library(pheatmap)

df <- read.csv("/media/md1200/public/TCFSH_project/CGC_data.csv")
fn <- read.table("/media/md1200/public/TCFSH_project/cancer_fullname.txt", header = T)
cgt <- read.table("/media/md1200/public/TCFSH_project/drivergene2muttool.txt",header = T)
nt <- read.table("/media/md1200/public/TCFSH_project/drivergene_ntool_mutation.txt",header = T)

#View(df) # all genes and details
#View(fn) # full name 
#View(cgt) # cancer type with gene and tool
#View(nt) # define by n_tools


# gene = Gym.Sympol
df <- df %>% 
  mutate(gene = Gene.Symbol) %>% 
  mutate(Gene.Symbol = NULL)

#-----statics of n_tools---------------
ggplot(nt) +
  geom_histogram(aes(x = n_tool), bins = 30, binwidth = 0.5)+
  stat_count(aes(x = n_tool),geom = "point")+
  stat_count(aes(x = n_tool),geom = "smooth", color  = "black")+
  geom_label(stat = "count", aes(x = n_tool, label = ..count..), position = position_dodge(0.9), vjust = 0)



nt_1 <- nt %>% 
  group_by(n_tool) %>% 
  summarise(count = n()) %>% 
  mutate(per = c(8.32, 12.33, 21.49, 47.98, 83.79, 98.61, 98.92, 100.00, 100.00, 100.00, 100.00, 100.00))

ggplot(nt_1)+
  geom_bar(aes(x = n_tool, y = log10(count)), stat = "identity")+
  geom_line(aes(x = n_tool, y = per*5/100))+
  geom_point(aes(x = n_tool, y = per*5/100), size = 3, shape = 21, fill = "white")+
  scale_y_continuous(name = "log10(count)", limits = c(0,5), sec.axis = sec_axis(~.*100/5, name = "percentage(%)"))+
  scale_x_continuous(breaks = seq(0, 12, by = 1)) 
#-----ratio of nts genes in CGC-----------------
df_all <- merge(df, nt, by = "gene")
ggplot(df_all)+
  geom_bar(aes((x = n_tool)))+
  geom_label(stat = "count", aes(x = n_tool, label =..count..))

# n_tools >= 3
nt_3 <- nt %>% 
  filter(n_tool > 2) # -> 3146 ob

df_all_1 <- merge(df, nt_3, by = "gene") # 676 ob
## -> 676/3146 =21.49%

# n_tools >= 4
nt_4 <- nt %>% 
  filter(n_tool > 3) # -> 719 ob

df_all_2 <- merge(df, nt_4, by = "gene") # -> 345 ob
## -> 345/719 = 47.98%

# n_tools >= 5
nt_5 <- nt %>% 
  filter(n_tool > 4) # -> 253 ob

df_all_3 <- merge(df, nt_5, by = "gene") # -> 212 ob
## -> 212/253 = 83.79%

# n_tools >= 6
nt_6 <- nt %>% 
  filter(n_tool > 5) # -> 144 ob

df_all_4 <- merge(df, nt_6, by = "gene") # -> 142 ob
## -> 142/144 = 98.61%

su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  all = c(21.49, 47.98, 83.79, 98.61),
  num = c("676/3146", "345/719", "212/253", "142/144")
)

ggplot(su)+
  geom_col(aes(x = n_tools, y = all))+
  geom_label(aes(x = n_tools, y = all, label = all))+
  ylab("cgc percentage(%)")

#-----the gene out of CGC dataframe-----------------------
# -> according the percentage, deciding to use n_tools >= 5 & n_tools >= 6
df_5_out_cgc <- anti_join(nt_5, df, by = "gene")

df_6_out_cgc <- anti_join(nt_6, df, by = "gene")

df_5_out_cgc <- df_5_out_cgc %>% 
  mutate(cancer_type_abbr = NULL,
         n_tool = NULL)

df_5_out_cgc <- unique(df_5_out_cgc)

df_5_out <- merge(df_5_out_cgc, nt, by = "gene")

df_5_out_ct <- df_5_out %>% 
  group_by(gene) %>% 
  summarise(count = length(gene)) %>% 
  mutate(gene = fct_reorder(gene, count))

ggplot(df_5_out_ct)+
  geom_col(aes(x = gene, y = count))+
  coord_flip()

#-----gene out of CGC (all cancer)------
#----------ACC--------
nt_acc <- nt %>% 
  filter(cancer_type_abbr == "ACC")
nt_acc_3 <- nt_acc %>% 
  filter(n_tool > 2)
nt_acc_3_in <- merge(nt_acc_3, df, by = "gene") 
# 18 obj / 102 obj = 17.65%
nt_acc_4 <- nt_acc %>% 
  filter(n_tool > 3)
nt_acc_4_in <- merge(nt_acc_4, df, by = "gene")
# 4 obj / 10 obj = 40.00%
nt_acc_5 <- nt_acc %>% 
  filter(n_tool > 4)
nt_acc_5_in <- merge(nt_acc_5, df, by = "gene")
# 2 obj / 3 obj = 66.67%
nt_acc_6 <- nt_acc %>% 
  filter(n_tool > 5)
nt_acc_6_in <- merge(nt_acc_6, df, by = "gene")
# 1 obj / 1 obj = 100%
nt_acc_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  acc = c(17.65, 40.00, 66.67, 100.00)
)

#----------BLCA----------
nt_blca <- nt %>% 
  filter(cancer_type_abbr == "BLCA")
nt_blca_3 <- nt_blca %>% 
  filter(n_tool > 2)
nt_blca_3_in <- merge(nt_blca_3, df, by = "gene")
# 41 obj / 139 obj = 29.50%
nt_blca_4 <- nt_blca %>% 
  filter(n_tool > 3)
nt_blca_4_in <- merge(nt_blca_4, df, by = "gene")
# 20 obj / 38 obj = 52.63%
nt_blca_5 <- nt_blca %>% 
  filter(n_tool > 4)
nt_blca_5_in <- merge(nt_blca_5, df, by = "gene")
# 16 obj / 19 obj = 84.21%
nt_blca_6 <- nt_blca %>% 
  filter(n_tool > 5)
nt_blca_6_in <- merge(nt_blca_6, df, by = "gene")
# 13 obj / 13 obj = 100.00%
nt_blca_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  blca = c(29.50, 52.63, 84.21, 100.00)
)

#----------BRCA--------
nt_brca <- nt %>% 
  filter(cancer_type_abbr == "BRCA")
nt_brca_3 <- nt_brca %>% 
  filter(n_tool > 2)
nt_brca_3_in <- merge(nt_brca_3, df, by = "gene")
# 36 obj / 122 obj = 29.51% 
nt_brca_4 <- nt_brca %>% 
  filter(n_tool > 3)
nt_brca_4_in <- merge(nt_brca_4, df, by = "gene")
# 22 obj / 43 obj = 51.16%
nt_brca_5 <- nt_brca %>% 
  filter(n_tool > 4)
nt_brca_5_in <- merge(nt_brca_5, df, by = "gene")
# 19 obj / 22 obj = 86.36%
nt_brca_6 <- nt_brca %>% 
  filter(n_tool > 5)
nt_brca_6_in <- merge(nt_brca_6, df, by = "gene")
# 12 obj / 12 obj = 100.00%
nt_brca_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  brca = c(29.51, 51.16, 86.36, 100.00)
)

#----------CESC------
nt_cesc <- nt %>% 
  filter(cancer_type_abbr == "CESC")
nt_cesc_3 <- nt_cesc %>% 
  filter(n_tool > 2)
nt_cesc_3_in <- merge(nt_cesc_3, df, by = "gene")
# 26 obj / 74 obj = 35.14%
nt_cesc_4 <- nt_cesc %>% 
  filter(n_tool > 3)
nt_cesc_4_in <- merge(nt_cesc_4, df, by = "gene")
# 18 obj / 21 obj = 85.71%
nt_cesc_5 <- nt_cesc %>% 
  filter(n_tool > 4)
nt_cesc_5_in <- merge(nt_cesc_5, df, by = "gene")
# 10 obj / 11 obj = 90.91%
nt_cesc_6 <- nt_cesc %>% 
  filter(n_tool > 5)
nt_cesc_6_in <- merge(nt_cesc_6, df, by = "gene")
# 7 obj / 7 obj = 100.00%
nt_cesc_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  cesc = c(35.14, 85.71, 90.91, 100.00)
)

#----------CHOL--------
nt_chol <- nt %>% 
  filter(cancer_type_abbr == "CHOL")
nt_chol_3 <- nt_chol %>% 
  filter(n_tool > 2)
nt_chol_3_in <- merge(nt_chol_3, df, by = "gene")
# 4 obj / 6 obj = 66.67%
nt_chol_4 <- nt_chol %>% 
  filter(n_tool > 3)
nt_chol_4_in <- merge(nt_chol_4, df, by = "gene")
# 1 obj / 1 obj = 100.00% 
nt_chol_5 <- nt_chol %>% 
  filter(n_tool > 4)
nt_chol_5_in <- merge(nt_chol_5, df, by = "gene")
# 1 obj / 1 obj = 100.00%
nt_chol_6 <- nt_chol %>% 
  filter(n_tool > 5)
nt_chol_6_in <- merge(nt_chol_6, df, by = "gene")
# 0 obj / 0 obj = NA
nt_chol_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  chol = c(66.67, 100.00, 100.00, NA)
)

#----------COAD-----
nt_coad <- nt %>% 
  filter(cancer_type_abbr == "COAD")
nt_coad_3 <- nt_coad %>%
  filter(n_tool > 2)
nt_coad_3_in <- merge(nt_coad_3, df, by = "gene")
# 66 obj / 312 obj = 21.15%
nt_coad_4 <- nt_coad %>%
  filter(n_tool > 3)
nt_coad_4_in <- merge(nt_coad_4, df, by = "gene")
# 35 obj / 98 obj = 35.71%
nt_coad_5 <- nt_coad %>%
  filter(n_tool > 4)
nt_coad_5_in <- merge(nt_coad_5, df, by = "gene")
# 20 obj / 35 obj = 57.14%
nt_coad_6 <- nt_coad %>%
  filter(n_tool > 5)
nt_coad_6_in <- merge(nt_coad_6, df, by = "gene")
# 13 obj / 13 obj = 100.00%
nt_coad_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  coad = c(21.15, 35.71, 57.14, 100.00)
)

#----------DLBC------
nt_dlbc <- nt %>% 
  filter(cancer_type_abbr == "DLBC")
nt_dlbc_3 <- nt_dlbc %>% 
  filter(n_tool > 2)
nt_dlbc_3_in <- merge(nt_dlbc_3, df, by = "gene")
# 13 obj / 25 obj = 52.00%
nt_dlbc_4 <- nt_dlbc %>% 
  filter(n_tool > 3)
nt_dlbc_4_in <- merge(nt_dlbc_4, df, by = "gene")
# 6 obj / 8 obj = 75.00%
nt_dlbc_5 <- nt_dlbc %>% 
  filter(n_tool > 4)
nt_dlbc_5_in <- merge(nt_dlbc_5, df, by = "gene")
# 2 obj / 3 obj = 66.67%
nt_dlbc_6 <- nt_dlbc %>% 
  filter(n_tool > 5)
nt_dlbc_6_in <- merge(nt_dlbc_6, df, by = "gene")
# 1 obj / 1 obj = 100.00%
nt_dlbc_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  dlbc = c(52.00, 75.00, 66.67, 100.00)
)

#----------ESCA----------
nt_esca <- nt %>% 
  filter(cancer_type_abbr == "ESCA")
nt_esca_3 <- nt_esca %>% 
  filter(n_tool > 2)
nt_esca_3_in <- merge(nt_esca_3, df, by = "gene")
# 9 obj / 22 obj = 40.91%
nt_esca_4 <- nt_esca %>% 
  filter(n_tool > 3)
nt_esca_4_in <- merge(nt_esca_4, df, by = "gene")
# 5 obj / 7 obj = 71.43%
nt_esca_5 <- nt_esca %>% 
  filter(n_tool > 4)
nt_esca_5_in <- merge(nt_esca_5, df, by = "gene")
# 3 obj / 3 obj = 100.00%
nt_esca_6 <- nt_esca %>% 
  filter(n_tool > 5)
nt_esca_6_in <- merge(nt_esca_6, df, by = "gene")
# 1 obj / 1 obj = 100.00%
nt_esca_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  esca = c(40.91, 71.43, 100.00, 100.00)
)

#-----------GBM----------
nt_gbm <- nt %>% 
  filter(cancer_type_abbr == "GBM")
nt_gbm_3 <- nt_gbm %>% 
  filter(n_tool > 2)
nt_gbm_3_in <- merge(nt_gbm_3, df, by = "gene")
# 26 obj / 77 obj = 33.77%
nt_gbm_4 <- nt_gbm %>% 
  filter(n_tool > 3)
nt_gbm_4_in <- merge(nt_gbm_4, df, by = "gene")
# 13 obj / 19 obj = 68.42%
nt_gbm_5 <- nt_gbm %>% 
  filter(n_tool > 4)
nt_gbm_5_in <- merge(nt_gbm_5, df, by = "gene")
# 7 obj / 7 obj = 100.00%
nt_gbm_6 <- nt_gbm %>% 
  filter(n_tool > 5)
nt_gbm_6_in <- merge(nt_gbm_6, df, by = "gene")
# 7 obj / 7 obj = 100.00%
nt_gbm_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  gbm = c(33.77, 68.42, 100.00, 100.00)
)

#----------HNSC----------
nt_hnsc <- nt %>% 
  filter(cancer_type_abbr == "HNSC")
nt_hnsc_3 <- nt_hnsc %>% 
  filter(n_tool > 2)
nt_hnsc_3_in <- merge(nt_hnsc_3, df, by = "gene")
# 24 obj / 87 obj = 27.59%
nt_hnsc_4 <- nt_hnsc %>% 
  filter(n_tool > 3)
nt_hnsc_4_in <- merge(nt_hnsc_4, df, by = "gene")
# 17 obj / 27 obj = 62.96%
nt_hnsc_5 <- nt_hnsc %>% 
  filter(n_tool > 4)
nt_hnsc_5_in <- merge(nt_hnsc_5, df, by = "gene")
# 13 obj / 14 obj = 92.86%
nt_hnsc_6 <- nt_hnsc %>% 
  filter(n_tool > 5)
nt_hnsc_6_in <- merge(nt_hnsc_6, df, by = "gene")
# 8 obj / 8 obj = 100.00%
nt_hnsc_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  hnsc = c(27.59, 62.96, 92.86, 100.00)
)

#----------KICH----------
nt_kich <- nt %>% 
  filter(cancer_type_abbr == "KICH")
nt_kich_3 <- nt_kich %>% 
  filter(n_tool > 2)
nt_kich_3_in <- merge(nt_kich_3, df, by = "gene")
# 1 obj / 1 obj = 100.00%
nt_kich_4 <- nt_kich %>% 
  filter(n_tool > 3)
nt_kich_4_in <- merge(nt_kich_4, df, by = "gene")
# 1 obj / 1 obj = 100.00% 
nt_kich_5 <- nt_kich %>% 
  filter(n_tool > 4)
nt_kich_5_in <- merge(nt_kich_5, df, by = "gene")
# 0 obj / 0 obj = NA
nt_kich_6 <- nt_kich %>% 
  filter(n_tool > 5)
nt_kich_6_in <- merge(nt_kich_6, df, by = "gene")
# 0 obj / 0 obj = NA
nt_kich_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  kich = c(100.00, 100.00, NA, NA)
)

#----------KIRC----------
nt_kirc <- nt %>% 
  filter(cancer_type_abbr == "KIRC")
nt_kirc_3 <- nt_kirc %>% 
  filter(n_tool > 2)
nt_kirc_3_in <- merge(nt_kirc_3, df, by = "gene")
# 8 obj / 24 obj = 33.33%
nt_kirc_4 <- nt_kirc %>% 
  filter(n_tool > 3)
nt_kirc_4_in <- merge(nt_kirc_4, df, by = "gene")
# 3 obj / 3 obj = 100.00%
nt_kirc_5 <- nt_kirc %>% 
  filter(n_tool > 4)
nt_kirc_5_in <- merge(nt_kirc_5, df, by = "gene")
# 1 obj / 1 obj = 100.00%
nt_kirc_6 <- nt_kirc %>% 
  filter(n_tool > 5)
nt_kirc_6_in <- merge(nt_kirc_6, df, by = "gene")
# 0 obj / 0 obj = NA
nt_kirc_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  kirc = c(33.33, 100.00, 100.00, NA)
)

#----------KIRP----------
nt_kirp <- nt %>% 
  filter(cancer_type_abbr == "KIRP")
nt_kirp_3 <- nt_kirp %>% 
  filter(n_tool > 2)
nt_kirp_3_in <- merge(nt_kirp_3, df, by = "gene")
# 11 obj / 39 obj = 28.21%
nt_kirp_4 <- nt_kirp %>% 
  filter(n_tool > 3)
nt_kirp_4_in <- merge(nt_kirp_4, df, by = "gene")
# 3 obj / 4 obj = 75.00%
nt_kirp_5 <- nt_kirp %>% 
  filter(n_tool > 4)
nt_kirp_5_in <- merge(nt_kirp_5, df, by = "gene")
# 1 obj / 1 obj = 100.00%
nt_kirp_6 <- nt_kirp %>% 
  filter(n_tool > 5)
nt_kirp_6_in <- merge(nt_kirp_6, df, by = "gene")
# 1 obj / 1 obj = 100.00%
nt_kirp_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  kirp = c(28.21, 75.00, 100.00, 100.00)
)

#----------LAML----------
nt_laml<- nt %>% 
  filter(cancer_type_abbr == "LAML")
nt_laml_3 <- nt_laml %>% 
  filter(n_tool > 2)
nt_laml_3_in <- merge(nt_laml_3, df, by = "gene")
# 16 obj / 21 obj = 76.19%
nt_laml_4 <- nt_laml %>% 
  filter(n_tool > 3)
nt_laml_4_in <- merge(nt_laml_4, df, by = "gene")
# 14 obj / 15 obj = 93.33%
nt_laml_5 <- nt_laml %>% 
  filter(n_tool > 4)
nt_laml_5_in <- merge(nt_laml_5, df, by = "gene")
# 13 obj / 13 obj = 100.00%
nt_laml_6 <- nt_laml %>% 
  filter(n_tool > 5)
nt_laml_6_in <- merge(nt_laml_6, df, by = "gene")
# 11 obj / 11 obj = 100.00%
nt_laml_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  laml = c(76.19, 93.33, 100.00, 100.00)
)

#-----------LGG----------
nt_lgg<- nt %>% 
  filter(cancer_type_abbr == "LGG")
nt_lgg_3 <- nt_lgg %>% 
  filter(n_tool > 2)
nt_lgg_3_in <- merge(nt_lgg_3, df, by = "gene")
# 14 obj / 32 obj = 43.75%
nt_lgg_4 <- nt_lgg %>% 
  filter(n_tool > 3)
nt_lgg_4_in <- merge(nt_lgg_4, df, by = "gene")
# 13 obj / 15 obj = 86.67%
nt_lgg_5 <- nt_lgg %>% 
  filter(n_tool > 4)
nt_lgg_5_in <- merge(nt_lgg_5, df, by = "gene")
# 7 obj / 7 obj = 100.00%
nt_lgg_6 <- nt_lgg %>% 
  filter(n_tool > 5)
nt_lgg_6_in <- merge(nt_lgg_6, df, by = "gene")
# 5 obj / 5 obj = 100.00%
nt_lgg_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  lgg = c(43.75, 86.67, 100.00, 100.00)
)

#----------LIHC----------
nt_lihc<- nt %>% 
  filter(cancer_type_abbr == "LIHC")
nt_lihc_3 <- nt_lihc %>% 
  filter(n_tool > 2)
nt_lihc_3_in <- merge(nt_lihc_3, df, by = "gene")
# 22 obj / 58 obj = 37.93%
nt_lihc_4 <- nt_lihc %>% 
  filter(n_tool > 3)
nt_lihc_4_in <- merge(nt_lihc_4, df, by = "gene")
# 7 obj / 13 obj = 53.87%
nt_lihc_5 <- nt_lihc %>% 
  filter(n_tool > 4)
nt_lihc_5_in <- merge(nt_lihc_5, df, by = "gene")
# 4 obj / 4 obj = 100.00%
nt_lihc_6 <- nt_lihc %>% 
  filter(n_tool > 5)
nt_lihc_6_in <- merge(nt_lihc_6, df, by = "gene")
# 4 obj / 4 obj = 100.00%
nt_lihc_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  lihc = c(37.93, 53.87, 100.00, 100.00)
)

#----------LUAD----------
nt_luad<- nt %>% 
  filter(cancer_type_abbr == "LUAD")
nt_luad_3 <- nt_luad %>% 
  filter(n_tool > 2)
nt_luad_3_in <- merge(nt_luad_3, df, by = "gene")
# 23 obj / 141 obj = 16.31%
nt_luad_4 <- nt_luad %>% 
  filter(n_tool > 3)
nt_luad_4_in <- merge(nt_luad_4, df, by = "gene")
# 12 obj / 29 obj = 41.38%
nt_luad_5 <- nt_luad %>% 
  filter(n_tool > 4)
nt_luad_5_in <- merge(nt_luad_5, df, by = "gene")
# 9 obj / 10 obj = 90.00%
nt_luad_6 <- nt_luad %>% 
  filter(n_tool > 5)
nt_luad_6_in <- merge(nt_luad_6, df, by = "gene")
# 7 obj / 7 obj = 100.00%
nt_luad_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  luad = c(16.31, 41.38, 90.00, 100.00)
)

#----------LUSC----------
nt_lusc<- nt %>% 
  filter(cancer_type_abbr == "LUSC")
nt_lusc_3 <- nt_lusc %>% 
  filter(n_tool > 2)
nt_lusc_3_in <- merge(nt_lusc_3, df, by = "gene")
# 20 obj / 112 obj = 17.86%
nt_lusc_4 <- nt_lusc %>% 
  filter(n_tool > 3)
nt_lusc_4_in <- merge(nt_lusc_4, df, by = "gene")
# 9 obj / 22 obj = 40.91%
nt_lusc_5 <- nt_lusc %>% 
  filter(n_tool > 4)
nt_lusc_5_in <- merge(nt_lusc_5, df, by = "gene")
# 4 obj / 4 obj = 100.00%
nt_lusc_6 <- nt_lusc %>% 
  filter(n_tool > 5)
nt_lusc_6_in <- merge(nt_lusc_6, df, by = "gene")
# 3 obj / 3 obj = 100.00%
nt_lusc_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  lusc = c(17.86, 40.91, 100.00, 100.00)
)

#----------MESO----------
nt_meso<- nt %>% 
  filter(cancer_type_abbr == "MESO")
nt_meso_3 <- nt_meso %>% 
  filter(n_tool > 2)
nt_meso_3_in <- merge(nt_meso_3, df, by = "gene")
# 2 obj / 3 obj = 66.67%
nt_meso_4 <- nt_meso %>% 
  filter(n_tool > 3)
nt_meso_4_in <- merge(nt_meso_4, df, by = "gene")
# 1 obj / 1 obj = 100.00%
nt_meso_5 <- nt_meso %>% 
  filter(n_tool > 4)
nt_meso_5_in <- merge(nt_meso_5, df, by = "gene")
# 1 obj / 1 obj = 100.00%
nt_meso_6 <- nt_meso %>% 
  filter(n_tool > 5)
nt_meso_6_in <- merge(nt_meso_6, df, by = "gene")
# 0 obj / 0 obj = NA
nt_meso_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  meso = c(66.67, 100.00, 100.00, NA)
)

#------------OV----------
nt_ov<- nt %>% 
  filter(cancer_type_abbr == "OV")
nt_ov_3 <- nt_ov %>% 
  filter(n_tool > 2)
nt_ov_3_in <- merge(nt_ov_3, df, by = "gene")
# 20 obj / 73 obj = 27.40%
nt_ov_4 <- nt_ov %>% 
  filter(n_tool > 3)
nt_ov_4_in <- merge(nt_ov_4, df, by = "gene")
# 6 obj / 15 obj = 40.00%
nt_ov_5 <- nt_ov %>% 
  filter(n_tool > 4)
nt_ov_5_in <- merge(nt_ov_5, df, by = "gene")
# 2 obj / 2 obj = 100.00%
nt_ov_6 <- nt_ov %>% 
  filter(n_tool > 5)
nt_ov_6_in <- merge(nt_ov_6, df, by = "gene")
# 1 obj / 1 obj = 100.00%
nt_ov_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  ov = c(27.40, 40.00, 100.00, 100.00)
)

#----------PAAD----------
nt_paad<- nt %>% 
  filter(cancer_type_abbr == "PAAD")
nt_paad_3 <- nt_paad %>% 
  filter(n_tool > 2)
nt_paad_3_in <- merge(nt_paad_3, df, by = "gene")
# 9 obj / 19 obj = 47.37%
nt_paad_4 <- nt_paad %>% 
  filter(n_tool > 3)
nt_paad_4_in <- merge(nt_paad_4, df, by = "gene")
# 4 obj / 4 obj = 100.00%
nt_paad_5 <- nt_paad %>% 
  filter(n_tool > 4)
nt_paad_5_in <- merge(nt_paad_5, df, by = "gene")
# 4 obj / 4 obj = 100.00%
nt_paad_6 <- nt_paad %>% 
  filter(n_tool > 5)
nt_paad_6_in <- merge(nt_paad_6, df, by = "gene")
# 1 obj / 1 obj = 100.00%
nt_paad_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  paad = c(47.37, 100.00, 100.00, 100.00)
)

#----------PCPG----------
nt_pcpg<- nt %>% 
  filter(cancer_type_abbr == "PCPG")
nt_pcpg_3 <- nt_pcpg %>% 
  filter(n_tool > 2)
nt_pcpg_3_in <- merge(nt_pcpg_3, df, by = "gene")
# 4 obj / 4 obj = 100.00%
nt_pcpg_4 <- nt_pcpg %>% 
  filter(n_tool > 3)
nt_pcpg_4_in <- merge(nt_pcpg_4, df, by = "gene")
# 4 obj / 4 obj = 100.00%
nt_pcpg_5 <- nt_pcpg %>% 
  filter(n_tool > 4)
nt_pcpg_5_in <- merge(nt_pcpg_5, df, by = "gene")
# 2 obj / 2 obj = 100.00%
nt_pcpg_6 <- nt_pcpg %>% 
  filter(n_tool > 5)
nt_pcpg_6_in <- merge(nt_pcpg_6, df, by = "gene")
# 1 obj / 1 obj = 100.00%
nt_pcpg_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  pcpg = c(100.00, 100.00, 100.00, 100.00)
)

#----------PRAD----------
nt_prad<- nt %>% 
  filter(cancer_type_abbr == "PRAD")
nt_prad_3 <- nt_prad %>% 
  filter(n_tool > 2)
nt_prad_3_in <- merge(nt_prad_3, df, by = "gene")
# 22 obj / 84 obj = 26.19%
nt_prad_4 <- nt_prad %>% 
  filter(n_tool > 3)
nt_prad_4_in <- merge(nt_prad_4, df, by = "gene")
# 12 obj / 20 obj = 60.00%
nt_prad_5 <- nt_prad %>% 
  filter(n_tool > 4)
nt_prad_5_in <- merge(nt_prad_5, df, by = "gene")
# 6 obj / 8 obj = 75.00%
nt_prad_6 <- nt_prad %>% 
  filter(n_tool > 5)
nt_prad_6_in <- merge(nt_prad_6, df, by = "gene")
# 5 obj / 5 obj = 100.00%
nt_prad_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  prad = c(26.19, 60.00, 75.00, 100.00)
)

#----------READ----------
nt_read<- nt %>% 
  filter(cancer_type_abbr == "READ")
nt_read_3 <- nt_read %>% 
  filter(n_tool > 2)
nt_read_3_in <- merge(nt_read_3, df, by = "gene")
# 27 obj / 94 obj = 28.72%
nt_read_4 <- nt_read %>% 
  filter(n_tool > 3)
nt_read_4_in <- merge(nt_read_4, df, by = "gene")
# 15 obj / 23 obj = 65.22%
nt_read_5 <- nt_read %>% 
  filter(n_tool > 4)
nt_read_5_in <- merge(nt_read_5, df, by = "gene")
# 8 obj / 8 obj = 100.00%
nt_read_6 <- nt_read %>% 
  filter(n_tool > 5)
nt_read_6_in <- merge(nt_read_6, df, by = "gene")
# 7 obj / 7 obj = 100.00%
nt_read_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  read = c(28.72, 65.22, 100.00, 100.00)
)

#----------SARC----------
nt_sarc<- nt %>% 
  filter(cancer_type_abbr == "SARC")
nt_sarc_3 <- nt_sarc %>% 
  filter(n_tool > 2)
nt_sarc_3_in <- merge(nt_sarc_3, df, by = "gene")
# 6 obj / 18 obj = 33.33%
nt_sarc_4 <- nt_sarc %>% 
  filter(n_tool > 3)
nt_sarc_4_in <- merge(nt_sarc_4, df, by = "gene")
# 3 obj / 4 obj = 75.00%
nt_sarc_5 <- nt_sarc %>% 
  filter(n_tool > 4)
nt_sarc_5_in <- merge(nt_sarc_5, df, by = "gene")
# 2 obj / 2 obj = 100.00% 
nt_sarc_6 <- nt_sarc %>% 
  filter(n_tool > 5)
nt_sarc_6_in <- merge(nt_sarc_6, df, by = "gene")
# 0 obj / 0 obj = NA
nt_sarc_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  sarc = c(33.33, 75.00, 100.00, NA)
)

#----------SKCM----------
nt_skcm<- nt %>% 
  filter(cancer_type_abbr == "SKCM")
nt_skcm_3 <- nt_skcm %>% 
  filter(n_tool > 2)
nt_skcm_3_in <- merge(nt_skcm_3, df, by = "gene")
# 48 obj / 366 obj = 13.11%
nt_skcm_4 <- nt_skcm %>% 
  filter(n_tool > 3)
nt_skcm_4_in <- merge(nt_skcm_4, df, by = "gene")
# 21 obj / 62 obj =  33.87%
nt_skcm_5 <- nt_skcm %>% 
  filter(n_tool > 4)
nt_skcm_5_in <- merge(nt_skcm_5, df, by = "gene")
# 7 obj / 7 obj = 100.00%
nt_skcm_6 <- nt_skcm %>% 
  filter(n_tool > 5)
nt_skcm_6_in <- merge(nt_skcm_6, df, by = "gene")
# 4 obj / 4 obj = 100.00%
nt_skcm_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  skcm = c(13.11, 33.87, 100.00, 100.00)
)

#----------STAD----------
nt_stad<- nt %>% 
  filter(cancer_type_abbr == "STAD")
nt_stad_3 <- nt_stad %>% 
  filter(n_tool > 2)
nt_stad_3_in <- merge(nt_stad_3, df, by = "gene")
# 37 obj / 324 obj = 11.42%
nt_stad_4 <- nt_stad %>% 
  filter(n_tool > 3)
nt_stad_4_in <- merge(nt_stad_4, df, by = "gene")
# 18 obj / 80 obj = 22.50%
nt_stad_5 <- nt_stad %>% 
  filter(n_tool > 4)
nt_stad_5_in <- merge(nt_stad_5, df, by = "gene")
# 15 obj / 23 obj = 65.22%
nt_stad_6 <- nt_stad %>% 
  filter(n_tool > 5)
nt_stad_6_in <- merge(nt_stad_6, df, by = "gene")
# 10 obj / 11 obj = 90.91%
nt_stad_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  stad = c(11.42, 22.50, 65.22, 90.91)
)

#----------TGCT----------
nt_tgct<- nt %>% 
  filter(cancer_type_abbr == "TGCT")
nt_tgct_3 <- nt_tgct %>% 
  filter(n_tool > 2)
nt_tgct_3_in <- merge(nt_tgct_3, df, by = "gene")
# 4 obj / 5 obj = 80.00%
nt_tgct_4 <- nt_tgct %>% 
  filter(n_tool > 3)
nt_tgct_4_in <- merge(nt_tgct_4, df, by = "gene")
# 3 obj / 3 obj = 100.00%
nt_tgct_5 <- nt_tgct %>% 
  filter(n_tool > 4)
nt_tgct_5_in <- merge(nt_tgct_5, df, by = "gene")
# 3 obj / 3 obj = 100.00%
nt_tgct_6 <- nt_tgct %>% 
  filter(n_tool > 5)
nt_tgct_6_in <- merge(nt_tgct_6, df, by = "gene")
# 2 obj / 2 obj = 100.00%
nt_tgct_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  tgct = c(80.00, 100.00, 100.00, 100.00)
)

#----------THCA----------
nt_thca<- nt %>% 
  filter(cancer_type_abbr == "THCA")
nt_thca_3 <- nt_thca %>% 
  filter(n_tool > 2)
nt_thca_3_in <- merge(nt_thca_3, df, by = "gene")
# 7 obj / 13 obj = 53.87%
nt_thca_4 <- nt_thca %>% 
  filter(n_tool > 3)
nt_thca_4_in <- merge(nt_thca_4, df, by = "gene")
# 3 obj / 4 obj = 75.00%
nt_thca_5 <- nt_thca %>% 
  filter(n_tool > 4)
nt_thca_5_in <- merge(nt_thca_5, df, by = "gene")
# 3 obj / 3 obj = 100.00%
nt_thca_6 <- nt_thca %>% 
  filter(n_tool > 5)
nt_thca_6_in <- merge(nt_thca_6, df, by = "gene")
# 3 obj / 3 obj = 100.00%
nt_thca_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  thca = c(53.87, 75.00, 100.00, 100.00)
)

#----------THYM----------
nt_thym<- nt %>% 
  filter(cancer_type_abbr == "THYM")
nt_thym_3 <- nt_thym %>% 
  filter(n_tool > 2)
nt_thym_3_in <- merge(nt_thym_3, df, by = "gene")
# 3 obj / 5 obj = 60.00%
nt_thym_4 <- nt_thym %>% 
  filter(n_tool > 3)
nt_thym_4_in <- merge(nt_thym_4, df, by = "gene")
# 1 obj / 3 obj = 33.33%
nt_thym_5 <- nt_thym %>% 
  filter(n_tool > 4)
nt_thym_5_in <- merge(nt_thym_5, df, by = "gene")
# 0 obj / 1 obj = 0.00%
nt_thym_6 <- nt_thym %>% 
  filter(n_tool > 5)
nt_thym_6_in <- merge(nt_thym_6, df, by = "gene")
# 0 obj / 1 obj = 0.00%
nt_thym_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  thym = c(60.00, 33.33, 0.00, 0.00)
)

#----------UCEC----------
nt_ucec <- nt %>% 
  filter(cancer_type_abbr == "UCEC")
nt_ucec_3 <- nt_ucec %>% 
  filter(n_tool > 2)
nt_ucec_3_in <- merge(nt_ucec_3, df, by = "gene")
# 90 obj / 719 obj = 12.52%
nt_ucec_4 <- nt_ucec %>% 
  filter(n_tool > 3)
nt_ucec_4_in <- merge(nt_ucec_4, df, by = "gene")
# 40 obj / 111 obj = 36.04%
nt_ucec_5 <- nt_ucec %>% 
  filter(n_tool > 4)
nt_ucec_5_in <- merge(nt_ucec_5, df, by = "gene")
# 18 obj / 22 obj = 81.82%
nt_ucec_6 <- nt_ucec %>% 
  filter(n_tool > 5)
nt_ucec_6_in <- merge(nt_ucec_6, df, by = "gene")
# 10 obj / 10 obj = 100.00%
nt_ucec_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  ucec = c(12.52, 36.04, 81.82, 100.00)
)

#-----------UCS----------
nt_ucs <- nt %>% 
  filter(cancer_type_abbr == "UCS")
nt_ucs_3 <- nt_ucs %>% 
  filter(n_tool > 2)
nt_ucs_3_in <- merge(nt_ucs_3, df, by = "gene")
# 14 obj / 20 obj = 70.00% 
nt_ucs_4 <- nt_ucs %>% 
  filter(n_tool > 3)
nt_ucs_4_in <- merge(nt_ucs_4, df, by = "gene")
# 6 obj / 6 obj = 100.00% 
nt_ucs_5 <- nt_ucs %>% 
  filter(n_tool > 4)
nt_ucs_5_in <- merge(nt_ucs_5, df, by = "gene")
# 5 obj / 5 obj = 100.00%
nt_ucs_6 <- nt_ucs %>% 
  filter(n_tool > 5)
nt_ucs_6_in <- merge(nt_ucs_6, df, by = "gene")
# 2 obj / 2 obj = 100.00% 
nt_ucs_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  ucs = c(70.00, 100.00, 100.00, 100.00)
)

#-----------UVM----------
nt_uvm <- nt %>% 
  filter(cancer_type_abbr == "UVM")
nt_uvm_3 <- nt_uvm %>% 
  filter(n_tool > 2)
nt_uvm_3_in <- merge(nt_uvm_3, df, by = "gene")
# 5 obj / 5 obj = 100.00%
nt_uvm_4 <- nt_uvm %>% 
  filter(n_tool > 3)
nt_uvm_4_in <- merge(nt_uvm_4, df, by = "gene")
# 5 obj / 5 obj = 100.00%
nt_uvm_5 <- nt_uvm %>% 
  filter(n_tool > 4)
nt_uvm_5_in <- merge(nt_uvm_5, df, by = "gene")
# 4 obj / 4 obj = 100.00% 
nt_uvm_6 <- nt_uvm %>% 
  filter(n_tool > 5)
nt_uvm_6_in <- merge(nt_uvm_6, df, by = "gene")
# 2 obj / 2 obj = 100.00% 
nt_uvm_su <- data.frame(
  n_tools = c(3, 4, 5, 6),
  uvm = c(100.00, 100.00, 100.00, 100.00)
)



#-----heatmap -------
nt_su <- merge(nt_acc_su, su, by = "n_tools")
nt_su <- merge(nt_blca_su, nt_su, by = "n_tools")
nt_su <- merge(nt_brca_su, nt_su, by = "n_tools")
nt_su <- merge(nt_cesc_su, nt_su, by = "n_tools")
nt_su <- merge(nt_chol_su, nt_su, by = "n_tools")
nt_su <- merge(nt_coad_su, nt_su, by = "n_tools")
nt_su <- merge(nt_dlbc_su, nt_su, by = "n_tools")
nt_su <- merge(nt_esca_su, nt_su, by = "n_tools")
nt_su <- merge(nt_gbm_su, nt_su, by = "n_tools")
nt_su <- merge(nt_hnsc_su, nt_su, by = "n_tools")
nt_su <- merge(nt_kich_su, nt_su, by = "n_tools")
nt_su <- merge(nt_kirc_su, nt_su, by = "n_tools")
nt_su <- merge(nt_kirp_su, nt_su, by = "n_tools")
nt_su <- merge(nt_laml_su, nt_su, by = "n_tools")
nt_su <- merge(nt_lgg_su, nt_su, by = "n_tools")
nt_su <- merge(nt_lihc_su, nt_su, by = "n_tools")
nt_su <- merge(nt_luad_su, nt_su, by = "n_tools")
nt_su <- merge(nt_lusc_su, nt_su, by = "n_tools")
nt_su <- merge(nt_meso_su, nt_su, by = "n_tools")
nt_su <- merge(nt_ov_su, nt_su, by = "n_tools")
nt_su <- merge(nt_paad_su, nt_su, by = "n_tools")
nt_su <- merge(nt_pcpg_su, nt_su, by = "n_tools")
nt_su <- merge(nt_prad_su, nt_su, by = "n_tools")
nt_su <- merge(nt_read_su, nt_su, by = "n_tools")
nt_su <- merge(nt_sarc_su, nt_su, by = "n_tools")
nt_su <- merge(nt_skcm_su, nt_su, by = "n_tools")
nt_su <- merge(nt_stad_su, nt_su, by = "n_tools")
nt_su <- merge(nt_tgct_su, nt_su, by = "n_tools")
nt_su <- merge(nt_thca_su, nt_su, by = "n_tools")
nt_su <- merge(nt_thym_su, nt_su, by = "n_tools")
nt_su <- merge(nt_ucec_su, nt_su, by = "n_tools")
nt_su <- merge(nt_ucs_su, nt_su, by = "n_tools")
nt_su <- merge(nt_uvm_su, nt_su, by = "n_tools")


nt_su <- nt_su %>% 
  mutate(n_tools = NULL,
         num = NULL)
nt_heatmap <- data.matrix(nt_su)
nt_heatmap <- t(nt_heatmap)
colnames(nt_heatmap) <- c(3, 4, 5, 6)
pheatmap(nt_heatmap, display_numbers = TRUE, cluster_rows = T, cluster_cols = F,
         show_rownames = T, show_colnames = T, cellheight = 10, cellwidth = 40)
#-----cancer type----------------------------

#ACC blca brca cesc chol coad dlbc esca gbm hnsc kich kirc kirp laml
#lgg lihc luad lusc meso ov paad pcpg prad read sarc skcm stad tgct
#thca thym ucec ucs uvm
#-----stastic of differ cancer number of genes define by n_tools------
nt_5_part1 <- nt_5[1:119, ]
nt_5_part2 <- nt_5[120:175, ]
nt_5_part3 <- nt_5[176:253, ]

ggplot(nt_5_part1) +
  geom_histogram(aes(x = n_tool))+
  facet_wrap(~cancer_type_abbr)+
  geom_text(stat = "count", aes(x = n_tool, label = ..count..))

ggplot(nt_5_part2) +
  geom_histogram(aes(x = n_tool))+
  facet_wrap(~cancer_type_abbr)+
  geom_text(stat = "count", aes(x = n_tool, label = ..count..))

ggplot(nt_5_part3) +
  geom_histogram(aes(x = n_tool))+
  facet_wrap(~cancer_type_abbr)+
  geom_text(stat = "count", aes(x = n_tool, label = ..count..))

# -> >= 5
#-----cancer type of -----------------------------------------------
nt_gtf2i <- nt %>% 
  filter(gene == "GTF2I")

nt_gng12 <- nt %>% 
  filter(gene == "GNG12")

#-----functional enrichment-------------------
df_vector <- df_5_out_ct %>% 
  mutate(count = NULL)

#### Configure 

rm(list=ls())

#comPATH='/media/md1200'

casePATH='/home/leotsai0127/functional'  #output path

outfolder1=file.path(casePATH,'fun.ann_35gene')  #output folder

environmentDir=file.path(comPATH,'/analysis/Script/function-annotation/environment/')

#### Choose your species 

source(file.path(environmentDir,'code','annotation_human.r'))  #Human (Cytoscape)

#source(file.path(environmentDir,'code','annotation_mouse.r')) #Mouse (Cytoscape)

#source(file.path(environmentDir,'code','annotation_rat.r')) #Rat (Cytoscape)

#### data input 

genelist = df_vector$gene
genelist = as.vector(genelist)

#### function annotation 

if (!file.exists(outfolder1)) dir.create(outfolder1)

setwd(outfolder1)

annotation(genelist,environmentDir)

#-----mutation of gng12----------------------------

mut_GNG12 <- read.delim("/home/leotsai0127/Mut.data_GNG12.txt")

gng12_chr <- mut_GNG12 %>%
  group_by(pos) %>% 
  summarise(pos_count = length(pos))

ggplot(gng12_chr)+
  geom_line(aes(x = pos, y = pos_count))+
  geom_point(aes(x = pos, y = pos_count))+
  geom_label(aes(x = pos, y = pos_count, label = pos))

gng12_cds <- mut_GNG12 %>% 
  group_by(cds_position) %>% 
  summarise(pos_count = length(cds_position)) %>% 
  filter(cds_position != "") %>% 
  mutate(cds_position = as.character(cds_position))

gng12_cds <- gng12_cds %>% 
  mutate(cds_position = ifelse(cds_position == "166-167", 166.5,
                               ifelse(cds_position == "186-187", 186.5,
                                      ifelse(cds_position == "202-203", 202.5, as.numeric(cds_position)))))

gng12_cds$cds_position <- as.numeric(gng12_cds$cds_position)

ggplot(gng12_cds)+
  geom_line(aes(x = cds_position, y = pos_count))+
  geom_point(aes(x = cds_position, y = pos_count))+
  geom_label(aes(x = cds_position, y = pos_count, label = cds_position), size = 3)

#View(mut_GNG12)

gng12_pro <- mut_GNG12 %>% 
  group_by(protein_position) %>% 
  summarise(pos_count = length(protein_position)) %>% 
  filter(protein_position != "") %>% 
  mutate(protein_position = as.character(protein_position))

gng12_pro <- gng12_pro %>% 
  mutate(protein_position = ifelse(protein_position == "62-63", 62.5, as.numeric(protein_position)))

ggplot(gng12_pro)+
  geom_col(aes(x = protein_position, y = pos_count)) 

gng12_pro2 <- mut_GNG12 %>% 
  select(cancer_project, protein_position) %>% 
  group_by(cancer_project, protein_position) %>% 
  summarise(pos_count = n()) %>% 
  filter(cancer_project != "",
         protein_position != "") %>% 
  mutate(protein_position = as.character(protein_position))

gng12_pro2 <- gng12_pro2 %>% 
  mutate(protein_position = ifelse(protein_position == "62-63", "62.5",
                                   ifelse(protein_position == "1", "01",
                                          ifelse(protein_position == "2", "02",
                                                 ifelse(protein_position == "7", "07",
                                                        as.character(protein_position))))))

ggplot(gng12_pro2)+
  geom_bar(aes(x = protein_position, y = pos_count), stat = "identity", fill = "navyblue")

gng12_pro2 <- merge(gng12_pro2, fn, by = "cancer_project")

gng12_pro2 <- gng12_pro2 %>% 
  filter(protein_position == "68")

gng12_pro2$cancer_type_abbr <- factor(gng12_pro2$cancer_type_abbr, levels = c("STAD", "COAD", "UCEC", "LGG", "KIRP", "CESC", "HNSC"))

ggplot(gng12_pro2)+
  geom_bar(aes(x = cancer_type_abbr, y = pos_count), stat = "identity", fill = "navyblue")


#-----mutation of gtf2i---------------------------
mut_GTF2I <- read.delim("/home/leotsai0127/Mut.data_GTF2I.txt")

gtf2i_chr <- mut_GTF2I %>% 
  group_by(pos) %>% 
  summarise(pos_count = length(pos))

ggplot(gtf2i_chr) +
  geom_point(aes(x = pos, y = pos_count)) +
  geom_line(aes(x = pos, y = pos_count))

gtf2i_cds <- mut_GTF2I %>% 
  group_by(cds_position) %>% 
  summarise(pos_count = length(cds_position)) %>% 
  filter(cds_position != "")

ggplot(gtf2i_cds) +
  geom_point(aes(x = cds_position, y = pos_count)) +
  geom_line(aes(x = cds_position, y = pos_count))

#-----mutation of epha2------------------
mut_EPHA2 <- read.delim("/home/leotsai0127/exp.data_EPHA2.txt")

epha2_chr <- mut_EPHA2 %>%
  group_by(pos) %>% 
  summarise(pos_count = length(pos))

ggplot(epha2_chr)+
  geom_line(aes(x = pos, y = pos_count))+
  geom_point(aes(x = pos, y = pos_count))

epha2_cds <- mut_EPHA2 %>% 
  group_by(cds_position) %>% 
  summarise(pos_count = length(cds_position)) %>% 
  filter(cds_position != "")

ggplot(epha2_cds)+
  geom_line(aes(x = cds_position, y = pos_count))+
  geom_point(aes(x = cds_position, y = pos_count))

epha2_pro <- mut_EPHA2 %>% 
  group_by(protein_position) %>% 
  summarise(pos_count = length(protein_position)) %>% 
  filter(protein_position != "") 

ggplot(epha2_pro)+
  geom_point(aes(x = protein_position, y = pos_count))+
  geom_line(aes(x = protein_position, y = pos_count))

#-----mutation of lama1------------------
mut_LAMA1 <- read.delim("/home/leotsai0127/exp.data_LAMA1.txt")

lama1_chr <- mut_LAMA1 %>%
  group_by(pos) %>% 
  summarise(pos_count = length(pos))

ggplot(lama1_chr)+
  geom_line(aes(x = pos, y = pos_count))+
  geom_point(aes(x = pos, y = pos_count))

lama1_cds <- mut_LAMA1 %>% 
  group_by(cds_position) %>% 
  summarise(pos_count = length(cds_position)) %>% 
  filter(cds_position != "")

ggplot(lama1_cds)+
  geom_line(aes(x = cds_position, y = pos_count))+
  geom_point(aes(x = cds_position, y = pos_count))

lama1_pro <- mut_LAMA1 %>% 
  group_by(protein_position) %>% 
  summarise(pos_count = length(protein_position)) %>% 
  filter(protein_position != "") 

ggplot(lama1_pro)+
  geom_point(aes(x = protein_position, y = pos_count))+
  geom_line(aes(x = protein_position, y = pos_count))

#-----mutation of sos1-------------------
mut_SOS1 <- read.delim("/home/leotsai0127/exp.data_SOS1.txt")

sos1_chr <- mut_SOS1 %>%
  group_by(pos) %>% 
  summarise(pos_count = length(pos))

ggplot(sos1_chr)+
  geom_line(aes(x = pos, y = pos_count))+
  geom_point(aes(x = pos, y = pos_count))

sos1_cds <- mut_SOS1 %>% 
  group_by(cds_position) %>% 
  summarise(pos_count = length(cds_position)) %>% 
  filter(cds_position != "")

ggplot(sos1_cds)+
  geom_line(aes(x = cds_position, y = pos_count))+
  geom_point(aes(x = cds_position, y = pos_count))

sos1_pro <- mut_SOS1 %>% 
  group_by(protein_position) %>% 
  summarise(pos_count = length(protein_position)) %>% 
  filter(protein_position != "") 

ggplot(sos1_pro)+
  geom_point(aes(x = protein_position, y = pos_count))+
  geom_line(aes(x = protein_position, y = pos_count))

#-----twelve tools name
ggplot(cgt) + 
  stat_count(aes(x = tool))
## activedriver comet dawnrank drivernet
## dendrix driverml e_driver ipac
## memo msea mutex mutsigcv 
## netbox oncodriveclust
#-----observation of GNG12-----
ob1 <- mut_GNG12 %>% 
  select(protein_position, type) %>% 
  filter(protein_position != "") %>% 
  mutate(protein_position = ifelse(protein_position == "62-63", "62.5",
                                   ifelse(protein_position == "1", "01",
                                          ifelse(protein_position == "2", "02",
                                                 ifelse(protein_position == "7", "07",
                                                        as.character(protein_position))))))

ggplot(ob1) +
  stat_count(aes(x = protein_position, fill = type))
#----------------------