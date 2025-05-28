library(tidyverse)

# Reads BD List
bd_list <- read_csv("/Users/mayukh/Desktop/Telomere genes update/CBC RRP Microarray in csv.csv", na = c("NA", "","N.A.", "N/A"))

# Renames RRP6-HL to RRP6_HL
bd_list <- rename(bd_list, "RRP6_HL" = "RRP6-HL")

# Reads Kupiec list
kupiec_list <- read_csv("/Users/mayukh/Desktop/Telomere genes update/Kupiec list in csv.csv", na = c("NA", "","N.A.", "N/A"))

#removes the blank line
kupiec_list <- kupiec_list[-1,]

#Fixing the incorrect systematic name against PRS3
kupiec_list <- kupiec_list %>% mutate(Sys_name = replace(Sys_name, Gene_name == "PRS3", "YHL011C"))

#Fixing MTC7 gene name for YEL033W
kupiec_list <- kupiec_list %>% mutate(Gene_name = replace(Gene_name, Sys_name == "YEL033W", "MTC7"))

# Removes the Puddu entry of YOL138C
kupiec_list <- kupiec_list %>% filter(!(kupiec_list$Sys_name == "YOL138C" & kupiec_list$Telomere_phenotype == "ss"))

# Copies bd_list to bd_folds
bd_folds <- bd_list

# Calculates and adds the fold changes
bd_folds <- bd_folds %>% mutate(CBC1_SSL_Fold = (CBC_SSL / WT_SSL), .after = CBC_SSL) %>% 
  mutate(CBC1_HL_Fold = (CBC1_HL / WT_HL), .after = CBC1_HL) %>% 
  mutate(RRP6_SSL_Fold = (RRP6_SSL / WT_SSL), .after = RRP6_SSL) %>% 
  mutate(RRP6_HL_Fold = (RRP6_HL / WT_HL), .after = RRP6_HL)

# Makes copy of Kupiec List
kupiec_merger <- kupiec_list

# Merges kupiec_merger & bd_folds
kupiec_merged <- kupiec_merger |> 
  left_join(bd_folds |> 
               select(sys_name, CBC1_SSL_Fold, CBC1_HL_Fold, RRP6_SSL_Fold, RRP6_HL_Fold), 
             join_by(Sys_name == sys_name))

# Relocates fold columns to beside telomere phenotype
kupiec_merged <- relocate(kupiec_merged, CBC1_SSL_Fold : RRP6_HL_Fold, .after = Telomere_phenotype)

# Creates CBC_SSL_2Fold
cbc_ssl_2fold <- filter(kupiec_merged, CBC1_SSL_Fold >= 1.9)

# Creates RRP6_SSL_2Fold
rrp6_ssl_2fold <- filter(kupiec_merged, RRP6_SSL_Fold >= 1.9)

# Merges cbc2fold with rrp62fold
cbc_rrp6_ssl_2F <- cbc_ssl_2fold |> inner_join(rrp6_ssl_2fold, join_by(Sys_name))


# Copying for temporary experimentation:
kupiec_temp <- kupiec_merged

# Kupiec long phenotype:
kupiec_long <- kupiec_merged |>
  filter(Telomere_phenotype == "sl" | Telomere_phenotype == "L" | Telomere_phenotype == "VL" | Telomere_phenotype == "DAmP Long")

# Kupiec short phenotype:
kupiec_short <- kupiec_merged |>
  filter(Telomere_phenotype == "ss" | Telomere_phenotype == "S" | Telomere_phenotype == "VS" | Telomere_phenotype == "DAmP Short")

# Filtering short by RRP6 2 fold increase
kupiec_short_RRP6_2f <- kupiec_short |>
  filter(RRP6_SSL_Fold > 1.7)

# Filtering long by RRP6 2 fold increase
kupiec_long_RRP6_2f <- kupiec_long |>
  filter(RRP6_SSL_Fold > 1.7)

# Filtering short by CBC1 2 fold increase
kupiec_short_CBC1_2f <- kupiec_short |>
  filter(CBC1_SSL_Fold > 1.7)

# Filtering long by CBC1 2 fold increase
kupiec_long_CBC1_2f <- kupiec_long |>
  filter(CBC1_SSL_Fold > 1.7)



####### ----- Experimentation ------ ########
# Grouping by:
kupiec_grouped <- kupiec_temp |>
  group_by(Telomere_phenotype == "sl" | Telomere_phenotype == "L" | Telomere_phenotype == "VL" | Telomere_phenotype == "DAmP Long")

# Rename:
kupiec_grouped <- rename(kupiec_grouped, "group" = "|...")

# Filtering CBC from Grouped:
kupiec_grouped_CBC <- kupiec_grouped |>
  select(-RRP6_SSL_Fold, -RRP6_HL_Fold)

# Visualise:
library(ggrepel)
genes_to_label <- c("EST2", "MTC7", "TEL1", "SRB5", "RIF1", "RIF2", "MLH1")

# LogFC transformation of SSL values and labelling increase and decrease
kupiec_grouped <- kupiec_grouped %>% mutate(Tel_length = factor(ifelse(group == "TRUE", "Decrease", "Increase")))
kupiec_grouped <- kupiec_grouped %>% mutate(LogFCcbc = log2(CBC1_SSL_Fold))
kupiec_grouped <- kupiec_grouped %>% mutate(LogFCrrp = log2(RRP6_SSL_Fold))

# Plot for CBC SSL:
kupiec_grouped |>
  ggplot(mapping = aes(y = CBC1_SSL_Fold, x = Gene_name))+
  geom_point(aes(color = group))+
  geom_hline(yintercept = c(1.5, 0.67), color = "blue", linetype = "dashed", size = 1)+
  scale_x_discrete(expand = expansion(add = 2))+
  scale_y_continuous(trans = "log10", breaks = c(0.1, 0.25, 0.5, 0.75, 1, 2.5, 5, 7.5, 10))+
  geom_text_repel(aes(label = ifelse(Gene_name %in% genes_to_label, Gene_name, "")),
            color = "purple", size = 4.2, box.padding = 0.5, point.padding = 0.0,
            segment.color = "black", segment.size = 0.5, nudge_y = 0.15)+
  scale_color_discrete(name = "tel phenotype", labels = c
                       ("short", "long"))

#Experimenting with CBC SSL above for logFC values
kupiec_grouped |>
  ggplot(mapping = aes(y = LogFCcbc, x = Gene_name))+
  geom_point(aes(color = Tel_length))+
  geom_hline(yintercept = c(0.5, -0.5), color = "blue", linetype = "dashed", size = 0.8)+
  scale_x_discrete(expand = expansion(add = 4))+
  scale_y_continuous(limits = c(-2.5, 2.5), breaks = c(-2, -1, 0, 1, 2))+
  geom_text_repel(aes(label = ifelse(Gene_name %in% genes_to_label, Gene_name, "")),
                  color = "purple", size = 4.2, box.padding = 0.8, point.padding = 0.0,
                  segment.color = "black", segment.size = 0.5, nudge_y = 0.15)+
  scale_color_discrete(name = "tel phenotype", labels = c
                       ("short", "long"))


# Plot for RRP6 SSL:
kupiec_grouped |>
  ggplot(mapping = aes(y = RRP6_SSL_Fold, x = Gene_name))+
  geom_point(aes(color = group))+
  geom_hline(yintercept = c(1.5, 0.67), color = "blue", linetype = "dashed", size = 1)+
  scale_x_discrete(expand = expansion(add = 2))+
  scale_y_continuous(trans = "log10", breaks = c(0.1, 0.25, 0.5, 0.75, 1, 2.5, 5, 7.5, 10))+
  geom_text_repel(aes(label = ifelse(Gene_name %in% genes_to_label, Gene_name, "")),
                  color = "purple", size = 4.2, box.padding = 0.5, point.padding = 0.0,
                  segment.color = "black", segment.size = 0.5, nudge_y = 0.15)+
  scale_color_discrete(name = "tel phenotype", labels = c
                       ("short", "long"))

#Experimenting with RRP SSL above for logFC values
kupiec_grouped |>
  ggplot(mapping = aes(y = LogFCrrp, x = Gene_name))+
  geom_point(aes(color = Tel_length))+
  geom_hline(yintercept = c(0.5, -0.5), color = "blue", linetype = "dashed", size = 0.8)+
  scale_x_discrete(expand = expansion(add = 4))+
  scale_y_continuous(limits = c(-2.5, 4), breaks = c(-3, -2, -1, 0, 1, 2, 3, 4))+
  geom_text_repel(aes(label = ifelse(Gene_name %in% genes_to_label, Gene_name, "")),
                  color = "purple", size = 4.2, box.padding = 0.8, point.padding = 0.0,
                  segment.color = "black", segment.size = 0.5, nudge_y = 0.15)+
  scale_color_discrete(name = "tel phenotype", labels = c
                       ("short", "long"))

# Experimental:

# kupiec_grouped |>
#   ggplot(aes(y = RRP6_SSL_Fold, x = Gene_name, label = ifelse(Gene_name %in% genes_to_label, Gene_name, "", color = group)))+
#   geom_point()+
#   scale_x_discrete(expand = expansion(add = 2))+
#   scale_y_continuous(trans = "log10")+
#   geom_text(aes(label = ifelse(Gene_name %in% genes_to_label, Gene_name, "")),
#                   size = 4.2, box.padding = 0.7, point.padding = 0.0,
#                   segment.color = "black", segment.size = 0.5, nudge_y = 0.1)+
#   scale_color_discrete(name = "tel phenotype", labels = c
#                        ("short", "long"))

specific_genes <- kupiec_grouped %>% filter(Gene_name %in% genes_to_label)

# Density plot for RRP6
kupiec_grouped |>
  ggplot(aes(x = LogFCrrp, fill = Tel_length))+
  geom_density(alpha = 0.4, bw = 0.4)+
  # geom_point(data = specific_genes, aes(y=0, color = group), size = 3)+
  # geom_text_repel(data = specific_genes, aes(y=0, label = Gene_name, colour = group),
  #                 box.padding = 3.2, point.padding = 0.5, segment_color = "gray")+
  scale_x_continuous(limits = c(-3,3))+
  scale_fill_manual(values = c("Increase" = "#ff0000", "Decrease" = "#0000FF"))+
  labs(title = "rrp6 density", x = "logFC", y = "Density")+
  theme_minimal()



# Density plot for CBC1
kupiec_grouped |>
  ggplot(aes(x = LogFCcbc, fill = Tel_length))+
  geom_density(alpha = 0.4, bw = 0.3)+
  # geom_point(data = specific_genes, aes(y=0, color = group), size = 3)+
  # geom_text_repel(data = specific_genes, aes(y=0, label = Gene_name, colour = group),
  #                 box.padding = 3.2, point.padding = 0.5, segment_color = "gray")+
  scale_x_continuous(limits = c(-3,3))+
  scale_fill_manual(values = c("Increase" = "#FF0000", "Decrease" = "#0000FF"))+
  labs(title = "cbc1 density", x = "logFC", y = "Density")+
  theme_minimal()


# kupiec_grouped |>
#   ggplot(mapping = aes(y = CBC1_SSL_Fold, x = Gene_name))+
#   geom_point(aes(color = group))+
#   scale_y_continuous(trans = "log10")+
#   geom_text(aes(label = ifelse(CBC1_SSL_Fold < 0.5, Gene_name, '')),
#             vjust = -1, hjust = 0.85, color = "purple")+
#   scale_color_discrete(name = "tel phenotype", labels = c
#                       ("short", "long"))

## filter(CBC1_SSL_Fold > 1.5) |>

# DBP2 experiment
# venn_list <- read_csv("/Users/mayukh/Downloads/cbc1-rrp6-dbp2 mRNAs.csv")
# venn_copy <- venn_list
# 
# venn_copy <- rename(venn_copy, "cbc1_rrp6_dbp2"="YPL164C...1")
# venn_copy$YPL164C...2 <- NULL
# venn_copy$YPL164C...3 <- NULL
# 
# venn_merged <- venn_copy |>
#   left_join(bd_folds |>
#               select(Symbol, sys_name, CBC1_SSL_Fold, CBC1_HL_Fold, RRP6_SSL_Fold, RRP6_HL_Fold), 
#             join_by(cbc1_rrp6_dbp2 == sys_name))
# 
# venn_merged_copy <- venn_merged
# venn_merged_copy_merged <- venn_merged_copy |>
#   left_join(kupiec_merger |>
#               select(Sys_name, Gene_name, Function),
#                      join_by(cbc1_rrp6_dbp2 == Sys_name))
# 
# venn_merged_copy_merged$Gene_name <- NULL
# 
# bd_func <- read_csv("/Users/mayukh/Desktop/Telomere genes update/CBC1RRP6 KO MicroarrayTable Latest.csv", na = c("NA", "","N.A.", "N/A"))
# 
# venn_merged_copy_merged <- venn_merged_copy_merged[-c(32, 33), ]
# 
# venn_merged_final <- venn_merged_copy_merged |>
#   left_join(bd_func|>
#               select(SGD, Descriptions),
#             join_by(cbc1_rrp6_dbp2 == SGD))
# 
# venn_merged_final <- venn_merged_final[-c(31,32), ]
# venn_merged_final$Function <- NULL
# 
# write_csv(venn_merged_final, "/Users/mayukh/Desktop/Telomere genes update/cbc-rrp-dbp mRNAs.csv")

# Joining BY4741 and W303 to rrp6 SSL 2-fold:
# Importing BY and W3 csv
BY4741EtOH <- read_csv("./BY4741 ethanol.csv", na = c("NA", "","N.A.", "N/A"))
W303EtOH <- read_csv("./w303 ethanol.csv", na = c("NA", "","N.A.", "N/A"))

# Add Linear Fold Change column
BY4741EtOH <- BY4741EtOH |> mutate(Linear_fold_BY1 = 2^Log2FC, .before = Log2FC)
W303EtOH <- W303EtOH |> mutate(Linear_fold_W3 = 2^Log2FC, .before = Log2FC)

# Merging linear fold columns to rrp6 ssl 2 fold
rrp6_final <- rrp6_ssl_2fold |> left_join(BY4741EtOH |> select(Linear_fold_BY1, Sys_name),
                                          join_by(Sys_name == Sys_name))

rrp6_final <- rrp6_final |> left_join(W303EtOH |> select(Linear_fold_W3, Sys_name),
                                      join_by(Sys_name == Sys_name))


rrp6_final <- rrp6_final |> relocate(Linear_fold_BY1:Linear_fold_W3, .before = Function)



## Further experiments:

special_messages_cbc_ssl_2f <- filter(bd_folds, CBC1_SSL_Fold >= 1.98)
special_messages_rrp_ssl_2f <- filter(bd_folds, RRP6_SSL_Fold >= 1.98)

special_messages_cbc_hl_2f <- filter(bd_folds, CBC1_HL_Fold >= 1.98)
special_messages_rrp_hl_2f <- filter(bd_folds, RRP6_HL_Fold >= 1.98)

special_cbc_rrp_ssl_2f <- filter(bd_folds, CBC1_SSL_Fold >= 1.98 & RRP6_SSL_Fold >=1.98)
special_cbc_rrp_hl_2f <- filter(bd_folds, CBC1_HL_Fold >=1.98 & RRP6_HL_Fold >=1.98)
