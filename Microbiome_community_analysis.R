# Import source data
source("PATH/TO/Function_source.R")
# Import end


##########
##########
##########
########## RAREFACTION
##########
##########
##########

#####
### OTU
#####
Sample_file <- "PATH/TO/Sample_metadata_1.csv"
Data_file <- "PATH/TO/Strain_1.csv"
Out_file <- "PATH/TO/Rarefy_std_Strain_1.txt"

Sample_meta <- read.csv(Sample_file, header = T, row.names = 1)
Data <- read.csv(Data_file, header = T, row.names = 1)

# Output Frequency data
Data_freq <- DF_FREQ(Data)
write.csv(Data_freq, "PATH/TO/Strain_freq_1.csv")

# rarefaction start
## Calculate max slope and rarefied (Coverage-based rarefaction)
max_list <- CALCMAXSLOPE(Data, Out_file)
Data_rarefied <- RAREFACTION(Data, max_list)
Data_rarefied <- DF_NONZERO(Data_rarefied)
write.csv(Data_rarefied, "PATH/TO/Strain_rarefied_1.csv")

Data_rarefied_freq <- DF_FREQ(Data_rarefied)
write.csv(Data_rarefied_freq, "PATH/TO/Strain_rarefied_freq_1.csv")

## Summary data for rarefaction
df_summary_raw <- DF_SUMMARY(Data)
df_summary_rarefied <- DF_SUMMARY(Data_rarefied)
colnames(df_summary_rarefied) <- c("Sum_values_rarefied", "Mean_values_rarefied", "Nonzero_values_rarefied")
df_summary <- data.frame(df_summary_raw, df_summary_rarefied)
write.csv(df_summary, "PATH/TO/Rarefy_summary_Strain_1.csv")


#####
### Genus
#####
Sample_file <- "PATH/TO/Sample_metadata_3.csv"
Data_file <- "PATH/TO/Genus_3.csv"
Out_file <- "PATH/TO/Rarefy_std_Strain_3.txt"

Sample_meta <- read.csv(Sample_file, header = T, row.names = 1)
Data <- read.csv(Data_file, header = T, row.names = 1)

# Output Frequency data
Data_freq <- DF_FREQ(Data)
write.csv(Data_freq, "PATH/TO/Genus_freq_3.csv")

# rarefaction start
## Calculate max slope and rarefied (Coverage-based rarefaction)
max_list <- CALCMAXSLOPE(Data, Out_file) #0.998430855953655
Data_rarefied <- RAREFACTION(Data, max_list)
Data_rarefied <- DF_NONZERO(Data_rarefied)
write.csv(Data_rarefied, "PATH/TO/Genus_rarefied_3.csv")

Data_rarefied_freq <- DF_FREQ(Data_rarefied)
write.csv(Data_rarefied_freq, "PATH/TO/Genus_rarefied_freq_3.csv")

## Summary data for rarefaction
df_summary_raw <- DF_SUMMARY(Data)
df_summary_rarefied <- DF_SUMMARY(Data_rarefied)
colnames(df_summary_rarefied) <- c("Sum_values_rarefied", "Mean_values_rarefied", "Nonzero_values_rarefied")
df_summary <- data.frame(df_summary_raw, df_summary_rarefied)
write.csv(df_summary, "PATH/TO/Rarefy_summary_Genus_3.csv")

#####
### Family
#####
Sample_file <- "PATH/TO/Sample_metadata_4.csv"
Data_file <- "PATH/TO/Family_4.csv"
Out_file <- "PATH/TO/Rarefy_std_Strain_4.txt"


Sample_meta <- read.csv(Sample_file, header = T, row.names = 1)
Data <- read.csv(Data_file, header = T, row.names = 1)

# Output Frequency data
Data_freq <- DF_FREQ(Data)
write.csv(Data_freq, "PATH/TO/Family_freq_4.csv")

# rarefaction start
## Calculate max slope and rarefied (Coverage-based rarefaction)
max_list <- CALCMAXSLOPE(Data, Out_file) #0.999487619204838
Data_rarefied <- RAREFACTION(Data, max_list)
Data_rarefied <- DF_NONZERO(Data_rarefied)
write.csv(Data_rarefied, "PATH/TO/Family_rarefied_4.csv")

Data_rarefied_freq <- DF_FREQ(Data_rarefied)
write.csv(Data_rarefied_freq, "PATH/TO/Family_rarefied_freq_4.csv")

## Summary data for rarefaction
df_summary_raw <- DF_SUMMARY(Data)
df_summary_rarefied <- DF_SUMMARY(Data_rarefied)
colnames(df_summary_rarefied) <- c("Sum_values_rarefied", "Mean_values_rarefied", "Nonzero_values_rarefied")
df_summary <- data.frame(df_summary_raw, df_summary_rarefied)
write.csv(df_summary, "PATH/TO/Rarefy_summary_Family_4.csv")

#####
### ASV
#####
Sample_file <- "PATH/TO/Sample_metadata_5.csv"
Data_file <- "PATH/TO/ASV_5.csv"
Out_file <- "PATH/TO/Rarefy_std_Strain_5.txt"

Sample_meta <- read.csv(Sample_file, header = T, row.names = 1)
Data <- read.csv(Data_file, header = T, row.names = 1)
ROWCHECK(Sample_meta, Data)
Sample_meta[Sample_meta$Type == "Soil", ]
Df_temp <- DF_ROWFILT(Data, Sample_meta, (Sample_meta$Sampling_day != 3 & !(Sample_meta$Sampling_day == 0 & Sample_meta$Condition == "Salt")))
Data2 <- Df_temp[1][[1]]
Sample_meta2 <- Df_temp[2][[1]]
Data2 <- DF_NONZERO(Data2)
nrow(Data2)

write.csv(Data2, "PATH/TO/ASV_5.csv")
write.csv(Sample_meta2, "PATH/TO/Sample_metadata_5.csv")

Sample_file <- "PATH/TO/Sample_metadata_5.csv"
Data_file <- "PATH/TO/ASV_5.csv"
Out_file <- "PATH/TO/Rarefy_std_Strain_5.txt"

Sample_meta <- read.csv(Sample_file, header = T, row.names = 1)
Data <- read.csv(Data_file, header = T, row.names = 1)

# Output Frequency data
Data_freq <- DF_FREQ(Data)
write.csv(Data_freq, "PATH/TO/ASV_freq_5.csv")

# rarefaction start
## Calculate max slope and rarefied (Coverage-based rarefaction)
max_list <- CALCMAXSLOPE(Data, Out_file) #0.833417583977742
Data_rarefied <- RAREFACTION(Data, max_list)
Data_rarefied <- DF_NONZERO(Data_rarefied)
write.csv(Data_rarefied, "PATH/TO/ASV_rarefied_5.csv")

Data_rarefied_freq <- DF_FREQ(Data_rarefied)
write.csv(Data_rarefied_freq, "PATH/TO/ASV_rarefied_freq_5.csv")

## Summary data for rarefaction
df_summary_raw <- DF_SUMMARY(Data)
df_summary_rarefied <- DF_SUMMARY(Data_rarefied)
colnames(df_summary_rarefied) <- c("Sum_values_rarefied", "Mean_values_rarefied", "Nonzero_values_rarefied")
df_summary <- data.frame(df_summary_raw, df_summary_rarefied)
write.csv(df_summary, "PATH/TO/Rarefy_summary_ASV_5.csv")

###
###
###
###
###
### Rarecurve plot
outdir <- "PATH/TO/"
#OTU
Sample_file <- "PATH/TO/Sample_metadata_1.csv"
Data_file <- "PATH/TO/Strain_1.csv"

Sample_meta <- read.csv(Sample_file, header = T, row.names = 1)
Data <- read.csv(Data_file, header = T, row.names = 1)
df_summary_raw <- DF_SUMMARY(Data)

col_list <- DF_CATE2NUM(Sample_meta$Type)
RARECURVE_PLOT(Data, df_summary_raw,
               out_name = Lc(outdir, "rarecurve_otu"),
               col_list = col_list, showonly = F, step = 1000)

# Genus
Sample_file <- "PATH/TO/Sample_metadata_3.csv"
Data_file <- "PATH/TO/Genus_3.csv"

Sample_meta <- read.csv(Sample_file, header = T, row.names = 1)
Data <- read.csv(Data_file, header = T, row.names = 1)
df_summary_raw <- DF_SUMMARY(Data)

col_list <- DF_CATE2NUM(Sample_meta$Type)
RARECURVE_PLOT(Data, df_summary_raw,
               out_name = Lc(outdir, "rarecurve_genus"),
               col_list = col_list, showonly = F, step = 1000)

# Family
Sample_file <- "PATH/TO/Sample_metadata_5.csv"
Data_file <- "PATH/TO/Family_4.csv"

Sample_meta <- read.csv(Sample_file, header = T, row.names = 1)
Data <- read.csv(Data_file, header = T, row.names = 1)
df_summary_raw <- DF_SUMMARY(Data)

col_list <- DF_CATE2NUM(Sample_meta$Type)
RARECURVE_PLOT(Data, df_summary_raw,
               out_name = Lc(outdir, "rarecurve_family"),
               col_list = col_list, showonly = F, step = 1000)


# ASV
Sample_file <- "PATH/TO/Sample_metadata_5.csv"
Data_file <- "PATH/TO/ASV_5.csv"

Sample_meta <- read.csv(Sample_file, header = T, row.names = 1)
Data <- read.csv(Data_file, header = T, row.names = 1)
df_summary_raw <- DF_SUMMARY(Data)

col_list <- DF_CATE2NUM(Sample_meta$Type)
RARECURVE_PLOT(Data, df_summary_raw,
               out_name = Lc(outdir, "rarecurve_asv"),
               col_list = col_list, showonly = F, step = 1000)












##########
##########
##########
########## ALPHA DIVERSITY
##########
##########
##########

Sample_file <- "PATH/TO/Sample_metadata_1.csv"
Data_file <- "PATH/TO/Strain_rarefied_freq_1.csv"
Out_file <- "PATH/TO/Alpha_out.txt"
Summary_file <- "PATH/TO/Rarefy_summary_Strain_1.csv"

Sample_meta <- read.csv(Sample_file, header = T, row.names = 1)
Data <- read.csv(Data_file, header = T, row.names = 1)
ROWCHECK(Sample_meta, Data)

### Alpha diversity
# Calculate diversity
Sample_meta <- CALC_ALPHA(Data, Sample_meta)

# For Plot
## Order change
Sample_meta <- ORDER_CHANGE(Sample_meta)
Sample_meta <- Sample_meta[Sample_meta$Inoculant != "CONT", ]

## Determine xlim and ylim
Sample_meta <- Sample_meta[Sample_meta$Inoculant != "CONT", ]
y_diff <- max(Sample_meta$Shannon_index) - min(Sample_meta$Shannon_index)
y_min <- min(Sample_meta$Shannon_index) - y_diff*0.02
y_max <- max(Sample_meta$Shannon_index) + y_diff*0.02

Sample_meta_root <- Sample_meta[Sample_meta$Type == "Root", ]
Sample_meta_soil <- Sample_meta[Sample_meta$Type == "Soil", ]

## Plot root samples
p_alpha1 <- Sample_meta_root %>%
  ggplot(aes(x = Inoculant, y =  Shannon_index, fill = Condition)) +#
  geom_boxplot(outlier.colour = NA, position = position_dodge2(preserve = "single", width = 0.5)) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_point(aes(x=Inoculant, y=Shannon_index, color=Condition),
             position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0), size=0.3) +
  theme(legend.position = "right", plot.title = element_text(size = 8))+
  ggtitle("") +
  labs(x = "", y = "A diversity [Shannon's H]")+
  ylim(y_min,y_max)

p_alpha2 <- Sample_meta_root %>%
  ggplot(aes(x = Host, y =  Shannon_index, fill = Condition)) +#
  geom_boxplot(outlier.colour = NA, position = position_dodge2(preserve = "single", width = 0.5)) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_point(aes(x=Host, y=Shannon_index, color=Condition),
             position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0), size=0.3) +
  theme(legend.position = "right", plot.title = element_text(size = 8))+
  ggtitle("") +
  labs(x = "", y = "A diversity [Shannon's H]")+
  ylim(y_min,y_max)

p_alpha3 <- Sample_meta_root %>%
  ggplot(aes(x = Host, y =  Shannon_index, fill = Inoculant)) +#
  geom_boxplot(outlier.colour = NA, position = position_dodge2(preserve = "single", width = 0.5)) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_point(aes(x=Host, y=Shannon_index, color=Inoculant),
             position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0), size=0.3) +
  theme(legend.position = "right", plot.title = element_text(size = 8))+
  ggtitle("") +
  labs(x = "", y = "A diversity [Shannon's H]")+
  ylim(y_min,y_max)


## Plot soil samples
Inoc_temp <- c()
for (i in 1:nrow(Sample_meta_soil)){
  if (Sample_meta_soil[i, "Sampling_day"] == "0"){
    Inoc_temp <- c(Inoc_temp, Lc(Sample_meta_soil[i, "Inoculant"], "_0"))
  }
  else{
    Inoc_temp <- c(Inoc_temp, as.character(Sample_meta_soil[i, "Inoculant"]))
  }
}
Sample_meta_soil$Inoculant2 <- Inoc_temp
Sample_meta_soil <- transform(Sample_meta_soil, Inoculant2 = factor(Inoculant2, levels = c("F5C_0", "F5S_0",
                                                                                           "F5C", "MIX", "F5S")))

p_alpha4 <- Sample_meta_soil %>%
  ggplot(aes(x = Inoculant2, y =  Shannon_index, fill = Condition)) +#
  geom_boxplot(outlier.colour = NA, position = position_dodge2(preserve = "single", width = 0.5)) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_point(aes(x=Inoculant2, y=Shannon_index, color=Condition),
             position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0), size=0.3) +
  theme(legend.position = "right", plot.title = element_text(size = 8))+
  ggtitle("") +
  labs(x = "", y = "A diversity [Shannon's H]")+
  ylim(y_min,y_max)

p_all <- gridExtra::grid.arrange(p_alpha2, p_alpha3, p_alpha1,p_alpha4, nrow = 2)
ggsave(file = "PATH/TO/Alpha_boxplot.pdf", plot = p_all, dpi = 300, width = 24, height = 24)
ggsave(file = "PATH/TO/Alpha_boxplot.png", plot = p_all, dpi = 300, width = 24, height = 24)


## No legends for papers
p_alpha1 <- Sample_meta_root %>%
  ggplot(aes(x = Inoculant, y =  Shannon_index, fill = Condition)) +#
  geom_boxplot(outlier.colour = NA, position = position_dodge2(preserve = "single", width = 0.5)) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_point(aes(x=Inoculant, y=Shannon_index, color=Condition),
             position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0), size=0.3) +
  theme(legend.position = "none", plot.title = element_text(size = 8))+
  ggtitle("") +
  labs(x = "", y = "A diversity [Shannon's H]")+
  ylim(y_min,y_max)

p_alpha2 <- Sample_meta_root %>%
  ggplot(aes(x = Host, y =  Shannon_index, fill = Condition)) +#
  geom_boxplot(outlier.colour = NA, position = position_dodge2(preserve = "single", width = 0.5)) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_point(aes(x=Host, y=Shannon_index, color=Condition),
             position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0), size=0.3) +
  theme(legend.position = "none", plot.title = element_text(size = 8))+
  ggtitle("") +
  labs(x = "", y = "A diversity [Shannon's H]")+
  ylim(y_min,y_max)

p_alpha3 <- Sample_meta_root %>%
  ggplot(aes(x = Host, y =  Shannon_index, fill = Inoculant)) +#
  geom_boxplot(outlier.colour = NA, position = position_dodge2(preserve = "single", width = 0.5)) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_point(aes(x=Host, y=Shannon_index, color=Inoculant),
             position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0), size=0.3) +
  theme(legend.position = "none", plot.title = element_text(size = 8))+
  ggtitle("") +
  labs(x = "", y = "A diversity [Shannon's H]")+
  ylim(y_min,y_max)

p_alpha4 <- Sample_meta_soil %>%
  ggplot(aes(x = Inoculant2, y =  Shannon_index, fill = Condition)) +#
  geom_boxplot(outlier.colour = NA, position = position_dodge2(preserve = "single", width = 0.5)) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_point(aes(x=Inoculant2, y=Shannon_index, color="Condition"),
             position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0), size=0.3) +
  theme(legend.position = "none", plot.title = element_text(size = 8))+
  ggtitle("") +
  labs(x = "", y = "A diversity [Shannon's H]")+
  ylim(y_min,y_max)


p_all <- gridExtra::grid.arrange(p_alpha2, p_alpha3, p_alpha1,p_alpha4, nrow = 2)
ggsave(file = "PATH/TO/Alpha_boxplot_nolegends.pdf", plot = p_all, dpi = 300, width = 24, height = 24)
ggsave(file = "PATH/TO/Alpha_boxplot_nolegends.png", plot = p_all, dpi = 300, width = 24, height = 24)

# ANOVA for alpha diversity
outputdir <- "PATH/TO/"
## For soil samples
## only 28 days
model1 <- lm(Shannon_index ~ Inoculant*Condition,
             data = Sample_meta_soil[Sample_meta_soil$Sampling_day != 0, ])
shapiro.test(resid(model1)) # p = 0.1192
bptest(model1) # P = 0.2817
model1_res <- data.frame(Anova(model2, type = 2, test = "F"))
write.csv(model1_res, Lc(outputdir, "alpha_anova_soil.csv"))

model1_tukey <- TukeyHSD(aov(model1))
WRITE_TUKEY_RESULT(model1_tukey, Lc(outputdir, "alpha_tukey_soil.xlsx"))

## only 0 day
model2 <- lm(Shannon_index ~ Inoculant,
             data = Sample_meta_soil[Sample_meta_soil$Sampling_day == 0, ])
shapiro.test(resid(model2)) # p = 0.4302
bptest(model2) # P = 0.2066
model2_res <- data.frame(Anova(model2, type = 2, test = "F")) # P = 0.0256

## For root samples
model3 <- lm(Shannon_index ~ Host*Inoculant*Condition,
             data = Sample_meta_root)
shapiro.test(resid(model3)) # p = 0.003648
bptest(model3) # P = 3.579e-09
model4 <- glm(Shannon_index ~ Host*Inoculant*Condition,
              data = Sample_meta_root,
              family = Gamma(log))
model4_res <- data.frame(Anova(model4, type = 2, test = "F"))
write.csv(model4_res, Lc(outputdir, "alpha_anova_root.csv"))
model3_tukey <- TukeyHSD(aov(model3))
WRITE_TUKEY_RESULT(model3_tukey, Lc(outputdir, "alpha_tukey_root.xlsx"))












##########
##########
##########
########## BETA DIVERSITY
##########
##########
##########




# MDS analysis



Sample_file <- "PATH/TO/Sample_metadata_1.csv"
Data_file <- "PATH/TO/Strain_rarefied_freq_1.csv"
outdir <- "PATH/TO/"
Summary_file <- "PATH/TO/Rarefy_summary_Strain_1.csv"

Sample_meta <- read.csv(Sample_file, header = T, row.names = 1)
Data <- read.csv(Data_file, header = T, row.names = 1)
ROWCHECK(Sample_meta, Data)
Sample_meta <- ORDER_CHANGE(Sample_meta)

#####
## All data
#####
Mds_data_1 <- metaMDS(Data, dist = "horn", autotransform = F, k = 2, try = 50, trymax = 100)
Sample_meta_1 <- data.frame(Sample_meta, Mds_data_1$points)


p_beta_all1 <- Sample_meta_1 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Inoculant, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "right", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta_all2 <- Sample_meta_1 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Host, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "right", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta_all3 <- Sample_meta_1 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Condition, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "right", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta_all4 <- Sample_meta_1 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Sampling_day, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "right", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta_all <- gridExtra::grid.arrange(p_beta_all1, p_beta_all2, p_beta_all3,p_beta_all4,  nrow = 2)
ggsave(file = Lc(outdir, "Beta_mds_all.pdf"), plot = p_beta_all, dpi = 300, width = 24, height = 24)
ggsave(file = Lc(outdir, "Beta_mds_all.png"), plot = p_beta_all, dpi = 300, width = 24, height = 24)

# No legends

p_beta_all1_no <- Sample_meta_1 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Inoculant, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "none", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta_all2_no <- Sample_meta_1 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Host, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "none", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta_all3_no <- Sample_meta_1 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Condition, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "none", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta_all4_no <- Sample_meta_1 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Sampling_day, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "none", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta_all_no <- gridExtra::grid.arrange(p_beta_all1_no, p_beta_all2_no, p_beta_all3_no,p_beta_all4_no,  nrow = 2)
ggsave(file = Lc(outdir, "Beta_mds_all_nolegends.pdf"), plot = p_beta_all_no, dpi = 300, width = 24, height = 24)
ggsave(file = Lc(outdir, "Beta_mds_all_nolegends.png"), plot = p_beta_all_no, dpi = 300, width = 24, height = 24)
# For supplemental figure S3

#####
# No non-inoculant & 28 days data
#####
Df_temp <- DF_ROWFILT(Data, Sample_meta, (Sample_meta$Inoculant != "CONT" & Sample_meta$Sampling_day == 28))
Data2 <- Df_temp[1][[1]]
Sample_meta2 <- Df_temp[2][[1]]
Data2 <- DF_NONZERO(Data2)

Mds_data_2 <- metaMDS(Data2, dist = "horn", autotransform = F, k = 2, try = 50, trymax = 100)
Sample_meta2_2 <- data.frame(Sample_meta2, Mds_data_2$points)


p_beta1 <- Sample_meta2_2 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Inoculant, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "right", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta2 <- Sample_meta2_2 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Condition, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "right", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta3 <- Sample_meta2_2 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Host, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "right", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))



# Only soils
Df_temp <- DF_ROWFILT(Data, Sample_meta, Sample_meta$Type == "Soil")
Data3 <- Df_temp[1][[1]]
Sample_meta3 <- Df_temp[2][[1]]
Data3 <- DF_NONZERO(Data3)

Mds_data_3 <- metaMDS(Data3, dist = "horn", autotransform = F, k = 2, try = 50, trymax = 100)
Sample_meta2_3 <- data.frame(Sample_meta3, Mds_data_3$points)
Sample_meta2_3$Condition <- factor(Sample_meta2_3$Condition, levels = c(levels(Sample_meta2_3$Condition), "Extract"))
Sample_meta2_3[Sample_meta2_3$Sampling_day == 0, ]$Condition <- "Extract"

p_beta4 <- Sample_meta2_3 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Inoculant, shape = Condition)) +
  geom_point(size = 3, alpha = 0.5) +
  theme(legend.position = "right", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(16, 15, 17))


p_beta_main <- gridExtra::grid.arrange(p_beta3, p_beta1, p_beta2, p_beta4, nrow = 2)
ggsave(file = Lc(outdir, "Beta_mds_main.pdf"), plot = p_beta_main, dpi = 300, width = 24, height = 24)
ggsave(file = Lc(outdir, "Beta_mds_main.png"), plot = p_beta_main, dpi = 300, width = 24, height = 24)

# no legends

p_beta1_no <- Sample_meta2_2 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Inoculant, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "none", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta2_no <- Sample_meta2_2 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Condition, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "none", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta3_no <- Sample_meta2_2 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Host, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "none", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta4_no <- Sample_meta2_3 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Inoculant, shape = Condition)) +
  geom_point(size = 3, alpha = 0.5) +
  theme(legend.position = "none", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(16, 15, 17))

p_beta_main_no <- gridExtra::grid.arrange(p_beta3_no, p_beta1_no, p_beta2_no, p_beta4_no, nrow = 2)
ggsave(file = Lc(outdir, "Beta_mds_main_nolegends.pdf"), plot = p_beta_main_no, dpi = 300, width = 24, height = 24)
ggsave(file = Lc(outdir, "Beta_mds_main_nolegends.png"), plot = p_beta_main_no, dpi = 300, width = 24, height = 24)
# main figure 2

#####
#####
# Genus
#####
Sample_file <- "PATH/TO/Sample_metadata_3.csv"
Data_file <- "PATH/TO/Genus_rarefied_freq_3.csv"

Sample_meta <- read.csv(Sample_file, header = T, row.names = 1)
Data <- read.csv(Data_file, header = T, row.names = 1)
ROWCHECK(Sample_meta, Data)
Sample_meta <- ORDER_CHANGE(Sample_meta)

Data <- Data[, colnames(Data) != "NotAssigned"]
Df_temp <- DF_ROWFILT(Data, Sample_meta, (Sample_meta$Inoculant != "CONT" & Sample_meta$Sampling_day == 28))
Data2 <- Df_temp[1][[1]]
Sample_meta2 <- Df_temp[2][[1]]
Data2 <- DF_NONZERO(Data2)

Mds_data_2 <- metaMDS(Data2, dist = "horn", autotransform = F, k = 2, try = 50, trymax = 100)
Sample_meta2_2 <- data.frame(Sample_meta2, Mds_data_2$points)

p_beta_genus1 <- Sample_meta2_2 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Inoculant, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "right", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta_genus2 <- Sample_meta2_2 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Condition, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "right", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta_genus3 <- Sample_meta2_2 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Host, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "right", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta_genus1_no <- Sample_meta2_2 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Inoculant, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "none", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta_genus2_no <- Sample_meta2_2 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Condition, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "none", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta_genus3_no <- Sample_meta2_2 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Host, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "none", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

# Family
#####
Sample_file <- "PATH/TO/Sample_metadata_4.csv"
Data_file <- "PATH/TO/Family_rarefied_freq_4.csv"

Sample_meta <- read.csv(Sample_file, header = T, row.names = 1)
Data <- read.csv(Data_file, header = T, row.names = 1)
ROWCHECK(Sample_meta, Data)
Sample_meta <- ORDER_CHANGE(Sample_meta)

Data <- Data[, colnames(Data) != "NotAssigned"]
Df_temp <- DF_ROWFILT(Data, Sample_meta, (Sample_meta$Inoculant != "CONT" & Sample_meta$Sampling_day == 28))
Data2 <- Df_temp[1][[1]]
Sample_meta2 <- Df_temp[2][[1]]
Data2 <- DF_NONZERO(Data2)

Mds_data_2 <- metaMDS(Data2, dist = "horn", autotransform = F, k = 2, try = 50, trymax = 100)
Sample_meta2_2 <- data.frame(Sample_meta2, Mds_data_2$points)

p_beta_family1 <- Sample_meta2_2 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Inoculant, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "right", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta_family2 <- Sample_meta2_2 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Condition, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "right", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta_family3 <- Sample_meta2_2 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Host, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "right", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta_family1_no <- Sample_meta2_2 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Inoculant, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "none", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta_family2_no <- Sample_meta2_2 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Condition, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "none", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta_family3_no <- Sample_meta2_2 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Host, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "none", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

# ASV
#####
Sample_file <- "PATH/TO/Sample_metadata_5.csv"
Data_file <- "PATH/TO/ASV_rarefied_freq_5.csv"

Sample_meta <- read.csv(Sample_file, header = T, row.names = 1)
Data <- read.csv(Data_file, header = T, row.names = 1)
ROWCHECK(Sample_meta, Data)
Sample_meta <- ORDER_CHANGE(Sample_meta)

Data <- Data[, colnames(Data) != "NotAssigned"]
Df_temp <- DF_ROWFILT(Data, Sample_meta, (Sample_meta$Inoculant != "CONT" & Sample_meta$Sampling_day == 28))
Data2 <- Df_temp[1][[1]]
Sample_meta2 <- Df_temp[2][[1]]
Data2 <- DF_NONZERO(Data2)

Mds_data_2 <- metaMDS(Data2, dist = "horn", autotransform = F, k = 2, try = 50, trymax = 100)
Sample_meta2_2 <- data.frame(Sample_meta2, Mds_data_2$points)

p_beta_ASV1 <- Sample_meta2_2 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Inoculant, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "right", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta_ASV2 <- Sample_meta2_2 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Condition, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "right", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta_ASV3 <- Sample_meta2_2 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Host, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "right", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))


p_beta_ASV1_no <- Sample_meta2_2 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Inoculant, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "none", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta_ASV2_no <- Sample_meta2_2 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Condition, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "none", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta_ASV3_no <- Sample_meta2_2 %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Host, shape = Type)) +
  geom_point(size = 4, alpha = 0.5) +
  theme(legend.position = "none", plot.title = element_text(size = 8),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))+
  scale_shape_manual(values = c(20, 18))

p_beta_taxa <- gridExtra::grid.arrange(p_beta_ASV3, p_beta_ASV1, p_beta_ASV2,
                                       p_beta_genus3,p_beta_genus1, p_beta_genus2,
                                       p_beta_family3, p_beta_family1, p_beta_family2,
                                       nrow = 3)
ggsave(file = Lc(outdir, "Beta_mds_othertaxa.pdf"), plot = p_beta_taxa, dpi = 300, width = 24, height = 24)
ggsave(file = Lc(outdir, "Beta_mds_othertaxa.png"), plot = p_beta_taxa, dpi = 300, width = 24, height = 24)

p_beta_taxa_no <- gridExtra::grid.arrange(p_beta_ASV3_no, p_beta_ASV1_no, p_beta_ASV2_no,
                                          p_beta_genus3_no,p_beta_genus1_no, p_beta_genus2_no,
                                          p_beta_family3_no, p_beta_family1_no, p_beta_family2_no,
                                          nrow = 3)
ggsave(file = Lc(outdir, "Beta_mds_othertaxa_nolegends.pdf"), plot = p_beta_taxa_no, dpi = 300, width = 24, height = 24)
ggsave(file = Lc(outdir, "Beta_mds_othertaxa_nolegends.png"), plot = p_beta_taxa_no, dpi = 300, width = 24, height = 24)




# PERMANOVA analysis



#####
# OTU
#####
Sample_file <- "PATH/TO/Sample_metadata_1.csv"
Data_file <- "PATH/TO/Strain_rarefied_freq_1.csv"
outdir <- "PATH/TO/"
Summary_file <- "PATH/TO/Rarefy_summary_Strain_1.csv"

Sample_meta <- read.csv(Sample_file, header = T, row.names = 1)
Data <- read.csv(Data_file, header = T, row.names = 1)
ROWCHECK(Sample_meta, Data)
Sample_meta <- ORDER_CHANGE(Sample_meta)

Df_temp <- DF_ROWFILT(Data, Sample_meta, (Sample_meta$Inoculant != "CONT" & Sample_meta$Type == "Root"))
Data2 <- Df_temp[1][[1]]
Sample_meta2 <- Df_temp[2][[1]]
Data2 <- DF_NONZERO(Data2)

Perm_otu <- adonis2(Data2 ~ Host*Inoculant*Condition, data = Sample_meta2, permutations = 9999, method = "horn")
write.csv(Perm_otu, Lc(outdir, "Beta_PERMANOVA_main.csv"))


#####
# OTU only Lotus japonicus
#####

Sample_file <- "PATH/TO/Sample_metadata_1.csv"
Data_file <- "PATH/TO/Strain_rarefied_freq_1.csv"
outdir <- "PATH/TO/"

Sample_meta <- read.csv(Sample_file, header = T, row.names = 1)
Data <- read.csv(Data_file, header = T, row.names = 1)
ROWCHECK(Sample_meta, Data)
Sample_meta <- ORDER_CHANGE(Sample_meta)

Df_temp <- DF_ROWFILT(Data, Sample_meta, (Sample_meta$Inoculant != "CONT" & Sample_meta$Type == "Root" & Sample_meta$Host != "Burttii"))
Data2 <- Df_temp[1][[1]]
Sample_meta2 <- Df_temp[2][[1]]
Data2 <- DF_NONZERO(Data2)

Perm_otu <- adonis2(Data2 ~ Host*Inoculant*Condition, data = Sample_meta2, permutations = 9999, method = "horn")
write.csv(Perm_otu, Lc(outdir, "Beta_PERMANOVA_onlyLja.csv"))

#####
# Other taxonomic levels
#####

Sample_file <- "PATH/TO/Sample_metadata_3.csv"
Data_file <- "PATH/TO/Genus_rarefied_freq_3.csv"
outdir <- "PATH/TO/"

Sample_meta <- read.csv(Sample_file, header = T, row.names = 1)
Data <- read.csv(Data_file, header = T, row.names = 1)
ROWCHECK(Sample_meta, Data)
Sample_meta <- ORDER_CHANGE(Sample_meta)
Data <- Data[, colnames(Data) != "NotAssigned"]


Df_temp <- DF_ROWFILT(Data, Sample_meta, (Sample_meta$Inoculant != "CONT" & Sample_meta$Type == "Root"))
Data2 <- Df_temp[1][[1]]
Sample_meta2 <- Df_temp[2][[1]]
Data2 <- DF_NONZERO(Data2)

Perm_genus <- adonis2(Data2 ~ Host*Inoculant*Condition, data = Sample_meta2, permutations = 9999, method = "horn")

#
Sample_file <- "PATH/TO/Sample_metadata_4.csv"
Data_file <- "PATH/TO/Family_rarefied_freq_4.csv"
outdir <- "PATH/TO/"

Sample_meta <- read.csv(Sample_file, header = T, row.names = 1)
Data <- read.csv(Data_file, header = T, row.names = 1)
ROWCHECK(Sample_meta, Data)
Sample_meta <- ORDER_CHANGE(Sample_meta)
Data <- Data[, colnames(Data) != "NotAssigned"]


Df_temp <- DF_ROWFILT(Data, Sample_meta, (Sample_meta$Inoculant != "CONT" & Sample_meta$Type == "Root"))
Data2 <- Df_temp[1][[1]]
Sample_meta2 <- Df_temp[2][[1]]
Data2 <- DF_NONZERO(Data2)

Perm_family <- adonis2(Data2 ~ Host*Inoculant*Condition, data = Sample_meta2, permutations = 9999, method = "horn")

Sample_file <- "PATH/TO/Sample_metadata_5.csv"
Data_file <- "PATH/TO/ASV_rarefied_freq_5.csv"
outdir <- "PATH/TO/"

Sample_meta <- read.csv(Sample_file, header = T, row.names = 1)
Data <- read.csv(Data_file, header = T, row.names = 1)
ROWCHECK(Sample_meta, Data)
Sample_meta <- ORDER_CHANGE(Sample_meta)
Data <- Data[, colnames(Data) != "NotAssigned"]


Df_temp <- DF_ROWFILT(Data, Sample_meta, (Sample_meta$Inoculant != "CONT" & Sample_meta$Type == "Root"))
Data2 <- Df_temp[1][[1]]
Sample_meta2 <- Df_temp[2][[1]]
Data2 <- DF_NONZERO(Data2)

Perm_ASV <- adonis2(Data2 ~ Host*Inoculant*Condition, data = Sample_meta2, permutations = 9999, method = "horn")

wb <- createWorkbook()
sheetname = "Genus"
addWorksheet(wb, sheetname)
writeData(wb, sheetname, Perm_genus, rowNames = TRUE)
sheetname = "Family"
addWorksheet(wb, sheetname)
writeData(wb, sheetname, Perm_family, rowNames = TRUE)
sheetname = "ASV"
addWorksheet(wb, sheetname)
writeData(wb, sheetname, Perm_ASV, rowNames = TRUE)
saveWorkbook(wb, Lc(outdir, "Beta_PERMANOVA_othertaxa.xlsx"), overwrite = TRUE)


#####
# Estimate only host effects 
#####
Sample_file <- "PATH/TO/Sample_metadata_1.csv"
Data_file <- "PATH/TO/Strain_rarefied_freq_1.csv"
outdir <- "PATH/TO/"
Summary_file <- "PATH/TO/Rarefy_summary_Strain_1.csv"

Sample_meta <- read.csv(Sample_file, header = T, row.names = 1)
Data <- read.csv(Data_file, header = T, row.names = 1)
ROWCHECK(Sample_meta, Data)
Sample_meta <- ORDER_CHANGE(Sample_meta)

Df_temp <- DF_ROWFILT(Data, Sample_meta, (Sample_meta$Inoculant != "CONT" & Sample_meta$Type == "Root"))
Data2 <- Df_temp[1][[1]]
Sample_meta2 <- Df_temp[2][[1]]
Data2 <- DF_NONZERO(Data2)

PERMANOVA_HOST(Data2, Sample_meta2, outdir)




# Comparison between distance of root microbiomes and plant genetic distance




Sample_file <- "PATH/TO/Sample_metadata_1.csv"
Data_file <- "PATH/TO/Strain_rarefied_freq_1.csv"
outdir <- "PATH/TO/"

Sample_meta <- read.csv(Sample_file, header = T, row.names = 1)
Data <- read.csv(Data_file, header = T, row.names = 1)
ROWCHECK(Sample_meta, Data)
Sample_meta <- ORDER_CHANGE(Sample_meta)

DF_temp <- DF_ROWFILT(Data, Sample_meta, (Sample_meta$Type == "Root" & Sample_meta$Inoculant != "CONT"))
Data2 <- DF_temp[1][[1]]
Sample_meta2 <- DF_temp[2][[1]]

Geno_data <- "PATH/TO/MG_line_SNPs_MG20mapping.csv"
Geno_df <- read.csv(Geno_data)

Geno_kinship <- GENOME_KINSHIP(Geno_df, Sample_meta2)
MDS_list <- AVE_BETA_HOST(Data2, Sample_meta2)
mantel_result <- MANTEL_MDS_KINSHIP(MDS_list, Geno_kinship, outdir)

kin_vector <- Geno_kinship[upper.tri(Geno_kinship, diag = FALSE)]
mds_vectors <- lapply(MDS_list[[1]], function(mds) 1 - mds[upper.tri(mds, diag = FALSE)])
Plot_df <- data.frame(
  InoculantCondition = rep(MDS_list[[2]], each = length(kin_vector)),
  Kinship = rep(kin_vector, times = length(MDS_list[[1]])),
  MDS_distance = unlist(mds_vectors))

p_distcomp <- Plot_df %>%
  ggplot(aes(x = Kinship, y = MDS_distance, color = InoculantCondition))+
  geom_point(size = 4, alpha = 0.25) +
  geom_smooth(method=lm, se=FALSE, alpha = 0.5) +
  theme(legend.position = "right", plot.title = element_text(size = 8)) +
  ggtitle("Compare of distances between genotype and microbiome") +
  xlab("IBS kinship") +
  ylab("1 - Distance with Horn index")

ggsave(file = "PATH/TO/Beta_mdswithkinship.pdf", plot = p_distcomp, dpi = 300, width = 8, height = 8)
ggsave(file = "PATH/TO/Beta_mdswithkinship.png", plot = p_distcomp, dpi = 300, width = 8, height = 8)

p_distcomp_no <- Plot_df %>%
  ggplot(aes(x = Kinship, y = MDS_distance, color = InoculantCondition))+
  geom_point(size = 4, alpha = 0.25) +
  geom_smooth(method=lm, se=FALSE, alpha = 0.5) +
  theme(legend.position = "none", plot.title = element_text(size = 8)) +
  ggtitle("Compare of distances between genotype and microbiome") +
  xlab("IBS kinship") +
  ylab("1 - Distance with Horn index")

ggsave(file = "PATH/TO/Beta_mdswithkinship_no.pdf", plot = p_distcomp_no, dpi = 300, width = 8, height = 8)
ggsave(file = "PATH/TO/Beta_mdswithkinship_no.png", plot = p_distcomp_no, dpi = 300, width = 8, height = 8)





##########
##########
##########
########## GSM ANALYSIS
##########
##########
##########

#####
# Taxa file making (DBid2Strain)
#####
Taxa_file <- "PATH/TO//Data/INROOT02_Rev4_taxa.csv"
taxa_data <- read.csv(Taxa_file, header = T, row.names = 1)
taxa_data <- taxa_data[order(rownames(taxa_data)), ]
Id_list <- c()
for (i in unique(taxa_data$Strain_information)){
  Id_list <- c(Id_list, rownames(taxa_data[taxa_data$Strain_information == i, ])[1])
}
taxa_data_edited <- taxa_data[Id_list, ]
rownames(taxa_data_edited) <- sapply(strsplit(taxa_data_edited$Strain_information, split = " "), "[", 1)

write.csv(taxa_data_edited, "PATH/TO/taxa_strain.csv")

#####
# GSM
#####
Sample_file <- "PATH/TO/Sample_metadata_1.csv"
Data_file <- "PATH/TO/Strain_freq_1.csv"
Taxa_file <- "PATH/TO/taxa_strain.csv"
outdir <- "PATH/TO/"

Sample_meta <- read.csv(Sample_file, header = T, row.names = 1)
Data <- read.csv(Data_file, header = T, row.names = 1)
ROWCHECK(Sample_meta, Data)
Sample_meta <- ORDER_CHANGE(Sample_meta)
Taxa_data <- read.csv(Taxa_file, header = T, row.names = 1)

DF_temp <- DF_ROWFILT(Data, Sample_meta, (Sample_meta$Type == "Root" & Sample_meta$Inoculant != "CONT"))
Data2 <- DF_temp[1][[1]]
Sample_meta2 <- DF_temp[2][[1]]
Data2 <- DF_NONZERO(Data2)

GSM_ANALYSIS(Data2, Sample_meta2, Taxa_data, outdir)

#####
# GSM Fisher
#####
GSM_result <- read.csv(Lc(outdir, "GSM_result.csv"), header = T, row.names = 1)
GSM_result <- GSM_result[!is.na(GSM_result$domain), ]

GSM_FISHER(GSM_result, Taxa_data, ourdir)

#####
# Relative ratio in enriched taxa
#####

Fish_result <- read.csv(Lc(outdir, "GSM_Fisher_family.csv"), header = T, row.names = 1)
Sample_meta <- read.csv("PATH/TO/Sample_metadata_4.csv", header = T, row.names = 1)
Data <- read.csv("PATH/TO/Family_freq_4.csv", header = T, row.names = 1)
ROWCHECK(Sample_meta, Data)
Sample_meta <- ORDER_CHANGE(Sample_meta)

DF_temp <- DF_ROWFILT(Data, Sample_meta, (Sample_meta$Type == "Root" & Sample_meta$Inoculant != "CONT"))
Data2 <- DF_temp[1][[1]]
Sample_meta2 <- DF_temp[2][[1]]

Fish_result <- Fish_result[, grep("BH", names(Fish_result), value = T)]
Fish_result_host <-  Fish_result[, grep("Host", names(Fish_result), value = T)]
Fish_result_host2 <- Fish_result_host < 0.05
Sig_family <- rownames(Fish_result_host2[rowSums(Fish_result_host2) != 0, ])
Data2_mean <- colMeans(Data2)
for (i in Sig_family){
  Sig_temp <- mean(Data2[, make.names(i)])
  ranks <- sum(Data2_mean > Sig_temp) + 1
  print(Lc(i, ":", Sig_temp, ", Rank : ", ranks))
}


#####
# Make venn diagram
#####
GSM_result <- read.csv(Lc(outdir, "GSM_result.csv"), header = T, row.names = 1)
GSM_result <- GSM_result[!is.na(GSM_result$domain), ]

GSM_VENN(GSM_result, outdir)



##########
##########
##########
########## 
##########
##########
##########


