# Import start
source("PATH/TO/Function_source.R")
# Import end

#####
# Preprocess
#####
Sample_file <- "PATH/TO/Sample_metadata.csv"
outdir <- "PATH/TO/Output/"

Sample_meta <- read.csv(Sample_file, header = T, row.names = 1)
Sample_meta <- Sample_meta[Sample_meta$Type != "Soil", ]
Sample_meta <- ORDER_CHANGE(Sample_meta)

#####
# Heatmap plot
#####

HEATMAP_PLOT(Sample_meta, "Length", outdir, "Heatmap_SL")
HEATMAP_PLOT(Sample_meta, "RootLength", outdir, "Heatmap_RL")
HEATMAP_PLOT(Sample_meta, "Number_of_leaf", outdir, "Heatmap_NOL")
HEATMAP_PLOT(Sample_meta, "Number_of_branch", outdir, "Heatmap_NOB")


#####
# Generalized linear model
#####
outdir <- "PATH/TO/Output/"
dir.create(outdir)

glm_model <- Length ~ Host*Inoculant*Condition
GLM_ANALYSIS(glm_model, Sample_meta, outdir, "SL_glm")
TUKEY_ANALYSIS_FULLMODEL(glm_model, Sample_meta, outdir, "SL_tuk")

glm_model <- RootLength ~ Host*Inoculant*Condition
GLM_ANALYSIS(glm_model, Sample_meta, outdir, "RL_glm")
TUKEY_ANALYSIS_FULLMODEL(glm_model, Sample_meta, outdir, "RL_tuk")

glm_model <- Number_of_branch ~ Host*Inoculant*Condition
GLM_ANALYSIS(glm_model, Sample_meta, outdir, "NOB_glm")
TUKEY_ANALYSIS_FULLMODEL(glm_model, Sample_meta, outdir, "NOB_tuk")

glm_model <- Number_of_leaf ~ Host*Inoculant*Condition
GLM_ANALYSIS(glm_model, Sample_meta, outdir, "NOL_glm")
TUKEY_ANALYSIS_FULLMODEL(glm_model, Sample_meta, outdir, "NOB_tuk")



# GLM without non-inoculant data
Sample_meta <- Sample_meta[Sample_meta$Inoculant != "CONT", ]

glm_model <- Length ~ Host*Inoculant*Condition
GLM_ANALYSIS(glm_model, Sample_meta, outdir, "SL_glm_noCont")
TUKEY_ANALYSIS_FULLMODEL(glm_model, Sample_meta, outdir, "SL_tuk_noCont")

glm_model <- RootLength ~ Host*Inoculant*Condition
GLM_ANALYSIS(glm_model, Sample_meta, outdir, "RL_glm_noCont")
TUKEY_ANALYSIS_FULLMODEL(glm_model, Sample_meta, outdir, "RL_tuk_noCont")

glm_model <- Number_of_branch ~ Host*Inoculant*Condition
GLM_ANALYSIS(glm_model, Sample_meta, outdir, "NOB_glm_noCont")
TUKEY_ANALYSIS_FULLMODEL(glm_model, Sample_meta, outdir, "NOB_tuk_noCont")

glm_model <- Number_of_leaf ~ Host*Inoculant*Condition
GLM_ANALYSIS(glm_model, Sample_meta, outdir, "NOL_glm_noCont")
TUKEY_ANALYSIS_FULLMODEL(glm_model, Sample_meta, outdir, "NOB_tuk_noCont")



#####
# GLM randomize pot
#####
Random_list <- DF_CHOSEPOTS(Sample_meta, rep = 10000)

glm_model <- Length ~ Host*Inoculant*Condition
glm_random <- GLM_RANDOM(glm_model, Sample_meta, Random_list)
glm_random_eta <- RONDOM_ETA_PLOT(glm_random, outdir, "SL_glm_violine")

glm_model <- RootLength ~ Host*Inoculant*Condition
glm_random <- GLM_RANDOM(glm_model, Sample_meta, Random_list)
glm_random_eta <- RONDOM_ETA_PLOT(glm_random, outdir, "RL_glm_violine")




#####
### Phenotype correlations
#####
PHENOTYPE_PPLOT(Sample_meta, "Host", outdir, "PairPlot_Host")
PHENOTYPE_PPLOT(Sample_meta, "Inoculant", outdir, "PairPlot_Inoculant")
PHENOTYPE_PPLOT(Sample_meta, "Condition", outdir, "PairPlot_Condition")



#####
### Phenotype ICC calculation
#####
ICC_CALC(Sample_meta, "Length", outdir, "ICC_SL")
ICC_CALC(Sample_meta, "RootLength", outdir, "ICC_RL")



#####
### Phenotype variance component
#####
Sample_file <- "PATH/TO/Sample_metadata_1.csv"
Data_file <- "PATH/TO/Strain_rarefied_freq_1.csv"
outdir <- "PATH/TO/Output/"

Sample_meta <- read.csv(Sample_file, header = T, row.names = 1)
Data <- read.csv(Data_file, header = T, row.names = 1)
ROWCHECK(Sample_meta, Data)
Sample_meta <- ORDER_CHANGE(Sample_meta)

Df_temp <- DF_ROWFILT(Data, Sample_meta, (Sample_meta$Inoculant != "CONT" & Sample_meta$Type == "Root"))
Data2 <- Df_temp[1][[1]]
Sample_meta2 <- Df_temp[2][[1]]
Data2 <- DF_NONZERO(Data2)

HERITABILITY(Data2, Sample_meta2)


#####
#####
#####

