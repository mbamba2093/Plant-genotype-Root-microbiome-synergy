library(tidyverse)
library(ggplot2)
library(viridis)
library(vegan)
library(reshape2)
library(car)
library(GGally)
library(lmtest)
library(openxlsx)
library(statgenGWAS)
library(ape)
library(npreg)
library(VennDiagram)
library(gaston)

###

### Common function
Lc <- function(..., f) {
  # Functions to easily connect strings
  if (missing(f)) 
    f <- paste(rep("%s", length(c(...))), collapse = "")
  sprintf(fmt = f, ...)
}

ROWCHECK = function(df1, df2){
  # Check if the line names of the two df's are the same.
  
  if (nrow(df1) == nrow(df2)){
    for (i in 1:nrow(df2)){
      if (rownames(df1)[i] != rownames(df2)[i]){
        stop("Warnings: Rownames were not matched between df1 and df2")
      }
    }
  }
  else{
    stop("Warnings: Number of df1 and df2 were different")
  }
  print("Ready for analysis")
  return(T)
}

DF_ROWFILT = function(df1, df2, filt){
  # Filter two DFs with the same boolean.
  # Do not break the structure of Microbiome data and Sample metadata.
  if (ROWCHECK(df1, df2)){
    df1_ <- df1[filt, ]
    df2_ <- df2[filt, ]
  }
  else{
    print("df1 and df2 did not match their rownames")
  }
  return(list(df1_, df2_))
}

DF_FREQ = function(df){
  # Converts the value of each row to a frequency and returns the replaced DF
  df_temp <- data.frame(df/rowSums(df))
  return(df_temp)
}

DF_NONZERO = function(df){
  # Remove columns not identified in all samples
  df <- df[, colSums(df) != 0]
  return(df)
}

DF_CATE2NUM = function(type_list){
  # Function to allow easy conversion of category data to color
  type_list_ref <- unique(type_list)
  out_list <- rep(-9, length(type_list))
  for (i in (1:length(type_list_ref))){
    for (j in (1:length(type_list))){
      if (type_list[j] == type_list_ref[i]){
        out_list[j] = i
      }
    }
  }
  return(out_list)
}

ORDER_CHANGE = function(Sample_data){
  # This function is used for re-ordering to plot easily.
  Sample_data <- transform(Sample_data, Inoculant= factor(Inoculant, levels = c("F5C", "MIX", "F5S", "CONT")))
  Sample_data <- transform(Sample_data, Host= factor(Host, levels = c("Gifu", "MG20", "Burttii",
                                                                      "MG11", "MG46", "MG56",
                                                                      "MG63", "MG67", "MG68", "Soil")))
  Sample_data <- transform(Sample_data, Condition= factor(Condition, levels = c("NonSalt", "Salt")))
  Sample_data <- transform(Sample_data, Type= factor(Type, levels = c("Root", "Soil")))
  Sample_data$Sampling_day <- as.character(Sample_data$Sampling_day)
  Sample_data <- transform(Sample_data, Sampling_day = factor(Sampling_day, levels = c("0", "28")))
  return(Sample_data)
}

###
#####
### Microbiome
#####
CALCMAXSLOPE = function(df, out){
  # Calculates the minimum slope for each sample and returns the value of the sample with the largest minimum slope
  df <- data.matrix(df)
  max_list <- c()
  for (i in 1:nrow(df)){
    max_list[i] <- rareslope(df[i,],sum(df[i,])-1)
  }
  print(paste("Max slope is", as.character(1 - max(max_list)), sep = " "))
  sink(out, append = T)
  print(paste("Max slope is", as.character(1 - max(max_list)), sep = " "))
  sink()
  return(max_list)
}

RAREFACTION = function(df, max_list){
  # rarefaction to match the maximum slope sample
  # Search for a value that does not exceed MAX SLOPE for each digit of the number.
  slope_list <- c()
  for (i in 1:nrow(df)){
    rareslope_start <- c()
    for (j in round(log10(sum(df[i,]))):0){
      if (j == round(log10(sum(df[i,])))){
        rareslope_tmp <- rareslope(df[i,],seq(1, sum(df[i,]), 10^j))
        rareslope_start[j+1] <- length(rareslope_tmp[rareslope_tmp >= max(max_list) ])-1
        rareslope_start[j+1] <- rareslope_start[j+1]*10^j
        rareslope_start[is.na(rareslope_start)] <- 0
      }
      else{
        rareslope_tmp <- rareslope(df[i,],seq(sum(rareslope_start) , sum(rareslope_start) + 10^(j+1), 10^j))
        rareslope_start[j+1] <- length(rareslope_tmp[rareslope_tmp >= max(max_list) ])-1
        rareslope_start[j+1] <- rareslope_start[j+1]*10^j
      }
    }
    slope_list[i] <- sum(rareslope_start)
  }
  df_rarefied <-rrarefy(df,slope_list)
  return(df_rarefied)
}

DF_SUMMARY = function(df){
  # Calculate summaries of Microbiome data.
  COUNT_NON_ZERO = function(Values){return(sum(Values > 0))}
  Sum_values <- rowSums(df)
  Mean_values <- rowMeans(df)
  Nonzero_values <- c()
  for (i in 1:nrow(df)){
    Nonzero_values <- c(Nonzero_values, COUNT_NON_ZERO(df[i, ]))
  }
  Out_df <- data.frame(Sum_values = Sum_values,
                       Mean_values = Mean_values,
                       Nonzero_values = Nonzero_values)
  rownames(Out_df) <- rownames(df)
  return(Out_df)
}

RARECURVE_PLOT = function(data, df_summary, out_name = "rarecurve", col_list = 4, showonly = T, step = 1000){
  # Plot rarefaction curve. df_summary is used for determin x and y lim.
  xlim_r <- max(df_summary$Sum_values ) + max(df_summary$Sum_values )*0.05
  ylim_r = max(df_summary$Nonzero_values ) + max(df_summary$Nonzero_values)*0.05
  rarecurve(data, step = step,
            col = col_list + 3,
            main = "Rarefaction curve",
            xlim = c(0,xlim_r), ylim = c(0, ylim_r),
            cex = 0.5, label = T,
            xlab = "Sample size", ylab = "Number of strains")
  if (showonly != T){
    pdf(paste(out_name, ".pdf", sep = ""), width = 14, height = 10)
    rarecurve(data, step = step,
              col = col_list + 3,
              main = "Rarefaction curve",
              xlim = c(0,xlim_r), ylim = c(0, ylim_r),
              cex = 0.5, label = TRUE,
              xlab = "Sample size", ylab = "Number of strains")
    dev.off()
    
    png(paste(out_name, ".png", sep = ""), width = 2400, height = 1697, res = 300)
    rarecurve(data, step = step,
              col = col_list + 3,
              main = "Rarefaction curve",
              xlim = c(0,xlim_r), ylim = c(0, ylim_r),
              cex = 0.5, label = TRUE,
              xlab = "Sample size", ylab = "Number of strains")
    dev.off()
  }
  
}

CALC_ALPHA = function(data, metadata){
  # Calculation alpha diversity
  Shannon_index <- diversity(data, index = "shannon", base = 2)
  Simpson_index <- diversity(data, index = "simpson")
  SimpsonInv_index <- diversity(data, index = "invsimpson")
  metadata <- data.frame(metadata, Shannon_index, Simpson_index, SimpsonInv_index)
  return(metadata)
}

WRITE_TUKEY_RESULT = function(tuk_result, output){
  wb <- createWorkbook()
  tuk_names <- names(tuk_result)
  for (temp_name in tuk_names){
    sheetname = make.names(temp_name)
    addWorksheet(wb, sheetname)
    writeData(wb, sheetname, as.data.frame(tuk_result[[temp_name]]), rowNames = TRUE)
  }
  saveWorkbook(wb, output, overwrite = TRUE)
}

PERMANOVA_HOST = function(data, sample, outdir){
  interaction_list <- interaction(sample$Inoculant, sample$Condition)
  num_temp <- 1
  out_df <- data.frame(matrix(nrow = 0, ncol = 5))
  for (i in unique(interaction_list)){
    Df_temp <- DF_ROWFILT(data, sample, interaction_list == i)
    data_temp <- Df_temp[1][[1]]
    sample_temp <- Df_temp[2][[1]]
    data_temp <- DF_NONZERO(data_temp)
    perm_temp <- adonis2(data_temp ~ Host, data = sample_temp, permutations = 9999, method = "horn")
    perm_temp <- as.data.frame(perm_temp["Host", ])
    rownames(perm_temp) <- i
    out_df <- rbind(out_df, perm_temp)
  }
  write.csv(out_df, Lc(outdir, "Beta_PERMANOVA_host.csv"))
}

GENOME_KINSHIP = function(Geno_df, Sample_meta){
  colnames(Geno_df) <- gsub("^MG00|^MG0", "MG", colnames(Geno_df))
  rownames(Geno_df) <- apply(Geno_df, 1, function(i){Lc(i[1], "_", i[2])})
  Geno_df <- Geno_df[, c(-1, -2)]
  
  host_list <- unique(Sample_meta$Host)
  host_list <- host_list[!(host_list == "Burttii" | host_list == "Soil")]
  
  Geno_df <- Geno_df[, as.vector(host_list)]
  Geno_df <- t(Geno_df)
  Geno_df <- DF_NONZERO(Geno_df)
  Geno_df <- Geno_df[order(rownames(Geno_df)), ]
  Geno_df <- Geno_df/2
  Geno_kinship <- kinship(Geno_df, method = "IBS")
  diag(Geno_kinship) <- 1
  return (Geno_kinship)
}

AVE_BETA_HOST = function(data, sample){
  dist_list <- list()
  Interaction_list <- interaction(sample$Inoculant, sample$Condition)
  host_list <- unique(sample$Host)
  host_list <- host_list[!(host_list == "Burttii" | host_list == "Soil")]
  host_list <- host_list[order(host_list)]
  for (i in unique(Interaction_list)){
    Df_temp <- DF_ROWFILT(data, sample, Interaction_list == i)
    data_temp <- Df_temp[1][[1]]
    sample_temp <- Df_temp[2][[1]]
    data_temp <- DF_NONZERO(data_temp)
    data_temp <- data.frame(Host = sample_temp$Host,
                            data_temp)
    data_ave <- data_temp %>%
      group_by(Host) %>%
      summarise(across(where(is.numeric), ~ mean(., na.rm = TRUE)))
    data_ave <- as.data.frame(data_ave)
    rownames(data_ave) <- data_ave$Host
    data_ave$Host = NULL
    data_ave <- data_ave[host_list, ]
    
    mds_dist <- as.matrix(vegdist(data_ave, method = "horn", diag = TRUE, upper = TRUE))
    mds_dist <- mds_dist[order(rownames(mds_dist)), order(colnames(mds_dist))]
    dist_list <- c(dist_list, list(mds_dist))
  }
  return(list(dist_list, unique(Interaction_list)))
}

MANTEL_MDS_KINSHIP = function(mds_list, geno_kinship, outdir){
  res_p <- c()
  res_z <- c()
  for (i in 1:length(mds_list[[2]])){
    mds_dist <- 1 - mds_list[[1]][[i]]
    ROWCHECK(mds_dist, geno_kinship)
    mantel_result <- mantel.test(mds_dist, geno_kinship, nperm = 9999)
    res_p <- c(res_p, mantel_result$p)
    res_z <- c(res_z, mantel_result$z.stat)
  }
  out_df <- data.frame(Inoculant_Condition = mds_list[[2]],
                       z_score = res_z,
                       p_value = res_p,
                       row.names = "Inoculant_Condition")
  out_df <- out_df[order(rownames(out_df)), ]
  write.csv(out_df, Lc(outdir, "Beta_mantel_kinship.csv"))
  return(out_df)
}

GSM_ANALYSIS = function(data, sample, taxa, outdir, mab = 6){
  data <- data[, as.vector(apply(data, 2, function(x) sum(x != 0) > 6))]
  colnames(data) <- sapply(strsplit(colnames(data), split = "\\."), "[", 1)
  taxa <- taxa[colnames(data), ]
  
  out_df <- data.frame(matrix(nrow = ncol(data), ncol = ncol(taxa) + 8))
  colnames(out_df) <- c(colnames(taxa),
                        "Host", "Inoculant", "Condition",
                        "HostxInoculant", "HostxCondition", "InoculantxCondition",
                        "HostxInoculantxCondition", "Residuals")
  
  rm_id <- c()
  for (i in 1:ncol(data)){#
    temp_df <- data.frame(Bact = data[, i],
                          sample) 
    temp_res <- tryCatch({
      summary(gsm(Bact ~ Host*Inoculant*Condition, data = temp_df))
    }, error = function(e){
      print(Lc(i, " Error"))
      rm_id <<- c(rm_id, i)
      return(NULL)
    })
    if (!is.null(temp_res)){
      temp_res_p <- temp_res$s.table$`Pr(>F)`
      taxa_info <- as.vector(t(taxa[colnames(data)[i],]))
      out_list <- c(taxa_info, temp_res_p)
      out_df[i, ] <- out_list
    }
    print(i/ncol(data))
  }
  write.csv(out_df, Lc(outdir, "GSM_result.csv"))
}

GSM_FISHER = function(gsm_res, taxa, outdir, target = "family"){
  GSM_pval <- gsm_res[, 8:14]
  GSM_taxa <- gsm_res[, 1:7]
  
  taxa_unique <- unique(gsm_res[, target])
  variable_names <- colnames(GSM_pval)
  out_df <- data.frame(matrix(nrow = length(taxa_unique), ncol = length(variable_names)))
  rownames(out_df) <- taxa_unique
  colnames(out_df) <- variable_names
  out_df <- out_df[order(rownames(out_df)), ]
  
  for (i in taxa_unique){
    temp_dfs <- DF_ROWFILT(GSM_pval, GSM_taxa, GSM_taxa[, target] == i)
    if (nrow(temp_dfs[1][[1]]) < 3){
      next
    }
    temp_pval <- temp_dfs[1][[1]]
    temp_taxa <- temp_dfs[2][[1]]
    pval_list <- c()
    for (j in variable_names){
      a1 <- sum(temp_pval[, j] < 0.05) # Significant in target taxa
      a2 <- nrow(temp_pval) - a1 # Not significant in target taxa
      a3 <- sum(GSM_pval[, j] < 0.05) - a1 # Significant not in target taxa
      a4 <- nrow(GSM_pval) - nrow(temp_pval) - a3 # Not significant not in target taxa
      conting_mat <- matrix(c(a1, a2, a3, a4), nrow = 2)
      fish_result <- fisher.test(conting_mat, alternative = "greater")$p.value
      pval_list <- c(pval_list, fish_result)
    }
    out_df[i, ] <- pval_list
  }
  
  out_df2 <- na.omit(out_df)
  BH_name <- c()
  out_df_BH <- out_df2
  for (j in variable_names){
    out_df_BH[, j] <- p.adjust(out_df2[, j], "BH")
    BH_name <- c(BH_name, Lc(j, "_BH"))
  }
  colnames(out_df_BH) <- BH_name
  out_df2 <- cbind(out_df2, out_df_BH)
  
  for (i in colnames(out_df_BH)){
    print(Lc(i, " ", sum(out_df_BH[, i] < 0.05)))
  }
  
  write.csv(out_df2, Lc(outdir, "GSM_Fisher_family.csv"))
}

GSM_VENN = function(gsm_result, outdir){
  GSM_pval <- gsm_result[, 8:14]
  aG <- rownames(GSM_pval[GSM_pval$Host < 0.05, ])
  aC <- rownames(GSM_pval[GSM_pval$Inoculant < 0.05, ])
  aE <- rownames(GSM_pval[GSM_pval$Condition < 0.05, ])
  aGC <- rownames(GSM_pval[GSM_pval$HostxInoculant < 0.05, ])
  aGE <- rownames(GSM_pval[GSM_pval$HostxInoculant < 0.05, ])
  aCE <- rownames(GSM_pval[GSM_pval$InoculantxCondition < 0.05, ])
  aGCE <- rownames(GSM_pval[GSM_pval$HostxInoculantxCondition < 0.05, ])
  
  venn_temp <- list(G = aG, M = aC, E = aE)
  venn.diagram(venn_temp, filename = Lc(outdir, "Venn_Sole.png"),
               imagetype = "png", height = 1500, width = 1500, fill = c(4,7,5), disable.logging = T)
  venn.diagram(venn_temp, filename = Lc(outdir, "Venn_Sole.svg"),
               imagetype = "svg", height = 5, width = 5, fill = c(4,7,5), disable.logging = T)
  
  venn_temp <- list(G = aG, GxM = aGC, GxE = aGE, GxMxE = aGCE)
  venn.diagram(venn_temp, filename = Lc(outdir, "Venn_Geffects.png"),
               imagetype = "png", height = 1500, width = 1500, fill = c(4,7,5,6), disable.logging = T)
  venn.diagram(venn_temp, filename = Lc(outdir, "Venn_Geffects.svg"),
               imagetype = "svg", height = 5, width = 5, fill = c(4,7,5,6), disable.logging = T)
  
  venn_temp <- list(M = aC, E = aE, MxE = aCE, GxMxE = aGCE)
  venn.diagram(venn_temp, filename = Lc(outdir, "Venn_notGeffects.png"),
               imagetype = "png", height = 1500, width = 1500, fill = c(4,7,5,6), disable.logging = T)
  venn.diagram(venn_temp, filename = Lc(outdir, "Venn_notGeffects.svg"),
               imagetype = "svg", height = 5, width = 5, fill = c(4,7,5,6), disable.logging = T)
}

#####
# Phenotype
#####

HEATMAP_PLOT <- function(Sample_meta, Factors, Out_dir, out_prefix){
  temp_df <- data.frame(Sample_meta,
                        Factors = Sample_meta[[Factors]])
  temp_df_agg <- aggregate.data.frame(temp_df$Factors,
                                      by = list(temp_df$Host, temp_df$Inoculant, temp_df$Condition), mean)
  temp_df_agg$MxE <- interaction(temp_df_agg[, "Group.3"], temp_df_agg[, "Group.2"])
  temp_df_agg <- data.frame(Host = temp_df_agg[, "Group.1"],
                            name = temp_df_agg$MxE,
                            value = temp_df_agg$x)
  temp_mat <- as.data.frame(dcast(temp_df_agg,
                                  Host ~ name))
  temp_mat <- as.matrix(as.data.frame(temp_mat, row.names = temp_mat$Host)[, 2:ncol(temp_mat)])
  temp_mat <- scale(t(temp_mat))
  temp_mat_melt <- melt(temp_mat)
  
  p1 <- temp_mat_melt %>%
    ggplot(aes(Var2, Var1, fill = value))+
    geom_tile()+
    scale_fill_distiller()+
    ggtitle("Heatmap")+
    xlab("Host")+
    ylab("Condition_Inoculant") +
    labs(fill = "Values")
  
  ggsave(file = Lc(Out_dir, out_prefix, ".png"), plot = p1, dpi = 300, width = 12, height = 12)
}

GLM_ANALYSIS = function(model, data, out_name = "GLM_test", family = "Gamma", out = NA){
  # family check start
  if (family == "Gamma"){
    GLM_model <- glm(model, data = data,
                     family = Gamma)
  }else if (family == "Gamma(log)"){
    GLM_model <- glm(model, data = data,
                     family = Gamma(log))
  }else if (family == "gaussian"){
    GLM_model <- glm(model, data = data,
                     family = gaussian)
  }else{
    print("Error: Gamma, Gamma(log), and gaussian were supported as a family")
    return()
  }
  # family check end
  # GLM analysis start
  GLM_anova <- Anova(GLM_model, type = 2, test.statistic = "F")
  GLM_eta <- GLM_anova$"Sum Sq" / sum(GLM_anova$"Sum Sq")
  GLM_eta_part <- c()
  for (i in GLM_anova$"Sum Sq"){
    GLM_eta_part <- c(GLM_eta_part, i / (i + GLM_anova["Residuals", "Sum Sq"]))
  }
  GLM_out <- data.frame(Eta_squared = GLM_eta,
                        Partial_eta_squared = GLM_eta_part,
                        row.names = rownames(GLM_anova))
  GLM_out <- data.frame(GLM_anova, GLM_out)
  # GLM analysis end
  # QQ plot start
  qq_data <- data.frame(value = resid(GLM_model))
  qq_plot <- qq_data %>% 
    ggplot(aes(sample = value)) +
    stat_qq_band(distribution = "norm", dparams = list(mean = 0, sd = 1)) +
    stat_qq_line(distribution = "norm", dparams = list(mean = 0, sd = 1)) +
    stat_qq_point(distribution = "norm", dparams = list(mean = 0, sd = 1)) +
    theme(
      legend.position="none",
      plot.title = element_text(size=8)
    )+
    ggtitle(paste("Likelihood: ", round(logLik(GLM_model), 4))) +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
  # QQ plot end
  # Output start
  ggsave(file = paste(out_name, "_qqplot.png", sep = ""), plot = qq_plot, dpi = 300, width = 12, height = 12)
  ggsave(file = paste(out_name, "_qqplot.pdf", sep = ""), plot = qq_plot, width = 12, height = 12)
  write.csv(GLM_out, paste(out_name, "_anova.csv", sep = ""))
  
  if (!is.na(out)){
    sink(out, append = T)
    print(paste(out_name,  ": model is", model, sep = " "))
    print(paste(out_name,  ": family is", family, sep = " "))
    sink()
  }
  # Output end
  return(GLM_out)
}

TUKEY_ANALYSIS_FULLMODEL = function(model, data, out_name = "Tukey_temp", out = NA){
  model_anova <- aov(lm(model, data = data))
  Tukey_data <- TukeyHSD(model_anova)
  col_revised <- c("Difference", "Lower", "Uper", "P-adjustment", "Type")
  type_list <- c("Host", "Inoculant", "Condition", "Host_Inoculant", "Host_Condition","Inoculant_Condition", "Host_Inoculant_Condition")
  out_df <- data.frame(matrix(nrow = 0, ncol = 5))
  colnames(out_df) <- col_revised
  j <- 1
  for (i in Tukey_data){
    df_1 <- data.frame(i)
    df_1$Type <- type_list[j]
    colnames(df_1) <- col_revised
    out_df <- rbind(out_df, df_1)
    j <- j + 1
  }
  write.csv(out_df, paste(out_name, ".csv"))
  if (!is.na(out)){
    sink(out, append = T)
    print(paste(out_name,  ": Tukey model is", model, sep = " "))
    print(paste(out_name,  ": Tukey output is", paste(out_name, ".csv"), sep = " "))
    sink()
  }
  return(out_df)
}

DF_CHOSEPOTS = function(Sample_meta, rep = 100){
  temp_df <- data.frame(Sample_meta,
                        GxMxE = interaction(Sample_meta$Host, Sample_meta$Inoculant, Sample_meta$Condition))
  out_bool_list <- vector("list", rep)
  for (j in 1:rep){
    comb_list <- c()
    for (i in unique(temp_df$GxMxE)){
      comb_list <- c(comb_list, sample(unique(temp_df[temp_df$GxMxE == i, "CombID"]), 1))
    }
    out_bool <- temp_df$CombID %in% comb_list
    out_bool_list[[j]] <- out_bool
  }
  return(out_bool_list)
}

GLM_RANDOM = function(model, Sample_meta, Random_list){
  Out_list <- vector("list", length(Random_list))
  for (j in 1:length(Random_list)){
    data_temp <- Sample_meta[Random_list[[j]], ]
    glm_model <- glm(model, data = data_temp,family = Gamma(log))
    glm_anova <- Anova(glm_model, type = 2, test.statistic = "F")
    glm_eta <- glm_anova$"Sum Sq" / sum(glm_anova$"Sum Sq")
    glm_eta_part <- c()
    for (i in glm_anova$"Sum Sq"){
      glm_eta_part <- c(glm_eta_part, i / (i + glm_anova["Residuals", "Sum Sq"]))
    }
    glm_out <- data.frame(Eta_squared = glm_eta,
                          Partial_eta_squared = glm_eta_part,
                          row.names = rownames(glm_anova))
    glm_out <- data.frame(glm_anova, glm_out)
    Out_list[[j]] <- glm_out
  }
  return(Out_list)
}

RONDOM_ETA_PLOT <- function(glm_random, out_dir, out_prefix){
  # Aggregate
  outmat <- matrix(ncol = nrow(glm_random[[1]]), nrow = length(glm_random))
  for (i in 1:length(glm_random)){
    outmat[i,] <- glm_random[[i]]$Eta_squared
  }
  outmat <- data.frame(outmat)
  colnames(outmat) <- rownames(glm_random[[1]])
  # Plot
  reps <- nrow(outmat)
  etalist <- as.data.frame(pivot_longer(outmat, cols = colnames(outmat), names_to = "name", values_to = "value"))
  etalist <- transform(etalist, name= factor(name, levels = unique(etalist$name)))
  p1 <- etalist %>%
    ggplot(aes(x = name, y = value, fill = name)) +
    geom_violin(alpha = 0.8) +
    theme(legend.position = "none") +
    ggtitle("") +
    xlab("Effects") +
    ylab("Eta_squared")
  ggsave(file = Lc(out_dir, out_prefix, ".png"), plot = p1, dpi = 300, width = 12, height = 12)
  return(etalist)
}

PHENOTYPE_PPLOT <- function(Sample_meta, factors, out_dir, out_prefix){
  Sample_temp <- data.frame(Factors = Sample_meta[[factors]],
                            SL = Sample_meta$Length,
                            RL = Sample_meta$RootLength,
                            NOL = Sample_meta$Number_of_leaf,
                            NOB = Sample_meta$Number_of_branch)
  p1 <- Sample_temp %>%
    ggpairs(mapping = aes(color = Factors),
            upper = list(combo  = "blank",
                         continuas = wrap("cor", size = 5)),
            diag = list(discrete  = "blankDiag",
                        continuous = wrap("densityDiag", alpha = 0.3, size = 0.3)),
            lower = list(
              combo = function(data, mapping){ggplot(data = data, mapping = mapping) + geom_violin(alpha = 0.5)},
              continuous = wrap("points", alpha = 0.3)))
  
  ggsave(file = Lc(out_dir, out_prefix, ".png"), plot = p1, dpi = 300, width = 12, height = 12)
}

ICC_CALC = function(Sample_meta, factors, out_dir, out_prefix){
  Sample_meta$GxMxE <- interaction(Sample_meta$Host, Sample_meta$Inoculant, Sample_meta$Condition)
  Sample_meta <- data.frame(Sample_meta,
                            Factors = Sample_meta[[factors]])
  ICC <- c()
  for (i in unique(Sample_meta$GxMxE)){
    temp_df <- Sample_meta[Sample_meta$GxMxE == i, ]
    if (length(unique(temp_df$CombID)) != 1){
      temp_glmm <- glmer(Factors ~  (1|CombID), data = temp_df, family = Gamma(log), nAGQ=0)
      temp_var <- as.data.frame(VarCorr(temp_glmm))
      ICC <- c(ICC, (temp_var$vcov[1] / sum(temp_var$vcov)))
    }
  }
  
  glmm_total <- glmer(Factors ~  (1|CombID), data = Sample_meta, family = Gamma(log), nAGQ=0)
  glmm_total_vardf <- as.data.frame(VarCorr(glmm_total))
  Total_ICC <- glmm_total_vardf$vcov[1] / sum(glmm_total_vardf$vcov)
  
  out_df <- data.frame(ICC = mean(ICC),
                       Total_ICC = Total_ICC)
  
  write.csv(out_df, Lc(out_dir, out_prefix, ".csv"))
}

HERITABILITY = function(data, sample, target="Length"){
  Beta <- 1 - as.matrix(vegdist(data, method = "horn", diag = T))
  Pheno <- as.matrix(sample[, target])
  rownames(Pheno) <- rownames(sample)
  
  if (!ROWCHECK(Beta, Pheno)){
    print("The structures of data and sample were not matched.")
    return()
  }
  
  n1 <- nrow(Beta)
  Beta_stand <- (n1-1)/sum((diag(n1)-matrix(1,n1,n1)/n1)*Beta)*Beta
  
  # Host standardized
  for (i in unique(sample$Host)){
    temp_pheno <- Pheno[sample$Host == i, 1]
    Pheno[sample$Host == i, 1] <- scale(temp_pheno)
  }
  
  nullmodel <-lmm.aireml(Pheno,K=Beta_stand, EMsteps = 10000, verbose = F)
  Heritability <-nullmodel$tau/(nullmodel$tau+nullmodel$sigma2)
  print(Lc("Heritability is ", Heritability))
  return()
}


