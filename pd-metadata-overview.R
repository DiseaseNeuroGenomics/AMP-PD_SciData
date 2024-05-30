library(ggplot2)
library(data.table)
library(ggpubr)
library(dplyr)
library(ggsci)
library(UpSetR)
library(reshape)


setwd("./") # FIXME: SET TO THE FOLDER WITH UNCOMPRESSED REPOSITORY
dry_lab_dir = "dry_lab/"
wet_lab_dir = "wet_lab/"
outDir = "output_dir"
dir.create(outDir)
mpdf = function(x, width=7,height=7, outDir=outDir, onefile=T)eval.parent(substitute({ pdf(paste0(outDir,"/plot_",make.names(x),".pdf"),useDingbats=F,width=width,height=height,onefile=onefile) })) #outDir must be defined as a global var

##############################################################
# Color functions
col_color_palette <- function(vals, grp_var, legend_show = TRUE){
  if (class(vals) != 'factor') {
    vals <- as.factor(vals)
  }
  # Remove the 'first' color (black) in this sequence - 'black' is reserved for text
  my_colors <- palette.colors(nlevels(vals)+1, palette = "Okabe-Ito")
  my_colors <- tail(my_colors, n=-1)
  names(my_colors) <-levels(vals)
  if (class(legend_show) == "character") {
    return(scale_colour_manual(legend_show, name = grp_var,values = my_colors))
  } else if (legend_show) {
    return(scale_colour_manual(name = grp_var,values = my_colors))
  } else {
    return(scale_colour_manual("", name = grp_var,values = my_colors))
  }
  
}

# Fill function
fill_color_palette <- function(vals, grp_var, legend_show = TRUE){
  if (class(vals) != 'factor') {
    vals <- as.factor(vals)
  }
  # Remove the 'first' color (black) in this sequence - 'black' is reserved for text
  my_colors <- palette.colors(nlevels(vals)+1, palette = "Okabe-Ito")
  my_colors <- tail(my_colors, n=-1)
  names(my_colors) <-levels(vals)
  
  if (legend_show) {
    return(scale_fill_manual(name = grp_var,values = my_colors))
  } else {
    return(scale_fill_manual("", name = grp_var,values = my_colors))
  }
  
}


################################################################
######## Create Color palettes##################################

{
  ######
  
  # Brain bank color palette
  br_bnk_color_pal <- col_color_palette(c("Harvard", "MSSM", "Udall", "UMBEB"), "brain_bank")
  br_bnk_color_pal_noLeg <- col_color_palette(c("Harvard", "MSSM", "Udall", "UMBEB"), "brain_bank", legend_show = FALSE)
  
  # Brain region color palette
  br_reg_color_pal <- col_color_palette(c("DMNX", "GPI", "PMC", "PFC", "PVC"), "brain_region")
  br_reg_color_pal_noLeg <- col_color_palette(c("DMNX", "GPI", "PMC", "PFC", "PVC"), "brain_region", legend_show = FALSE)
  
  # Diagnosis color palette
  diag_color_pal <- col_color_palette(c("Control", "PD or Parkinsonism"), "diagnosis_final")
  diag_color_pal_noLeg <- col_color_palette(c("Control", "PD or Parkinsonism"), "diagnosis_final", legend_show = FALSE)
  
  # Diagnosis color palette
  sex_color_pal <- col_color_palette(c("Male", "Female"), "sex")
  sex_color_pal_noLeg <- col_color_palette(c("Male", "Female"), "sex", legend_show = FALSE)
  
  # Race color palette
  race_color_pal <- col_color_palette(c("AFR", "AMR", "EUR", "SAS"), "race")
  race_color_pal_noLeg <- col_color_palette(c("AFR", "AMR", "EUR", "SAS"), "race", legend_show = FALSE)
  
  # Dissection status palette
  dissect_color_pal <- col_color_palette(c("kept", "filtered", "removed"), "removed")
  dissect_color_pal_noLeg <- col_color_palette(c("kept", "filtered", "removed"), "removed", legend_show = FALSE)
  
  # Dissection status palette (v2)
  dissect2_color_pal <- col_color_palette(c("kept", "filtered", "removed"), "dissections")
  dissect2_color_pal_noLeg <- col_color_palette(c("kept", "filtered", "removed"), "dissections", legend_show = FALSE)
  
  ######
  
  ######
  
  # Brain bank color palette
  br_bnk_fill_pal <- fill_color_palette(c("Harvard", "MSSM", "Udall", "UMBEB"), "brain_bank")
  br_bnk_fill_pal_noLeg <- fill_color_palette(c("Harvard", "MSSM", "Udall", "UMBEB"), "brain_bank", legend_show = FALSE)
  
  # Brain region color palette
  br_reg_fill_pal <- fill_color_palette(c("DMNX", "GPI", "PMC", "PFC", "PVC"), "brain_region")
  br_reg_fill_pal_noLeg <- fill_color_palette(c("DMNX", "GPI", "PMC", "PFC", "PVC"), "brain_region", legend_show = FALSE)
  
  # Diagnosis color palette
  diag_fill_pal <- fill_color_palette(c("Control", "PD or Parkinsonism"), "diagnosis_final")
  diag_fill_pal_noLeg <- fill_color_palette(c("Control", "PD or Parkinsonism"), "diagnosis_final", legend_show = FALSE)
  
  # Diagnosis color palette
  sex_fill_pal <- fill_color_palette(c("Male", "Female"), "sex")
  sex_fill_pal_noLeg <- fill_color_palette(c("Male", "Female"), "sex", legend_show = FALSE)
  
  # Race color palette
  race_fill_pal <- fill_color_palette(c("AFR", "AMR", "EUR", "SAS"), "race")
  race_fill_pal_noLeg <- fill_color_palette(c("AFR", "AMR", "EUR", "SAS"), "race", legend_show = FALSE)
  
  # Dissection status palette
  dissect_fill_pal <- fill_color_palette(c("kept", "filtered", "removed"), "removed")
  dissect_fill_pal_noLeg <- fill_color_palette(c("kept", "filtered", "removed"), "removed", legend_show = FALSE)
  
  # Dissection status palette (v2)
  dissect2_fill_pal <- fill_color_palette(c("kept", "filtered", "removed"), "dissections")
  dissect2_fill_pal_noLeg <- fill_color_palette(c("kept", "filtered", "removed"), "dissections", legend_show = FALSE)
  
  ######
}
##############################################################


##############################################################
### LOAD & FIX METADATA (CLINICAL, ALIGNMENT) ################

allInfo = read.csv(file.path(wet_lab_dir, "R_plt_allInfo.csv"))
allInfoDf = read.csv(file.path(wet_lab_dir, "R_plt_allInfoDf.csv"))
upset_input = read.csv(file.path(dry_lab_dir, "R_UpSetPlot_input_NoPersIdent.csv"))
upset_input_wgs = read.csv(file.path(dry_lab_dir, "SciData_fig4.csv"))
notRmDf = read.csv(file.path(dry_lab_dir, "R_plt_notRmdf_NoPersIdent.csv"))
rboundDf = read.csv(file.path(dry_lab_dir, "R_plt_3rbinddf_NoPersIdent.csv"))

##############################################################
### PLOT DATASET CHARACTERISTICS / CLINICAL METADATA #########

{
  # Plot: Samples by brain bank
  # To include the 'removed' 3 donors
  # used allInfoDf rather than allInfo
  df = allInfoDf %>% count(brain_bank)
  df$labs = paste0(df$brain_bank, " (", df$n/sum(df$n)*100, "%)")
  df$brain_bank <- ordered(df$brain_bank, levels=c("Harvard", "MSSM", "Udall", "UMBEB"))
  myPlot = ggpie(df, "n", fill="brain_bank", lab.pos="in") + ggtitle("Brain bank")
  myPlot + br_bnk_fill_pal
  mpdf("brainbank", width=6, height=4); print(myPlot + br_bnk_fill_pal); dev.off()
  
  # Plot: Samples by diagnosis
  df = allInfoDf %>% count(diagnosis_final)
  myPlot = ggpie(df, "n", fill="diagnosis_final", lab.pos="in") + ggtitle("Diagnosis")
  myPlot + diag_fill_pal
  mpdf("diagnosis", width=4, height=4); print(myPlot + diag_fill_pal); dev.off()

  
  # Plot: Sex
  allInfoDf$sex = as.factor(allInfoDf$sex)
  df = allInfoDf %>% count(sex)
  myPlot = ggpie(df, "n", fill="sex", lab.pos="in") + ggtitle("Sex")
  myPlot + sex_fill_pal
  mpdf("sex", width=4, height=4); print(myPlot + sex_fill_pal); dev.off()
  
  # Plot: Sex by brainbank
  df = data.frame(allInfoDf %>% group_by(brain_bank, sex) %>% count("sex"))
  myPlot = ggplot(df, aes(fill=sex, y=n, x=brain_bank)) + geom_bar(position="fill", stat="identity") + theme(axis.text.x=element_blank()) + theme_classic() + ylab("Male -vs- Female ratio") + xlab("Brain bank") + ggtitle("Sex by brain bank")
  myPlot + sex_fill_pal
  mpdf("sex_by_brainbank", width=4, height=4); print(myPlot + sex_fill_pal); dev.off()
  
  # Plot: Race
  allInfoDf$race = as.factor(allInfoDf$race)
  # df = allInfoDf[, .N, by = race]
  df = allInfoDf %>% count(race)
  myPlot = ggpie(df, "n", fill="race", lab.pos="out") + ggtitle("Ethnicity")
  myPlot + race_fill_pal
  mpdf("race", width=4, height=4); print(myPlot + race_fill_pal); dev.off()

  # Plot: Braak by brainbank
  mu = plyr::ddply(allInfoDf, "brain_bank", summarise, grp.mean=mean(path_braak_lb))
  myPlot = ggplot(allInfoDf, aes(x=path_braak_lb, color=brain_bank, fill=brain_bank)) + geom_density(alpha=0.3) + geom_vline(data=mu, aes(xintercept=grp.mean, color=brain_bank), linetype="dashed") + labs(title="Braak by brain bank", x="Braak (LBD)", y = "Density") + theme_classic()
  myPlot + br_bnk_color_pal + br_bnk_fill_pal
  myPlot = ggplot(allInfoDf, aes(x=path_braak_lb, color=brain_bank, fill=brain_bank)) + geom_bar(position="stack") + labs(title="Braak by brain bank", x="Braak (LBD)", y = "Density") + theme_classic() # + geom_vline(data=mu, aes(xintercept=grp.mean, color=brain_bank), linetype="dashed") + labs(title="Braak by brain bank", x="Braak (LBD)", y = "Density") + theme_classic()
  myPlot + br_bnk_color_pal + br_bnk_fill_pal
  mpdf("braak_by_brainbank", width=4, height=4); print(myPlot + br_bnk_color_pal + br_bnk_fill_pal); dev.off()

  # Plot: Age by sex
  mu = plyr::ddply(allInfoDf, "sex", summarise, grp.mean=mean(age_at_baseline))
  myPlot = ggplot(allInfoDf, aes(x=age_at_baseline, color=sex, fill=sex)) + geom_density(alpha=0.3) + geom_vline(data=mu, aes(xintercept=grp.mean, color=sex), linetype="dashed") + labs(title="Age by sex", x="Age (baseline)", y = "Density") + theme_classic()
  myPlot + sex_color_pal + sex_fill_pal
  mpdf("age_by_sex", width=4, height=4); print(myPlot + sex_color_pal + sex_fill_pal); dev.off()
  
  # Plot: Age by diagnosis
  mu = plyr::ddply(allInfoDf, "diagnosis_final", summarise, grp.mean=mean(age_at_baseline))
  myPlot = ggplot(allInfoDf, aes(x=age_at_baseline, color=diagnosis_final, fill=diagnosis_final)) + geom_density(alpha=0.3) + geom_vline(data=mu, aes(xintercept=grp.mean, color=diagnosis_final), linetype="dashed") + labs(title="Age by diagnosis", x="Age (baseline)", y = "Density") + theme_classic()
  myPlot + diag_color_pal + diag_fill_pal
  mpdf("age_by_diagnosis", width=5, height=4); print(myPlot + diag_color_pal + diag_fill_pal); dev.off()
  
  # Plot: Diagnosis by brainbank
  df = data.frame(allInfoDf %>% group_by(brain_bank, diagnosis_final) %>% count("diagnosis_final"))
  myPlot = ggplot(df, aes(fill=diagnosis_final, y=n, x=brain_bank)) + geom_bar(position="fill", stat="identity") + theme(axis.text.x=element_blank()) + theme_classic() + ylab("Case -vs- Control ratio") + xlab("Brain bank") + ggtitle("Diagnosis by brain bank")
  myPlot + diag_fill_pal
  mpdf("diagnosis_by_brainbank", width=5, height=4); print(myPlot + diag_fill_pal); dev.off()
}

##############################################################
### PLOT DATASET CHARACTERISTICS / CLINICAL METADATA #########

{
  # Plot cell count per dissection stratified by the information whether the dissection was removed or kept
  df = rboundDf
  df$cell_count = as.numeric(df$cell_count)
  df$removed = ordered(gsub("TRUE", "removed", gsub("FALSE", "kept", gsub("filtered", "removed", df$removed))))
  df$removed = ordered(gsub("TRUE", "removed", gsub("FALSE", "kept", df$removed)), levels=c("kept", "removed"))
  
  colnames(df)[which(names(df) == "removed")] <- "dissections"
  mu = plyr::ddply(df, "dissections", summarise, grp.mean=mean(cell_count))
  myPlot = ggplot(df, aes(x=cell_count, color=dissections, fill=dissections)) + 
    geom_density(alpha=0.3) + dissect_color_pal + dissect_fill_pal + 
    geom_vline(data=mu, aes(xintercept=grp.mean, color=dissections), linetype="dashed") + 
    labs(title="Cell count per dissection", x="Cell count", y = "Density", color = "Was dissection removed?") + 
    theme_classic()
  myPlot
  mpdf("cellCount2", width=10, height=4); print(myPlot); dev.off()

  # Plot samples per brain region
  df = data.frame(notRmDf %>% group_by(brain_region) %>% count("cell_count"))
  df$brain_region = ordered(df$brain_region, levels=c("DMNX", "GPI", "PMC", "PFC", "PVC"))
  myPlot = ggplot(df, aes(x=brain_region, y=n, fill=brain_region)) + geom_bar(stat="identity", width=0.5) + geom_text(aes(label=n), vjust=-0.5) + 
    labs(title="Sample count per brain region", x="brain region", y = "sample count", fill = "brain region") + theme_classic()
  myPlot + br_reg_fill_pal
  mpdf("sampleCount_byBrainRegion", width=5, height=5); print(myPlot + br_reg_fill_pal); dev.off()
  
  # Plot avg cells per brain region
  df = plyr::ddply(notRmDf, "brain_region", summarise, grp.mean=mean(cell_count))
  df$grp.mean = floor(df$grp.mean)
  df$brain_region = ordered(df$brain_region, levels=c("DMNX", "GPI", "PMC", "PFC", "PVC"))
  myPlot = ggplot(df, aes(x=brain_region, y=grp.mean, fill=brain_region)) + geom_bar(stat="identity", width=0.5) + geom_text(aes(label=grp.mean), vjust=-0.5) + 
    labs(title="Average cell count per dissection", x="brain region", y = "cell count", fill = "brain region") + theme_classic()
  myPlot + br_reg_fill_pal
  mpdf("avgCellCountPerDiss_byBrainRegion", width=5, height=5); print(myPlot + br_reg_fill_pal); dev.off()
  
  # Plot total cells per brain region
  df = plyr::ddply(notRmDf, "brain_region", summarise, grp.mean=sum(cell_count))
  df$grp.mean = floor(df$grp.mean)
  df$brain_region = ordered(df$brain_region, levels=c("DMNX", "GPI", "PMC", "PFC", "PVC"))
  myPlot = ggplot(df, aes(x=brain_region, y=grp.mean, fill=brain_region)) + geom_bar(stat="identity", width=0.5) + geom_text(aes(label=grp.mean), vjust=-0.5) + 
    labs(title="Total cell count per brain region", x="brain region", y = "cell count", fill = "brain region") + theme_classic()
  myPlot + br_reg_fill_pal
  mpdf("totalCellCount_byBrainRegion", width=5, height=5); print(myPlot + br_reg_fill_pal); dev.off()
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat) {
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

{
  allInfo$diagnosis_binary = ifelse(as.character(allInfo$diagnosis_final) == "Control", 0, 1)
  corrMatrix = data.frame(allInfo[,c("diagnosis_binary", "path_braak_lb", "path_braak_nft")])
  df = corrMatrix
  corrMatrixOut = data.frame(sapply(1:ncol(corrMatrix), function(x) sapply(1:ncol(corrMatrix), function(y) {
    df = corrMatrix[!is.na(df[,x]) & !is.na(df[,y]),]
    cor.test(as.numeric(corrMatrix[,x]), as.numeric(corrMatrix[,y]), method="spearman")$estimate 
  } )))
  colnames(corrMatrixOut) = c("PD case/control","PD Braak", "AD Braak")
  rownames(corrMatrixOut) = colnames(corrMatrixOut)
  corrMatrixOut = get_upper_tri(corrMatrixOut)
  corrMatrixOut$ID = rownames(corrMatrixOut)
  corrMatrixOut = melt(corrMatrixOut, id.vars="ID", factorsAsStrings=T)
  corrMatrixOut$variable = ordered(corrMatrixOut$variable,levels=c("PD case/control","PD Braak", "AD Braak"))
  corrMatrixOut$ID = ordered(corrMatrixOut$ID,levels=c("PD case/control","PD Braak", "AD Braak"))
  
  corr = ggplot(corrMatrixOut, aes(variable, ID, fill = value)) + geom_tile() + scale_fill_material("red") +
    coord_equal() + theme_bw() + 
    theme(axis.title.x=element_blank(), axis.text.x=element_text(colour = "black", 
         angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(colour = "black"), 
         axis.ticks.x=element_blank(), axis.title.y=element_blank(),  axis.ticks.y=element_blank()
         )
  corr = corr + geom_text(aes(variable, ID, label = round(value, 2)), color = "black", size = 3)
  corr
  mpdf("pd_phenotypes_correlation_spearman", width=4, height=4); print(corr); dev.off() 
}

##############################################################
### PLOT DATASET CHARACTERISTICS / FIG. 4: QC OF WGS DATA ####

{
  # Fig. 4a: Distribution of mean coverage indicating the average number of high-quality sequencing reads per base after applying all QC steps.
  print(paste0("The mean mapped depth was ", round(mean(upset_input_wgs$wgsMetricsMean.Coverage)), " (sd: ", (round(sd(upset_input_wgs$wgsMetricsMean.Coverage),2)), "; ", round(min(upset_input_wgs$wgsMetricsMean.Coverage)), "-", round(max(upset_input_wgs$wgsMetricsMean.Coverage)), ")"))
  coverageDistrPlot = ggplot(upset_input_wgs, aes(x = wgsMetricsMean.Coverage)) + geom_density(fill = "#E69F00", alpha = 0.5) + labs(x = "Mean Coverage", y = "Density") + theme_classic() + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
  mpdf("fig_4a_coverageDistrPlot", width=4, height=3); print(coverageDistrPlot); dev.off()
  
  # Fig. 4b: The fraction of the genome sequenced at different depths. The center line (black) indicates the median, the box shows the interquartile range, and the whiskers indicate the highest/lowest values within 1.5× the interquartile range.
  print(paste0("The paired-end short reads mapped to around ", round(mean(upset_input_wgs$wgsMetricsPCT_1X),3), " (sd:", round(sd(upset_input_wgs$wgsMetricsPCT_1X),3), ") of human reference genome; 30x: ", round(mean(upset_input_wgs$wgsMetricsPCT_10X),3), " (", round(sd(upset_input_wgs$wgsMetricsPCT_10X),3) , ") "))
  tempDf = upset_input_wgs[grepl(pattern="wgsMetricsPCT_[0-9]", colnames(upset_input_wgs))]
  tempDf = tempDf[,1:7]
  colnames(tempDf) = gsub("wgsMetricsPCT_", "", colnames(tempDf))
  tempDf = reshape::melt(tempDf)
  coverageBoxPlot = ggplot(tempDf, aes(x=variable, y=value)) + labs(x="Coverage", y="Genome fraction") + 
    geom_boxplot() + theme_classic() + scale_y_continuous(limits = c(0.25, 1), expand = c(0, 0))
  coverageBoxPlot
  mpdf("fig_4b_coverageBoxPlot", width=4, height=3); print(coverageBoxPlot); dev.off()
  
  # Fig. 4c-d: Number of per-sample SNPs and indels.
  paste0("> Total SNPs per donor = ", round(mean(upset_input_wgs$wgsMetricsVariant_TOTAL_SNPS)), " (sd ", round(sd(upset_input_wgs$wgsMetricsVariant_TOTAL_SNPS)), ")") 
  paste0("> Total indels per donor = ", round(mean(upset_input_wgs$wgsMetricsVariant_TOTAL_INDELS)), " (sd ", round(sd(upset_input_wgs$wgsMetricsVariant_TOTAL_INDELS)), ")") 
  snpsPlot = ggplot(upset_input_wgs, aes(x = wgsMetricsVariant_TOTAL_SNPS)) + geom_density(fill = "#E69F00", alpha = 0.5) + labs(x = "Per-sample SNPs", y = "Density") + theme_classic() + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
  mpdf("fig_4c_snpsPlot", width=4, height=3); print(snpsPlot); dev.off()
  indelPlot = ggplot(upset_input_wgs, aes(x = wgsMetricsVariant_TOTAL_INDELS)) + geom_density(fill = "#E69F00", alpha = 0.5) + labs(x = "Per-sample indels", y = "Density") + theme_classic() + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
  mpdf("fig_4d_indelsPlot", width=4, height=3); print(indelPlot); dev.off()
  
  # Fig. 4f The first two PCs of genetic ancestry.
  genoAncestryMergedPlot = ggplot(upset_input_wgs, aes(x=ancestry_ALL.PC1, y=ancestry_ALL.PC2, color=Demographics.ethnicity)) + theme_bw() + geom_point() + race_color_pal
  genoAncestryMergedPlot
  mpdf("geno_ancestry_merged", width=6, height=4); print(genoAncestryMergedPlot); dev.off()
}

##############################################################
### UPSET PLOTS #########

{
  mylist <- list()
  for (i in unique(upset_input$brain_region)) {
    mylist[i] <- list(unique(upset_input[upset_input$brain_region == i, c('final_vcfID')]))
  }
  myPlot <- upset(fromList(mylist), order.by = 'freq', number.angles = 30, point.size = 3.5, line.size = 2, 
        mainbar.y.label = "Intersection of donors", sets.x.label = "Donors per brain region", 
        text.scale = c(1.3, 1.3, 1, 1, 2, 0.75))
  mpdf("sampleCount_byBrainRegion_Upset", width=7, height=5); print(myPlot); dev.off()
}

