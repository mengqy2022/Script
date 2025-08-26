# OPENGWAS_JWT=''
usethis::edit_r_environ()
ieugwasr::get_opengwas_jwt()
ieugwasr::user()


library(clusterProfiler)
library(TwoSampleMR)
library(ggplot2)

DECRGs <- read.csv('../data/02_DECRGs/DECRGs.csv')
DECRGs_enterzID <- bitr(DECRGs$DECRGs,fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = "org.Hs.eg.db")  
DECRGs_enterzID$ENSEMBL <- paste0('eqtl-a-',DECRGs_enterzID$ENSEMBL)
write.csv(DECRGs_enterzID,'../data/02_DECRGs/DECRGs_enterzID.csv',row.names = F)
#提取暴露BMI的工具变量 ####
bmi_exp_dat <- data.frame()
for (i in DECRGs_enterzID$ENSEMBL) {
  skip_to_next <- FALSE
  tryCatch({
    bmi <- extract_instruments(outcomes = i,clump = TRUE, p1 = 5e-8, p2 = 5e-8,kb = 100)
    bmi_exp_dat <- rbind(bmi_exp_dat, bmi)
    # write.csv(bmi_exp_dat,'../data/03_MR/bmi_exp_dat.csv',row.names = F)
  }, error = function(e) {
    message("Error: ", i)
    skip_to_next <<- TRUE
  })
  
  if(skip_to_next) { next }
  
}

##补充F和R
# Calculate R2 and F stat for exposure data
# Liberal hypertension F stat
exposure_dat$r2 <- (2 * (exposure_dat$beta.exposure^2) * exposure_dat$eaf.exposure * (1 - exposure_dat$eaf.exposure)) /
  (2 * (exposure_dat$beta.exposure^2) * exposure_dat$eaf.exposure * (1 - exposure_dat$eaf.exposure) +
     2 * exposure_dat$samplesize.exposure * exposure_dat$eaf.exposure * (1 - exposure_dat$eaf.exposure) * exposure_dat$se.exposure^2)

exposure_dat$F <- exposure_dat$r2 * (exposure_dat$samplesize.exposure - 2) / (1 - exposure_dat$r2)
write.csv(exposure_dat,'../data/03_MR/exposure_dat.csv',row.names = F)

# bmi_exp_dat <- read.csv('../data/03_MR/bmi_exp_dat.csv',row.names = 1)
#结局变量 ####
exposure_dat <- bmi_exp_dat
chd_out_dat <- extract_outcome_data(snps = exposure_dat$SNP, outcomes = 'finn-b-N14_CALCUKIDUR')
write.csv(chd_out_dat,'../data/03_MR/outcome.csv',quote = F)
#将暴露和结局的数据进行合并，产生用于进行MR分析的数据

dat <- harmonise_data(
  exposure_dat=exposure_dat,
  outcome_dat=chd_out_dat,
  action= 2)

write.csv(dat,'../data/03_MR/dat.csv',quote = F)


#MR分析的主要结果:默认用5种方法进行MR分析 并计算OR值####
OR <- generate_odds_ratios(mr(dat))
write.csv(OR,'../data/03_MR/MR_res.csv')

#提取IVW中p小于0.05的变量 ####
IVW_OR <- subset(OR,OR$pval < 0.05 & OR$method == 'Inverse variance weighted')
write.csv(IVW_OR,'OR_p0.05_or1.csv',quote = F)

#数据敏感性分析
#异质性检验
het <- mr_heterogeneity(dat)
write.csv(het,'异质性检验.csv',quote = F)
#多效性检验
pleio <- mr_pleiotropy_test(dat)
write.csv(pleio,'多效性检验.csv',quote = F)

#Steiger方向检验
Steiger_results <- data.frame()
for(i in IVW_OR$id.exposure){
  steiger_data <- subset(dat,id.exposure == i)
  out <- directionality_test(steiger_data)
  Steiger_results <- rbind(Steiger_results,out)
  
}

write.csv(Steiger_results,'Steiger方向检验.csv',quote = F)

path <- '../data/04_MR/plot/'
for(i in mr_res$id.exposure){
  if (!file.exists(paste0(path, i))) {
    dir.create(paste0(path, i))
  }
  
  gene_res <- subset(mr_res,mr_res$id.exposure == i)
  gene_dat <- subset(mr_dat,mr_dat$id.exposure == i)
  
  single <- mr_leaveoneout(gene_dat)
  p1 <- mr_leaveoneout_plot(single)[[1]]
  p1 <- p1+scale_color_brewer(palette = 'Set1')+
    theme_classic(base_size = 18)+
    theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
          axis.ticks = element_line(size = 1),
          legend.position = 'none')+
    update_geom_defaults("line", list(size = 5))+
    labs(x='AD')
  
  pdf(paste0(path,i,'/leave-one-out forest plot.pdf'),width = 6,height = 8)
  print(p1)
  dev.off()
  
  png(file = paste0(path,i,'/leave-one-out forest plot.png'),w=6,h=8,units = 'in',res = 600)
  print(p1)
  dev.off()
  
  p2 <- mr_scatter_plot(gene_res,gene_dat)[[1]]
  p2$labels$x <- unlist(strsplit(p2$labels$x,split = ' \\|\\| '))[1]
  # p2$labels$y <- 'SNP effect on AD'
  p2 <- p2+scale_color_brewer(palette = 'Set1')+
    theme_classic(base_size = 18)+
    theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
          axis.ticks = element_line(size = 1),
          legend.position = 'top')+
    update_geom_defaults("line", list(size = 5))
  
  pdf(paste0(path,i,'/scatter plot.pdf'),width = 8,height = 6)
  print(p2)
  dev.off()
  
  png(file = paste0(path,i,'/scatter plot.png'),w=6,h=8,units = 'in',res = 600)
  print(p2)
  dev.off()
  
  res_single <- mr_singlesnp(gene_dat)
  res_single <- res_single[-nrow(res_single),]
  p3 <- mr_forest_plot(res_single)[[1]]
  p3 <- p3+scale_color_brewer(palette = 'Set1')+
    theme_classic(base_size = 18)+
    theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
          axis.ticks = element_line(size = 1),
          legend.position = 'none')+
    update_geom_defaults("line", list(size = 5))
  
  pdf(paste0(path,i,'/forest plot.pdf'),width = 8,height = 6)
  print(p3)
  dev.off()
  
  png(file = paste0(path,i,'/forest plot.png'),w=8,h=6,units = 'in',res = 600)
  print(p3)
  dev.off()
  
  
  p4 <- mr_funnel_plot(res_single)[[1]]
  p4 <- p4+
    scale_color_brewer(palette = 'Set1')+
    theme_classic(base_size = 18)+
    theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
          axis.ticks = element_line(size = 1),
          legend.position = 'top')+
    update_geom_defaults("line", list(size = 5))
  
  pdf(paste0(path,i,'/funnel plot.pdf'),width = 8,height = 6)
  print(p4)
  dev.off()
  
  png(file = paste0(path,i,'/funnel plot.png'),w=8,h=6,units = 'in',res = 600)
  print(p4)
  dev.off()
}
