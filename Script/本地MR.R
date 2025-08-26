library(TwoSampleMR)
library(ieugwasr)
library(plinkbinr)
library(dplyr)
library(data.table)
library(gwasglue)
# library(clusterProfiler)
# library(stringr)
# library(org.Hs.eg.db)

setwd('/data/nas2/MR_results/pvalue5e-8_kb10000/')
if (! dir.exists("./Type.2.diabetes")){
  dir.create("./Type.2.diabetes")
}
setwd("./Type.2.diabetes")
##读取IEU openGWAS数据库中全部的基因里面

ids <- read.table("/data/nas2/database/MR/eQTL_vcf/eQTL_list") %>%
  select(V1) %>%      # 选择第一列（默认列名为 V1）
  pull()              # 提取为向量

# #### GWAScatlog 下载的结局文件处整理
# outcome<-readr::read_tsv('/data/nas2/database/MR/outcome_files/GWAScatlog/GCST90018917_buildGRCh37.tsv.gz')
# outcome$n_cases <- 647
# outcome$n_controls <- 482264
# 
# #如果SNP没有信息，可以转换。根据chr:pos查找RS号  https://blog.csdn.net/wendy_milk/article/details/121642843
# outcome$`chromosome:start` <- paste(outcome$chromosome,outcome$base_pair_location,sep=':')
# ### tes = read.table("test.txt",header=T,check.names=F,sep="\t") #注意这里我设置的是制表符分隔符，如果你的文件不是制表符的话，需要修改成对应的分隔符
# match <- fread("/data/nas2/database/MR/snp150_hg19.txt.gz",header=T,check.names=F,sep="\t")
# # system.time(save(match, file = "/data/nas2/database/MR/outcome_files/match.Rdata"))
# need <- dplyr::left_join(outcome,match,by="chromosome:start") #如果snp150_hg19.txt文件中有对应的RS号，则比对到test.txt文件中，如果没有的话，就变为NA
# ##使用format_data()函数将该数据框转化成TwoSampleMR的格式
# save(need, file='need.RData')
# outcome_dat <- format_data(
#   dat=need,
#   type = "outcome",
#   snps = NULL,
#   snp_col = "name",
#   beta_col = "beta",
#   se_col = "standard_error",
#   effect_allele_col = "effect_allele",
#   other_allele_col = "other_allele",
#   pval_col = "p_value",
#   ncase_col = "n_cases",
#   ncontrol_col = "n_controls",
#   chr_col = "chromosome",
#   pos_col = "base_pair_location"
# )
# 
# outcome_all <- outcome_dat

#### IEU openGWAS直接下载的vcf文件处理 
# outcome_vcf <- 'ukb-e-250_CSA.vcf'
# vcf2 <- readVcf(outcome_vcf)
# outcome <- gwasvcf_to_TwoSampleMR(vcf2,type = "outcome")
# outcome_all <- outcome

### finngen数据库下载的结局文件整理
outcome<-fread('/data/nas2/database/MR/outcome_files/FinnGen/finngen_R12_T2D.gz') %>% as.data.frame()
ALCOLIV_outcome <- outcome
#colnames(outcome) <- c("chrom","pos","ref","alt","rsids","nearest_genes","pval","mlogp","beta","sebeta","af_alt","af_alt_cases","af_alt_controls")
#outcome$pval = 10^(-outcome$pval)
outcome$n = 82878	+	403489
outcome$n_cases = 82878
outcome$n_controls = 403489

outcome_all <- format_data(
  dat=outcome,
  type = "outcome",
  snps = NULL,
  #header = TRUE,
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  pval_col = "pval",
  samplesize_col = 'n',
  ncase_col = "ncase",
  ncontrol_col = "ncontrol",
  chr_col = "#chrom",
  pos_col = "pos"
)
table(outcome_all$mr_keep.outcome)
range(outcome_all$pval.outcome)



outcome_all$id.outcome <- "Type.2.diabetes"
outcome_all$outcome <- "Type.2.diabetes || id:finngen_R12_T2D"
system.time(save(outcome_all, file = "outcome_dat.Rdata"))

dat.res <- list()
expo.res <- list()
steiger.res <- list()

results_all<-data.frame()
hetpvalue_all <- data.frame()
pelpvalue_all <- data.frame()

final_res <- data.frame(exposure = NA, method = NA, pvalue = NA, 
                        het.qvalue = NA, pel.pvalue = NA, steiger.direction = NA, steiger.pvalue = NA)
i <- 1
ids <- ids[-c(1638,2597,4976,6095,8317,9257,12955)]
for (i in 1:length(ids)){
  print(i)
  ####暴露因素
  expos_vcf = paste('/data/nas2/database/MR/eQTL_vcf/',ids[i],'.vcf.gz',sep='')
  vcf1 <- VariantAnnotation::readVcf(expos_vcf,"hg19")
  exposure_dat <- gwasvcf_to_TwoSampleMR(vcf1,type = "exposure")
  exposure_dat$id.exposure = ids[i]
  exposure_dat$exposure = paste(' || id:',ids[i],sep='')
  exp_data <- exposure_dat
  temp_dat <- exp_data[exp_data$pval.exposure<5e-08,]
  if (nrow(temp_dat) < 3) {
    next
  }
  # 补充ld_clump需要的三列
  temp_dat$id <- temp_dat$id.exposure
  temp_dat$rsid <- temp_dat$SNP
  temp_dat$pval <- temp_dat$pval.exposure
  # 执行ld_clump  bfile指定参考文件的路径
  # exp_dat <- ld_clump(temp_dat,
  #                     plink_bin = get_plink_exe(),
  #                     bfile = '/data/nas2/database/MR/1000g/EUR/EUR' ,
  #                     clump_kb =10000, clump_r2 = 0.001) 
  
  skip_to_next <- FALSE
  
  tryCatch({
    exp_dat <- ld_clump(temp_dat,
                        plink_bin = get_plink_exe(),
                        bfile = '/data/nas2/database/MR/1000g/EUR/EUR',
                        clump_kb = 10000, 
                        clump_r2 = 0.001)
  }, error = function(e) {
    message("Error: ", ids[i],",去连锁不平衡后没SNP了")
    skip_to_next <<- TRUE
  })
  
  if(skip_to_next) { next }
  
  ## 结局变量
  outcome_dat<- subset(outcome_all,outcome_all$SNP %in% exp_dat$SNP)
  if (nrow(outcome_dat) == 0) {
    next
  }
  dat <- harmonise_data(exposure_dat =exp_dat,outcome_dat = outcome_dat)
  
  dat <- subset(dat,mr_keep == TRUE)
  ###判断F值---
  clumping = mutate(dat,R=get_r_from_bsen(dat$beta.exposure, dat$se.exposure, dat$samplesize.exposure))
  clumping = mutate(clumping,F=(samplesize.exposure-2)*((R*R)/(1-R*R)))
  dat = clumping[clumping$F >10,]
  
  #判断snp数是否少于3个
  if (dim(dat)[1] < 3) {
    next
  }
  
  #异质性检验
  het <- mr_heterogeneity(dat)
  
  
  #多效性检验
  ple <- mr_pleiotropy_test(dat)
  pelpvalue_all <- rbind(pelpvalue_all,ple)
  
  #方向性检验
  steiger_sl<-directionality_test(dat)
  
  expo.res[[i]] <- exp_dat
  dat.res[[i]] <- dat
  steiger.res[[i]] <- steiger_sl
  
  #判断方向性检验的p值、水平多效性检验的p值、异质性检验的p值是否为NA
  if(is.na(steiger_sl$steiger_pval) | is.na(ple$pval) | is.na(het$Q_pval[2])){
    next
  }
  
  
  #判断水平多效性检验和方向性检验是否能过
  if (ple$pval < 0.05 | steiger_sl$correct_causal_direction == FALSE | steiger_sl$steiger_pval > 0.05) {
    next
  }
  
  #判断是否有异质性
  if(het$Q_pval[2] < 0.05){
    result <- generate_odds_ratios(mr_res = mr(dat,method_list = 'mr_ivw_mre'))
    #重新进行异质性检验
    het <- mr_heterogeneity(dat,method_list = 'mr_ivw_mre')
    hetpvalue_all <- rbind(hetpvalue_all,het)
  }else{
    result <- generate_odds_ratios(mr_res = mr(dat,method_list = 'mr_ivw'))
    hetpvalue_all <- rbind(hetpvalue_all,het)
  }
  
  results_all<-rbind(results_all,result)
  
  #判断基因MR结果中的 IVW算法对应的P值是否显著
  if(result$method %in% c('Inverse variance weighted','Inverse variance weighted (multiplicative random effects)') & result$pval < 0.05){
    res <- data.frame(exposure = ids[i], method = result$method, pvalue = result$pval, 
                      het.qvalue = het$Q_pval[2], pel.pvalue = ple$pval, 
                      steiger.direction = steiger_sl$correct_causal_direction, 
                      steiger.pvalue = steiger_sl$steiger_pval)
    final_res <- rbind(final_res,res)
  }
  write.csv(final_res,'02.results_sig.csv')
}

save.image(file = 'Type.2.diabetes.Rdata')
#合并数据并保存

expo.all <- do.call(rbind,expo.res)
dat.all <- do.call(rbind,dat.res)
steiger.all <- do.call(rbind,steiger.res)

write.csv(expo.all,file="expo.all.csv")
write.csv(dat.all,file="dat.all.csv")
write.csv(steiger.all,file = 'steiger.all.csv')

results_all <- results_all[results_all$nsnp > 2,]
write.csv(results_all,'01.results_all.csv')

final_res <- na.omit(final_res)
write.csv(final_res,'02.results_sig.csv')


# pelpvalue <- pelpvalue_all[pelpvalue_all$pval>0.05,]
# pelpvalue$ENSEMBL <- pelpvalue$id.exposure
# pelpvalue <- merge(pelpvalue,id,by='ENSEMBL')
# 
# gene <- intersect(results$SYMBOL,pelpvalue$SYMBOL)
# print(gene)
# write.csv(gene,file="gene.csv")





# ####反向--------------------------------
# ###table exposure
# exposure<-read.table('finngen_R9_T2D.gz',sep = "\t")
# ALCOLIV_exposure <- exposure
# # exposure <- ALCOLIV_exposure
# colnames(exposure) <- c("chrom","pos","ref","alt","rsids","nearest_genes","pval","mlogp","beta","sebeta","af_alt","af_alt_cases","af_alt_controls")
# exposure_all <- format_data(
#   dat=exposure,
#   type = "exposure",
#   snps = NULL,
#   header = TRUE,
#   snp_col = "rsids",
#   beta_col = "beta",
#   se_col = "sebeta",
#   eaf_col = "af_alt",
#   effect_allele_col = "alt",
#   other_allele_col = "ref",
#   pval_col = "pval",
#   # ncase_col = "n_cases",
#   # ncontrol_col = "n_controls",
#   chr_col = "chrom",
#   pos_col = "pos"
# )
# 
# exposure_all$id.exposure <- "T2DM"
# exposure_all$exposure <- " || id:finngen_R9_T2D"
# system.time(save(exposure_all, file = "exposure_all.Rdata"))
# 
# 
# temp_dat <- exposure_all[exposure_all$pval.exposure<5e-05,]
# temp_dat$id <- temp_dat$id.exposure
# temp_dat$rsid <- temp_dat$SNP
# temp_dat$pval <- temp_dat$pval.exposure
# # 执行ld_clump  bfile指定参考文件的路径
# exp_dat <- ld_clump(temp_dat,
#                     plink_bin = get_plink_exe(),
#                     bfile = '/data/nas1/yuanyt/pipeline/MR/g1000_eur/g1000_eur' ,
#                     clump_kb =10, clump_r2 = 0.001) 
# 
# 
# 
# 
# results_all<-data.frame()
# 
# ids = id$ENSEMBL[which(id$SYMBOL %in% gene)]
# 
# for (i in 1:length(ids))  {
#   ####暴露因素
#   out_vcf = paste('/data/nas1/luchunlin/epQTL_coloc/eQTL_vcf/',ids[i],'.vcf.gz',sep='')
#   vcf1 <- VariantAnnotation::readVcf(out_vcf,"hg19")
#   outcome_dat <- gwasvcf_to_TwoSampleMR(vcf1,type = "outcome")
#   outcome_dat$id.outcome = ids[i]
#   outcome_dat$outcome = paste(' || id:',ids[i],sep='')
#   outcome_dat<- subset(outcome_dat,outcome_dat$SNP %in% exp_dat$SNP)
#   if (nrow(outcome_dat) == 0) {
#     next
#   }
#   dat <- harmonise_data(exposure_dat =exp_dat,outcome_dat = outcome_dat)
#   result <- generate_odds_ratios(mr_res = mr(dat))
#   if (result$nsnp[1] < 3) {
#     next
#   }
#   results_all<-rbind(results_all,result)
#   # clumping = mutate(dat,R=get_r_from_bsen(dat$beta.exposure, dat$se.exposure, dat$samplesize.exposure))
#   # clumping = mutate(clumping,F=(samplesize.exposure-2)*((R*R)/(1-R*R)))
#   # clumping
# }
# 
# 
# write.csv(results_all,'01.reverse_mr.csv')
