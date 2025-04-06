getwd()
rm(list = ls())
setwd("./PX-LASSO-SVM-BROUTA-XG")
if (! dir.exists('./machine')){
  dir.create('./machine')
}
setwd('./machine')

library(magrittr)
library(ggplot2)
library(lance)

# intersect <- read.csv('./01.Intersection_gene20.csv') %>% as.data.frame()
# intersect <- data.frame(x=intersect$symbol)
# dat_expr<-read.csv('./01.count_mRNA_normalized.csv',row.names = 1,check.names = F) %>% lc.tableToNum()
# group <- read.csv('./01.group.csv')
# colnames(group) <- c('sample', 'type')
# 
# dat<-dat_expr[intersect$x,group$sample]%>%t%>%as.data.frame()
# dat<-merge(dat,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
# table(group$type)
# dat$type<-factor(dat$type,levels = c('PZ','TZ'))
# #dat$type <- ifelse(dat$type=='PZ',0,1)

save(dat, file = "dat.Rdata")
load("dat.Rdata")

#LASSO----------
library(glmnet)
?cv.glmnet
set.seed(1)
library(glmnet)
cv_fit <- cv.glmnet(as.matrix(dat[-ncol(dat)]),
                    dat$type, family = "binomial",
                    type.measure = "deviance",
                    alpha=1,
                    nfolds = 5)
?cv.glmnet
#alpha 对回归的类型进行判别,1代表Lasso回归，0代表岭回归，0-1之间代表弹性网络

# 提取最优lambda下的系数
coef.min <- coef(cv_fit, s = "lambda.min")
active.min <- which(coef.min@i != 0)
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i + 1]
lasso_geneids <- lasso_geneids[-1]
lasso_geneids
lasso_geneids.write <- data.frame(symbol=lasso_geneids)
write.csv(lasso_geneids.write,file = '01.lasso.gene.csv',row.names = F,quote = F)

cv_model  <- cv_fit
best_lambda <- cv_model$lambda.min
lambda_lse <- cv_model$lambda.1se
lambda_values <- cv_model$lambda
coef_list <- lapply(lambda_values, function(lambda) coef(cv_model, s = lambda))
coef_matrices <- lapply(coef_list, function(coef) as.matrix(coef))
coef_df_list <- lapply(coef_matrices, function(matrix) as.data.frame(t(matrix)))
coef_df <- do.call(rbind, coef_df_list)
coef_df$lambda <- rep(lambda_values, each = nrow(coef_df_list[[1]]))
coef_df <- coef_df[, -1]
coef_df <- coef_df[!rownames(coef_df) %in% "(Intercept)", ]
colnames(coef_df) <- c(colnames(coef_df)[-ncol(coef_df)], "lambda")
library(reshape2)
coef_long <- melt(coef_df, id.vars = "lambda", variable.name = "Variable", value.name = "Coefficient")
coef_long$log_lambda <- log(coef_long$lambda)
# 加载必要的包
library(ggplot2)
library(ggsci)
table(coef_long$Variable)
p1 <- ggplot(coef_long, aes(x = log_lambda, y = Coefficient, color = Variable)) + 
  geom_vline(xintercept = log(best_lambda), size = 0.8, color = 'grey60', alpha = 0.8, linetype = 2) +
  # geom_vline(xintercept = log(lambda_lse), size = 0.8, color = 'grey60', alpha = 0.8, linetype = 2) +
  geom_line(size = 1) +
  xlab("Log Lambda") + 
  ylab("Coefficients") + 
  theme_bw(base_rect_size = 2) +
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 15, color = 'black'), 
        axis.text = element_text(size = 12, color = 'black'), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12, color = 'black'), 
        legend.position = 'right') +
  annotate('text', x = log(best_lambda),
           y = 5, 
           label = paste0('log(lambda.min)\n', round(log(best_lambda), 4)), 
           color = 'black', size = 5, family = 'serif') +
  # annotate('text', x = log(lambda_lse),
  #          y = 10, 
  #          label = paste0('log(lambda.lse)\n', round(log(lambda_lse), 4)), 
  #          color = 'black', size = 5, family = 'serif') +
  guides(col = guide_legend(ncol = 2))
p1
ggsave('02.lasso.Coef.pdf',p1,w=11,h=7)
ggsave('02.lasso.Coef.png',p1,w=11,h=7)

xx <- data.frame(
  lambda = cv_fit$lambda,
  cvm = cv_fit$cvm,
  cvsd = cv_fit$cvsd,
  cvup = cv_fit$cvup,
  cvlo = cv_fit$cvlo,
  nozezo = cv_fit$nzero
)
xx$ll <- log(xx$lambda)
xx$NZERO <- paste0(xx$nozezo, ' vars')
dev.off()
p2 <- ggplot(xx, aes(ll, cvm, color = NZERO)) +
  geom_errorbar(aes(ymin = cvlo, ymax = cvup), width = 0.05, size = 1) +
  geom_vline(xintercept = log(best_lambda), size = 0.8, color = 'grey60', alpha = 0.8, linetype = 2) +
  # geom_vline(xintercept = log(lambda_lse), size = 0.8, color = 'grey60', alpha = 0.8, linetype = 2) +
  geom_point(size = 2) +
  xlab("Log Lambda") +
  ylab('Partial Likelihood Deviance') +
  theme_bw(base_rect_size = 1.5) +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 15, color = 'black'),
    axis.text = element_text(size = 12, color = 'black'),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, color = 'black'),
    legend.position = 'bottom'
  ) +
  annotate('text', x = log(best_lambda), y = 1.6, 
           label = paste0('log(lambda.min)\n', round(log(best_lambda), 4)), 
           color = 'black', size = 4, family = 'serif') +
  # annotate('text', x = log(lambda_lse), y =1.8, 
  #          label = paste0('log(lambda.1se)\n', round(log(lambda_lse), 4)), 
  #          color = 'black', size = 4, family = 'serif') +
  # guides(col = guide_legend(ncol = 9))
  guides(col = guide_legend(ncol = 8))
p2
ggsave('03.lasso.CV.pdf',p2,w=11,h=7)
ggsave('03.lasso.CV.png',p2,w=11,h=7)

###调整结果
result <- data.frame()
j <- 1
for (j in c(1:50)) {
  set.seed(j)
  res.lasso <- cv.glmnet(as.matrix(dat[-ncol(dat)]), dat$type, family = "binomial",
                         type.measure = "deviance",alpha=1)
  l.coef<-coef(res.lasso$glmnet.fit,s=res.lasso$lambda.min,exact= F)
  l.coef
  coef.min = coef(res.lasso, s = "lambda.min")  ## lambda.min & lambda.1se 取一个
  res.lasso$lambda.min
  active.min = which(coef.min@i != 0)
  lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1] %>% as.data.frame()
  print(j)
  ifelse (nrow(lasso_geneids) == 0,
          result[1,j] <- NA,
          result[c(1:nrow(lasso_geneids)),j] <- lasso_geneids)
  colnames(result)[j] <- paste('seed', j)
  
  re <- unique(t(result)) %>% as.data.frame()
  write.csv(re,file = '04.lasso2000.csv',quote = F,row.names = T)
}

re2 <- read.csv('./04.lasso2000.csv',row.names = 1)
re2


#svm-rfe----------
rm(list = ls())
setwd("./PX-LASSO-SVM-BROUTA-XG")
if (! dir.exists('./machine')){
  dir.create('./machine')
}
setwd('./machine')

library(magrittr)
library(tibble)
# 加载库和包
library(mlbench)
library(caret)
library(lance)

load("dat.Rdata")

?rfe
set.seed(2)
dat1 <- dat
##除目标变量外的特征数量
num <- ncol(dat1)-1
##递归特征消除（RFE）的控制参数
control <- rfeControl(functions = caretFuncs,method = "cv", number = 5)

# 执行SVM-RFE算法
results <- rfe(dat1[,1:num],            #数据框格式
               as.factor(dat1[,num+1]),   #因子格式
               sizes = c(1:num),
               rfeControl = control,
               method = "svmRadial"
)

save.image(file = 'svmrfe.RData')
load(file = 'svmrfe.RData')

# # 结果分析
svmrfe_result <- predictors(results)
svmrfe_result
length(svmrfe_result)
svmrfe_result.write <- data.frame(symbol=svmrfe_result)
write.csv(svmrfe_result.write,"05.svmrfe_gene.csv",row.names = F,quote = F)

# 绘制结果
pdf(file = paste0("06.SVM_RFE_Accuracy.pdf"),width = 5,height = 4,family='Times')
a <- dev.cur()  
png(file = paste0("06.SVM_RFE_Accuracy.png"),width = 5, height=4, units="in", res=600,family='Times')
dev.control("enable")
par(mar = c(2,2,2,2));
plot(results, type=c("o"),
     xgap.axis = 1
)

dev.copy(which = a)  
dev.off()
dev.off()

#循环运行很慢，挂后台
result <- data.frame()
i <- 1
for (i in c(1:500)) {
  print(i)
  set.seed(i)
  dat1 <- dat
  num <- ncol(dat1)-1
  control <- rfeControl(functions = caretFuncs,method = "cv", number = 5)
  results <- rfe(dat1[,1:num],            
                 as.factor(dat1[,num+1]),  
                 sizes = c(1:num),
                 rfeControl = control,
                 method = "svmRadial"
  )
  print(results)
  svmrfe_result <- predictors(results)
  svmrfe_result
  ifelse (length(table(svmrfe_result)) == 0,
          result[1,i] <- NA,
          result[c(1:length(table(svmrfe_result))),i] <- svmrfe_result)
  colnames(result)[i] <- paste('seedd', i)
  re <- unique(t(result)) %>% as.data.frame()
  write.csv(re,file = '07.svm500.csv',quote = F,row.names = T)
}
re2 <- read.csv('./07.svm500.csv',row.names = 1)


# # Boruta ------------------------------------------------------------------
rm(list = ls())
setwd("./PX-LASSO-SVM-BROUTA-XG")
if (! dir.exists('./machine')){
  dir.create('./machine')
}
setwd('./machine')


load("dat.Rdata")

library("Boruta")
?Boruta
set.seed(44242424)
res.Boruta<-Boruta(x=dat[,1: ncol(dat)-1], 
                   y=as.factor(dat[,ncol(dat)]),
                   pValue=0.01, #置信水平（p 值）
                   maxRuns=100 #特征重要性计算的最大运行次数
                   )
Boruta<-attStats(res.Boruta) #给出Boruta算法的结果
#通过访问finalDecision属性,
#可以获取每个特征是否被认为是对模型有贡献的("Confirmed"),
#还是被认为是噪声("Rejected"),
#或者是未确定的("Tentative")
table(res.Boruta$finalDecision)
boruta_geneids<-Boruta[Boruta$decision=='Confirmed',]%>%rownames(.)
boruta_geneids
write.csv(boruta_geneids,"08.boruta_gene.csv",row.names = F,quote = F)

##定义一个函数提取每个变量对应的重要性值。
# x=res.Boruta
library(dplyr)
boruta.imp <- function(x){
  imp <- reshape2::melt(x$ImpHistory, na.rm=T)[,-1]
  colnames(imp) <- c("Variable","Importance")
  imp <- imp[is.finite(imp$Importance),]
  variableGrp <- data.frame(Variable=names(x$finalDecision),
                            finalDecision=x$finalDecision)
  showGrp <- data.frame(Variable=c("shadowMax", "shadowMean", "shadowMin"),
                        finalDecision=c("shadowMax", "shadowMean", "shadowMin"))
  variableGrp <- rbind(variableGrp, showGrp)
  boruta.variable.imp <- merge(imp, variableGrp, all.x=T)
  sortedVariable <- boruta.variable.imp %>% group_by(Variable) %>%
    summarise(median=median(Importance)) %>% arrange(median)
  sortedVariable <- as.vector(sortedVariable$Variable)
  boruta.variable.imp$Variable <- factor(boruta.variable.imp$Variable, levels=sortedVariable)
  invisible(boruta.variable.imp)
}

boruta.variable.imp <- boruta.imp(res.Boruta)

head(boruta.variable.imp)
library(YSX)
library(ImageGP)

sp_boxplot(boruta.variable.imp, melted=T, xvariable = "Variable", yvariable = "Importance",
           legend_variable = "finalDecision", legend_variable_order = c("Tentative", "Confirmed", "Rejected", "shadowMax", "shadowMean", "shadowMin"),
           xtics_angle = 90)

pdf("09.Boruta.pdf",w = 7, h = 6,family='Times')
sp_boxplot(boruta.variable.imp, melted=T, xvariable = "Variable", yvariable = "Importance",
           legend_variable = "finalDecision", legend_variable_order = c("Tentative", "Confirmed", "Rejected", "shadowMax", "shadowMean", "shadowMin"),
           xtics_angle = 90)
dev.off()

png("09.Boruta.png",w = 7, h = 6,family='Times',units='in',res=600)
sp_boxplot(boruta.variable.imp, melted=T, xvariable = "Variable", yvariable = "Importance",
           legend_variable = "finalDecision", legend_variable_order = c("Tentative", "Confirmed", "Rejected", "shadowMax", "shadowMean", "shadowMin"),
           xtics_angle = 90)
dev.off()


#xgboost -----------------------------------------------------------------
rm(list = ls())
setwd("./PX-LASSO-SVM-BROUTA-XG")
if (! dir.exists('./machine')){
  dir.create('./machine')
}
setwd('./machine')

load("dat.Rdata")

library('xgboost')
library("Matrix")
library('PRROC')
library('ggplot2')

set.seed(1765727)
colnames(dat) <- gsub("-", "_", colnames(dat))
train_matrix <- sparse.model.matrix(type ~.-1, data = dat)
traindata1<- data.matrix(dat[,-ncol(dat)])
traindata2 <- Matrix(traindata1,sparse = T)
train_y <- as.numeric(dat[,ncol(dat)])-1 # 将因变量转换为numeric类型，-1是为了从0开始计数

traindata <- list(data=traindata2,label=train_y)
dtrain <- xgb.DMatrix(data = traindata$data, label = traindata$label)
set.seed(989655)
res.xgb <- xgboost(data = dtrain,max_depth=2, eta=0.3,
                   objective='binary:logistic', nround=25)
# eta：range: [0,1]，参数值越大，越可能无法收敛。把学习率 eta 设置的小一些，小学习率可以使得后面的学习更加仔细。
# max_depth：每颗树的最大深度，树高越深，越容易过拟合
# nround：迭代次数，RMSE(均方根误差)开始增加，是过度拟合训练数据的迹象。
# objective：定义最小化损失函数类型，常用参数binary:logistic,multi:softmax,multi:softprob

xgb_importance <- xgb.importance(train_matrix@Dimnames[[2]], model = res.xgb)    ##特征重要度
xgb_importance $Feature<-gsub('.','-',xgb_importance $Feature,fixed = T )

ggplot(xgb_importance, aes(x= reorder( Feature,Gain), y=Gain,fill='blue')) +
  geom_bar(stat="identity") +
  theme_classic() +
  guides(fill=FALSE)+
  #theme(legend.position = )+
  #geom_hline(yintercept = mean(xgb_importance$Gain),lty = 4,col = "darkred",lwd = 0.8) +
  coord_flip()+
  theme_bw()+
  ggtitle('XGBoost')+
  theme(plot.title = element_text(size=24,color='black', face = "bold",family='Times'),
        axis.title.x =element_text(size=18,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=16, color='black', face = "bold",family='Times'),
        axis.title.y =element_blank(),
        axis.text.y=element_text(size=16,   color='black',face = "bold",family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18, color='black', face = "bold",family='Times'),
        title=element_text(size=20, color='black', face = "bold",family='Times'),
        strip.text = element_text(size = 14,family = "Times", face = "bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x="gene",y="Gain",fill="")
ggsave('10.XGBoost_importance.pdf',w=8,h=8)
ggsave('10.XGBoost_importance.png',w=8,h=8)
xgb.write <- data.frame(symbol=xgb_importance$Feature)
xgb.write$symbol
write.csv(xgb.write,'11.XGBoost.gene.csv',quote=F,row.names=F)


