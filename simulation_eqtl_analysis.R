##################################################
######## eQTL analysis of simulated data #########
##################################################

library(snpStats)
library(dplyr)
library(tidyr)
library(stringr)
library(plyr)
library(foreach)  
library(doParallel)
library(stringr)
library(reshape2)
library(reshape)
library(qvalue)

# Read simulated snp genotypes
snp.sim <- read.plink(file.path(data_dir,"simulated_genotypes"))
snp.mat <- snp.sim$genotypes
snp.df <- as.data.frame(snp.mat)

# Replace genotypes from 0, 1, 2 to 03, 02, 01
gt <- as.data.frame(sapply(snp.df, as.numeric))
rownames(gt) <- rownames(snp.df)
gt[gt == 3] <- 0
gt[gt == 2] <- 1
gt[gt == 1] <- 2
germ.geno <- t(gt)

# Simulate TF expression
set.seed(123) 
nsamples <- 670
ntfs <- 20000
TFexpr <- matrix(rnorm(nsamples*ntfs, mean = 0, sd = 1), nrow = nsamples, ncol = ntfs)
rownames(TFexpr) <- rownames(gt)
colnames(TFexpr) <- paste0("TF", 1:ntfs)
all_tfs <- colnames(TFexpr)


# Generate grid of effect sizes
tf_eff_all <- c(seq(-0.6, -0.1, 0.1), 0, seq(0.1, 0.6, 0.1))
snp_eff_all <- c(seq(-0.6, -0.1, 0.1), 0, seq(0.1, 0.6, 0.1))
tf_eff_add_all <- c(seq(-0.3, -0.1, 0.1), 0, seq(0.1, 0.3, 0.1))
tf_snp_eff_grid <- data.frame(expand.grid(tf_eff = tf_eff_all, snp_eff = snp_eff_all, tf_eff_add = tf_eff_add_all))

# Generate TG-TF-SNP trios for each grid row
ntgs <- 800
all_snps <- colnames(gt)
all_tfs <- colnames(TFexpr)
all_tgs <- paste0("TG",1:ntgs)
s_classes <- c("A", "B", "C", "D")
tg_snp_pairs <- data.frame(TG=all_tgs,SNP=all_snps)

grid.trios <- list()
for(r in 1:nrow(tf_snp_eff_grid)) {
	print(r)
	grid <- tf_snp_eff_grid[r, ]
	tf_eff <- grid$tf_eff
	snp_eff <- grid$snp_eff
	tf_eff_add <- grid$tf_eff_add
	sclass.trios <- list()
	for (s in 1:length(s_classes)) {
		print(s)
		sclass <- s_classes[[s]]
		snp.sclass <- grep(paste0(sclass,"_"),all_snps,value=TRUE)
		sclass_tg_snp_pairs <- tg_snp_pairs[which(tg_snp_pairs$SNP %in% snp.sclass), ]
		tg.sclass <- sclass_tg_snp_pairs$TG
		tf.sclass <- sample(all_tfs, 200)
		trio <- data.frame(TG = tg.sclass, SNP = snp.sclass,TF = tf.sclass, tf_eff = tf_eff, snp_eff = snp_eff, tf_eff_add = tf_eff_add, snp_class = sclass)
		sclass.trios[[s]] <- trio
	}
	grid.trios[[r]] <- do.call(rbind,sclass.trios)
}
all_trios <- do.call(rbind,grid.trios)

# Simulate intercept and error terms for the regression equation
nint=800
nerr=800
nsamples=670
intercept_df <- matrix(rnorm(nsamples*nint, mean = 0, sd = 1), nrow = nsamples, ncol = nint)
rownames(intercept_df) <- rownames(TFexpr)
colnames(intercept_df) <- c(paste0("IntA_", 0:199),paste0("IntB_", 0:199), paste0("IntC_", 0:199), paste0("IntD_", 0:199))

error_df <- matrix(rnorm(nsamples*nerr, mean = 0, sd = 0.1), nrow = nsamples, ncol = nerr)
rownames(error_df) <- rownames(TFexpr)
colnames(error_df) <- c(paste0("ErrA_", 0:199),paste0("ErrB_", 0:199), paste0("ErrC_", 0:199), paste0("ErrD_", 0:199))


# Function to run reg-eQTL for trios
regGLM <- function(trios,rnaExpr,germ.geno) {
	glm_res <- list()
	for(p in 1:nrow(trios)) {  
	    cur.trio=trios[p,]
	    # print(paste0("Trio", p))
	    TG=tg=as.character(cur.trio$TG)
	    # print(tg)
	    tf=TF=as.character(cur.trio$TF)
	    # print(tf)
	    SNP=snp=as.character(cur.trio$SNP)
	    # print(snp)
		maj.count <- 0
		min.count <- 0
		tot.count <- 0	
		maf <- 0

		# Subset germ.geno for rsID in snp_gene_tf_filter
	    germ.geno_sub <- germ.geno[which(rownames(germ.geno) == snp), ]
	    germ.geno_sub <- as.data.frame(germ.geno_sub)
	    names(germ.geno_sub) <- snp
		germ.geno_sub$sample_id <- rownames(germ.geno_sub)

	    # Get MAF and proceed if maf > 0.01
	    geno.alleles.tab <- as.data.frame(table(germ.geno_sub[,snp]))
	    if(nrow(geno.alleles.tab) > 1) {
		    geno.alleles <- as.character(geno.alleles.tab$Var1)
		    maj.count = sum(2*geno.alleles.tab[which(geno.alleles.tab$Var1 == 0), ]$Freq,geno.alleles.tab[which(geno.alleles.tab$Var1 == 1), ]$Freq)
		    min.count <- sum(2*geno.alleles.tab[which(geno.alleles.tab$Var1 == 2), ]$Freq,geno.alleles.tab[which(geno.alleles.tab$Var1 == 1), ]$Freq)
		    tot.count <- maj.count+min.count
		    maf <- min.count/tot.count
		    # print(maf)
			}
	    if(maf > 0.01 & maf < 1) {
			rnaExpr_sub <- rnaExpr[which(rownames(rnaExpr) %in% c(tg,tf)), ]
			rnaExpr_sub <- t(rnaExpr_sub)
			rnaExpr_sub <- as.data.frame(rnaExpr_sub)
			rnaExpr_sub$sample_id <- rownames(rnaExpr_sub)

			# Merge rna,genotype & clinical 
			glm_data_2 <- merge(rnaExpr_sub, germ.geno_sub, by="sample_id")

		    # subset glm_data
		    cur_glm_data <- glm_data_2[ ,colnames(glm_data_2) %in% c(tg,tf,snp)]
		    names(cur_glm_data)[names(cur_glm_data) == tg] <- "TG"
		    names(cur_glm_data)[names(cur_glm_data) == tf] <- "TF"
		    names(cur_glm_data)[names(cur_glm_data) == snp] <- "SNP"
		    cur_glm_data <- na.omit(cur_glm_data)

		    # Run glm 
		    glm_model <- glm(TG~TF+SNP+TF:SNP,data=cur_glm_data)
		    f2 <- as.data.frame(effectsize(glm_model, type = "f2"))
			x <- summary(glm_model)
		    vars <- rownames(x$coef)
		    # if(c("mrna") %in% vars & c("gene") %in% vars) {
		    x_res <- x$coefficients
		    x_res <- as.data.frame(x_res)
		    x_res <- x_res[rownames(x_res) != "(Intercept)",c("Estimate","Pr(>|t|)")]
		    names(x_res) <- c("coef","pval")
		    row.names(x_res)[row.names(x_res) %in% c("TF:SNP")] <- "TF_SNP"
		    x_res$cov <- rownames(x_res)
		    row.names(x_res) <- 1:nrow(x_res)
		    # Get lmer estimates
		    x_coeff <- x_res[,c("cov","coef")]
		    x_coeff_wide <- spread(x_coeff, key = cov, value = coef)
		    # Add missing columns
		    all_covs <- c("SNP", "TF", "TF_SNP")
		    covs.miss <- setdiff(all_covs, names(x_coeff_wide))
		    x_coeff_wide[covs.miss] <- NA 
		    x_coeff_wide <- x_coeff_wide[all_covs]
		    names(x_coeff_wide) <- paste0(names(x_coeff_wide),".coef")

		    # Get lmer pvalues
		    x_pval <- x_res[,c("cov","pval")]
		    x_pval_wide <- spread(x_pval, key = cov, value = pval)
		    # Add missing columns
		    x_pval_wide[covs.miss] <- NA 
		    x_pval_wide <- x_pval_wide[all_covs]
		    names(x_pval_wide) <- paste0(names(x_pval_wide),".pval")

		    # Cbind and Return
		    glm_all_res <- cbind(TG,TF,SNP,x_coeff_wide,x_pval_wide)
		    glm_all_res$MAF <- maf
		    glm_all_res$TF_f2 <- f2[which(f2$Parameter == "TF"), ]$Std_Coefficient
		    glm_all_res$SNP_f2 <- f2[which(f2$Parameter == "SNP"), ]$Std_Coefficient
		    glm_all_res$TF_SNP_f2 <- f2[which(f2$Parameter == "TF:SNP"), ]$Std_Coefficient
		    glm_res[[p]] <- glm_all_res
		}
	}
	glm_final_res <- do.call(rbind,glm_res)
	# qvalue
	snp_qobj <- qvalue(p=glm_final_res$SNP.pval, fdr.level=0.05, lambda=0, pi0.method="smoother")
	glm_final_res$SNP.qval <- snp_qobj$qvalues
	tf_qobj <- qvalue(p=glm_final_res$TF.pval, fdr.level=0.05, lambda=0, pi0.method="smoother")
	glm_final_res$TF.qval <- tf_qobj$qvalues
	tf_snp_qobj <- qvalue(p=glm_final_res$TF_SNP.pval, fdr.level=0.05, lambda=0, pi0.method="smoother")
	glm_final_res$TF_SNP.qval <- tf_snp_qobj$qvalues
	return(glm_final_res)
}

# Function to run s-eQTL for TG-SNP pairs
simGLM <- function(pairs,rnaExpr,germ.geno) {
	glm_simple_res <- list()
	for(p in 1:nrow(pairs)) {  # nrow(trios)
	    cur.pair=pairs[p,]
	    # print(paste0("Trio", p))
	    TG=tg=as.character(cur.pair$TG)
	    # print(tg)
	    SNP=snp=as.character(cur.pair$SNP)
	    # print(snp)
		maj.count <- 0
		min.count <- 0
		tot.count <- 0	
		maf <- 0

		# Subset germ.geno for rsID in snp_gene_tf_filter
	    germ.geno_sub <- germ.geno[which(rownames(germ.geno) == snp), ]
	    germ.geno_sub <- as.data.frame(germ.geno_sub)
	    names(germ.geno_sub) <- snp
		germ.geno_sub$sample_id <- rownames(germ.geno_sub)

	    # Get MAF and proceed if maf > 0.01
	    geno.alleles.tab <- as.data.frame(table(germ.geno_sub[,snp]))
	    if(nrow(geno.alleles.tab) > 1) {
		    geno.alleles <- as.character(geno.alleles.tab$Var1)
		    maj.count = sum(2*geno.alleles.tab[which(geno.alleles.tab$Var1 == 0), ]$Freq,geno.alleles.tab[which(geno.alleles.tab$Var1 == 1), ]$Freq)
		    min.count <- sum(2*geno.alleles.tab[which(geno.alleles.tab$Var1 == 2), ]$Freq,geno.alleles.tab[which(geno.alleles.tab$Var1 == 1), ]$Freq)
		    tot.count <- maj.count+min.count
		    maf <- min.count/tot.count
		    # print(maf)
		}
	    if(maf > 0.01 & maf < 1) {
	    	# Filter rnaExpr for TG
			rnaExpr_sub <- rnaExpr[which(rownames(rnaExpr) %in% c(tg)), ]
			# rnaExpr_sub <- t(rnaExpr_sub)
			rnaExpr_sub <- as.data.frame(rnaExpr_sub)
			names(rnaExpr_sub) <- tg
			rnaExpr_sub$sample_id <- rownames(rnaExpr_sub)

			# Merge rna,genotype & clinical 
			glm_data_2 <- merge(rnaExpr_sub, germ.geno_sub, by="sample_id")

		    # subset glm_data
		    cur_glm_data <- glm_data_2[ ,colnames(glm_data_2) %in% c(tg,snp)]
		    names(cur_glm_data)[names(cur_glm_data) == tg] <- "TG"
		    names(cur_glm_data)[names(cur_glm_data) == snp] <- "SNP"
		    cur_glm_data <- na.omit(cur_glm_data)


		    # Run glm models
		    glm_model <- glm(TG~SNP,data=cur_glm_data)
		    f2 <- as.data.frame(effectsize(glm_model, type = "f2"))
			x <- summary(glm_model)
		    vars <- rownames(x$coef)
		    # if(c("mrna") %in% vars & c("gene") %in% vars) {
		    x_res <- x$coefficients
		    x_res <- as.data.frame(x_res)
		    x_res <- x_res[rownames(x_res) != "(Intercept)",c("Estimate","Pr(>|t|)")]
		    names(x_res) <- c("coef","pval")
		    # snp_class <- row.names(x_res)[row.names(x_res) %in% c("SNP_Het1","SNP_Dom1","SNP_Rec1","SNP_Add")]
		    # row.names(x_res)[row.names(x_res) %in% c("SNP_Het1","SNP_Dom1","SNP_Rec1","SNP_Add")] <- "SNP"
		    # row.names(x_res)[row.names(x_res) %in% c("TF:SNP")] <- "TF_SNP"
		    x_res$cov <- rownames(x_res)
		    row.names(x_res) <- 1:nrow(x_res)
		    # Get lmer estimates
		    x_coeff <- x_res[,c("cov","coef")]
		    x_coeff_wide <- spread(x_coeff, key = cov, value = coef)
		    # Add missing columns
		    all_covs <- c("SNP")
		    covs.miss <- setdiff(all_covs, names(x_coeff_wide))
		    x_coeff_wide[covs.miss] <- NA 
		    x_coeff_wide <- x_coeff_wide[all_covs]
		    names(x_coeff_wide) <- paste0(names(x_coeff_wide),".coef")

		    # Get lmer pvalues
		    x_pval <- x_res[,c("cov","pval")]
		    x_pval_wide <- spread(x_pval, key = cov, value = pval)
		    # Add missing columns
		    x_pval_wide[covs.miss] <- NA 
		    x_pval_wide <- x_pval_wide[all_covs]
		    names(x_pval_wide) <- paste0(names(x_pval_wide),".pval")

		    # Cbind and Return
		    glm_all_res <- cbind(TG,SNP,x_coeff_wide,x_pval_wide)
		    glm_all_res$MAF <- maf
		    glm_all_res$SNP_f2 <- f2[which(f2$Parameter == "SNP"), ]$Std_Coefficient
		    glm_simple_res[[p]] <- glm_all_res
		}
	}
	glm_sim_final_res <- do.call(rbind,glm_simple_res)
	# qvalue
	snp_qobj <- qvalue(p=glm_sim_final_res$SNP.pval, fdr.level=0.05, lambda=0, pi0.method="smoother")
	glm_sim_final_res$SNP.qval <- snp_qobj$qvalues
	return(glm_sim_final_res)
}

# Run reg-eQTl and s-eQTL for each grid row gridi
#gridi=1
tf_eff = tf_snp_eff_grid[gridi, ]$tf_eff
snp_eff = tf_snp_eff_grid[gridi, ]$snp_eff
tf_eff_add = tf_snp_eff_grid[gridi, ]$tf_eff_add
snp_classes <- c("A", "B", "C", "D")

tg_file <- file.path(data_dir,"TG", paste0("TGexpr_TF",tf_eff,"_SNP",snp_eff,"_TFadd",tf_eff_add,".RDS"))
reg_file <- file.path(data_dir,"regGLM",paste0("regGLM_TF",tf_eff,"_SNP",snp_eff,"_TFadd",tf_eff_add,".RDS"))
sim_file <- file.path(data_dir,"simGLM",paste0("simGLM_TF",tf_eff,"_SNP",snp_eff,"_TFadd",tf_eff_add,".RDS"))
trio_file <- file.path(data_dir,"trios", paste0("Trios_TF",tf_eff,"_SNP",snp_eff,"_TFadd",tf_eff_add,".RDS"))

if(!(file.exists(reg_file))) {
	tg_expr_list <- list()
	reg_glm_list <- list()
	sim_glm_list <- list()
	trio_list <- list()
	for (s in 1:length(snp_classes)) {  
		sclass <- snp_classes[s]
		c_trios <- all_trios[which(all_trios$tf_eff == tf_eff & all_trios$snp_eff == snp_eff & all_trios$tf_eff_add == tf_eff_add & all_trios$snp_class == sclass), ]
		tf <- as.character(c_trios$TF)
		snp <- as.character(c_trios$SNP)
		tg <- as.character(c_trios$TG)

		# TG expr simulation
		intrcpt <- paste0("Int",snp)
		err <- paste0("Err",snp)
		c_intercept <- int_df[ ,intrcpt]
		c_error <- err_df[ ,err]
		tf_expr_mtx <- as.matrix(TFexpr[ ,tf])
		snp_gt_mtx <- as.matrix(gt[ ,snp])
		tf_mtx <- tf_expr_mtx*tf_eff
		snp_mtx <- snp_gt_mtx*snp_eff
		tf_snp_mtx <- tf_expr_mtx * snp_gt_mtx * tf_eff_add
		TGexpr <- c_intercept + tf_mtx + snp_mtx + tf_snp_mtx + c_error
		colnames(TGexpr) <- tg
		tg_expr_list[[s]] <- TGexpr

		# regGLM
		rnaExpr <- cbind(TGexpr,TFexpr)
		rnaExpr <- t(rnaExpr)
		reg_glm_res <- regGLM(c_trios,rnaExpr,germ.geno)
		reg_glm_res$TF_SNP <- paste0(reg_glm_res$TF,"_",reg_glm_res$SNP)
		reg_glm_res$TG_SNP <- paste0(reg_glm_res$TG,"_",reg_glm_res$SNP)
		reg_glm_res$snp_class <- gsub("_.*","",reg_glm_res$SNP)
		reg_glm_list[[s]] <- reg_glm_res

		# simGLM
		pairs <- distinct(c_trios[ ,c("TG", "SNP")])
		rnaExpr <- t(TGexpr)
		simple_glm_res <- simGLM(pairs,rnaExpr,germ.geno)
		simple_glm_res$TG_SNP <- paste0(simple_glm_res$TG,"_",simple_glm_res$SNP)
		simple_glm_res$snp_class <- gsub("_.*","",simple_glm_res$SNP)
		sim_glm_list[[s]] <- simple_glm_res
	}
	tg_expr_df <- do.call(cbind,tg_expr_list)
	saveRDS(tg_expr_df,tg_file)
	reg_glm_df <- do.call(rbind,reg_glm_list)
	saveRDS(reg_glm_df,reg_file)
	sim_glm_df <- do.call(rbind,sim_glm_list)
	saveRDS(sim_glm_df,sim_file)
}



# Combine eQTL results
grid.reg.all <- list()
grid.sim.all <- list()

for(g in 1:nrow(eff_grid)) {
	cur.grid <- eff_grid[g, ]
	tf_eff <- cur.grid$tf_eff
	snp_eff <- cur.grid$snp_eff
	tf_eff_add <- cur.grid$tf_eff_add

	reg_file <- file.path(data_dir,"regGLM",paste0("regGLM_TF",tf_eff,"_SNP",snp_eff,"_TFadd",tf_eff_add,".RDS"))
	sim_file <- file.path(data_dir,"simGLM",paste0("simGLM_TF",tf_eff,"_SNP",snp_eff,"_TFadd",tf_eff_add,".RDS"))

	reg_glm <- readRDS(reg_file)
	reg_glm$SNP_EFF <- snp_eff
	reg_glm$TF_EFF <- tf_eff
	reg_glm$TF_EFF_ADD <- tf_eff_add

	sim_glm <- readRDS(sim_file)
	sim_glm$SNP_EFF <- snp_eff
	sim_glm$TF_EFF <- tf_eff
	sim_glm$TF_EFF_ADD <- tf_eff_add

	grid.reg.all[[g]] <- reg_glm
	grid.sim.all[[g]] <- sim_glm
}

grid.reg.df <- do.call(rbind,grid.reg.all)
grid.sim.df <- do.call(rbind,grid.sim.all)

grid.reg.df <- grid.reg.df %>%
  mutate(
    case = if_else(SNP_EFF != 0, "positive", "negative"))
grid.reg.df <- grid.reg.df %>%
  mutate(
    truth = if_else(case == "positive", 1, 0)
  )
grid.reg.df <- grid.reg.df %>%
  mutate(
    predict = if_else(SNP.qval < 0.05, 1, 0))


grid.sim.df <- grid.sim.df %>%
  mutate(
    case = if_else(SNP_EFF != 0, "positive", "negative"))
grid.sim.df <- grid.sim.df %>%
  mutate(
    truth = if_else(case == "positive", 1, 0)
  )
grid.sim.df <- grid.sim.df %>%
  mutate(
    predict = if_else(SNP.qval < 0.05, 1, 0))


### ROC - AUC ###
list.by.maf.sim = split(grid.sim.df, grid.sim.df$snp_class)
list.by.maf.reg = split(grid.reg.df, grid.reg.df$snp_class)
aucs = c()
pdf(file.path(data_dir,'plot.rocs.pdf'), width=20, height=15)
par(mfrow=c(4, 5))
for(m in 1:length(list.by.maf.sim)) {
	df.sim = list.by.maf.sim[[m]]
	df.sim$tf.snp.eff = paste(df.sim$TF_EFF, df.sim$TF_EFF_ADD, df.sim$SNP_EFF, sep=':')
	df.sim.list = split(df.sim, df.sim$tf.snp.eff)

	df.reg = list.by.maf.reg[[m]]
	df.reg$tf.snp.eff = paste(df.reg$TF_EFF,  df.reg$TF_EFF_ADD, df.reg$SNP_EFF,  sep=':')
	df.reg.list = split(df.reg, df.reg$tf.snp.eff)

	nms = names(df.reg.list)
	for(nm in nms) {
		nm.split <- strsplit(nm, ":")[[1]]
		if(!(nm.split[[3]] == 0)) {
			nm.split[[3]] <- 0
			nm.0 = paste(nm.split, collapse = ":")
			cat(nm, '...')
			pos.sim = df.sim.list[[nm]]
			neg.sim = df.sim.list[[nm.0]]
			df.sim = rbind(pos.sim, neg.sim)
			roc_curve <- roc(predictor=1-df.sim$SNP.qval, response=df.sim$truth)
			plot(roc_curve, main = paste0("ROC Curve [", names(list.by.maf.sim)[m], ' ', nm, ']'), col = "blue", lwd = 2)
			pos.reg = df.reg.list[[nm]]
			neg.reg = df.reg.list[[nm.0]]
			df.reg = rbind(pos.reg, neg.reg)
			roc_curve_reg <- roc(predictor=1-df.reg$SNP.qval, response=df.reg$truth)
			plot(roc_curve_reg, col = "red",  lwd = 2,add = TRUE) 
			roc_test <- roc.test(roc_curve, roc_curve_reg, method = "delong")
			aucs = rbind(aucs, c(roc_curve$auc, roc_curve_reg$auc, roc_test$p.value, names(list.by.maf.sim)[m], nm))
		}
	}
}
dev.off()
ss = do.call(rbind, strsplit(aucs[, 5], ':'))
auc.df = data.frame(auc.sim=as.numeric(aucs[,1]), auc.reg=as.numeric(aucs[,2]), delong.pval = as.numeric(aucs[,3]), maf=aucs[,4], tf.eff=as.numeric(ss[, 1]), tf.eff.add=as.numeric(ss[, 2]), snp.eff=as.numeric(ss[, 3]))
auc.df$diff = auc.df$auc.reg - auc.df$auc.sim
auc.df$color = ifelse(auc.df$diff == 0, '0', ifelse(auc.df$diff > 0, '+', '-'))
saveRDS(auc.df,file.path(data_dir,"REG_SIM_AUC_COMPARE.RDS"))


#### Precision - Recall - F1 ####
q.cutoff <- seq(0,1,0.01)
list.by.maf.sim = split(grid.sim.df, grid.sim.df$snp_class)
list.by.maf.reg = split(grid.reg.df, grid.reg.df$snp_class)
metrics = c()
for(m in 1:length(list.by.maf.sim)) {
	df.sim = list.by.maf.sim[[m]]
	df.sim$tf.snp.eff = paste(df.sim$TF_EFF, df.sim$TF_EFF_ADD, df.sim$SNP_EFF, sep=':')
	df.sim.list = split(df.sim, df.sim$tf.snp.eff)

	df.reg = list.by.maf.reg[[m]]
	df.reg$tf.snp.eff = paste(df.reg$TF_EFF,  df.reg$TF_EFF_ADD, df.reg$SNP_EFF,  sep=':')
	df.reg.list = split(df.reg, df.reg$tf.snp.eff)

	nms = names(df.reg.list)
	for(nm in nms) {
		nm.split <- strsplit(nm, ":")[[1]]
		if(!(nm.split[[3]] == 0)) {
			nm.split[[3]] <- 0
			nm.0 = paste(nm.split, collapse = ":")
			cat(nm, nm.0, '...')

			pos.nm.sim = df.sim.list[[nm]]
			neg.nm.sim = df.sim.list[[nm.0]]
			df.nm.sim = rbind(pos.nm.sim, neg.nm.sim)
			pos.nm.reg = df.reg.list[[nm]]
			neg.nm.reg = df.reg.list[[nm.0]]
			df.nm.reg = rbind(pos.nm.reg, neg.nm.reg)

			for(q in 1:length(q.cutoff)) {
				q_thr <- q.cutoff[[q]]

				# Simple
				df.nm.sim <- df.nm.sim %>%
								mutate(predict = if_else(SNP.qval < q_thr, 1, 0))
				df.nm.sim$predict <- factor(df.nm.sim$predict, levels = c(0,1))
				sim.cf <- table(Actual = df.nm.sim$truth, Predicted = df.nm.sim$predict)
				TP <- sim.cf[2, 2]
				TN <- sim.cf[1, 1]
				FP <- sim.cf[1, 2]
				FN <- sim.cf[2, 1]
				sim.precision <- TP / (TP + FP)
				sim.precision
				sim.recall <- TP / (TP + FN)
				sim.specificity <- TN / (TN + FP)
				sim.accuracy <- (TP + TN) / (TP + TN + FP + FN)
				sim.roc_curve <- roc(predictor=1-df.sim$SNP.qval, response=df.sim$truth)
				sim.auc <- sim.roc_curve$auc

				# Reg
				df.nm.reg <- df.nm.reg %>%
								mutate(predict = if_else(SNP.qval < q_thr, 1, 0))
				df.nm.reg$predict <- factor(df.nm.reg$predict, levels = c(0,1))
				reg.cf <- table(Actual = df.nm.reg$truth, Predicted = df.nm.reg$predict)
				TP <- reg.cf[2, 2]
				TN <- reg.cf[1, 1]
				FP <- reg.cf[1, 2]
				FN <- reg.cf[2, 1]
				reg.precision <- TP / (TP + FP)
				reg.recall <- TP / (TP + FN)
				reg.specificity <- TN / (TN + FP)
				reg.accuracy <- (TP + TN) / (TP + TN + FP + FN)
				reg.roc_curve <- roc(predictor=1-df.reg$SNP.qval, response=df.reg$truth)
				reg.auc <- reg.roc_curve$auc

				metrics = rbind(metrics, c(q_thr, nm, names(list.by.maf.sim)[m], sim.precision, reg.precision, 
										sim.recall, reg.recall, 
										sim.specificity, reg.specificity, 
										sim.accuracy, reg.accuracy, 
										sim.auc, reg.auc))
			}
		}
	}
}
ss = do.call(rbind, strsplit(metrics[, 2], ':'))
metrics.df = data.frame(q_thr=as.numeric(metrics[,1]), eff.grid=metrics[,2],  
								tf.eff=as.numeric(ss[, 1]), tf.eff.add=as.numeric(ss[, 2]), snp.eff=as.numeric(ss[, 3]), maf=metrics[,3], 
								sim.prec=as.numeric(metrics[,4]), reg.prec=as.numeric(metrics[,5]), 
								sim.sens=as.numeric(metrics[,6]), reg.sens=as.numeric(metrics[,7]), 
								sim.spec=as.numeric(metrics[,8]), reg.spec=as.numeric(metrics[,9]), 
								sim.acc=as.numeric(metrics[,10]), reg.acc=as.numeric(metrics[,11]),
								sim.auc=as.numeric(metrics[,12]), reg.auc=as.numeric(metrics[,13]))

maf.classes <- c("A","B","C","D")
metrics_res <- list()
for(m in 1:length(maf.classes)) {
	maf <- maf.classes[[m]]
	m.file <- paste0("MAF_",maf,"REG_SIM_PERF_COMPARE.RDS")
	print(m.file)
	metric.file <- readRDS(file.path(data_dir,m.file))
	metrics_res[[m]] <- metric.file
}
metrics.df <- do.call(rbind,metrics_res) # 4368 unique dataset - eff.grid+maf


#### calculate highest F1 score by cohens effect sizes ####
metrics.df.2 <- metrics.df[which(!(metrics.df$q_thr %in% c(0,1))), ]

max_reg_f1 <- metrics.df.2 %>%
  group_by(eff.grid, maf) %>%
  filter(reg.f1 == max(reg.f1)) %>%
  arrange(q_thr) %>%
  slice(1) %>%
  select(q_thr, eff.grid, tf.eff, tf.eff.add, snp.eff, maf, reg.f1)
max_reg_f1 <- as.data.frame(max_reg_f1)
names(max_reg_f1) <- c("reg_q_thr", "eff.grid", "tf.eff", "tf.eff.add", "snp.eff", "maf", "max.reg.f1")

max_sim_f1 <- metrics.df.2 %>%
  group_by(eff.grid, maf) %>%
  filter(sim.f1 == max(sim.f1)) %>%
  arrange(q_thr) %>%
  slice(1) %>%
  select(q_thr, eff.grid, tf.eff, tf.eff.add, snp.eff, maf, sim.f1)
max_sim_f1 <- as.data.frame(max_sim_f1)
names(max_sim_f1) <- c("sim_q_thr", "eff.grid", "tf.eff", "tf.eff.add", "snp.eff", "maf", "max.sim.f1")

# Merge sim and reg f1 scores
max_f1 <- merge(max_reg_f1,max_sim_f1)

# Merge f1 scores to cohen's
mean_f2_df <- readRDS(file.path(data_dir,"Mean_Cohens_effect_size.RDS"))
pos_mean_f2 <- mean_f2_df[which(mean_f2_df$snp_eff != 0), c("snp_class", "eff.grid", "tf_f2",
													"snp_f2.reg", "snp_f2.sim", "tf_snp_f2")]
names(pos_mean_f2)[names(pos_mean_f2) == "snp_class"] <- "maf"
cohens_f1_df_1 <- merge(max_f1,pos_mean_f2)
cohens_f1_df_1 <- cohens_f1_df_1 %>%
  mutate(cohens.reg = case_when(
  	abs(snp_f2.reg) < 0.15 ~ "small",
    abs(snp_f2.reg) >= 0.15 & abs(snp_f2.reg) < 0.35 ~ "medium",
    abs(snp_f2.reg) >= 0.35 ~ "large"
  ))

