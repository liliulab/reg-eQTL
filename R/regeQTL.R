#' Process reg-eQTL Analysis
#'
#' This function runs reg-eQTL analysis on expression data, genotype data, and TG-TF-SNV trios.
#' @param expr.data Expression data matrix.
#' @param cov.data Covariates data matrix.
#' @param trio.data TG-TF-SNV trio information.
#' @param gt.data Genotype data matrix.
#' @param out.dir Directory to save the output.
#' @return The function outputs the reg-eQTL results into the specified directory.
#' @export


process.regeqtl = function(expr.data, cov.data, trio.data, gt.data, out.dir) {
    tryCatch(
        {
            process.regeqtl.internal(expr.data, cov.data, trio.data, gt.data, out.dir)
        }, 
        error = function(err) {
            message(paste0("Error running regeQTL: ", err$message))  # Print the specific error
        },
        warning = function(warn) {
            message(paste0("Warning in regeQTL: ", warn$message))
        },
        finally = {
            message("regeQTL completed")
        }
    )
}


regGLM <- function(cur_glm_data,covs,all_covs,TG,TF,SNP,maf) {
    glm_model <- glm(reformulate(paste0("TF + SNP + TF:SNP + ",covs), response='TG'), data=cur_glm_data)
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
    # Get glm estimates
    x_coeff <- x_res[,c("cov","coef")]
    x_coeff_wide <- spread(x_coeff, key = cov, value = coef)
    # Add missing columns
    all_vars <- c("SNP", "TF", "TF_SNP",all_covs)
    covs.miss <- setdiff(all_vars, names(x_coeff_wide))
    x_coeff_wide[covs.miss] <- NA 
    x_coeff_wide <- x_coeff_wide[all_vars]
    names(x_coeff_wide) <- paste0(names(x_coeff_wide),".coef")

    # Get glm pvalues
    x_pval <- x_res[,c("cov","pval")]
    x_pval_wide <- spread(x_pval, key = cov, value = pval)
    # Add missing columns
    x_pval_wide[covs.miss] <- NA 
    x_pval_wide <- x_pval_wide[all_vars]
    names(x_pval_wide) <- paste0(names(x_pval_wide),".pval")

    # Standardize & GLM
    std_cur_glm_data <- standardize(cur_glm_data)
    std.glm <-  glm(reformulate(paste0("TF + SNP + TF:SNP + ",covs), response='TG'), data=std_cur_glm_data)
    x.std <- summary(std.glm)
    x_res_std <- x.std$coefficients
    x_res_std <- as.data.frame(x_res_std)
    x_res_std <- x_res_std[rownames(x_res_std) != "(Intercept)",c("Estimate","Pr(>|t|)")]
    names(x_res_std) <- c("std.coef","pval")
    row.names(x_res_std)[row.names(x_res_std) %in% c("TF:SNP")] <- "TF_SNP"
    x_res_std$cov <- rownames(x_res_std)
    row.names(x_res_std) <- 1:nrow(x_res_std)
    if(length(covs.miss)>0) {
        rows.miss <- data.frame(cov = covs.miss, std.coef = NA, pval = NA, 
                    stringsAsFactors = FALSE)
        x_res_std <- rbind(x_res_std, rows.miss)     
    }

    # Cbind and Return
    glm_all_res <- cbind(TG,TF,SNP,x_coeff_wide,x_pval_wide)
    glm_all_res$MAF <- maf
    glm_all_res$TF.coef.std <- x_res_std[which(x_res_std$cov == "TF"), ]$std.coef
    glm_all_res$SNP.coef.std <- x_res_std[which(x_res_std$cov == "SNP"), ]$std.coef
    glm_all_res$TF_SNP.coef.std <- x_res_std[which(x_res_std$cov == "TF_SNP"), ]$std.coef
    return(glm_all_res)
}


process.regeqtl.internal = function(expr.data, cov.data, trio.data, gt.data, out.dir) {

    # Expression data, genotypes and TG-TF-SNV trios
    rnaExpr <- expr.data
    gt <- gt.data
    germ.geno <- t(gt)
    trios_1 <- trio.data
    trios <-  trios_1[trios_1$gene != trios_1$TF, ]

    # Covariates
    covariates <- cov.data
    all_covs <- setdiff(colnames(covariates),"sample_id")
    all_peer_covs <- grep("Inferred",all_covs,value=TRUE)
    other_covs <- setdiff(all_covs,all_peer_covs)

    # Output files
    reg_oufile_1 <- file.path(out.dir,"regeQTL_All.RDS") 
    reg_oufile_2 <- file.path(out.dir,"regeQTL_Select.RDS") 

    # Run reg-eQTL for each trio
    reg_glm_res <- list()
    for(p in 1:nrow(trios)) {    
        print(p)
        cur.trio=trios[p,]
        TG=tg=as.character(cur.trio$gene)
        tf=TF=as.character(cur.trio$TF)
        SNP=snp=as.character(cur.trio$SNP)
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
            if(maf > 0.5) {
                maf <- 1-maf
            }
        }
        if(maf > 0.01 & maf < 0.5) {
            tg_ens <- grep(paste0("^",tg,"_"),rownames(rnaExpr),value=TRUE)
            tf_ens <- grep(paste0("^",tf,"_"),rownames(rnaExpr),value=TRUE)
            rnaExpr_sub <- rnaExpr[which(rownames(rnaExpr) %in% c(tg_ens,tf_ens)), ]
            rnaExpr_sub <- t(rnaExpr_sub)
            rnaExpr_sub <- as.data.frame(rnaExpr_sub)
            rnaExpr_sub$sample_id <- rownames(rnaExpr_sub)

            # Merge rna,genotype & clinical 
            glm_data_1 <- merge(rnaExpr_sub, covariates, by="sample_id")
            glm_data_2 <- merge(glm_data_1, germ.geno_sub, by="sample_id")

            # subset glm_data
            cur_glm_data <- glm_data_2[ ,c(tg_ens,tf_ens,snp,all_covs)]
            names(cur_glm_data)[names(cur_glm_data) == tg_ens] <- "TG"
            names(cur_glm_data)[names(cur_glm_data) == tf_ens] <- "TF"
            names(cur_glm_data)[names(cur_glm_data) == snp] <- "SNP"
            cur_glm_data <- na.omit(cur_glm_data)
            cur_glm_data$SNP <- as.numeric(cur_glm_data$SNP)

            if(nrow(cur_glm_data) < 150) {
                peer_covs <- all_peer_covs[1:15]
            } else if (nrow(cur_glm_data) >= 150 & nrow(cur_glm_data) < 250) {
                peer_covs <- all_peer_covs[1:30]
            } else if (nrow(cur_glm_data) >= 250 & nrow(cur_glm_data) < 350) {
                peer_covs <- all_peer_covs[1:45]
            } else {
                peer_covs <- all_peer_covs
            }
            covs <- c(other_covs, peer_covs)
            cur_glm_data <- cur_glm_data[ ,c("TG", "TF", "SNP", covs)]
            if(nrow(cur_glm_data) >= 20) {
                reg_glm_res[[p]] <- regGLM(cur_glm_data,covs, all_covs,TG,TF,SNP,maf)
            } 
        }
    }
    if(length(reg_glm_res) > 0) {
        reg_glm_final_res <- do.call(rbind,reg_glm_res)
        reg_glm_final_res_1 <- reg_glm_final_res[ ,c("TG", "TF", "SNP", 
                                "TF.coef", "TF.pval", "SNP.coef", "SNP.pval", 
                                "TF_SNP.coef", "TF_SNP.pval",
                                "TF.coef.std", "SNP.coef.std", "TF_SNP.coef.std", "MAF")]
        saveRDS(reg_glm_final_res, file = reg_oufile_1)
        saveRDS(reg_glm_final_res_1, file = reg_oufile_2)   
    }
}




#' Process s-eQTL Analysis
#'
#' This function runs s-eQTL analysis on expression data, genotype data, and TG-SNV pairs.
#' @param expr.data Expression data matrix.
#' @param cov.data Covariates data matrix.
#' @param pair.data TG-SNV pair information.
#' @param gt.data Genotype data matrix.
#' @param out.dir Directory to save the output.
#' @return The function outputs the s-eQTL results into the specified directory.
#' @export

process.seqtl = function(expr.data, cov.data, pair.data, gt.data, out.dir) {
    tryCatch(
        {
            process.seqtl.internal(expr.data, cov.data, pair.data, gt.data, out.dir)
        }, 
        error = function(err) {
            message(paste0("Error running seQTL: ", err$message))  
        },
        warning = function(warn) {
            message(paste0("Warning in seQTL: ", warn$message))
        },
        finally = {
            message("seQTL completed")
        }
    )
}


simGLM <- function(cur_glm_data,covs,all_covs,TG,SNP,maf) {
    glm_model <- glm(reformulate(paste0("SNP + ",covs), response='TG'), data=cur_glm_data)

    x <- summary(glm_model)
    vars <- rownames(x$coef)
    x_res <- x$coefficients
    x_res <- as.data.frame(x_res)
    x_res <- x_res[rownames(x_res) != "(Intercept)",c("Estimate","Pr(>|t|)")]
    names(x_res) <- c("coef","pval")
    x_res$cov <- rownames(x_res)
    row.names(x_res) <- 1:nrow(x_res)

    # Get glm estimates
    x_coeff <- x_res[,c("cov","coef")]
    x_coeff_wide <- spread(x_coeff, key = cov, value = coef)
    # Add missing columns
    all_vars <- c("SNP",all_covs)
    covs.miss <- setdiff(all_vars, names(x_coeff_wide))
    x_coeff_wide[covs.miss] <- NA 
    x_coeff_wide <- x_coeff_wide[all_vars]
    names(x_coeff_wide) <- paste0(names(x_coeff_wide),".coef")

    # Get glm pvalues
    x_pval <- x_res[,c("cov","pval")]
    x_pval_wide <- spread(x_pval, key = cov, value = pval)
    # Add missing columns
    x_pval_wide[covs.miss] <- NA 
    x_pval_wide <- x_pval_wide[all_vars]
    names(x_pval_wide) <- paste0(names(x_pval_wide),".pval")

    # Standardize & GLM
    std_cur_glm_data <- standardize(cur_glm_data)
    std.glm <-  glm(reformulate(paste0("SNP + ",covs), response='TG'), data=std_cur_glm_data)
    x.std <- summary(std.glm)
    x_res_std <- x.std$coefficients
    x_res_std <- as.data.frame(x_res_std)
    x_res_std <- x_res_std[rownames(x_res_std) != "(Intercept)",c("Estimate","Pr(>|t|)")]
    names(x_res_std) <- c("std.coef","pval")
    x_res_std$cov <- rownames(x_res_std)
    row.names(x_res_std) <- 1:nrow(x_res_std)
    if(length(covs.miss)>0) {
        rows.miss <- data.frame(cov = covs.miss, std.coef = NA, pval = NA, 
                    stringsAsFactors = FALSE)
        x_res_std <- rbind(x_res_std, rows.miss)     
    }

    # Cbind and Return
    glm_all_res <- cbind(TG,SNP,x_coeff_wide,x_pval_wide)
    glm_all_res$MAF <- maf
    glm_all_res$SNP.coef.std <- x_res_std[which(x_res_std$cov == "SNP"), ]$std.coef
    return(glm_all_res)
}


process.seqtl.internal = function(expr.data, cov.data, pair.data, gt.data, out.dir) {
    # Expression data, genotypes and TG-SNV pairs
    rnaExpr <- expr.data
    gt <- gt.data
    germ.geno <- t(gt)
    pairs <- pair.data

    # Covariates
    covariates <- cov.data
    all_covs <- setdiff(colnames(covariates),"sample_id")
    all_peer_covs <- grep("Inferred",all_covs,value=TRUE)
    other_covs <- setdiff(all_covs,all_peer_covs)

    # Output files
    s_oufile_1 <- file.path(out.dir,"seQTL_All.RDS") 
    s_oufile_2 <- file.path(out.dir,"seQTL_Select.RDS") 

    # Run seQTL for each pair
    sim_glm_res <- list()
    for(p in 1:nrow(pairs)) {    
        print(p)
        cur.pair=pairs[p,]
        TG=tg=as.character(cur.pair$gene)
        SNP=snp=as.character(cur.pair$SNP)
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
            if(maf > 0.5) {
                maf <- 1-maf
            }
        }
        if(maf > 0.01 & maf < 0.5) {
            tg_ens <- grep(paste0("^",tg,"_"),rownames(rnaExpr),value=TRUE)
            rnaExpr_sub <- rnaExpr[which(rownames(rnaExpr) %in% c(tg_ens)), ]
            rnaExpr_sub <- t(rnaExpr_sub)
            rnaExpr_sub <- as.data.frame(rnaExpr_sub)
            rnaExpr_sub$sample_id <- rownames(rnaExpr_sub)

            # Merge rna,genotype & clinical 
            glm_data_1 <- merge(rnaExpr_sub, covariates, by="sample_id")
            glm_data_2 <- merge(glm_data_1, germ.geno_sub, by="sample_id")

            # subset glm_data
            cur_glm_data <- glm_data_2[ ,c(tg_ens,snp,all_covs)]
            names(cur_glm_data)[names(cur_glm_data) == tg_ens] <- "TG"
            names(cur_glm_data)[names(cur_glm_data) == snp] <- "SNP"
            cur_glm_data <- na.omit(cur_glm_data)
            cur_glm_data$SNP <- as.numeric(cur_glm_data$SNP)

            if(nrow(cur_glm_data) < 150) {
                peer_covs <- all_peer_covs[1:15]
            } else if (nrow(cur_glm_data) >= 150 & nrow(cur_glm_data) < 250) {
                peer_covs <- all_peer_covs[1:30]
            } else if (nrow(cur_glm_data) >= 250 & nrow(cur_glm_data) < 350) {
                peer_covs <- all_peer_covs[1:45]
            } else {
                peer_covs <- all_peer_covs
            }
            covs <- c(other_covs, peer_covs)
            cur_glm_data <- cur_glm_data[ ,c("TG", "SNP", covs)]
            if(nrow(cur_glm_data) >= 20) {
                sim_glm_res[[p]] <- simGLM(cur_glm_data,covs, all_covs,TG,SNP,maf)
            } 
        }
    }
    if(length(sim_glm_res) > 0) {
        sim_glm_final_res <- do.call(rbind,sim_glm_res)
        sim_glm_final_res_1 <- sim_glm_final_res[ ,c("TG", "SNP", "SNP.coef", "SNP.pval", 
                                        "SNP.coef.std","MAF")]
        saveRDS(sim_glm_final_res, file = s_oufile_1)
        saveRDS(sim_glm_final_res_1, file = s_oufile_2)
    }
}



