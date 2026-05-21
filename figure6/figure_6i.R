################################################################################
###                              Figure 6C:                                  ###
####             Individual taxa at M6 association with                     ####
####                    Rotavirus-IgA titers at M12                         ####

### Load libraries needed ##
library(tidyverse); library(data.table)

###################### Functions needed   ################################

###################### Function to Reorganize ASV ###############
reorganize.asv.profile <- function(asv.table, taxonomy.table,
                                   cutoff.proportion, highest.rank=CLASS){
  tmp.asv.table <- asv.table
  tax.tab <- taxonomy.table
  cutoff <- nrow(tmp.asv.table)*cutoff.proportion
  fdat <- NULL
  for(j in 7:highest.rank){
    tmp.dat <- apply(tmp.asv.table, 1,
                     function(x,y) tapply(x,y,sum), y=tax.tab[,j])
    if(j==highest.rank){
      fdat <- cbind(fdat, t(tmp.dat))
    } else{
      sel.taxa <- names(which(apply(tmp.dat, 1, function(x) sum(x>0)) > cutoff))
      sel.taxa <- sel.taxa[!(sel.taxa %in% "__")]
      fdat <- cbind(fdat, t(tmp.dat[rownames(tmp.dat) %in% sel.taxa,,drop=FALSE]))
      rm.taxa <- which(tax.tab[,j] %in% sel.taxa)
      tax.tab <- tax.tab[-rm.taxa,]
      tmp.asv.table <- tmp.asv.table[,-rm.taxa]
    }
  }
  sel.taxa <- names(which(apply(fdat, 2, function(x) sum(x>0)) > (cutoff*0.01/cutoff.proportion)))
  sel.taxa <- sel.taxa[!(sel.taxa %in% "__")]
  return(fdat[,sel.taxa])
}


################### Log Contrast Modeling  ######################

################# Functions for the Log Contrast ###############

### Estimate the initial value of lambda for scaled lasso
get.init.lambda <- function(n.sample, n.feature, tol){
  f <- function(x, p) {(qnorm(1-x/p))^4 + 2*((qnorm(1-x/p))^2) - x}
  k <- uniroot(f, lower=0, upper=n.feature-1, tol=tol, p=n.feature)$root
  lambda_0 <- sqrt(2/n.sample)*qnorm(1-k/n.feature)
  return(lambda_0)
}

### Lasso with a linear constraint
lasso_constr <- function(y, x, contr2, lam, tol, max.iter){
  n <- nrow(x); p <- ncol(x); k <- dim(contr2)[1]
  gramC <- crossprod(contr2); gramX <- crossprod(x)
  diagC <- diag(gramC); diagX <- diag(gramX)
  dediagC <- gramC - diag(diagC); dediagX <- gramX - diag(diagX)
  covXY <- crossprod(x, y)
  
  mu <- 1; bet <- rep(1,p)/p; bet0 <- rep(0,p); iter <- 0
  if (sum(abs(contr2))==0){
    term0 <- (covXY-dediagX%*%bet)/n
    term2 <- diagX/n
    while (sum(abs(bet-bet0))>tol & iter<max.iter){
      bet0 <- bet
      for(j in 1:p){
        term1 <- sign(term0[j])*max(0, abs(term0[j])-lam)
        bet[j] <- term1/term2[j]
        dif <- bet[j] - bet0[j]
        term0 <- term0 - dediagX[,j]*dif/n
      }
      iter <- iter+1
    }
  } else{
    ksi <- rep(0,k); ksi0 <- rep(1,k)
    term0 <- (covXY-dediagX%*%bet)/n - mu*(t(contr2)%*%ksi + dediagC%*%bet)
    term2 <- diagX/n + mu*diagC
    while (sum(abs(ksi-ksi0))>tol && iter<max.iter){
      ksi0 <- ksi; iter2 <- 0; bet0 <- bet0 + 1
      while (sum(abs(bet-bet0))>tol && iter2<1000){
        bet0 <- bet
        for(j in 1:p){
          term1 <- sign(term0[j])*max(0, abs(term0[j])-lam)
          bet[j] <- term1/term2[j]
          dif <- bet[j] - bet0[j]
          term0 <- term0 - dediagX[,j]*dif/n - dediagC[,j]*dif*mu
        }
        iter2 <- iter2 + 1
      }
      dif2 <- contr2%*%bet
      ksi <- ksi + dif2
      term0 <- term0 - mu*t(contr2)%*%dif2
      iter <- iter + 1
    }
  }
  return(bet)
}

### Scaled lasso for compositional data
cmm_slr <- function(y, Q, contr, tol, max.iter){
  n <- nrow(Q); p <- ncol(Q)
  mn.y <- mean(y); mn.Q <- colMeans(Q)
  cent.Q <- (Q-tcrossprod(rep(1, n), mn.Q))
  contr2 <- t(as.matrix(contr))
  cent.y <- y-mn.y
  lam0 <- get.init.lambda(n, p, tol)
  sigma <- 1; sigma_2 <- 1; sigma_s <- 0.5; iter <- 1
  while (abs(sigma-sigma_s)>0.01 & iter<100){
    iter <- iter + 1
    sigma <- (sigma_s + sigma_2)/2
    lam <- sigma*lam0
    bet2 <- lasso_constr(cent.y, cent.Q, contr2, lam, tol, max.iter)
    s <- sum(abs(bet2)>0.001)
    s <- min(s, n-1)
    sigma_s <- base::norm(cent.y-cent.Q%*%bet2, type="2")/sqrt(n-s-1)
    sigma_2 <- sigma
  }
  if(iter==100) print("Not converge!")
  sigma <- sigma_s
  bet <- bet2
  intercp <- mn.y - mn.Q%*%bet
  return(list(beta=bet, intercept=intercp, sigma=sigma, lambda0=lam0))
}

slcmc <- function(y, M, x, tol=1e-6, max.iter=1000, padj.method="BH"){
  if(!is.vector(y)) y <- as.vector(y)
  n <- length(y); k <- ncol(M)
  if(is.null(x)){
    Z <- log(M)
    n.x <- 0
  } else{
    Z <- as.matrix(cbind(log(M), x)) 
    if(!is.matrix(x)) x <- as.matrix(x)
    n.x <- ncol(x)
  }
  n.vrs <- k + n.x
  contr <- c(rep(1/sqrt(k), k), rep(0, n.x))
  Z.til <- Z %*% (diag(n.vrs)-tcrossprod(contr))
  est.param <- cmm_slr(y, Z.til, contr, tol, max.iter)
  cent.y <- scale(y, scale=FALSE)
  cent.Z.til <- scale(Z.til, scale=FALSE)
  gam0 <- est.param$lambda0/3; Sig <- crossprod(cent.Z.til)/n
  Sig2 <- Sig - diag(diag(Sig))
  Q <- diag(n.vrs) - tcrossprod(contr)
  M.til <- matrix(0, n.vrs, n.vrs)
  for(i in 1:n.vrs){
    gam <- gam0/2
    while(gam<0.5){
      gam <- gam*2
      mi <- rep(1,n.vrs)
      mi0 <- rep(0,n.vrs)
      iter <- 1
      while(sum(abs(mi-mi0))>tol & iter<max.iter){
        mi0 <- mi
        for(j in 1:n.vrs){
          v <- -Sig2[j,]%*%mi+Q[j,i]
          mi[j] <- sign(v)*max(0, abs(v)-gam)/Sig[j,j]
        }
        iter <- iter + 1
      }
      if(iter<max.iter) break
    }
    M.til[i,] <- mi
  }
  M.til <- Q%*%M.til
  debias.B <- est.param$beta + M.til%*%t(cent.Z.til)%*%(cent.y-cent.Z.til%*%est.param$beta-rep(est.param$intercept, n))/n
  cov.debias.B <- est.param$sigma^2*M.til%*%Sig%*%t(M.til)/n
  pval <- 2*pnorm(abs(debias.B)/sqrt(diag(cov.debias.B)), lower.tail=FALSE)
  adjp <- c(p.adjust(pval[1:k], method=padj.method), pval[(k+1):n.vrs])
  rslt <- data.frame(beta=debias.B, se=sqrt(diag(cov.debias.B)), p=pval, adjp=adjp)
  rownames(rslt) <- c(colnames(M), colnames(x))
  return(rslt)
}

###########################################
in.dat <- read.csv("figure6ghi.csv") %>% as.data.table()

mapping<- in.dat %>% select(sid,pid,month,type,source,feeding_method,rec_pca_weeks,
                            bld_pca_weeks_12,c_diarrhea_12,
                            rvv_status_12,log_titer_12)

indv.dat <- in.dat %>%
  mutate(sdepth = rowSums(select(., -c(sid:log_titer_12)))) %>%
  filter(sdepth > 100) %>%
  select(sid:log_titer_12, sdepth, everything()) %>%
  dplyr::select(sid:log.titer.m6,sdepth,
                names(which(apply(dplyr::select(., -(sid:log_titer_12), -sdepth), 
                                  2, function(x) sum(x) > 0)))) %>% as.data.table()

tmp.indv <- tmp.indv %>%
  mutate(other_taxa = rowSums(select(., all_of(taxa_cols)))) %>%
  select(-all_of(taxa_cols),-s__uncultured_bacterium, -s__uncultured_organism,-s__human_gut)

### We assign the pseudocount value by determining the lowest #####

psudocount<- min(tmp.indv[tmp.indv > 0])*0.1

## Pseudocount added to 0 values. 
tmp.indv[tmp.indv == 0] <- tmp.indv[tmp.indv == 0] + psudocount


tmp.p.indv <- sweep(tmp.indv, 1, rowSums(tmp.indv), "/") %>%
  bind_cols(in.dat[, sid:sdepth] %>% as.data.frame()) %>%
  filter(!is.na(type)) %>%
  column_to_rownames(var = "sid") %>%
  select(pid:sdepth, everything())

# Analysis objects
iga.level<- as.numeric(tmp.p.indv$log_titer_12)

matrix<- tmp.p.indv %>% select(-c(pid:sdepth)) %>% as.matrix()

# Set up metadata needed for log contrast. In this cse we adjust for age and history of gastroenteritis.
# Categorical variables like history of gastroenteritis has to be converted to binary.

tmp.meta<- tmp.p.indv %>%
  select(pid:sdepth) %>%
  mutate(c_diarrhea_num = case_when(
    c_diarrhea_12 == "TRUE" ~ 1,
    c_diarrhea_12 == "FALSE" ~ 0,
    TRUE ~ NA_real_))%>%
  select(bld_pca_weeks_12,rec_pca_weeks,c_diarrhea_num) %>%
  column_to_rownames(var = "sid") %>%
  as.matrix()

# Run log contrast test 
test <- slcmc(y=iga.level, M= as.matrix(matrix), x=tmp.meta,tol=1e-6, max.iter=50, padj.method="BH")


###########   Figure Forest Plot     ##############

#### Set up data for plot

test <- test %>% 
  rownames_to_column(var="taxa") %>%
  select(taxa,everything()) %>%
  mutate( CI_lower = beta - 1.96 * se, 
          CI_upper = beta + 1.96 * se) %>%
  mutate(CI = paste0("(", round(CI_lower, 3), ", ", round(CI_upper, 3), ")")) %>%
  select(-c(CI_lower,CI_upper)) %>%
  mutate(across(where(is.numeric), ~ round(., 5)))

## Forest plot
test %>%
  filter(adjp < 0.2) %>%
  filter(!taxa=="c_diarrhea_num_12") %>%
  select(taxa, beta, se, CI, p, adjp) %>%
  mutate(
    taxa = gsub("g__", "", taxa),
    taxa = gsub("s__", "", taxa),
    taxa = gsub("_", " ", taxa),
    taxa = factor(taxa, levels=c("Streptococcus","Anaerococcus","Peptococcus","Peptoniphilus",
                                 "Prevotella buccalis","Campylobacter ureolyticus",
                                 "Murdochiella asaccharolytica","Dialister propionicifaciens"),
                  labels =c("Streptococcus","Anaerococcus","Peptococcus","Peptoniphilus",
                            "Prevotella\nbuccalis","Campylobacter\nureolyticus",
                            "Murdochiella\nasaccharolytica","Dialister\npropionicifaciens")))%>%
  separate(CI, into = c("ci_lower", "ci_upper"), sep = ", ", remove = FALSE) %>%
  mutate(
    ci_lower = as.numeric(gsub("[\\(]", "", ci_lower)),
    ci_upper = as.numeric(gsub("[\\)]", "", ci_upper)),
    val_string = as.character(signif(adjp, 2)),
    label_base = ifelse(adjp < 0.0001, 
                        "q < 0.0001", 
                        paste0("q = ", val_string)),
    label_final = gsub("\\.", "\u00B7", label_base)) %>%
  arrange(beta, taxa) %>%
  ggplot(aes(x = beta, y = reorder(taxa, beta))) +
  geom_point(size = 4, color = "black") +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, size = 1, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_text(aes(label = label_final),
            nudge_y = 0.3,
            size = 5,
            color = "black") +
  theme_bw() +
  scale_x_continuous(limits = c(-0.5, 0.5),labels = scales::label_number(decimal.mark = "\u00B7")) +
  labs(x = "Beta Estimate (Effect Size)", y = NULL) +
  theme(legend.position = "none", 
        text = element_text(size = 16), 
        axis.title = element_text(size = 16),
        axis.text.y = element_text(face = "italic"),
        axis.title.x = element_text(margin = margin(t = 10), size = 16))
