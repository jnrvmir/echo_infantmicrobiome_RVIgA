################################################################################
###                              Figure 6G:                                  ###
####   Association between alpha diversity at M6 and Rotavirus-IgA at M12   ####


######## Load libraries ######## 
library(tidyverse); library(vegan); library(data.table);library(emmeans)


##############################################################################
######################### Function to Reorganize ASV
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

##############################################################################
########################   Data Set Up  ###########################
alpha.dat <- read.csv("figure6ghi.csv") %>% as.data.table()

#Extract metadata
mapping<- alpha.dat %>% select(sid,pid,month,type,source,feeding_method,rec_pca_weeks,
                               bld_pca_weeks_12,c_diarrhea_12,
                               rvv_status_12,log_titer_12)

alpha.dat$sdepth <- rowSums(alpha.dat[,-(sid:log_titer_12)])

#Remove samples with less than 100 reads total. 

alpha.dat<- alpha.dat %>% dplyr::filter(sdepth>100)

alpha.dat <- alpha.dat %>% 
  dplyr::select(sid:log_titer_12,sdepth,everything())

alpha.dat <- alpha.dat %>% 
  dplyr::select(sid:log_titer_12,sdepth,
                names(which(apply(dplyr::select(., -(sid:log_titer_12), -sdepth), 
                                  2, function(x) sum(x) > 0)))) %>% as.data.table()

colnames(alpha.dat) <- str_replace_all(colnames(alpha.dat), "Peptostreptococcales.Tissierellales", 
                                       "Peptostreptococcales_Tissierellales") %>%
  str_replace_all(., "Veillonellales.Selenomonadales", "Veillonellales_Selenomonadales") %>%
  str_replace_all(., "WPS.2", "WPS2") %>%
  str_replace_all(., "Escherichia.Shigella", "Escherichia_Shigella") %>%
  str_replace_all(., ".Ruminococcus._", "Ruminococcus_") %>%
  str_replace_all(., ".Clostridium._", "Clostridium_") %>%
  str_replace_all(., ".Eubacterium._", "Eubacterium_") %>%
  str_replace_all(., "S5.A14a", "S5A14a") %>%
  str_replace_all(., "UCG.014", "UCG014") %>%
  str_replace_all(., "UCG.004", "UCG004") %>%
  str_replace_all(., "R.7", "R7") %>%
  str_replace_all(., "Allorhizobium.Neorhizobium.Pararhizobium.Rhizobium",
                  "Allorhizobium_Neorhizobium_Pararhizobium_Rhizobium") %>%
  str_replace_all(., "Methylobacterium.Methylorubrum", "Methylobacterium_Methylorubrum") %>%
  str_replace_all(., "Hafnia.Obesumbacterium", "Hafnia_Obesumbacterium") %>%
  str_replace_all(., "g_Eubacterium__", "g_Eubacterium.__") %>%
  str_replace_all(., "g__Erysipelotrichaceae_UCG.003.__", "g__Erysipelotrichaceae_UCG003.__") %>%
  str_replace_all(., "g__UCG.003.__", "g__UCG003.__")

### Create taxonomic rank table

txnames <- colnames(subset(alpha.dat, select = -(sid:sdepth)))
KINGDOM=1; PHYLUM=2; CLASS=3; ORDER=4; FAMILY=5; GENUS = 6; SPECIES=7
taxonomy <- str_split_fixed(txnames, "\\.", 7)

# Aggregate the taxa using the reorganize.asv.profile function.
# For this analysis we will do aggregation at 20% Phylum level. 

tmp.alpha <- reorganize.asv.profile(as.matrix(alpha.dat[,-(sid:sdepth)]), taxonomy, 0.2, PHYLUM)
rownames(tmp.alpha) <- alpha.dat$sid

# For each taxa determine in what proportion of samples the taxa is present. Then,
# we remove the taxa not present in more than 20% of samples. Lastly, we remove
# unclassified taxa. 

# Filter, clean, and normalize. We remove unclassified/uncultured taxa
tmp.p.alpha <- tmp.alpha[, apply(tmp.alpha, 2, function(x) sum(x>0)/length(x)) > 0.2] %>%
  as.data.frame() %>%
  dplyr::select(-s__uncultured_bacterium, -s__uncultured_organism) %>%
  {sweep(., 1, rowSums(.), "/")}

# Build final analysis table directly

meta.sh <- alpha.dat[, sid:sdepth] %>%
  as.data.frame() %>%
  mutate(shannon = vegan::diversity(tmp.p.alpha, index = "shannon")) %>% #calculate Shannon diversity
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.factor), ~ fct_relevel(., levels(.)[which.max(table(.))]))) %>%
  as.data.table()

rownames(meta.sh) <- meta.sh$sid

### Linear regression model adjusting for feeding method and postconception age, and history of gastroenteritis
model <- lm(log_titer_12~shannon+feeding_method+rec_pca_weeks+bld_pca_weeks_12+c_diarrhea_12, data=meta.sh) 
confint(model, level = 0.95)

B1 <- summary(model)$coefficients[2,1]
pval <- round(summary(model)$coefficients[2,4],3)

m6m12_shan <- meta.sh%>%
  ggplot(., aes(x=shannon, y=log_titer_12))+
  geom_point(size=3)+
  geom_smooth(method="lm", se=FALSE)+
  theme_bw() +
  scale_y_continuous(limits=c(-1,13))+
  scale_x_continuous(labels = scales::label_number(decimal.mark = "\u00B7")) +
  labs(x="Shannon Diversity (M6)", y=expression(""*Log[2]*" Rotavirus-IgA Titer (M12)"))+
  annotate("label", x = 2.5, y = 12, 
           label = paste0("ß: ",gsub("\\.", "\u00B7", round(B1, 2)), "\np: ", gsub("\\.", "\u00B7", pval)), 
           hjust =0 , vjust = 1, size=5)+
  theme(text=element_text(size=16))
