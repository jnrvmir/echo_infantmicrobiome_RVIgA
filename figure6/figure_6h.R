################################################################################
###                              Figure 6B:                                  ###
####        Overall microbiome composition (beta-diversity) at M6           ####
####             by Rotavirus-IgA titers at M12                             ####



################### Load libraries  #########################
library(tidyverse); library(vegan); library(data.table);
library(MiRKAT); library(scales)

#####################   Beta Diversity ################################

###################  Data Set Up  ############################
in.dat <- read.csv("figure6ghi.csv") %>% as.data.table()

mapping<- in.dat %>% select(sid,pid,month,type,source,feeding_method,rec_pca_weeks,
                            bld_pca_weeks_12,c_diarrhea_12,
                            rvv_status_12,log_titer_12)


in.dat$sdepth <- rowSums(in.dat[,-(sid:log_titer_12)])

#remove samples with less than a 100 total reads

in.dat<- in.dat %>% dplyr::filter(sdepth>100)

in.dat <- in.dat %>% 
  dplyr::select(sid:log_titer_12,sdepth,everything())

##############################
# Build relative abundance matrix, removing taxa present in <2 samples
m12_tax_ra <- in.dat %>%
  select(-one_of(names(which(colSums(in.dat[,-(sid:sdepth)] > 0) < 2)))) %>%
  column_to_rownames(var = "sid") %>%
  select(-c(pid:sdepth)) %>%
  {t(apply(., 1, function(x) x/sum(x)))}

# Filter metadata to matched samples and extract components
tmp.meta <- mdat %>% filter(sample_id %in% rownames(m12_tax_ra))

cov.ma <- mapping %>% filter(sid %in% rownames(m12_tax_ra)) %>%
  select(sample_id,bld_pca_weeks_12,rec_pca_week) %>%
  column_to_rownames(var = "sid") %>%
  as.matrix()

iga.level <- mapping %>% 
  filter(sid %in% rownames(m12_tax_ra)) %>% pull(log_titer_12)

#create a distance matrix with bray curtis values
set.seed(164)

distm12.iga <- vegdist(m12_tax_ra, method = "bray", sample=100, iterations = 100) %>% as.dist()
## create distance kernel needed for MiRKAT
dist.kernel <- MiRKAT::D2K(as.matrix(distm12.iga))

beta.test <- MiRKAT::MiRKAT(y = iga.level,X=cov.ma,Ks = dist.kernel, nperm = 100000,
                            out_type = "C", omnibus = "permutation", 
                            returnKRV = TRUE, returnR2 = TRUE)

############## Beta diversity Figure ###############

PCoAs<- ape::pcoa(distm12.iga)
barplot(PCoAs$values$Relative_eig[1:10])
# Dimension (i.e., Axis 1 (PCOA1))
Axis1.percent <- PCoAs$values$Relative_eig[[1]] * 100
# Dimension (i.e., Axis 2 (PCOA2))
Axis2.percent <- PCoAs$values$Relative_eig[[2]] * 100
## convert PCoA result into a data frame
PCoAs.data <- data.frame(sid = rownames(PCoAs$vectors),
                         X = PCoAs$vectors[, 1],
                         Y = PCoAs$vectors[, 2])
PCoAs.met.dat <-left_join(PCoAs.data, mapping, by = "sid")

R2<- beta.test$R2
pval <- beta.test$p_values

PCoAs.met.dat %>%
  ggplot(., aes(x = X, y = Y)) +
  geom_point(aes(color = log_titer_12), size=5)  +
  xlab(paste("PCoA1 - ", gsub("\\.", "\u00B7", round(Axis1.percent, 2)), "%", sep = "")) +
  ylab(paste("PCoA2 - ", gsub("\\.", "\u00B7", round(Axis2.percent, 2)), "%", sep = "")) +
  scale_x_continuous(labels = scales::label_number(decimal.mark = "\u00B7")) +
  scale_y_continuous(labels = scales::label_number(decimal.mark = "\u00B7")) +
  theme_bw() +
  theme(text = element_text(size = 16),
        legend.title = element_text(margin = margin(b = 20)),
        legend.key.height = unit(2.5, "cm"), size=14) +
  labs(color = expression(atop(Log[2] ~ "Rotavirus-IgA", "Titer (M6)")))+
  scale_color_gradientn(
    colors = paletteer::paletteer_c("grDevices::Dark Mint", 30),
    labels = scales::label_number(decimal.mark = "\u00B7"),  # replaces "." with "·"
    guide = guide_colorbar(
      barheight = unit(10, "cm"),
      barwidth = unit(1, "cm"),
      title.position = "top"))+ 
  annotate("label", x = 0.0, y = 0.5, 
           label = paste0("R2: ",gsub("\\.", "\u00B7", round(R2, 4)), "\np: ", gsub("\\.", "\u00B7", round(pval,3))), 
           hjust =0 , vjust = 1, size=5)