################################################################################
#################               Figure 4                       ################# 
################################################################################

######Load in Required Packages #######

library(data.table); library(tidyverse); library(glue); 
library(ggpubr); library(lme4);
library(emmeans)


######### Load Data ###########
full.dat <- read.csv("data_figure4.csv")

################################################################################
########           Figure 4A: Rotavirus-IgA titers overtime           ##########

### Analysis include a linear mixed effects model, adjusted for gestational age and repeated measured.

### Lienar mixed effects model
clean_model <- lmer(log.titer ~ as.factor(month) + ga_best + (1 | participant_id), 
                    data = full.dat)
### Estimated marginal means with adjustment for multiple comparisons
emm_res <- emmeans(clean_model, pairwise ~ month, adjust="BH")
p_data  <- as.data.frame(emm_res$contrasts)

raw_p1 <- p_data$p.value[p_data$contrast == "Birth - Month 6"]
label_0_6 <- if(raw_p1 < 0.0001) "p < 0\u00B70001" else paste0("p <", gsub("\\.", "\u00B7", format(round(raw_p1, 4), nsmall=4)))

raw_p2 <- p_data$p.value[p_data$contrast == "Month 6 - Month 12"]
label_6_12 <- if(raw_p2 < 0.0001) "p < 0\u00B70001" else paste0("p < ", gsub("\\.", "\u00B7", format(round(raw_p2, 4), nsmall=4)))



visit.cols = c("Birth"= "#999999",
               "Month 6"= "#56B4E9",
               "Month 12"= "darkblue")

### Figure plot 
full.dat %>%
  group_by(month) %>%
  summarize(
    mean_titer = mean(log.titer, na.rm = TRUE), 
    sd_titer = sd(log.titer, na.rm = TRUE), 
    .groups = 'drop') %>%
  ggplot(aes(x = month, y = mean_titer, colour = month, group=1)) +
  geom_point(size = 7) +
  geom_line(size = 0.95) +
  geom_errorbar(aes(ymin = mean_titer - sd_titer, ymax = mean_titer + sd_titer),
                width = 0.2) +
  scale_y_continuous(limits = c(-1, 11),
                     labels = scales::label_number(decimal.mark = "\u00B7")) +
  theme_bw() +
  scale_color_manual(values = visit.cols) +
  labs(y = expression(Log[2] * "Rotavirus-IgA Titer")) +
  theme(
    text = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.title.x = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 18, face = "bold"))+
  annotate("segment", x = 1, xend = 1.9, y = 9.5, yend = 9.5, colour = "black") +
  annotate("text", x = 1.5, y = 10, label = label_0_6, size = 6, colour = "black") +
  annotate("segment", x = 2.1, xend = 3, y = 9.5, yend = 9.5, colour = "black") +
  annotate("text", x = 2.5, y = 10, label = label_6_12, size = 6, colour = "black")+
  scale_x_discrete(labels = c("Month 6" = "M6", "Month 12" = "M12")) 


################################################################################
########  Figure 4b: Rotavirus IgA Serostatus Propotions Over Time   ##########

full.dat %>%
  group_by(month, seropositivity) %>%
  summarise(n = n(), .groups = "drop_last") %>% 
  mutate(N_total = sum(n), 
         prop = n / N_total,
         pct_text = gsub("\\.", "\u00B7", round(prop * 100, 1)),
         bar_label = paste0(n, " (", pct_text, "%)"),
         axis_label = paste0(month, "\n(n=", N_total, ")")) %>%
  ungroup() %>%
  ggplot(., aes(x = month, y = prop, fill = seropositivity)) +
  geom_col(position = "fill", width = 0.7,color = "black") +
  geom_text(aes(label = bar_label), 
            position = position_stack(vjust = 0.5), 
            color = "black", fontface="bold",size = 4) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("seronegative" = "#DCD6B2FF", "seropositive" = "#80944EFF")) +
  labs(x = NULL, y = "Proportion", fill = "Rotavirus-IgA\nSerostatus") +
  theme_bw() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size=10),
    legend.title = element_text(size=10))+
  scale_x_discrete(labels = c("Month 6" = "M6", "Month 12" = "M12")) 


################################################################################
########  Figure 4c: Rotavirus IgA by Vaccine Schedule Status   ##########

m6mat<- full.dat %>% filter(month=="Month 6")

fit <- lm(log.titer~schedule_status+bld_pca_weeks, data=m6mat)
p_val <- summary(fit)$coefficients[[2,4]]
p_label_text <- paste0("p = ", sprintf("%.3f", p_val))
p_label_custom <- gsub("\\.", "\u00B7", p_label_text)
y_bracket <- 11 + 0.2
y_text <- y_bracket + 0.5

m6_schedu<- m6mat %>%
  filter(month=="Month 6") %>%
  mutate(schedule_status = factor(
    schedule_status,
    levels = c("Full schedule", "Partial schedule"), 
    labels = c("Full\nschedule", "Partial\nSchedule")),
    month= factor(month, levels="Month 6", labels="M6")) %>%
  ggplot(aes(x = schedule_status, y = log.titer, fill = schedule_status)) +
  geom_boxplot(outliers=F) +
  geom_jitter(width = 0.2) +
  theme_bw() +
  scale_fill_manual(values = c("Full\nschedule" = "#CABEE9FF", 
                               "Partial\nSchedule" = "#FAE093FF")) +
  labs(x = "Vaccine Schedule Status",
       y = expression(Log[2] * " Rotavirus-IgA Titer at M6"))+
  annotate("segment", x = 1, xend = 2, y = y_bracket, yend = y_bracket) + # Horizontal
  annotate("segment", x = 1, xend = 1, y = y_bracket, yend = y_bracket - 0.05) + # Left Tick
  annotate("segment", x = 2, xend = 2, y = y_bracket, yend = y_bracket - 0.05) + # Right Tick
  annotate("text", x = 1.5, y = y_text, label = p_label_custom, size = 5) +
  theme(
    text = element_text(size = 16),
    axis.title.y = element_text(size=14),
    strip.background = element_rect(fill = "#EEE9E9"),
    strip.text = element_text(face = "bold", size = 16),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none")

m12mat<- full.dat %>% filter(month=="Month 12")

fit <- lm(log.titer~schedule_status+bld_pca_weeks+c_diarrhea, data=m12mat)
p_val <- summary(fit)$coefficients[[2,4]]
p_label_text <- paste0("p = ", sprintf("%.3f", p_val))
p_label_custom <- gsub("\\.", "\u00B7", p_label_text)
y_bracket <- 9 + 0.2
y_text <- y_bracket + 0.5

m12_schedu<- m12mat %>%
  filter(month=="Month 12") %>%
  mutate(schedule_status = factor(
    schedule_status,
    levels = c("Full schedule", "Partial schedule"), 
    labels = c("Full\nschedule", "Partial\nSchedule")),
    month= factor(month, levels="Month 12", labels="M12")) %>%
  ggplot(aes(x = schedule_status, y = log.titer, fill = schedule_status)) +
  geom_boxplot(outliers=F) +
  geom_jitter(width = 0.2) +
  theme_bw() +
  scale_fill_manual(values = c("Full\nschedule" = "#CABEE9FF", 
                               "Partial\nSchedule" = "#FAE093FF")) +
  labs(x = "Vaccine Schedule Status",
       y = expression(Log[2] * " Rotavirus-IgA Titer at M12"))+
  annotate("segment", x = 1, xend = 2, y = y_bracket, yend = y_bracket) + # Horizontal
  annotate("segment", x = 1, xend = 1, y = y_bracket, yend = y_bracket - 0.05) + # Left Tick
  annotate("segment", x = 2, xend = 2, y = y_bracket, yend = y_bracket - 0.05) + # Right Tick
  annotate("text", x = 1.5, y = y_text, label = p_label_custom, size = 5) +
  theme(
    text = element_text(size = 16),
    axis.title.y = element_text(size=14),
    strip.background = element_rect(fill = "#EEE9E9"),
    strip.text = element_text(face = "bold", size = 16),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none")

################################################################################
########  Figure 4d: Proportions of vaccine seroconversion   ##########

 full.dat %>%
  filter(!month=="Birth") %>%
  drop_na(rv5_seroconverter) %>%
  group_by(month, rv5_seroconverter) %>%
  mutate(rv5_seroconverter= factor(rv5_seroconverter, levels = c("yes", "no"), labels=c("Seroconversion", "No Seroconversion"))) %>%
  summarise(n = n(), .groups = "drop_last") %>% 
  mutate(N_total = sum(n), 
         prop = n / N_total,
         pct_text = gsub("\\.", "\u00B7", round(prop * 100, 1)),
         bar_label = paste0(n, " (", pct_text, "%)"),
         axis_label = paste0(month, "\n(n=", N_total, ")")) %>%
  ungroup() %>%
  ggplot(., aes(x = month, y = prop, fill = rv5_seroconverter)) +
  geom_col(position = "fill", width = 0.7, color="black") +
  geom_text(aes(label = bar_label), 
            position = position_stack(vjust = 0.5), 
            color = "black", fontface="bold",size = 4) +
  scale_y_continuous(labels = scales::percent) +
  ggsci::scale_fill_lancet()+
  scale_fill_manual(values = c("Seroconversion" = "#DA95AAFF", "No Seroconversion" = "#E5D8BD")) +
  labs(x = NULL, y = "Proportion", fill = "Rotavirus-IgA\nSeroconversion") +
  theme_bw() +
  theme(text = element_text(size = 16),
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size=10),
        legend.title = element_text(size=10))+
  scale_x_discrete(labels = c("Month 6" = "M6", "Month 12" = "M12")) 

##################################################

### Month 6 visit 
m6mat<- full.dat %>% filter(month=="Month 6")

fit <- lm(log.titer~schedule_status+bld_pca_weeks, data=m6mat)
p_val <- summary(fit)$coefficients[[2,4]]
p_label_text <- paste0("p = ", sprintf("%.3f", p_val))
p_label_custom <- gsub("\\.", "\u00B7", p_label_text)
y_bracket <- 11 + 0.2
y_text <- y_bracket + 0.5

 m6mat %>%
  filter(month=="Month 6") %>%
  mutate(schedule_status = factor(
    schedule_status,
    levels = c("Full schedule", "Partial schedule"), 
    labels = c("Full\nschedule", "Partial\nSchedule")),
    month= factor(month, levels="Month 6", labels="M6")) %>%
  ggplot(aes(x = schedule_status, y = log.titer, fill = schedule_status)) +
  geom_boxplot(outliers=F) +
  geom_jitter(width = 0.2) +
  theme_bw() +
  scale_fill_manual(values = c("Full\nschedule" = "#CABEE9FF", 
                               "Partial\nSchedule" = "#FAE093FF")) +
  labs(x = "Vaccine Schedule Status",
       y = expression(Log[2] * " Rotavirus-IgA Titer at M6"))+
  annotate("segment", x = 1, xend = 2, y = y_bracket, yend = y_bracket) + # Horizontal
  annotate("segment", x = 1, xend = 1, y = y_bracket, yend = y_bracket - 0.05) + # Left Tick
  annotate("segment", x = 2, xend = 2, y = y_bracket, yend = y_bracket - 0.05) + # Right Tick
  annotate("text", x = 1.5, y = y_text, label = p_label_custom, size = 5) +
  theme(
    text = element_text(size = 16),
    axis.title.y = element_text(size=14),
    strip.background = element_rect(fill = "#EEE9E9"),
    strip.text = element_text(face = "bold", size = 16),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none")

### Month 12 visit
 
m12mat<- full.dat %>% filter(month=="Month 12")

fit <- lm(log.titer~schedule_status+bld_pca_weeks+c_diarrhea, data=m12mat)
p_val <- summary(fit)$coefficients[[2,4]]
p_label_text <- paste0("p = ", sprintf("%.3f", p_val))
p_label_custom <- gsub("\\.", "\u00B7", p_label_text)
y_bracket <- 9 + 0.2
y_text <- y_bracket + 0.5

m12mat %>%
  filter(month=="Month 12") %>%
  mutate(schedule_status = factor(
    schedule_status,
    levels = c("Full schedule", "Partial schedule"), 
    labels = c("Full\nschedule", "Partial\nSchedule")),
    month= factor(month, levels="Month 12", labels="M12")) %>%
  ggplot(aes(x = schedule_status, y = log.titer, fill = schedule_status)) +
  geom_boxplot(outliers=F) +
  geom_jitter(width = 0.2) +
  theme_bw() +
  scale_fill_manual(values = c("Full\nschedule" = "#CABEE9FF", 
                               "Partial\nSchedule" = "#FAE093FF")) +
  labs(x = "Vaccine Schedule Status",
       y = expression(Log[2] * " Rotavirus-IgA Titer at M12"))+
  annotate("segment", x = 1, xend = 2, y = y_bracket, yend = y_bracket) + # Horizontal
  annotate("segment", x = 1, xend = 1, y = y_bracket, yend = y_bracket - 0.05) + # Left Tick
  annotate("segment", x = 2, xend = 2, y = y_bracket, yend = y_bracket - 0.05) + # Right Tick
  annotate("text", x = 1.5, y = y_text, label = p_label_custom, size = 5) +
  theme(
    text = element_text(size = 16),
    axis.title.y = element_text(size=14),
    strip.background = element_rect(fill = "#EEE9E9"),
    strip.text = element_text(face = "bold", size = 16),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none")
