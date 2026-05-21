################################################################################
#################               Figure 5                       ################# 
################################################################################

########### Load libraries ###########


library(tidyverse);library(data.table);library(gt);library(rstatix);
library(emmeans);library(sjPlot);library(markdown)

########### Load Data ###########

dat <- read.csv("data_figure5.csv")

################################################################################
################# Figure 5A: Month 6 Visit Forest Plot
################################################################################

######## Month 6 Data setup ###########

m6.dat<- dat %>% filter(monht==6) %>% 
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.factor), ~ fct_relevel(., levels(.)[which.max(table(.))]))) #This line is to order the categorical variables so the level with the most observation is first and used as reference in linear models
  

####### Month 6 Linear Models ########

# Variables to test
vars <- c("bld_pca_weeks", "sam_season_2", "infant_sex","z_score", "pre.bmi", 
          "parity_cat", "use_childcare","c_diarrhea","medicaid_in_preg","infant_ethrace_2",
          "mother_ethrace_2","feeding_method_m1", "feeding_method_m6",
          "exposed_bm_m1","exposed_bm_m6","bfdur_rec_m1_pct","bfdur_bld_m6_pct",
          "schedule_status","weeks_diff_pl")

adjusted_list <- lapply(formulas, function(fmla) {
  # Fit the model
  mod <- lm(fmla, data = m6.dat)
  # Get the independent variable (first variable in the formula)
  independent_var <- all.vars(fmla)[2]
  tidy(mod, conf.int = TRUE) %>%
    mutate(model = "adjusted", independent_variable = independent_var)  # Add independent variable info
})

# Combine results into a single data table
adjusted_dt <- rbindlist(adjusted_list, fill = TRUE)

filter_by_independent_variable <- function(df) {
  df %>%
    filter(startsWith(term, independent_variable)) 
}

adjusted_df <- filter_by_independent_variable(adjusted_dt ) %>% select(-independent_variable)

###############################################
#### Month 6 Linear Model Results Visualization 
#### This next code is to format the results into a figure. Markdown was used to format the variable names.

vars_order <- c("infant_sex", "infant_ethrace_2", "z_score", "c_diarrhea", 
                "bld_pca_weeks", "sam_season_2", "mother_ethrace_2", 
                "pre.bmi", "medicaid_in_preg", "parity_cat", "use_childcare", 
                "feeding_method_m1", "feeding_method_m6", "exposed_bm_m1", 
                "exposed_bm_m6", "bfdur_rec_m1_pct", "bfdur_bld_m6_pct", 

plot_data <- adjusted_df %>%
  mutate(
    # Create the "Pretty" Label with HTML Formatting
    nice_label = case_when(
      term == "bld_pca_weeks"            ~ "Age at Sample<br>(weeks)",
      term == "z_score"                  ~ "Weight<br>(Z-score)",
      term == "pre.bmi"                  ~ "Maternal BMI<br>(Pre-pregnancy)",
      term == "bfdur_rec_m1_pct"         ~ "BF Duration M1<br>(%)",
      term == "bfdur_bld_m6_pct"         ~ "BF Duration M6<br>(%)",
      term == "weeks_diff_pl"            ~ "Time from dosee<br>(weeks)",
      
      term == "feeding_method_m1Formula"      ~ "Feeding Type M1<br><i>(Formula vs. Excl. BF)</i><sup>1</sup>",
      term == "feeding_method_m1Combination"  ~ "Feeding Type M1<br><i>(Combination vs. Excl. BF)</i><sup>1</sup>",
      term == "feeding_method_m6Formula"      ~ "Feeding Type M6<br><i>(Formula vs. Excl. BF)</i><sup>1</sup>",
      term == "feeding_method_m6Combination"  ~ "Feeding Type M6<br><i>(Combination vs. Excl. BF)</i><sup>1</sup>",
      
      str_detect(term, "mother_ethrace_2")    ~ "Maternal Ethnicity-Race<br><i>(Non-White vs. White)</i><sup>1</sup>",
      str_detect(term, "infant_ethrace_2")    ~ "Infant Race<br><i>(Non-White vs. White)</i><sup>1</sup>",
      term == "infant_sexFemale"              ~ "Infant Sex<br><i>(Female vs. Male)</i><sup>1</sup>", 
      
      str_detect(term, "use_childcare")       ~ "Childcare Use<br><i>(Yes vs. No)</i><sup>1</sup>",
      str_detect(term, "medicaid_in_preg")    ~ "Medicaid<br><i>(Yes vs. No)</i><sup>1</sup>",
      str_detect(term, "exposed_bm_m1")       ~ "BM Exposure M1<br><i>(No vs. Yes)</i><sup>1</sup>",
      str_detect(term, "exposed_bm_m6")       ~ "BM Exposure M6<br><i>(No vs. Yes)</i><sup>1</sup>",
      str_detect(term, "sam_season_2")        ~ "Season<br><i>(Cold vs. Warm)</i><sup>1</sup>",
      str_detect(term, "c_diarrhea")          ~ "History of Gastroenteritis<br><i>(Yes vs. No)</i><sup>1</sup>",
      str_detect(term, "parity_cat")          ~ "Maternal Parity<br><i>(Yes vs. No)</i><sup>1</sup>",
      str_detect(term, "schedule_status")     ~ "Vaccine Schedule Status<br><i>(Full vs. Partial)</i><sup>1</sup>",
      TRUE ~ term),
    
    # Create a helper column to match 'term' to your 'vars_order' list
raw_var_match = case_when(
      str_detect(term, "infant_sex") ~ "infant_sex",
      str_detect(term, "infant_ethrace_2") ~ "infant_ethrace_2",
      str_detect(term, "z_score") ~ "z_score",
      str_detect(term, "c_diarrhea") ~ "c_diarrhea",
      str_detect(term, "bld_pca_weeks") ~ "bld_pca_weeks",
      str_detect(term, "sam_season_2") ~ "sam_season_2",
      str_detect(term, "mother_ethrace_2") ~ "mother_ethrace_2",
      str_detect(term, "pre.bmi") ~ "pre.bmi",
      str_detect(term, "medicaid_in_preg") ~ "medicaid_in_preg",
      str_detect(term, "parity_cat") ~ "parity_cat",
      str_detect(term, "use_childcare") ~ "use_childcare",
      str_detect(term, "feeding_method_m1") ~ "feeding_method_m1",
      str_detect(term, "feeding_method_m6") ~ "feeding_method_m6",
      str_detect(term, "exposed_bm_m1") ~ "exposed_bm_m1",
      str_detect(term, "exposed_bm_m6") ~ "exposed_bm_m6",
      str_detect(term, "bfdur_rec_m1_pct") ~ "bfdur_rec_m1_pct",
      str_detect(term, "bfdur_bld_m6_pct") ~ "bfdur_bld_m6_pct",
      str_detect(term, "schedule_status") ~ "schedule_status",
      str_detect(term, "weeks_diff_pl") ~ "weeks_diff_pl",
      TRUE ~ term),
    significance = ifelse(p.value < 0.05, "Significant\n(p<0.05)", "Non-Significant\n(p>0.05)"),
    p_label = ifelse(p.value < 0.05, paste0("p=", signif(p.value, 2)), NA_character_)) %>%
  mutate(order_index = match(raw_var_match, vars_order)) %>%
  arrange(desc(order_index)) %>% 
  mutate(nice_label = factor(nice_label, levels = unique(nice_label)))

### Forest Plot Month 6 
ggplot(plot_data, aes(x = estimate, y = nice_label, color = significance)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray20", linewidth = 1.0) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.3, color="grey20") +
  geom_point(aes(shape = significance), size = 4) +
  geom_text(aes(label = p_label), vjust = -0.8, size = 5, show.legend = FALSE, color = "black") + 
  scale_color_manual(values = c("Non-Significant\n(p>0.05)" = "black", "Significant\n(p<0.05)" = "#A9011BFF")) +
  scale_shape_manual(values = c("Non-Significant\n(p>0.05)" = 1, "Significant\n(p<0.05)" = 19)) +
  scale_x_continuous(limits = c(-6, 6), breaks = c(-6,-4, -2, 0, 2, 4, 6)) +
  labs(
    title = "Predictors of M6 Rotavirus-IgA Titers",
    x = "Estimate (95% CI)",
    y = NULL) +
  theme_bw()+
  theme(text=element_text(size = 16),
        title = element_text(size=10),
        legend.position = "bottom",
        axis.text.y = ggtext::element_markdown(size = 10, color = "black"), #this has to be markdown bc I used the HTML for fancy labls 
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_line(color = "gray95"))


################################################################################
################# Figure 5B: Month 12 Visit Forest Plot
################################################################################
 
#### Month 12 data setup
m12.dat<- dat %>% filter(monht==12) %>% 
   mutate(across(where(is.character), as.factor)) %>%
   mutate(across(where(is.factor), ~ fct_relevel(., levels(.)[which.max(table(.))]))) #This line is to order the categorical variables so the level with the most observation is first and used as reference in linear models
 
####### Month 12 Linear Models ########
 
 # Variables to test
vars <- c("bld_pca_weeks", "sam_season_2", "infant_sex","z_score", "pre.bmi", 
           "parity_cat", "c_diarrhea","medicaid_in_preg","infant_ethrace_2",
           "mother_ethrace_2","feeding_method_m6","feeding_method","exposed_bm_m6","bfdur_bld_m6_pct",
           "bf_duration_thru12mvisit","schedule_status","weeks_diff_pl")
 
adjusted_list <- lapply(formulas, function(fmla) {
   # Fit the model
   mod <- lm(fmla, data = m12.dat)
   independent_var <- all.vars(fmla)[2]
   tidy(mod, conf.int = TRUE) %>%
     mutate(model = "adjusted", independent_variable = independent_var)
 })
 
 # Combine results into a single data table
adjusted_dt <- rbindlist(adjusted_list, fill = TRUE)
 
filter_by_independent_variable <- function(df) {
   df %>%
     filter(startsWith(term, independent_variable)) 
 }
 
adjusted_df <- filter_by_independent_variable(adjusted_dt ) %>% select(-independent_variable)
 
vars_order <- c("infant_sex", "infant_ethrace_2", "z_score", "c_diarrhea", 
                 "bld_pca_weeks", "sam_season_2", "mother_ethrace_2", 
                 "pre.bmi", "medicaid_in_preg", "parity_cat", "use_childcare", 
                 "feeding_method", "feeding_method_m6", "exposed_bm_m1", 
                 "exposed_bm_m6", "bfdur_bld_m6_pct",  "bf_duration_thru12mvisit",
                 "schedule_status","weeks_diff_pl")

###############################################
#### Month 12 Linear Model Results Visualization 
#### This next code is to format the results into a figure. Markdown was used to format the variable names.
 
 plot_data <- adjusted_df %>%
   mutate(
     nice_label = case_when(
       term == "bld_pca_weeks"              ~ "Age at Sample<br>(weeks)",
       term == "z_score"                    ~ "Weight<br>(Z-score)",
       term == "pre.bmi"                    ~ "Maternal BMI<br>(Pre-pregnancy)",
       term == "weeks_diff_pl"              ~ "Time from dose<br>(weeks)",
       term == "bfdur_bld_m6_pct"           ~ "BF Duration M6<br>(% of month)",
       term == "bf_duration_thru12mvisit"   ~ "BF Duration M12<br>(weeks)",
       
       str_detect(term, "feeding_methodStill BF")       ~ "Feeding Type M12<br><i>(Still BF vs. Not BF)</i><sup>1</sup>",
       
       str_detect(term, "feeding_method_m6Formula")     ~ "Feeding Type M6<br><i>(Formula vs. Excl. BF)</i><sup>1</sup>",
       str_detect(term, "feeding_method_m6Combination") ~ "Feeding Type M6<br><i>(Combination vs. Excl. BF)</i><sup>1</sup>",
       
       str_detect(term, "mother_ethrace_2")    ~ "Maternal Race<br><i>(Non-White vs. White)</i><sup>1</sup>",
       str_detect(term, "infant_ethrace_2")    ~ "Infant Race<br><i>(Non-White vs. White)</i><sup>1</sup>",
       term == "infant_sexFemale"              ~ "Infant Sex<br><i>(Female vs. Male)</i><sup>1</sup>", 
       
       str_detect(term, "medicaid_in_preg")    ~ "Medicaid<br><i>(Yes vs. No)</i><sup>1</sup>",
       str_detect(term, "exposed_bm_m6")       ~ "BM Exposure M6<br><i>(No vs. Yes)</i><sup>1</sup>",
       str_detect(term, "sam_season_2")        ~ "Season<br><i>(Cold vs. Warm)</i><sup>1</sup>",
       str_detect(term, "c_diarrhea")          ~ "History of Gastroenteritis<br><i>(Yes vs. No)</i><sup>1</sup>",
       str_detect(term, "parity_cat")          ~ "Maternal Parity<br><i>(Yes vs. No)</i><sup>1</sup>",
       str_detect(term, "schedule_status")     ~ "Vaccine Schedule Status<br><i>(Full vs. Partial)</i><sup>1</sup>",
       TRUE ~ term ),
     raw_var_match = case_when(
       str_detect(term, "infant_sex") ~ "infant_sex",
       str_detect(term, "infant_ethrace_2") ~ "infant_ethrace_2",
       str_detect(term, "z_score") ~ "z_score",
       str_detect(term, "c_diarrhea") ~ "c_diarrhea",
       str_detect(term, "bld_pca_weeks") ~ "bld_pca_weeks",
       str_detect(term, "sam_season_2") ~ "sam_season_2",
       str_detect(term, "mother_ethrace_2") ~ "mother_ethrace_2",
       str_detect(term, "pre.bmi") ~ "pre.bmi",
       str_detect(term, "medicaid_in_preg") ~ "medicaid_in_preg",
       str_detect(term, "parity_cat") ~ "parity_cat",
       str_detect(term, "feeding_method_m6") ~ "feeding_method_m6", 
       str_detect(term, "feeding_method") ~ "feeding_method",
       str_detect(term, "exposed_bm_m6") ~ "exposed_bm_m6",
       str_detect(term, "bf_duration_thru12mvisit") ~ "bf_duration_thru12mvisit",
       str_detect(term, "bfdur_bld_m6_pct") ~ "bfdur_bld_m6_pct",
       str_detect(term, "schedule_status") ~ "schedule_status",
       str_detect(term, "weeks_diff_pl") ~ "weeks_diff_pl",
       TRUE ~ term),
     significance = ifelse(p.value < 0.05, "Significant\n(p<0.05)", "Non-Significant\n(p>0.05)"),
     p_label = ifelse(p.value < 0.05, paste0("p=", signif(p.value, 2)), NA_character_)) %>%
   mutate(order_index = match(raw_var_match, vars_order)) %>%
   arrange(desc(order_index)) %>% 
   mutate(nice_label = factor(nice_label, levels = unique(nice_label)))
 
### Forest Plot 
ggplot(plot_data, aes(x = estimate, y = nice_label, color = significance)) +
   geom_vline(xintercept = 0, linetype = "dashed", color = "gray30", linewidth = 1.0) +
   geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.3, color="grey40") +
   geom_point(aes(shape = significance), size = 4) +
   
   geom_text(aes(label = p_label), vjust = -0.8, size = 5, show.legend = FALSE, color = "black") + 
   
   scale_color_manual(values = c("Non-Significant\n(p>0.05)" = "black", "Significant\n(p<0.05)" = "#A9011BFF")) +
   scale_shape_manual(values = c("Non-Significant\n(p>0.05)" = 1, "Significant\n(p<0.05)" = 19)) +
   scale_x_continuous(limits = c(-7, 7), breaks = c(-6,-4, -2, 0, 2, 4,6)) +
   labs(
     title = "Predictors of M12 Rotavirus-IgA Titers",
     x = "Estimate (95% CI)",
     y = NULL) +
   theme_bw() +
   theme(text=element_text(size = 16),
         title = element_text(size=10),
         legend.position = "bottom",
         axis.text.y = ggtext::element_markdown(size = 10, color = "black"), 
         axis.ticks.y = element_blank(),
         panel.grid.major.y = element_line(color = "gray95"))
 