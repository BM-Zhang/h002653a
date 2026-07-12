
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, psych, lattice, corrplot, lavaan, semPlot, semptools, readr, stringr)



iri_files <- list.files(pattern = ".*_IRI_Scale_.*\\.csv")
iri_raw <- map_dfr(iri_files, read_csv, show_col_types = FALSE)

iri_wide <- iri_raw %>%
  select(SubjectID, Item, Score) %>%
  pivot_wider(names_from = Item, values_from = Score) %>%
  mutate(SubjectID = as.numeric(SubjectID)) %>%
  mutate(
    Q2  = 4 - Q2,
    Q5  = 4 - Q5,
    Q10 = 4 - Q10,
    Q11 = 4 - Q11,
    Q14 = 4 - Q14
  )



exp_files <- list.files(pattern = ".*_multimodalemotion_.*\\.csv")
exp_raw <- map_dfr(exp_files, read_csv, show_col_types = FALSE)

exp_clean <- exp_raw %>%
  mutate(
    RT = parse_number(as.character(mouse_resp.time)),
    ACC = parse_number(as.character(mouse_resp.corr)),
    clicked_emotion = str_extract(str_to_lower(mouse_resp.clicked_name), "happy|sad|anger|fear|disgust|surprise"),
    trial_type = str_to_lower(trial_type)
  )

exp_clean <- exp_clean %>%
  group_by(participant) %>%
  mutate(
    rt_mean = mean(RT, na.rm = TRUE),
    rt_sd = sd(RT, na.rm = TRUE),
    is_outlier = abs(RT - rt_mean) > (3 * rt_sd)
  ) %>%
  filter(!is_outlier) %>% 
  ungroup()

valid_participants <- exp_clean %>%
  group_by(participant) %>%
  summarise(overall_mean_RT = mean(RT, na.rm = TRUE)) %>%
  filter(overall_mean_RT > 3.0) %>%  
  pull(participant)

exp_clean <- exp_clean %>%
  filter(participant %in% valid_participants)

print(paste("过滤后，满足平均RT>3秒的有效实验被试数为:", length(valid_participants)))

exp_summary <- exp_clean %>%
  group_by(participant, trial_type) %>%
  summarise(
    Mean_RT = mean(RT[ACC == 1], na.rm = TRUE), 
    N_actual = n(),
    N_hits = sum(ACC == 1, na.rm = TRUE),
    .groups = "drop"
  )

response_counts <- exp_clean %>%
  group_by(participant, clicked_emotion) %>%
  summarise(N_responses = n(), .groups = "drop") %>%
  rename(trial_type = clicked_emotion) %>%
  filter(!is.na(trial_type))

exp_final <- exp_summary %>%
  left_join(response_counts, by = c("participant", "trial_type")) %>%
  mutate(
    N_responses = replace_na(N_responses, 0),
    Wagner_Hu = ifelse(N_actual > 0 & N_responses > 0, 
                       (N_hits / N_actual) * (N_hits / N_responses), 0)
  ) %>%
  select(participant, trial_type, Mean_RT, Wagner_Hu) %>%
  pivot_wider(
    names_from = trial_type, 
    values_from = c(Mean_RT, Wagner_Hu),
    names_glue = "{trial_type}_{.value}"
  ) %>%
  rename(SubjectID = participant) %>%
  mutate(SubjectID = as.numeric(SubjectID))



full_data <- inner_join(iri_wide, exp_final, by = "SubjectID")
mydata <- full_data %>% drop_na(Q1)
write_csv(mydata, "1_Merged_Data.csv")
print(paste("最终合并成功的配套被试人数为:", nrow(mydata)))

iri_cols <- mydata %>% select(paste0("Q", 1:22))

d2 <- outlier(iri_cols, plot = FALSE) 
alpha_level <- 0.001
p <- ncol(iri_cols) 
cutoff <- qchisq(1 - alpha_level, df = p)
is_outlier_subj <- d2 > cutoff
num_outliers <- sum(is_outlier_subj, na.rm = TRUE)

cat("Statistical Cutoff (D2):", round(cutoff, 2), "\n")
cat("Number of participants flagged (Mahalanobis):", num_outliers, "\n")

if(num_outliers > 0) {
  worst_indices <- order(d2, decreasing = TRUE)[1:num_outliers]
  clean_means <- colMeans(iri_cols[-worst_indices, ], na.rm = TRUE)
  person_corrs <- apply(iri_cols[worst_indices, ], 1, function(x) {
    cor(as.numeric(x), clean_means, use = "complete.obs")
  })
  results <- data.frame(
    SubjectID = mydata$SubjectID[worst_indices],
    Mahalanobis_D2 = round(d2[worst_indices], 2),
    Person_Rest_Cor = round(person_corrs, 3)
  )
  print(results)
}

clean_data <- mydata[!is_outlier_subj, ]
cat("排雷后，准备进入 CFA 的最终被试数:", nrow(clean_data), "\n")



model_iri <- ' 
  PT =~ Q6 + Q9 + Q15 + Q19 + Q22
  FS =~ Q3 + Q5 + Q10 + Q12 + Q17 + Q20
  EC =~ Q1 + Q2 + Q7 + Q11 + Q14 + Q16
  PD =~ Q4 + Q8 + Q13 + Q18 + Q21
'

fit_iri <- cfa(model_iri, data = clean_data, estimator = "MLR")

summary(fit_iri, fit.measures = TRUE, standardized = TRUE)

p_plot <- semPaths(fit_iri, 
                   style="lisrel", 
                   residuals = FALSE,
                   what = "std",           
                   whatLabels = "std", 
                   minimum = 0.2,          
                   edge.label.cex = 0.6, 
                   layout ='tree2',
                   edge.label.color = "black", 
                   rotation = 2, 
                   label.scale = FALSE, 
                   label.cex = 0.7,
                   equalizeManifests = FALSE, 
                   optimizeLatRes = TRUE, 
                   node.width = 1.0, 
                   edge.width = 0.5, 
                   shapeMan = "rectangle", 
                   shapeLat = "ellipse", 
                   shapeInt = "triangle", 
                   sizeMan = 3.5, 
                   sizeInt = 2, 
                   sizeLat = 5, 
                   curve=2, 
                   unCol = "#070b8c")

plot(p_plot)

print(modindices(fit_iri, sort = TRUE, maximum.number = 20))

model_iri_revised <- ' 
  PT =~ Q6 + Q9 + Q15 + Q19 + Q22
  FS =~ Q3 + Q5 + Q10 + Q12 + Q17 + Q20
  EC =~ Q1 + Q2 + Q7 + Q11 + Q14 + Q16
  PD =~ Q4 + Q8 + Q13 + Q18 + Q21

  Q5 ~~ Q10
  Q2 ~~ Q14
'
fit_iri_revised <- cfa(model_iri_revised, data = clean_data, estimator = "MLR")
summary(fit_iri_revised, fit.measures = TRUE, standardized = TRUE)

model_emotion <- '

Emotion =~
happy_Wagner_Hu +
sad_Wagner_Hu +
anger_Wagner_Hu +
fear_Wagner_Hu +
disgust_Wagner_Hu +
surprise_Wagner_Hu

'

fit_emotion <- cfa(
  model_emotion,
  data = clean_data,
  estimator = "MLR"
)

summary(fit_emotion,
        fit.measures = TRUE,
        standardized = TRUE)

semPaths(fit_emotion,
         style="lisrel",
         what="std",
         whatLabels="std",
         layout="tree2",
         residuals=FALSE)

modindices(fit_emotion,
           sort=TRUE,
           maximum.number=20)

describe(clean_data %>%
           select(ends_with("Wagner_Hu")))

emotion_data <- clean_data %>%
  select(
    happy_Wagner_Hu,
    sad_Wagner_Hu,
    anger_Wagner_Hu,
    fear_Wagner_Hu,
    disgust_Wagner_Hu,
    surprise_Wagner_Hu
  )

alpha(emotion_data)
omega(emotion_data)
splitHalf(emotion_data)

model_sem <- '

PT =~ Q6 + Q9 + Q15 + Q19 + Q22
FS =~ Q3 + Q5 + Q10 + Q12 + Q17 + Q20
EC =~ Q1 + Q2 + Q7 + Q11 + Q14 + Q16
PD =~ Q4 + Q8 + Q13 + Q18 + Q21

Emotion =~
happy_Wagner_Hu +
sad_Wagner_Hu +
anger_Wagner_Hu +
fear_Wagner_Hu +
disgust_Wagner_Hu +
surprise_Wagner_Hu

Emotion ~ PT + FS + EC + PD

Q5 ~~ Q10
Q2 ~~ Q14

'

fit_sem <- sem(
  model_sem,
  data = clean_data,
  estimator = "MLR"
)

summary(
  fit_sem,
  fit.measures = TRUE,
  standardized = TRUE,
  rsquare = TRUE
)

fitMeasures(
  fit_sem,
  c("chisq","df","pvalue",
    "cfi","tli",
    "rmsea","srmr",
    "aic","bic")
)

standardizedSolution(fit_sem)

modindices(
  fit_sem,
  sort = TRUE,
  maximum.number = 20
)

semPaths(
  fit_sem,
  style="lisrel",
  what="std",
  whatLabels="std",
  layout="tree2",
  residuals=FALSE,
  edge.label.cex=0.7
)

model_EC <- '

EC =~
Q1 + Q2 + Q7 + Q11 + Q14 + Q16

Emotion =~
happy_Wagner_Hu +
sad_Wagner_Hu +
anger_Wagner_Hu +
fear_Wagner_Hu +
disgust_Wagner_Hu +
surprise_Wagner_Hu

Emotion ~ EC

'

fit_EC <- sem(
  model_EC,
  data = clean_data,
  estimator = "MLR"
)

summary(
  fit_EC,
  fit.measures = TRUE,
  standardized = TRUE,
  rsquare = TRUE
)

fitMeasures(fit_EC,
            c("cfi","tli","rmsea","srmr","aic","bic"))
