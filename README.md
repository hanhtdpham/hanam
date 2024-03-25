## Overview

This package provides functions to fit homophily-adjusted network autocorrelation models.

## Installation

``` r
devtools::install_github('hanhtdpham/hanam')
```

## Example fitting HANE and HAND models to peer influence data

``` r
if(!require(hanam)) devtools::install_github('hanhtdpham/hanam')


library(readxl)
library(igraph)
library(tidygraph)
library(ggraph)
library(tidyverse)
library(magrittr)
library(janitor)
library(Matrix)
library(network)
library(latentnet)
library(sna)
library(hanam)

set.seed(1)

# Import and clean data ---------------------------------------------

# Download data from Reality commons and set working directory via setwd().

# Import and clean call log data
callLog =
  read_csv("CallLog.csv") %>% 
  clean_names()
callLog %<>%
  select(contains("participant")) %>% 
  drop_na() %>% 
  rename(sender = participant_id_a,
         receiver = participant_id_b) %>% 
  filter(sender != receiver) %>% 
  distinct()

# Create tidygraph object
callLog %<>% 
  graph_from_data_frame() %>% 
  as_tbl_graph()
  
# Visualize the call log network
callLog %>%
  as_tbl_graph() %>%
  activate(nodes) %>% 
  mutate(name = as.integer(factor(name))) %>%
  ggraph(layout = layout_nicely(callLog)) +
  geom_edge_arc(edge_colour = gray(0.5,0.25),
                strength = 0.1) +
  geom_node_point(size = 2,
                  pch = 16) + 
  theme_graph()
  
# Calculate FunFit average daily activity
funfit =
  read_excel("FunFit.xlsx", sheet = "Daily Stats") %>%
  clean_names()
funfit %<>%
  rename(name = person_participant_id) %>% 
  group_by(name) %>% 
  summarize(activity_sum = mean(activity_sum))
  
#Import and clean gender variable
initSurvey1  = 
  read_excel("SurveyMonthly.xlsx", 
             sheet = "Init_survey.csv") %>% 
  clean_names() %>% 
  select(participant_id, gender, sleep_actual)
initSurvey2  = 
  read_excel("SurveyMonthly.xlsx", 
             sheet = "Init_survey_2010_10.csv") %>% 
  clean_names() %>% 
  select(participant_id, gender, sleep_actual)
initSurvey3  = 
  read_excel("SurveyMonthly.xlsx", 
             sheet = "Init_survey_2010_10_continuinue") %>% 
  clean_names() %>% 
  select(participant_id, gender, sleep_actual)

initSurvey = 
  as.data.frame(rbind(initSurvey1, initSurvey2, initSurvey3))

initSurvey %<>%
  rename(name = participant_id) %>%
  mutate(gender = tolower(substr(gender,1,1)),
         sleep_actual = fct_collapse(sleep_actual,
                                     `5 or less` = c("5", "less than 5"),
                                     `6` = "6",
                                     `7` = "7",
                                     `8 or more` = c("8", "9", "more than 9")))
initSurvey %<>% 
  group_by(name) %>% 
  mutate(n_NA = is.na(gender) + is.na(sleep_actual)) %>% 
  arrange(n_NA) %>% 
  mutate(keep = 1:length(gender)) %>% 
  ungroup() %>% 
  filter(keep == 1) %>% 
  select(-n_NA,-keep)
  
# Import and clean quality of life variable
{
ql = list()
ql[[1]] = 
  read_excel("SurveyMonthly.xlsx", 
             sheet = "2010_05.csv") %>%
  clean_names()
ql[[2]] = 
  read_excel("SurveyMonthly.xlsx", 
             sheet = "2010_11.csv") %>%
  clean_names()
}

for(j in 1:length(ql)){
  ql[[j]] %<>%
    rename(name = participant_id) %>% 
    select(name,contains("more_agree_disagree")) %>% 
    select(-more_agree_disagree_i_am_generally_not_very_happy_although_i_am_not_depressed_i_never_seem_as_happy_as_i_might_be) %>% 
    mutate(across(contains("more_agree_disagree"),
                  function(x){
                    ifelse(x == "Strongly Agree", 7,
                           ifelse(x == "Agree", 6,
                                  ifelse("Slightly Agree", 5,
                                         ifelse("Neutral",4,
                                                ifelse("Slightly Disagree",3,
                                                       ifelse("Disagree",2,
                                                              ifelse("Strongly Disagree",1,NA)))))))
                  })) %>%
    mutate(denominator =  
             8 - 
             is.na(more_agree_disagree_in_most_ways_my_life_is_close_to_my_ideal) - 
             is.na(more_agree_disagree_the_conditions_in_my_life_are_excellent) -
             is.na(more_agree_disagree_i_am_satisfied_with_my_life) -
             is.na(more_agree_disagree_so_far_i_have_gotten_the_important_things_i_want_in_life) -
             is.na(more_agree_disagree_if_i_could_live_my_life_over_i_would_change_almost_nothing) -
             is.na(more_agree_disagree_in_general_i_consider_myself_a_happy_person) -
             is.na(more_agree_disagree_compared_to_most_of_my_peers_i_consider_myself_more_happy) -
             is.na(more_agree_disagree_i_am_generally_very_happy_i_enjoy_life_regardless_of_what_is_going_on_getting_the_most_out_of_everything)) %>% 
    # mutate(across(where(is.numeric), ))
    mutate(across(contains("more_agree_disagree"), ~ ifelse(is.na(.x),0,.x/denominator))) %>% 
    select(-denominator) %>% 
    mutate(quality_of_life = 
             rowSums(across(where(is.numeric)))) %>% 
    select(name, quality_of_life)
}

ql = do.call(bind_rows,ql)
ql %<>% 
  group_by(name) %>% 
  summarize(quality_of_life = mean(quality_of_life))
  
# Add nodal attribute to network and keep complete cases
callLog_complete = callLog %>%
  activate(nodes) %>% 
  left_join(funfit,
            by = "name") %>% 
  left_join(initSurvey,
            by = "name") %>% 
  left_join(ql,
            by = "name") %>% 
  filter(!is.na(name),
         !is.na(gender),
         !is.na(sleep_actual),
         !is.na(quality_of_life),
         !is.na(activity_sum))
         
# Fit homophily adjusted network models for the giant component -----

# Fit latent space model

# Create network object for ergmm input:
D = 3
A = callLog_complete %>% 
  as_adjacency_matrix() %>% 
  as.matrix()

set.seed(1)
K = which.min(randnet::ECV.block(A, max.K = 15, B=100)$l2)
g = network(A, direct = T)

set.vertex.attribute(g, "gender", V(callLog_complete)$gender)
set.vertex.attribute(g, "sleep_actual", as.character(V(callLog_complete)$sleep_actual))
set.vertex.attribute(g, "quality_of_life" , V(callLog_complete)$quality_of_life)

# Row normalize A:
A = diag(1 / ifelse(rowSums(A) == 0,1,rowSums(A))) %*% A

# Obtain posterior sample of U from ergmm() (takes ~10 minutes):
set.seed(1)
ergmmFit  = ergmm(g ~ euclidean(d=D, G = K) + rreceiver + rsender + 
                    nodematch("gender") + nodematch("sleep_actual") + 
                    absdiff("quality_of_life"),
                  seed = 1, 
                  family = "Bernoulli.logit",
                  control = ergmm.control(sample.size = 5e3, burnin = 5e3))
                  
# Prepare inputs to fit models
y = scale(V(callLog_complete)$activity_sum)[,1]
X = model.matrix(~ gender + sleep_actual + quality_of_life,
                 data = 
                   callLog_complete %>% 
                   activate(nodes) %>% 
                   as_tibble() )
X[,ncol(X)] = scale(X[,ncol(X)])
Udraw = ergmmFit$sample$Z

# Effects model
haneFit      = HANE(y = y, X = X, A = A, Usample = Udraw)
nemBayesFit  = NEM( y = y, X = X, A = A)
nemMLEFit    = lnam(y = y, x = X, W1= A)

# Disturbances model
handFit      = HAND(y = y, X = X, A = A, Usample = Udraw)
ndmBayesFit  = NDM( y = y, X = X, A = A)
ndmMLEFit    = lnam(y = y, x = X, W2= A)

# Output 
p = ncol(X)
out = data.frame(
  Model = rep(c("Effects", "Disturbances"), each = p+1),
  Variable = rep( c(paste0("beta", 0:(p-1)), "rho"), 2), 
  hanam  = sprintf("%.2f (%.2f, %.2f)",
                   c(haneFit$coefs$estimates[c(1:p,p+D+1)], 
                     handFit$coefs$estimates[c(1:p,p+D+1)]), 
                   c(haneFit$coefs$LB[c(1:p,p+D+1)], 
                     handFit$coefs$LB[c(1:p,p+D+1)]), 
                   c(haneFit$coefs$UB[c(1:p,p+D+1)], 
                     handFit$coefs$UB[c(1:p,p+D+1)])), 
  Bayes  =  sprintf("%.2f (%.2f, %.2f)",
                    c(nemBayesFit$coefs$estimates, 
                      ndmBayesFit$coefs$estimates), 
                    c(nemBayesFit$coefs$LB, 
                      ndmBayesFit$coefs$LB), 
                    c(nemBayesFit$coefs$UB, 
                      ndmBayesFit$coefs$UB)), 
  MLE  = sprintf("%.2f (%.2f, %.2f)",
                 c(nemMLEFit$beta, nemMLEFit$rho1, 
                   ndmMLEFit$beta, ndmMLEFit$rho2),
                 c(nemMLEFit$beta, nemMLEFit$rho1, 
                   ndmMLEFit$beta, ndmMLEFit$rho2) - 
                   1.96*c(nemMLEFit$beta.se, nemMLEFit$rho1.se, 
                          ndmMLEFit$beta.se, ndmMLEFit$rho2.se), 
                 c(nemMLEFit$beta, nemMLEFit$rho1, 
                   ndmMLEFit$beta, ndmMLEFit$rho2) +
                   1.96*c(nemMLEFit$beta.se, nemMLEFit$rho1.se, 
                          ndmMLEFit$beta.se, ndmMLEFit$rho2.se))
)

out
```
