## Overview

This package provides functions to fit homophily-adjusted network autocorrelation models.

## Installation

``` r
devtools::install_github('hanhtdpham/hanam')
```

## Example fitting HANE and HAND models to peer influence data

library(igraph)
library(tidygraph)
library(ggraph)
library(tidyverse)
library(magrittr)
library(hanam)
library(latentnet)
library(network)
library(sna)
library(Matrix)
# if(!require(NetData)) install.packages("https://cran.r-project.org/src/contrib/Archive/NetData/NetData_0.3.tar.gz",
#                                        type="source",repos=NULL)

set.seed(1)

# Import and clean data of students in 2nd semester -------------------
data('studentnets.peerinfl',package = "NetData")

# Get list of IDs with both network and attitudes data in 2nd sem
ids = intersect(unique(c(sem2$std_id,sem2$alter_id)), 
                attitudes$std_id)

# Clean attitudes data & calculate measure for classroom engagement
attitudes_sem2 =
  attitudes %>% 
  filter(sem_id == 2, 
         std_id %in% ids) %>% 
  group_by(std_id) %>% 
  summarize_if(is.integer,mean,na.rm=T) %>%
  mutate(engage = tlks + like_c + sub + tot + voln - misb) %>%
  select(std_id, sem_id, engage, call, imp, egrd) %>%
  na.omit()

# Filter network to focus on students with attitudes info
sem2 %<>%
  filter(std_id %in% attitudes_sem2$std_id,
         alter_id %in% attitudes_sem2$std_id)

# Create tidygraph object
peerinfl = 
  sem2 %>% 
  filter(relationship > 0) %>% 
  select(std_id,alter_id) %>% 
  graph_from_data_frame() %>%
  igraph::simplify() 

# Get giant component
peerinfl_gc = 
  peerinfl %>% 
  igraph::largest_component()

# Look at giant component
lo_gc = layout_nicely(peerinfl_gc)

peerinfl_gc %>%
  as_tbl_graph() %>%
  activate(nodes) %>% 
  mutate(name = as.integer(factor(name))) %>%
  ggraph(layout = lo_gc) +
  geom_edge_arc(edge_colour = gray(0.5,0.25),
                strength = 0.1) +
  geom_node_point(size = 2,
                  pch = 16) + 
  theme_graph()

# Get attitudes data of the giant component
attitudes_sem2_gc <-
  attitudes_sem2 %>%
  mutate(std_id = as.character(std_id)) %>%
  filter(std_id %in% V(peerinfl_gc)$name)


# Fit homophily adjusted network models for the giant component--------
# Create network object for ergmm input
D      <- 3
A_raw  <- as_adjacency_matrix(peerinfl_gc)
K      <- which.min(randnet::ECV.block(A_raw, max.K = 30, B=100)$l2)
g      <- network(A_raw, direct = T)

# organize rows of attitudes data to match adjacency matrix
attitudes_sem2_gc <- attitudes_sem2_gc[match(network.vertex.names(g), 
                                             paste(attitudes_sem2_gc$std_id)), ]
# identical(network.vertex.names(g), paste(attitudes_sem2_gc$std_id))

set.vertex.attribute(g, "call", attitudes_sem2_gc$call)
set.vertex.attribute(g, "imp" , attitudes_sem2_gc$imp)
set.vertex.attribute(g, "egrd", attitudes_sem2_gc$egrd)

# Obtain posterior sample of U from ergmm() 
ergmmFit  <- ergmm(g ~ euclidean(d=D, G = K) + rreceiver + rsender +
                     absdiff("call") + absdiff("imp") + absdiff("egrd"),
                    seed = 1, family = "Bernoulli.logit", verbose = 2,
                   tofit = c("mcmc", "procrustes"),
                   control = ergmm.control(sample.size = 2e3, burnin = 5e4,
                                           interval = 50))
# mcmc.diagnostics(ergmmFit)

# Row normalized A
A <- as.matrix(diag(1/sapply(rowSums(A_raw),function(x)ifelse(x==0,1,x)))%*%A_raw)

# Fit models
y <- attitudes_sem2_gc$engage
X <- cbind(1, 
           attitudes_sem2_gc$call, 
           attitudes_sem2_gc$imp, 
           attitudes_sem2_gc$egrd)
Udraw <- ergmmFit$sample$Z

haneFit      <- HANE(y = y, X = X, A = A, Usample = Udraw, verbose = T)
nemBayesFit  <- NEM( y = y, X = X, A = A)
nemMLEFit    <- lnam(y = y, x = X, W1= A)

round(cbind(HANE     = haneFit$coefs$estimates[c(1:4,8)], 
            nemBayes = nemBayesFit$coefs$estimates, 
            nemMLE   = c(nemMLEFit$beta, nemMLEFit$rho1)), 2)


handFit   <- HAND(y = y, X = X, A = A, Usample = Udraw)
ndmBayes  <- NDM( y = y, X = X, A = A)
ndmMLEFit <- lnam(y = y, x = X, W2= A)
round(cbind(HAND    = handFit$coefs$estimates[c(1:4,8)], 
            ndmBayes= ndmBayes$coefs$estimates,
            ndmMLE  = c(ndmMLEFit$beta, ndmMLEFit$rho2)), 2)
