# Clean session
rm(list=ls())
pacman::p_unload(all)
gc()
remove.packages("hanam")

# Document, including making NAMESPACE.  Safest to delete old one first.
file.remove("C:/Users/dksewell/Documents/hanam/NAMESPACE")
devtools::document("C:/Users/dksewell/Documents/hanam")

# Build tarball
devtools::build(pkg = "C:/Users/dksewell/Documents/hanam",
                path = "C:/Users/dksewell/Downloads",
                vignettes = FALSE)

# Install from tarball
pacman::p_unload(hanam)
install.packages("C:/Users/dksewell/Downloads/hanam_1.0.0.tar.gz",
                 repos=NULL,type='source')
pacman::p_load(hanam)
beepr::beep(4)

# Install from github
remotes::install_github("dksewell/hanam")

{
  library(network)
  library(igraph)
  library(intergraph)
  library(ggraph)
  library(tidygraph)
  
  library(lubridate)
  library(janitor)
  library(tidyverse)
  library(magrittr)
  library(Matrix)
  library(RColorBrewer)
  library(kableExtra)
  
  library(sna)
  library(latentnet)
  if(!require(hanam)){remotes::install_github("hanhtdpham/hanam");library(hanam)}
  
  library(NetData)
  data(studentnets.peerinfl,package = "NetData")
  
  id_key = unique(c(sem1[,1],
                    sem1[,2],
                    sem2[,1],
                    sem2[,2],
                    attitudes$std_id))
  id_key = tibble(std_id = id_key,
                  new_id = 1:NROW(id_key))
  sem1 %<>%
    mutate(ego = id_key$new_id[match(sem1[,"std_id"],unlist(id_key[,"std_id"]))],
           alter = id_key$new_id[match(sem1[,"alter_id"],unlist(id_key[,"std_id"]))])
  sem2 %<>%
    mutate(ego = id_key$new_id[match(sem2[,"std_id"],unlist(id_key[,"std_id"]))],
           alter = id_key$new_id[match(sem2[,"alter_id"],unlist(id_key[,"std_id"]))])
  attitudes %<>%
    mutate(new_id = id_key$new_id[match(attitudes[,"std_id"],unlist(id_key[,"std_id"]))])
  
  N = nrow(id_key)
  
  friendship_length = sem1_A = sem2_A = Matrix(0L,N,N)
  sem1_A[as.matrix(sem1[,c("ego","alter")])] = 1L
  sem2_A[as.matrix(sem2[,c("ego","alter")])] = 1L
  friendship_length[as.matrix(sem1[,c("ego","alter")])] = sem1$timea
  friendship_length[which(is.na(friendship_length),arr.ind=T)] = 0L
  
  attitudes %<>%
    mutate(socialization = tlks/4 + frn/4 + cmt/4,
           class_attitude = 0.0)
  for(i in 1:nrow(attitudes)){
    attitudes$class_attitude[i] =
      mean(c(attitudes$ike_c[i]/4, attitudes$tlkt[i]/4,
             (1 + attitudes$egrd[i])/5, attitudes$jgrd[i]/3,
             attitudes$sub[i]/4, attitudes$tot[i]/4,
             attitudes$tch[i]/4, attitudes$voln[i]/5,
             (6-attitudes$misb[i])/5, attitudes$sftch[i]/5,
             attitudes$tfstd[i]/5, attitudes$chal[i]/4), na.rm=T)
  }
  
  actor_attributes =
    left_join(id_key,
              attitudes %>%
                group_by(new_id) %>%
                summarize_all(mean, na.rm=T),
              by = "new_id") %>%
    select(std_id.x, imp, socialization, class_attitude) %>%
    rename(std_id = std_id.x)
  
  ind_2_rm = unique(which(is.na(actor_attributes),arr.ind=T)[,1])
  
  sem1_A = sem1_A[-ind_2_rm,-ind_2_rm]
  sem2_A = sem2_A[-ind_2_rm,-ind_2_rm]
  friendship_length = friendship_length[-ind_2_rm,-ind_2_rm]
  actor_attributes = actor_attributes[-ind_2_rm,]
  
  sem1_A %<>% as.matrix()
  sem2_A %<>% as.matrix()
  friendship_length %<>% as.matrix()
  
  diag(sem1_A) = 0
  diag(sem2_A) = 0
  diag(friendship_length) = 0
  
  
  row_norm = function(A){
    rns = rowSums(A);
    rns[which(rns == 0.0)] = 1.0
    
    Diagonal(x = 1.0 / rns) %*% A
  }
  
  y = actor_attributes$class_attitude
  X = model.matrix(class_attitude ~ imp + socialization,
                   data = actor_attributes)
  
  sem1_A_rn = row_norm(sem1_A)
  friendship_length_rn = row_norm(friendship_length)
  
  lsm_sem1 = 
    readRDS(url("https://myweb.uiowa.edu/dksewell/teaching/SAND/lms_fit-sem1.RDS"))
  Udraw = lsm_sem1$sample$Z
  
  
  sem1_hane = 
    HANE(y,
         X,
         sem1_A_rn,
         Udraw,
         verbose = TRUE,
         Usample.eps = 1e-1)
  
  sem1_hand = 
    HAND(y,
         X,
         sem1_A_rn,
         Udraw,
         verbose = TRUE,
         Usample.eps = 1e-1)
  
  
  beepr::beep(4)
}

{
  sem1_hane
  sem1_hand
  
  summary(sem1_hane)
  summary(sem1_hane,
          CI_level = 0.8)
  summary(sem1_hand)
  summary(sem1_hand,
          CI_level = 0.8)
  
}