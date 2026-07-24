#' Author: Samuel Lucas da S. Delgado Mendes, Paulo Cesar de Paiva,
#' Rodolfo L. Nascimento
#' Subject: Functional dispersion framework reveals that
#' depth and upwelling occurrence drive the assembly of 
#' soft-bottom polychaete communities
#' Journal: Marine Ecology Progress Series

# Packages ---------------------------------------------------------------------

{
  lib <- .libPaths()[1] 
  required.packages <- c("vegan", "writexl", "openxlsx",
                         "readxl",  "here","tidyverse", 
                         "ade4", "hclust", "gawdis") 
  i1 <- !(required.packages %in% row.names(installed.packages())) 
  if(any(i1)) { 
    install.packages(required.packages[i1], dependencies = TRUE, lib = lib) 
  } 
  lapply(required.packages, require, character.only = TRUE)
  
} # required packages


# Import -----------------------------------------------------------------------

{
  traits<-read.xlsx(here("tables", "data.xlsx"), sheet = "selected_traits") 
  
  abund<-read.xlsx(here("tables", "data.xlsx"), sheet = "abund_betadiv") 
  
  env<-read.xlsx(here("tables", "data.xlsx"), sheet = "env_betadiv")
  
}


# Fuzzy scores calculation -----------------------------------------------------

# Fuzzy trait calculation ------------------------------------------------------

#1. fuzzy scores and distance matrix
{
  
  #size
  {
    length<-prep.fuzzy(traits[2:4], col.blocks = 3)
  }
  
  #feeding strategy
  {
    feeding<-prep.fuzzy(traits[5:9], col.blocks = 5)
  }
  
  #reproduction
  {
    larval_development<-prep.fuzzy(traits[10:12], col.blocks = 3)
    brooding_strategy<-prep.fuzzy(traits[13:18], col.blocks = 6)
    assexual_reproduction<-prep.fuzzy(traits[19:20], col.blocks = 2)
  }
  
  #motility
  {
    motility<-prep.fuzzy(traits[21:23], col.blocks = 3)
  }
  
  #habitat occupation
  {
    substrate_occupation<-prep.fuzzy(traits[24:26], col.blocks = 3)
  }
  #morphology
  {
    
    head_apd <- prep.fuzzy(traits[27:29], col.blocks = 3)
    palps<- prep.fuzzy(traits[30:33], col.blocks = 4)
    teeth<-prep.fuzzy(traits[34:35], col.blocks = 2)
  }
  
  #fuzzy scores dataframe
  {
    db.trait = data.frame(length, feeding, larval_development, 
                          brooding_strategy, assexual_reproduction,motility,
                          substrate_occupation, head_apd, palps, teeth)
    rownames(db.trait) = traits$Genus
  }
  
  #Gower distance matrix
  {
    gawdist<- gawdis(db.trait, w.type = "equal", 
                     groups = c(1,1,1,
                                2,2,2,2,2,
                                3,3,3,
                                4,4,4,4,4,4,
                                5,5,
                                6,6,6,
                                7,7,7,
                                8,8,8,
                                9,9,9,9,
                                10,10), fuzzy = TRUE) #gower distance among genera
    
    cluster = hclust(gawdist) #checking genera functional dendrogram
    plot(cluster)
    
  }
  
}




# RLQ and fourthcorner analysis ------------------------------------------------

{
  #Data preparation
  {
    {
      R = data.frame(env[c(2:10,13)]) #for all continuous variables
      rownames(R) <- env$id 
      str(R) 
      
      L = data.frame(abund[4:44]) 
      rownames(L) <- abund$comm
      str(L)
      
      L_pa <- ifelse(L > 0, 1, 0)
      
      
      # Combining all fuzzy objects into a single data frame
      Q <- data.frame(db.trait)
      
      # Defining col.blocks (number of trait modalities, i.e. columns, for each fuzzy trait)
      
      col.blocks<- c(1,1,1,
                     2,2,2,2,2,
                     3,3,3,
                     4,4,4,4,4,4,
                     5,5,
                     6,6,6,
                     7,7,7,
                     8,8,8,
                     9,9,9,9,
                     10,10)
      
      
      str(col.blocks)
      
      attr(Q, "col.blocks") <- col.blocks
      
      row.w <- rep(1, nrow(Q)) #equal weights for fuzzy coded traits' FCA
      
      attr(Q, "row.w") <- row.w
      
    }#for RLQ and fourthcorner (preparation of R, L and Q objects)
  }
  
  #1. ordinations
  {
    {
      L_dudi <- dudi.coa(L, scannf = FALSE, nf = 2) #COA for L matrix
      R_dudi <- dudi.pca(R, row.w = L_dudi$lw, scale = T, scannf = FALSE, nf = 2) #PCA for R matrix
      Q_dudi <- dudi.fca(Q, scannf = FALSE, nf = 2) #FCA for Q matrix
      
      #Correcting line weights after the ordinations for the RLQ analysis
      R_dudi$lw <- L_dudi$lw  
      Q_dudi$lw <- L_dudi$cw  
    }
  }
  
  #2. RLQ analysis and significance
  {
    RLQ <- rlq(dudiR = R_dudi, dudiL = L_dudi, dudiQ = Q_dudi, scannf = FALSE, nf = 2) #RLQ analysis
    plot(RLQ) #plotting the RLQ 
    summary(RLQ) 
    
    rand_test_abiotic = randtest(RLQ, nrepet = 49999, modeltype = 6) #Significance test for the multivariate pattern
    
    write.xlsx(as.data.frame( summary(rand_test_abiotic)), "rand_test_abiotic.xlsx") 
    
    adjusted_p_value <- p.adjust(rand_test_abiotic$pvalue, method = "fdr") #Adjusting p values according to False Discovery Rate (FDR) method
    
    p_value_SRLQ = fourthcorner2(R, L, Q,
                                 modeltype = 6, p.adjust.method.G = "fdr", nrepet = 49999)#Fourthcorner statistic 
    
  }
  
  #3. Fourthcorner analysis
  {
    fc_adj = fourthcorner(tabR = R, tabL = L, tabQ = Q, 
                          modeltype = 6, p.adjust.method.G = "fdr", 
                          p.adjust.method.D = "fdr", nrepet = 49999) #Fourthcorner analysis (trait x enviroment correlations and their p values FDR corrected)
    
    write.xlsx(as.data.frame( summary(fc_adj)), "fc_adj.xlsx") 
    
    plot(fc_adj)
    plot(fc_adj, x.rlq = RLQ, stat = "D2", type = "biplot")
    
    summary(fc_adj)
    
    fc_Q = fourthcorner.rlq(RLQ, modeltype = 6,
                            typetest = "Q.axes", nrepet = 49999, p.adjust.method.G = "fdr",
                            p.adjust.method.D = "fdr") #Significance of each trait correlation with RLQ axes (FDR corrected)
    
    write.xlsx(as.data.frame( summary(fc_Q)), "fc_Q.xlsx") 
    
    fc_R = fourthcorner.rlq(RLQ, modeltype = 6,
                            typetest = "R.axes", nrepet = 49999, p.adjust.method.G = "fdr",
                            p.adjust.method.D = "fdr") #Significance of each enviromental variable correlation with the RLQ axes (FDR corrected)
    write.xlsx(as.data.frame( summary(fc_R)), "fc_R.xlsx") 
    
  }
}


r_weights <- RLQ$l1 
l_weights <- RLQ$lQ 
q_weights <- RLQ$c1


traits_abund <- RLQ$mQ


library(pvclust)

dist_sp = vegdist(x = l_weights, method = "euclidean")
dist_traits = vegdist(x = q_weights, method = "euclidean")

plot(hclust(dist_sp))
plot(hclust(dist_traits))


site_scores <- as.data.frame(RLQ$lR)
plot_scores = data.frame(site_scores, env)

# Adicionar setas e rótulos
env_r = s.arrow(r_weights, clabel = 1, boxes = F)
species = s.arrow(l_weights, clabel = 1, boxes = F)
traits = s.arrow(q_weights, clabel = 1, boxes = F) #juntos

tamanho = data.frame(q_weights[1:3,])
alimentação = data.frame(q_weights[4:8,])
reprodução = data.frame(q_weights[9:19,])
mobilidade = data.frame(q_weights[20:22,])
habitat = data.frame(q_weights[23:25,])
morfologia = data.frame(q_weights[26:34,])

s.arrow(morfologia)

r_PLOT = ggplot(plot_scores, aes(x = AxcR1, y = AxcR2, color = regiao, shape = month)) + 
  geom_point(size = 5, alpha = 0.9) + 
  theme_bw() + 
  labs(x = "Axis 1", y = "Axis 2", color = "local", shape = "month") + 
  geom_hline(yintercept = 0, color = "black", linetype = "solid") + 
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  scale_color_manual(values = c("blue", "red"))


r_PLOT = ggplot(plot_scores, aes(x = AxcR1, y = AxcR2, color = region, shape = month)) + 
  geom_point(size = 5, alpha = 0.9) + 
  theme_bw() + 
  labs(x = "Axis 1", y = "Axis 2", color = "Região", shape = "Mês") + 
  geom_hline(yintercept = 0, color = "black", linetype = "solid") + 
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  scale_color_manual(
    values = c("red", "blue"),
    labels = c("protegida", "exposta")
  )

ggsave("arraial_rlq.tiff", r_PLOT, dpi = 300)


