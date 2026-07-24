
#' Author: Samuel Lucas da S. Delgado Mendes, Paulo Cesar de Paiva,
#' Rodolfo L. Nascimento
#' Subject: Functional dispersion framework reveals that
#' depth and upwelling occurrence drive the assembly of 
#' soft-bottom polychaete communities
#' Journal: Deciding

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
  
  trait_groups <- data.frame(
    modality = rownames(RLQ$c1),
    group = c("Size", "Size", "Size",
              "Feeding", "Feeding","Feeding", "Feeding","Feeding",
              "Reproduction","Reproduction","Reproduction", "Reproduction", "Reproduction",
              "Reproduction", "Reproduction", "Reproduction", "Reproduction", "Reproduction", "Reproduction",
              "Motility","Motility","Motility","Habitat", "Habitat", "Habitat",
              "Morphology", "Morphology", "Morphology", "Morphology", "Morphology", "Morphology", "Morphology",
              "Morphology", "Morphology")
  )
  
  
  
  library(dplyr)
  library(tidyr)
  
  # species scores (RLQ2)
  sp_scores <- data.frame(
    species = rownames(RLQ$lQ),
    RLQ2 = RLQ$lQ[, 2]
  )
  
  # trait modalities
  mod_scores <- data.frame(
    modality = rownames(RLQ$c1)
  )
  
  # combine everything (cartesian product)
  
  # 1. Convert Q to long format and FILTER for presence (value == 1)
  Q_long <- Q %>%
    as.data.frame() %>%
    tibble::rownames_to_column("species") %>%
    pivot_longer(-species, names_to = "modality", values_to = "presence") %>%
    filter(presence > 0) # THIS IS THE KEY: Only keep species that have the trait
  
  # 2. Join with the scores and groups
  plot_df <- Q_long %>%
    left_join(sp_scores, by = "species") %>%
    left_join(trait_groups, by = "modality")
  
  # 3. Calculate summary stats on the filtered data
  summary_df <- plot_df %>%
    group_by(group, modality) %>%
    summarise(
      median = median(RLQ2, na.rm = TRUE),
      q25 = quantile(RLQ2, 0.25, na.rm = TRUE),
      q75 = quantile(RLQ2, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  
  ggplot(plot_df, aes(x = RLQ2, y = modality)) +
    
    # species points
    geom_point(
      size = 0.8,
      alpha = 0.4,
      position = position_jitter(height = 0.1)
    ) +
    
    # IQR bars
    geom_segment(
      data = summary_df,
      aes(x = q25, xend = q75, y = modality, yend = modality),
      linewidth = 1
    ) +
    
    # median points
    geom_point(
      data = summary_df,
      aes(x = median, y = modality),
      size = 3
    ) +
    
    facet_wrap(~ group, scales = "free_y") +
    
    theme_bw() +
    labs(
      x = "Species position on RLQ axis 2",
      y = "Trait modalities"
    )
  
  
  