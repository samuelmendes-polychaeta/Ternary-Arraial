#' TERNARY PLOT ARRAIAL JANUARY

# Packages ---------------------------------------------------------------------

{
  lib <- .libPaths()[1] 
  required.packages <- c("vegan", "FD", "dplyr",
                         "writexl", "openxlsx",
                         "readxl",  "here", "tidyr",
                         "tidyverse", "car", "tidyselect",
                         "ade4", "pvclust", "tidygraph", 
                         "gawdis", "ggeffects", "patchwork",
                         "MuMIn", "DHARMa", "adiv", 
                         "adegraphics", "PERMANOVA",
                         "gridExtra", "FactoMineR", 
                         "factoextra", "jtools", "picante") 
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
  
  abund_july<-abund %>% filter(month=="jul_17")
  
  abund_january<-abund %>% filter(month=="jan_18")
  
  env_july<-data.frame(abund_july[2])
  env_january<-data.frame(abund_january[2])
  
  abund_july<-abund_july[4:44]
  abund_january<-abund_january[4:44]
  
}


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


#fundis <- dist_obj/max(dist_obj)
fundis <- gawdist

plot(hclust(fundis))

####JULY####
prop <- sweep(abund_january, 1, rowSums(abund_january), "/")

groups <- env_january$regiao

propsplitted <- split(prop, as.factor(groups))
prop_Exposed <- propsplitted$Exposed
prop_Protected <- propsplitted$Protected



frameDKG_Exposed<- betaUniqueness(prop_Exposed, fundis)
#frameDKG_out_17 <- betaUniqueness(prop_out_17, fundis)
frameDKG_Protected <- betaUniqueness(prop_Protected, fundis)
#frameDKG_abr_18 <- betaUniqueness(prop_abr_18, fundis)

####################
#Exposed############
###################

D_KG_Exposed<- frameDKG_Exposed$DKG  # Pairwise functional dissimilarities between plots of early successional stage. 
S_BC_Exposed<- 1-frameDKG_Exposed$DR # Pairwise species similarity between plots of early successional stage.
R_beta_Exposed<- frameDKG_Exposed$DR-frameDKG_Exposed$DKG

dim(D_KG_Exposed)

D_KG_bar_Exposed <- sapply(1:14, function(i) mean(D_KG_Exposed[i, -i]))
S_BC_bar_Exposed<- sapply(1:14, function(i) mean(S_BC_Exposed[i, -i]))
R_beta_bar_Exposed <- sapply(1:14, function(i) mean(R_beta_Exposed[i, -i]))


####################
###Protected########
####################

D_KG_Protected <- frameDKG_Protected$DKG # Pairwise functional dissimilarities between plots of late successional stage.
S_BC_Protected<- 1-frameDKG_Protected$DR # Pairwise species similarity between plots of late successional stage.
R_beta_Protected <- frameDKG_Protected$DR- frameDKG_Protected$DKG # Pairwise beta redundancy between plots of late successional stage.

dim(D_KG_Protected)

D_KG_bar_Protected <- sapply(1:15, function(i) mean(D_KG_Protected[i, -i]))
S_BC_bar_Protected<- sapply(1:15, function(i) mean(S_BC_Protected[i, -i]))
R_beta_bar_Protected<- sapply(1:15, function(i) mean(R_beta_Protected[i, -i]))


TAB_Exposed <- cbind.data.frame(D_KG_bar_Exposed, S_BC_bar_Exposed, R_beta_bar_Exposed)
TAB_Protected <- cbind.data.frame(D_KG_bar_Protected, S_BC_bar_Protected, R_beta_bar_Protected)


names(TAB_Exposed) <- names(TAB_Protected) <- c("D_KG", "S_BC", "R_beta")
TAB2 <- rbind.data.frame(TAB_Exposed, TAB_Protected)


library(ggtern)
groups_TAB2 <- factor(
  c(rep("Exposed", nrow(TAB_Exposed)), rep("Protected", nrow(TAB_Protected))),
  levels = c("Exposed", "Protected")
)

TAB_ggtern = data.frame(TAB2, groups_TAB2)

ptri_jan <-  ggtern(TAB_ggtern, aes(D_KG, R_beta, S_BC, color = groups_TAB2)) +
  geom_point(size = 3, alpha = 0.70) +
  scale_color_manual(values = c(
    "Protected"     = "red",
    "Exposed"     = "blue"
  )) +  theme_bw();ptri_jan

#good themes: custom, 


#PERMANOVA TEST
TAB4Ptest <- DistContinuous(TAB2, , "Bray_Curtis")

Ptest <- PERMANOVA(TAB4Ptest, as.factor(groups), nperm=10000, PostHoc = "fdr")

Ptest
summary(Ptest)

#Anovas

df <- data.frame(TAB2,  groups = groups_TAB2)


resultado_univariado_D_KG <- adonis2(
  df$D_KG  ~ groups,
  data = df,
  permutations = 9999,
  method = "euclidean" # Usa a distância euclidiana, que é equivalente à soma dos quadrados na univariada
)



resultado_univariado_S_BC <- adonis2(
  df$S_BC  ~ groups,
  data = df,
  permutations = 9999,
  method = "euclidean" # Usa a distância euclidiana, que é equivalente à soma dos quadrados na univariada
)



resultado_univariado_R_beta <- adonis2(
  df$R_beta  ~ groups,
  data = df,
  permutations = 9999,
  method = "euclidean" # Usa a distância euclidiana, que é equivalente à soma dos quadrados na univariada
)

table <- data.frame(
  Partition = c(
    "D_KG_bar_Protected", "D_KG_bar_Exposed",
    "R_beta_bar_Protected", "R_beta_bar_Exposed",
    "S_BC_bar_Protected", "S_BC_bar_Exposed"),
  Mean = c(
    mean(D_KG_bar_Protected), mean(D_KG_bar_Exposed),
    mean(R_beta_bar_Protected), mean(R_beta_bar_Exposed), 
    mean(S_BC_bar_Protected), mean(S_BC_bar_Exposed)),
  SD = c(
    sd(D_KG_bar_Protected), sd(D_KG_bar_Exposed),
    sd(R_beta_bar_Protected), sd(R_beta_bar_Exposed), 
    sd(S_BC_bar_Protected), sd(S_BC_bar_Exposed)))


# 1. Extract the p-values from your adonis objects
p_values <- c(
  D_KG = resultado_univariado_D_KG$`Pr(>F)`[1],
  S_BC = resultado_univariado_S_BC$`Pr(>F)`[1],
  R_beta = resultado_univariado_R_beta$`Pr(>F)`[1]
)
   
# 2. Apply Bonferroni correction
p_corrected <- p.adjust(p_values, method = "bonferroni")

# 3. View the results
print(p_corrected)

# ================================
# Stacked barplot – mean partition contribution per month
# ================================

# 1. Calcular médias por mês
mean_region <- data.frame(
  Region = c("Exposed", "Protected"),
  D_KG   = c(mean(D_KG_bar_Exposed),
             mean(D_KG_bar_Protected)),
  S_BC   = c(mean(S_BC_bar_Exposed),
             mean(S_BC_bar_Protected)),
  R_beta = c(mean(R_beta_bar_Exposed),
             mean(R_beta_bar_Protected))
)

t(mean_region)

# 2. Converter para formato longo
mean_region_long <- pivot_longer(
  mean_region,
  cols = c(D_KG, S_BC, R_beta),
  names_to = "Partition",
  values_to = "Value"
)

# (opcional, mas recomendável para interpretação):
# padronizar para proporções (cada barra soma 1)
mean_region_long <- mean_region_long %>%
  group_by(Region) %>%
  mutate(Value = Value / sum(Value))

mean_region_long$Region <- factor(
  mean_region_long$Region,
  levels = c("Exposed", "Protected")
)

p_bar_jan <- ggplot(mean_region_long,
                aes(x = Region, y = Value, fill = Partition)) +
  geom_bar(stat = "identity", width = 0.65, color = "black") +
  scale_fill_manual(
    values = c(
      "D_KG"   = "black",
      "S_BC"   = "grey",
      "R_beta" = "white"
    ),
    labels = c(
      "D_KG"   = "Functional dissimilarity (DKG)",
      "S_BC"   = "Species similarity (SBC)",
      "R_beta" = "Beta redundancy (RBETA)"
    )
  ) +
  labs(
    x = "Regions",
    y = "Relative contribution (%)",
    fill = "Partition"
  ) +
  theme_classic(base_size = 14)

p_bar_jan

final_plot_jan <- ptri_jan + p_bar_jan + 
  plot_layout(guides = 'collect') &  # This collects and shares the legend
  theme(legend.position = 'bottom'); final_plot_jan   # Position the shared legend 



