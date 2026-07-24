#' Author: Samuel Lucas da S. Delgado Mendes, Paulo Cesar de Paiva,
#' Rodolfo L. Nascimento
#' Subject: Deciding
#' Journal: Deciding
#' 
#' # Packages ---------------------------------------------------------------------

{
  lib <- .libPaths()[1] 
  required.packages <- c("vegan", "FD",
                         "writexl", "openxlsx",
                         "readxl",  "here",
                         "tidyverse", "car",
                         "ade4", "hclust", 
                         "gawdis", "ggeffects",
                         "MuMIn", "DHARMa", "corrplot",
                         "gridExtra", "FactoMineR", 
                         "factoextra", "jtools", "devtools") 
  i1 <- !(required.packages %in% row.names(installed.packages())) 
  if(any(i1)) { 
    install.packages(required.packages[i1], dependencies = TRUE, lib = lib) 
  } 
  lapply(required.packages, require, character.only = TRUE)
  
} # required packages

devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis) #adonis2

# Import -----------------------------------------------------------------------

{
  
  env<-read.xlsx(here("tables", "data.xlsx"), sheet = "env_betadiv")

}

# PCA -----------------------------------------------------------------------
#1. exploring collinearity 
{

  {
    pca.res1<-prcomp(env[c(2:10, 13)], scale. = TRUE)
    biplot(pca.res1)
    
    scores <- as.data.frame(pca.res1$x)
    loadings <- as.data.frame(pca.res1$rotation)
    
    scores$Month <- env$month
    scores$Region <- env$regiao
    
    
    ggplot(data = scores, aes(x = PC1, y = PC2, color = Region)) +
      geom_point(aes(size = 3, shape = Month)) +
      geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1*5, yend = PC2*5),
                   arrow = arrow(length = unit(0.2, "cm")), color = "black") +
      geom_text(data = loadings, aes(x = PC1*5, y = PC2*5, label = rownames(loadings)),
                color = "black", vjust = -0.5) +
      scale_shape_manual( values = c(16,17,18,15))+
      scale_color_manual(values = c("red", "royalblue"))+
      theme_bw() +
      xlab(paste("PC1 (", round(pca.res1$sdev[1]^2/sum(pca.res1$sdev^2)*100, 1), "%)",
                 sep = "")) +
      ylab(paste("PC2 (", round(pca.res1$sdev[2]^2/sum(pca.res1$sdev^2)*100, 1), "%)", 
                 sep = "")) +
      ggtitle("PCA Biplot with Calibrated Axes")
    
  } # PCA piblos com ggplot
  
  
  cor_matrix <- cor(env[,c(10:22)], use = "complete.obs")  
  corrplot(cor_matrix, method = "color", type = "upper", order = "hclust", 
           addCoef.col = "black", tl.cex = 0.8) # no collinearity
  
  
  
}

#' # PERMANOVA ---------------------------------------------------------------------



# Padronizando os dados (opcional, se necessário)
env_std = decostand(env[c(2:10, 13)], method = "standardize")

# Criando a matriz de distância (Euclidiana)
dist_env_all <- dist(env_std, method = "euclidean")

# Testando a dispersão dos grupos (homogeneidade de variâncias)
dispersion_test <- betadisper(dist_env_all, env$regiao)
permtest<- permutest(dispersion_test, pairwise = T)

dispersion_test2 <- betadisper(dist_env_all, env$month)
permtest2<- permutest(dispersion_test2, pairwise = T)


plot(dispersion_test)
plot(dispersion_test2)

# Rodando a PERMANOVA
permanova <- adonis2(dist_env_all ~ month * regiao, data = env, permutations = 999)
print(permanova)

write.xlsx(as.data.frame(permanova), "permanova.xlsx") 

# Pairwise permanova

paiwise_months = pairwise.adonis(dist_env_all, env$month, p.adjust.m = "bonferroni")
paiwise_regions = pairwise.adonis(dist_env_all, env$regiao, p.adjust.m = "bonferroni")

write.xlsx(as.data.frame(paiwise_period), "paiwise_period.xlsx") 
write.xlsx(as.data.frame(paiwise_groups), "paiwise_groups.xlsx") 

#plotting dispersion tests
{
  #dispersion test 1
  {
    # Extrair os scores e colocar em um data.frame
    scores_df <- as.data.frame(scores(dispersion_test, display = "sites"))  # só os pontos (sites)
    scores_df$month <- env$month
    scores_df$region <- env$regiao
    
    # Para dispersion_test (por Period)
    eigvals1 <- dispersion_test$eig
    var_exp1 <- eigvals1 / sum(eigvals1) * 100
    round(var_exp1[1:2], 1)
    
    # Centroides
    centroids <- as.data.frame(dispersion_test$centroids)
    centroids$Period <- rownames(centroids)
    
    # Plot com ggplot
    bt1 = ggplot(scores_df, aes(x = PCoA1, y = PCoA2, color = region)) +
      geom_point(size = 3, alpha = 0.8) +
      stat_ellipse(aes(fill = region), geom = "polygon", alpha = 0.2, color = NA) +
      geom_point(data = centroids, aes(x = PCoA1, y = PCoA2), 
                 shape = 4, size = 4, color = "black", stroke = 1.2) +
      theme_bw() +
      scale_color_manual(values = c( "dodgerblue", "tomato")) +
      scale_fill_manual(values = c("dodgerblue", "tomato")) +
      labs(x = "PCoA 1 (36.5%)", y = "PCoA 2 (21.2%)") 
  }
  
  #dispersion test 2
  {
    # Extrair as coordenadas dos pontos no espaço
    scores_df2 <- as.data.frame(scores(dispersion_test2, display = "sites"))
    scores_df2$months <- env$month  # adiciona a informação de grupo
    # Para dispersion_test2 (por Group)
    eigvals2 <- dispersion_test2$eig
    var_exp2 <- eigvals2 / sum(eigvals2) * 100
    round(var_exp2[1:2], 1)
    
    # Extrair os centroides (para o "X" preto)
    centroids2 <- as.data.frame(dispersion_test2$centroids)
    centroids2$Group <- rownames(centroids2)
    
    # Plot com ggplot2
    bt2 = ggplot(scores_df2, aes(x = PCoA1, y = PCoA2, color = months)) +
      geom_point(size = 3, alpha = 0.8) +
      stat_ellipse(aes(fill = months), geom = "polygon", alpha = 0.2, color = NA) +
      geom_point(data = centroids2, aes(x = PCoA1, y = PCoA2), 
                 shape = 4, size = 4, color = "black", stroke = 1.2) +
      theme_bw() +
      scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")) +
      scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"))+
      labs(x = "PCoA 1 (35.6%)", y = "PCoA 2 (21.2%)")
  }
  
  library(patchwork)
  prancha_bt = bt1 + bt2 + plot_layout(ncol = 2, nrow = 1)
  ggsave("S3.pdf", prancha_bt, width = 16, height = 8, dpi = 300)
}

str(env)

