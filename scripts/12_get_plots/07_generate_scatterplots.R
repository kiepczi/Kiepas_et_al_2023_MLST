library('ggplot2')
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggplot2)
library('ggExtra')

#Generate Genome coverage pyANI plots for connected components
data <- read.csv("../../data/pyani_analysis/pyani_coevarge_identity_df.csv")

data2 <- data[!data$unique_taxa_names_per_connected_component %in% c('6', '5','4', '15', '14', '8','7', '3', '9'), ]
data3 <- data[!data$unique_taxa_names_per_connected_component %in% c('1', '2'), ]
pdf(file='../../figures/pyani_analysis_connected_components/scatter_plots/coverage_connected_components.pdf', width=30, height=15)


p1<-ggplot(data2, aes(x=cluster_id, y=coverage, col=comparision_type_species)) +geom_jitter(alpha=0.2)  + facet_wrap(unique_taxa_names_per_connected_component ~., scales = 'free_x', nrow=1, shrink=TRUE, dir='h') + scale_y_continuous(name="Genome Coverage (%)") +
  scale_x_continuous(breaks = seq(0, 50, 5), name ="") + theme(legend.position="none") 

p2<-ggplot(data3, aes(x=cluster_id, y=coverage, col=comparision_type_species)) +geom_jitter(alpha=0.2) + 

  facet_wrap(unique_taxa_names_per_connected_component ~., scales = 'free_x', nrow=2, shrink=TRUE, dir='h') + scale_y_continuous(name="") +
  scale_x_continuous(breaks = seq(0, 10, 1), name ="Cluster ID") + theme(legend.position="none")

p4<-ggarrange(p1, p2,
              ncol = 1, nrow = 2, heights = c(1, 1, 1)) 
plot(p4)
dev.off()

#Generate Genome identity pyANI plots for connected components
data <- read.csv("../../data/pyani_analysis/pyani_coevarge_identity_df.csv")

data2 <- data[!data$unique_taxa_names_per_connected_component %in% c('6', '5','4', '15', '14', '8','7', '3', '9'), ]
data3 <- data[!data$unique_taxa_names_per_connected_component %in% c('1', '2'), ]
pdf(file='../../figures/pyani_analysis_connected_components/scatter_plots/identity_connected_components.pdf', width=30, height=15)


p1<-ggplot(data2, aes(x=cluster_id, y=identity, col=comparision_type_genus)) +geom_jitter(alpha=0.2)  + facet_wrap(unique_taxa_names_per_connected_component ~., scales = 'free_x', nrow=1, shrink=TRUE, dir='h') + scale_y_continuous(name="Genome Identity (%)") +
  scale_x_continuous(breaks = seq(0, 50, 5), name ="") + theme(legend.position="none") 

p2<-ggplot(data3, aes(x=cluster_id, y=identity, col=comparision_type_genus)) +geom_jitter(alpha=0.2) + 
  
  facet_wrap(unique_taxa_names_per_connected_component~., scales = 'free_x', nrow=2, shrink=TRUE, dir='h') + scale_y_continuous(name="") +
  scale_x_continuous(breaks = seq(0, 10, 1), name ="Cluster ID") + theme(legend.position="none")

p4<-ggarrange(p1, p2,
              ncol = 1, nrow = 2, heights = c(1, 1, 1)) 
plot(p4)
dev.off()



#Generate Genome identity pyANI plots for ST comparisions
data <- read.csv("~/Desktop/Kiepas_et_al_2023_MLST/data/pyani_analysis/pyani_coevarge_identity_ST_comparisions_df.csv")
data2 <- data[!data$unique_taxa_names_per_ST %in% c('3', '4'), ]
data3 <- data[!data$unique_taxa_names_per_ST %in% c('1', '2'), ]

pdf(file='../../figures/pyani_analysis_connected_components/scatter_plots/identity_STs.pdf', width=30, height=15)

p1<-ggplot(data2, aes(x=cluster_id, y=identity, col=comparision_type_genus)) +geom_jitter(alpha=0.2)  + facet_wrap(unique_taxa_names_per_ST ~., scales = 'free_x', nrow=4, shrink=TRUE, dir='h') + scale_y_continuous(name="Genome Identity (%)") +
  scale_x_continuous(breaks = seq(0, 100, 2), name ="") + theme(legend.position="none") 

p2<-ggplot(data3, aes(x=cluster_id, y=identity, col=comparision_type_genus)) +geom_jitter(alpha=0.2) + 
  
  facet_wrap(unique_taxa_names_per_ST~., scales = 'free_x', nrow=1, shrink=TRUE, dir='h') + scale_y_continuous(name="") +
  scale_x_continuous(breaks = seq(0, 5, 1), name ="Cluster ID") + theme(legend.position="none")

p5<-ggarrange(p1, p2,
              ncol = 1, nrow = 2, heights = c(1, 1)) 
plot(p5)
dev.off()



#Generate Genome coverage pyANI plots for ST comparisions
data <- read.csv("~/Desktop/Kiepas_et_al_2023_MLST/data/pyani_analysis/pyani_coevarge_identity_ST_comparisions_df.csv")
data2 <- data[!data$unique_taxa_names_per_ST %in% c('3', '4'), ]
data3 <- data[!data$unique_taxa_names_per_ST %in% c('1', '2'), ]

pdf(file='../../figures/pyani_analysis_connected_components/scatter_plots/coverage_STs.pdf', width=30, height=15)

p1<-ggplot(data2, aes(x=cluster_id, y=coverage, col=comparision_type_genus)) +geom_jitter(alpha=0.2)  + facet_wrap(unique_taxa_names_per_ST ~., scales = 'free_x', nrow=4, shrink=TRUE, dir='h') + scale_y_continuous(name="Genome Coverage (%)") +
  scale_x_continuous(breaks = seq(0, 100, 2), name ="") + theme(legend.position="none") 

p2<-ggplot(data3, aes(x=cluster_id, y=coverage, col=comparision_type_genus)) +geom_jitter(alpha=0.2) + 
  
  facet_wrap(unique_taxa_names_per_ST~., scales = 'free_x', nrow=1, shrink=TRUE, dir='h') + scale_y_continuous(name="") +
  scale_x_continuous(breaks = seq(0, 5, 1), name ="Cluster ID") + theme(legend.position="none")

p5<-ggarrange(p1, p2,
              ncol = 1, nrow = 2, heights = c(1, 1)) 
plot(p5)
dev.off()

