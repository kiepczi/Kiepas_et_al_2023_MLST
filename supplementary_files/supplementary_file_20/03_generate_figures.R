# This R script was used to generate graphs for MLST analysis of Streptomyces
#Set Up
library('ggplot2')
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggplot2)
library('ggExtra')


#========Generate Suplementary Figure 28. Genome count per ST bs Unique names per ST=====
data <- read.csv("ST_ancillary_info.csv") #Loading data

print(data)
pdf(file='../supplementary_file_28.pdf', width=10, height=5)

p <- ggplot(data, aes(x=Genome_Count, y=Unique_Species_Names_Count)) + 
  scale_x_continuous(breaks = seq(0, 30, 1), name ="Genome count per ST") +
  scale_y_continuous(breaks = seq(1, 4, 1), name="Unique NCBI species names per ST") + 
  geom_count(aes(colour=after_stat(n)), alpha=0.4)  + 
  scale_size_area(max_size=20) + 
  theme(legend.position = "none") +
  scale_color_viridis_c() +
  theme(legend.position = "none") + 
  coord_cartesian(clip = "off")
p + geom_text(data = ggplot_build(p)$data[[1]], 
              aes(x, y, label = n), color = "black")
dev.off()



#======Generate Genome identity pyANI plots for ST comparisions =====
data <- read.csv("pyani_coverage_identity_ST_comparisions_df.csv")
data2 <- data[!data$unique_taxa_names_per_ST %in% c('3', '2'), ]
data3 <- data[!data$unique_taxa_names_per_ST %in% c('1'), ]

pdf(file='main_text_figures/identity_STs.pdf', width=20, height=10)

p1<-ggplot(data2, aes(x=cluster_id, y=identity, col=comparision_type_genus)) +geom_jitter(alpha=0.4)  + facet_wrap(unique_taxa_names_per_ST ~., scales = 'free_x', nrow=4, shrink=TRUE, dir='h') + scale_y_continuous(name="Genome Identity (%)") +
  scale_x_continuous(breaks = seq(0, 100, 2), name ="") + theme(legend.position="none") +
  scale_color_manual(values=c( "#76B947", "#FF0000")) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  theme(legend.position="none",
        axis.title.x = element_text(size = 18),  # Adjust x-axis label size
        axis.title.y = element_text(size = 18),  # Adjust y-axis label size
        axis.text.x = element_text(size = 18),   # Adjust x-axis tick label size
        axis.text.y = element_text(size = 18),
        strip.text = element_text(size = 18)) 

p2<-ggplot(data3, aes(x=cluster_id, y=identity, col=comparision_type_genus)) +geom_jitter(alpha=0.4) + 
  
  facet_wrap(unique_taxa_names_per_ST~., scales = 'free_x', nrow=1, shrink=TRUE, dir='h') + scale_y_continuous(name="") +
  scale_x_continuous(breaks = seq(0, 30, 1), name ="Cluster ID") + theme(legend.position="none") +
  scale_color_manual(values=c( "#76B947", "#FF0000"))+
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  theme(legend.position="none",
        axis.title.x = element_text(size = 18),  # Adjust x-axis label size
        axis.title.y = element_text(size = 18),  # Adjust y-axis label size
        axis.text.x = element_text(size = 18),   # Adjust x-axis tick label size
        axis.text.y = element_text(size = 18),
        strip.text = element_text(size = 18)) 

p5<-ggarrange(p1, p2,
              ncol = 1, nrow = 2, heights = c(1, 1)) 
plot(p5)
dev.off()



#=====Generate Genome coverage pyANI plots for ST comparisions=====
data <- read.csv("pyani_coverage_identity_ST_comparisions_df.csv")
data2 <- data[!data$unique_taxa_names_per_ST %in% c('3', '2'), ]
data3 <- data[!data$unique_taxa_names_per_ST %in% c('1'), ]

pdf(file='main_text_figures/coverage_STs.pdf', width=20, height=10)

p1<-ggplot(data2, aes(x=cluster_id, y=coverage, col=comparision_type_genus)) +geom_jitter(alpha=0.4)  + facet_wrap(unique_taxa_names_per_ST ~., scales = 'free_x', nrow=4, shrink=TRUE, dir='h') + scale_y_continuous(name="Genome Coverage (%)") +
  scale_x_continuous(breaks = seq(0, 100, 2), name ="") + theme(legend.position="none") +
  scale_color_manual(values=c("#f01e2c", "#023E8A"),
                     labels=c("between-species", "within-species"),
                     name="Comparisons") +  # Add legend title here
  scale_y_continuous(name="Genome Coverage (%)")+
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  theme(legend.position="none",
        axis.title.x = element_text(size = 18),  # Adjust x-axis label size
        axis.title.y = element_text(size = 18),  # Adjust y-axis label size
        axis.text.x = element_text(size = 18),   # Adjust x-axis tick label size
        axis.text.y = element_text(size = 18),
        strip.text = element_text(size = 18)) 

p2<-ggplot(data3, aes(x=cluster_id, y=coverage, col=comparision_type_genus)) +geom_jitter(alpha=0.4) + 
  
  facet_wrap(unique_taxa_names_per_ST~., scales = 'free_x', nrow=1, shrink=TRUE, dir='h') + scale_y_continuous(name="") +
  scale_x_continuous(breaks = seq(0, 30, 1), name ="Cluster ID") + theme(legend.position="none") +
  scale_color_manual(values=c("#f01e2c", "#023E8A"),
                     labels=c("between-species", "within-species"),
                     name="Comparisons") +  # Add legend title here
  scale_y_continuous(name="Genome Coverage (%)")+
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  theme(legend.position="none",
        axis.title.x = element_text(size = 18),  # Adjust x-axis label size
        axis.title.y = element_text(size = 18),  # Adjust y-axis label size
        axis.text.x = element_text(size = 18),   # Adjust x-axis tick label size
        axis.text.y = element_text(size = 18),
        strip.text = element_text(size = 18)) 

p5<-ggarrange(p1, p2,
              ncol = 1, nrow = 2, heights = c(1, 1)) 
plot(p5)
dev.off()

#=====Generate Genome coverage pyANI plots for connected components=====
data <- read.csv("pyani_coverage_identity_df.csv")

data2 <- data[!data$unique_taxa_names_per_connected_component %in% c('13', '11','7', '5', '4', '3', '6'), ]
data3 <- data[!data$unique_taxa_names_per_connected_component %in% c('1', '2'), ]
pdf(file='../supplementary_file_29.pdf', width=20, height=10)


p1<-ggplot(data2, aes(x=cluster_id, y=coverage, col=comparision_type_species)) +geom_jitter(alpha=0.2)  + 
  facet_wrap(unique_taxa_names_per_connected_component ~., scales = 'free_x', nrow=2, shrink=TRUE, dir='h') + 
  scale_y_continuous(name="Genome Coverage (%)") +
  scale_x_continuous(breaks = seq(0, 100, 5), name ="") +
  scale_color_manual(values=c("#f01e2c", "#023E8A"),
                     labels=c("between-species", "within-species"),
                     name="Comparisons") +  # Add legend title here
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  theme(legend.position="none",
        axis.title.x = element_text(size = 18),  # Adjust x-axis label size
        axis.title.y = element_text(size = 18),  # Adjust y-axis label size
        axis.text.x = element_text(size = 18),   # Adjust x-axis tick label size
        axis.text.y = element_text(size = 18),
        strip.text = element_text(size = 18)) 

p2<-ggplot(data3, aes(x=cluster_id, y=coverage, col=comparision_type_species)) +geom_jitter(alpha=0.2) + 
  facet_wrap(unique_taxa_names_per_connected_component ~., scales = 'free_x', nrow=2, shrink=TRUE, dir='h') + 
  scale_y_continuous(name="") +
  scale_color_manual(values=c("#f01e2c", "#023E8A"),
                     labels=c("between-species", "within-species"),
                     name="Comparisons") +  # Add legend title here
  scale_x_continuous(breaks = seq(0, 10, 1), name ="Cluster ID") + theme(legend.position="none")+
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  theme(legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20),  # Adjust legend title size
        legend.position = c(1, 0.4), 
        legend.justification = c(1, 1.4), 
        legend.box.just = "right",
        legend.background = element_rect(fill = "white", color = "black"), 
        legend.key = element_rect(fill = "white"),
        axis.title.x = element_text(size = 18),  # Adjust x-axis label size
        axis.title.y = element_text(size = 18),  # Adjust y-axis label size
        axis.text.x = element_text(size = 18),   # Adjust x-axis tick label size
        axis.text.y = element_text(size = 18), strip.text = element_text(size = 18)) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) 

p4<-ggarrange(p1, p2, ncol = 1, nrow = 2, heights = c(1, 1, 1)) 
plot(p4)
dev.off()

#======Generate Genome identity pyANI plots for connected components========
data <- read.csv("pyani_coverage_identity_df.csv")

data2 <- data[!data$unique_taxa_names_per_connected_component %in% c('13', '11','7', '5', '4', '3', '6'), ]
data3 <- data[!data$unique_taxa_names_per_connected_component %in% c('1', '2'), ]
pdf(file='../supplementary_file_30.pdf', width=20, height=10)


p1<-ggplot(data2, aes(x=cluster_id, y=identity, col=comparision_type_genus)) +geom_jitter(alpha=0.2)  + 
  facet_wrap(unique_taxa_names_per_connected_component ~., scales = 'free_x', nrow=2, shrink=TRUE, dir='h') + 
  scale_y_continuous(name="Genome Identity (%)") +
  scale_x_continuous(breaks = seq(0, 50, 5), name ="") + theme(legend.position="none") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  scale_color_manual(values=c("purple", "orange"),
                     labels=c("between-genus", "within-genus"),
                     name="Comparisons") +  # Add legend title here
  theme(legend.position="none",
        axis.title.x = element_text(size = 18),  # Adjust x-axis label size
        axis.title.y = element_text(size = 18),  # Adjust y-axis label size
        axis.text.x = element_text(size = 18),   # Adjust x-axis tick label size
        axis.text.y = element_text(size = 18),
        strip.text = element_text(size = 18)) 

p2<-ggplot(data3, aes(x=cluster_id, y=identity, col=comparision_type_genus)) +geom_jitter(alpha=0.2) + 
  
  facet_wrap(unique_taxa_names_per_connected_component~., scales = 'free_x', nrow=2, shrink=TRUE, dir='h') + 
  scale_y_continuous(name="") +
  scale_x_continuous(breaks = seq(0, 10, 1), name ="Cluster ID") + theme(legend.position="none") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  scale_color_manual(values=c("purple", "orange"),
                     labels=c("between-genus", "within-genus"),
                     name="Comparisons") +  # Add legend title here
  theme(legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20),  # Adjust legend title size
        legend.position = c(1, 0.4), 
        legend.justification = c(1, 1.4), 
        legend.box.just = "right",
        legend.background = element_rect(fill = "white", color = "black"), 
        legend.key = element_rect(fill = "white"),
        axis.title.x = element_text(size = 18),  # Adjust x-axis label size
        axis.title.y = element_text(size = 18),  # Adjust y-axis label size
        axis.text.x = element_text(size = 18),   # Adjust x-axis tick label size
        axis.text.y = element_text(size = 18), strip.text = element_text(size = 18)) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) 

p4<-ggarrange(p1, p2,
              ncol = 1, nrow = 2, heights = c(1, 1, 1)) 
plot(p4)
dev.off()

#========Generate Genome coverage pyANI plots for shared name comparisions=======
data <- read.csv("pyani_coverage_identity_shared_NCBI_names_df.csv")
data2 <- data[!data$connected_components_count %in% c('2', '3', '4'), ]
data3 <- data[!data$connected_components_count %in% c('1'), ]


pdf(file='main_text_figures/identity_shared_name.pdf', width=20, height=10)
# Modify the ggplot code
p1 <- ggplot(data2, aes(x = organism, y = identity, col = comparision_type_genus)) +
  geom_jitter(alpha = 0.4) +
  facet_wrap(connected_components_count ~ ., scales = 'free_x', nrow = 4, shrink = TRUE, dir = 'h') +
  scale_y_continuous(name = "Genome Identity (%)",  limits = c(0.8, 1.0)) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values=c("#023E8A","#f01e2c"),
                     labels=c("between-genus", "within-genus"),
                     name="Comparisons") +  # Add legend title here
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  theme(legend.position="none",
        axis.title.x = element_text(size = 10),  # Adjust x-axis label size
        axis.title.y = element_text(size = 10),  # Adjust y-axis label size
        axis.text.x = element_text(size = 10),   # Adjust x-axis tick label size
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 18)) 

p2 <- ggplot(data3, aes(x = organism, y = identity, col = comparision_type_genus)) +
  geom_jitter(alpha = 0.4) +
  facet_wrap(connected_components_count ~ ., scales = 'free_x', nrow = 1, shrink = TRUE, dir = 'h') +
  scale_y_continuous(name = "Genome Identity (%)",  limits = c(0.8, 1.0)) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values=c("#f01e2c", "#023E8A"),
                     labels=c("between-genus", "within-genus"),
                     name="Comparisons") +  # Add legend title here
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  theme(legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20),  # Adjust legend title size
        legend.position = c(0.99, 0.1), 
        legend.justification = c(1, 1.4), 
        legend.box.just = "right",
        legend.background = element_rect(fill = "white", color = "black"), 
        legend.key = element_rect(fill = "white"),
        axis.title.x = element_text(size = 10),  # Adjust x-axis label size
        axis.title.y = element_text(size = 10),  # Adjust y-axis label size
        axis.text.x = element_text(size = 10),   # Adjust x-axis tick label size
        axis.text.y = element_text(size = 10), strip.text = element_text(size = 18)) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) 

p5 <- ggarrange(p1, p2, ncol = 1, nrow = 2, heights = c(1, 1))
plot(p5)
dev.off()



#========Generate Genome coverage pyANI plots for shared name comparisions======
data <- read.csv("pyani_coverage_identity_shared_NCBI_names_df.csv")
data2 <- data[!data$connected_components_count %in% c('2', '3', '4', '5'), ]
data3 <- data[!data$connected_components_count %in% c('1'), ]



pdf(file='main_text_figures/coverage_shared_name.pdf', width=20, height=10)

p1<-ggplot(data2, aes(x=organism, y=coverage, col=comparision_type_species)) +
  geom_jitter(alpha=0.4)  + 
  facet_wrap(connected_components_count ~., scales = 'free_x', nrow=4, shrink=TRUE, dir='h') + 
  scale_y_continuous(name="Genome Coverage (%)",  limits = c(0.1, 1.0)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  scale_color_manual(values=c("purple", "orange"),
                     labels=c("between-species", "within-species"),
                     name="Comparisons") +  # Add legend title here
  theme(legend.position="none",
        axis.title.x = element_text(size = 10),  # Adjust x-axis label size
        axis.title.y = element_text(size = 10),  # Adjust y-axis label size
        axis.text.x = element_text(size = 10),   # Adjust x-axis tick label size
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 18)) 

p2<-ggplot(data3, aes(x=organism, y=coverage, col=comparision_type_species)) +geom_jitter(alpha=0.4) + 
  
  facet_wrap(connected_components_count~., scales = 'free_x', nrow=1, shrink=TRUE, dir='h') + scale_y_continuous(name="") +
  theme(legend.position="none") +
  scale_y_continuous(name="Genome Coverage (%)",  limits = c(0.1, 1.0)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  scale_color_manual(values=c("purple", "orange"),
                     labels=c("between-species", "within-species"),
                     name="Comparisons") +  # Add legend title here
  theme(legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20),  # Adjust legend title size
        legend.position = c(0.99, 0.1), 
        legend.justification = c(1, 1.4), 
        legend.box.just = "right",
        legend.background = element_rect(fill = "white", color = "black"), 
        legend.key = element_rect(fill = "white"),
        axis.title.x = element_text(size = 10),  # Adjust x-axis label size
        axis.title.y = element_text(size = 10),  # Adjust y-axis label size
        axis.text.x = element_text(size = 10),   # Adjust x-axis tick label size
        axis.text.y = element_text(size = 10), strip.text = element_text(size = 18)) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) 


p5<-ggarrange(p1, p2,
              ncol = 1, nrow = 2, heights = c(1, 1)) 
plot(p5)
dev.off()



#==Generate Suplementary Figure 31. Number of genomes per pyANI species vs number of unique NCBU names per pyANI species
data <- read.csv("pyANI_species_count_info.csv") #Loading data

print(data)
pdf(file='../supplementary_file_31.pdf', width=15, height=8)

p <- ggplot(data, aes(x=Genome_Count, y=unique_ST_count)) + 
  scale_x_continuous(breaks = seq(1, 61, 2), name ="Number of genomes per ANIm species") +
  scale_y_continuous(breaks = seq(1, 46, 2), name="Number of unique STs per ANIm species") + 
  geom_count(aes(colour=after_stat(n)), alpha=0.4)  + 
  scale_size_area(max_size=20) + 
  theme(legend.position = "none") +
  scale_color_viridis_c() +
  theme(legend.position = "none") + 
  coord_cartesian(clip = "off") +
  theme(legend.position="none",
        axis.title.x = element_text(size = 12),  # Adjust x-axis label size
        axis.title.y = element_text(size = 12),  # Adjust y-axis label size
        axis.text.x = element_text(size = 12),   # Adjust x-axis tick label size
        axis.text.y = element_text(size = 12))
p + geom_text(data = ggplot_build(p)$data[[1]], 
              aes(x, y, label = n), color = "black")
dev.off()



#=====Generate Suplementary Figure 8. Number of genomes per pyANI species vs nnumber of unique NCBU names per pyANI species
data <- read.csv("pyANI_species_count_info.csv") #Loading data

print(data)
pdf(file='../supplementary_file_32.pdf', width=15, height=8)

p <- ggplot(data, aes(x=Genome_Count, y=Unique_names_Count)) + 
  scale_x_continuous(breaks = seq(0, 100, 2), name ="Number of genomes per ANIm species") +
  scale_y_continuous(name="Number of unique NCBI species names per ANIm species") + 
  geom_count(aes(colour=after_stat(n)), alpha=0.4)  + 
  scale_size_area(max_size=10) + 
  theme(legend.position = "none") +
  scale_color_viridis_c() +
  theme(legend.position = "none") + 
  coord_cartesian(clip = "off")+
  theme(legend.position="none",
        axis.title.x = element_text(size = 12),  # Adjust x-axis label size
        axis.title.y = element_text(size = 12),  # Adjust y-axis label size
        axis.text.x = element_text(size = 12),   # Adjust x-axis tick label size
        axis.text.y = element_text(size = 12))
p + geom_text(data = ggplot_build(p)$data[[1]], 
              aes(x, y, label = n), color = "black")
dev.off()

# Load necessary packages
library(ggplot2)
library(dplyr)
library(gridExtra) # Load gridExtra for arranging multiple plots

# Read and filter the dataframe
df <- read.csv("species_st_distibution.csv")

df <- df %>%
  filter(organism != "Streptomyces sp.")

# Identify unique organisms and split them into two groups
unique_organisms <- unique(df$organism)
midpoint <- ceiling(length(unique_organisms) / 2)

# Create two groups of organisms
organisms_group1 <- unique_organisms[1:midpoint]
organisms_group2 <- unique_organisms[(midpoint + 1):length(unique_organisms)]

# Split the dataframe into two subsets based on the unique organism groups
df1 <- df %>% filter(organism %in% organisms_group1)
df2 <- df %>% filter(organism %in% organisms_group2)

# Define a colorblind-friendly palette
colorblind_palette <- c("#E69F00", "#56B4E9", "#009E73")  # Orange, Sky Blue, and Bluish Green

# Generate the first scatter plot for the first subset of data
plot1 <- ggplot(df1, aes(x = organism, y = count, color = type)) +
  geom_point(size = 3, alpha = 0.6) +  # Increased alpha for better visibility
  labs(x = "Organism", y = "Count", color = "Type") +
  scale_color_manual(values = colorblind_palette) +  # Apply color palette
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45, 
      hjust = 1,
      size = 14,
      margin = margin(t = 10)  # Increase space above the x-axis text
    ),
    axis.title.x = element_text(size = 14),  # Adjust x-axis label size
    axis.title.y = element_text(size = 14),  # Adjust y-axis label size
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 14),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 100)  # Adjust the plot margins if needed
  ) +
  facet_grid(rows = vars(type), scales = "fixed")

# Generate the second scatter plot for the second subset of data
plot2 <- ggplot(df2, aes(x = organism, y = count, color = type)) +
  geom_point(size = 3, alpha = 0.6) +  # Increased alpha for better visibility
  labs(x = "Organism", y = "Count", color = "Type") +
  scale_color_manual(values = colorblind_palette) +  # Apply color palette
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45, 
      hjust = 1,
      size = 14,
      margin = margin(t = 10)  # Increase space above the x-axis text
    ),
    axis.title.x = element_text(size = 14),  # Adjust x-axis label size
    axis.title.y = element_text(size = 14),  # Adjust y-axis label size
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 14),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 100)  # Adjust the plot margins if needed
  ) +
  facet_grid(rows = vars(type), scales = "fixed")

# Create the PDF file and set the dimensions
pdf(file='../supplementary_file_35.pdf', width=23, height=15)  # Increased height to fit both plots on one page

# Arrange both plots on a single page using grid.arrange
grid.arrange(plot1, plot2, ncol = 1)  # Arrange in one column

# Close the PDF device
dev.off()

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(viridis)  # Load viridis package for color scales

# Step 1: Read the DataFrame from the CSV file
df <- read.csv("subgraph_sizes.csv")

# Check if the CSV was loaded correctly (optional)
print(head(df))

# Ensure that the dataframe has the correct column name. 
# Assuming that the column with the sizes is named "Size".
# If the column is named differently, replace "Size" with the correct column name.

# Step 2: Calculate the count of each subgraph size
df_count <- df %>%
  group_by(Size) %>%
  summarise(Count = n())

# Display the count DataFrame (optional)
print(df_count)

# Step 3: Open a PDF device to save the plot
pdf(file = "connected_subgraphs_sizes.pdf", width = 12, height = 8)  # Set the filename and dimensions

# Step 4: Generate the scatter plot with viridis color scale
ggplot(df_count, aes(x = Size, y = Count, color = Count)) +  # Add color aesthetic mapped to Count
  geom_point(size = 5, alpha=0.5) +  # Adjust point size if needed
  scale_color_viridis_c(option = "D") +  # Use viridis color scale
  coord_cartesian(clip = "off") + 
  scale_x_continuous(
    limits = c(1, 85),  # Set x-axis limits
    breaks = seq(1, 85, by = 10)  # Set x-axis breaks
  ) +
  scale_y_log10(  # Change y-axis to logarithmic scale
    limits = c(1, 200),  # Set y-axis limits
    breaks = c(1, 10, 100, 200)  # Set log-scale y-axis breaks
  ) +
  theme_minimal() +  # Use minimal theme
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 12),  # Adjust x-axis label size
    axis.title.y = element_text(size = 12),  # Adjust y-axis label size
    axis.text.x = element_text(size = 12),   # Adjust x-axis tick label size
    axis.text.y = element_text(size = 12)    # Adjust y-axis tick label size
  ) +
  labs(x = "Size of Connected Subgraph", y = "Count (log scale)")

# Step 5: Close the PDF device
dev.off()

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(viridis)

# Step 1: Read the DataFrame from the CSV file
df <- read.csv("ST_degrees.csv")

# Ensure that the dataframe has the correct column name. 
# Assuming that the column with the degrees is named "Degree".

# Step 2: Calculate the count of each degree
df_count <- df %>%
  group_by(Degree) %>%
  summarise(Count = n())

# Step 3: Open a PDF device to save the plot
pdf(file = "degree_count.pdf", width = 12, height = 8)  # Set the filename and dimensions

# Step 4: Generate the scatter plot with a custom color scale
ggplot(df_count, aes(x = Degree, y = Count, color = Count)) +  # Add color aesthetic mapped to Count
  geom_point(size = 5, alpha = 0.7) +  # Adjust point size if needed
  scale_color_viridis_c(option = "D") +  # Use viridis color scale
  coord_cartesian(clip = "off") + 
  scale_x_log10(  # Use log scale for x-axis
    limits = c(1, 15),  # Set x-axis limits
    breaks = c(1, 2, 3, 5, 10, 15)  # Set log-scale x-axis breaks
  ) +  
  scale_y_log10(  # Use log scale for y-axis
    limits = c(1, 350),  # Set y-axis limits
    breaks = c(1, 10, 100, 350)  # Set log-scale y-axis breaks
  ) +  
  theme_minimal() +  # Use minimal theme
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 12),  # Adjust x-axis label size
    axis.title.y = element_text(size = 12),  # Adjust y-axis label size
    axis.text.x = element_text(size = 12),   # Adjust x-axis tick label size
    axis.text.y = element_text(size = 12)    # Adjust y-axis tick label size
  ) +  
  labs(x = "STs connections (degree, log scale)", y = "Count (log scale)")

# Step 5: Close the PDF device
dev.off()
