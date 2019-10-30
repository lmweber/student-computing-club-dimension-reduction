###################################################
# Student computing club - dimension reduction demo
# Lukas Weber, 2019-10-29
###################################################


##################
# Install packages
##################

# from CRAN
install.packages(c("Rtsne", "uwot", "ggplot2"))

# from Bioconductor
# see website: http://bioconductor.org/
install.packages("BiocManager")
BiocManager::install(c("SumarizedExperiment", "HDCytoData"))



###########
# Load data
###########

library(HDCytoData)
d_SE <- Levine_32dim_SE()

# show data
d_SE
dim(d_SE)
head(assay(d_SE))
rowData(d_SE)
colData(d_SE)



###############
# Preprocessing
###############

# extract data
d_sub <- assay(d_SE[, colData(d_SE)$marker_class == "type"])
population <- rowData(d_SE)$population_id
dim(d_sub)
stopifnot(nrow(d_sub) == length(population))

# transform data
cofactor <- 5
d_sub <- asinh(d_sub / cofactor)
dim(d_sub)
summary(d_sub)

# subset
n <- 2000
set.seed(123)
ix <- sample(seq_len(nrow(d_sub)), n)
d_sub <- d_sub[ix, ]
population <- population[ix]
dim(d_sub)
stopifnot(nrow(d_sub) == length(population))

# remove any near-duplicate rows (required by Rtsne)
dups <- duplicated(d_sub)
d_sub <- d_sub[!dups, ]
population <- population[!dups]
dim(d_sub)
stopifnot(nrow(d_sub) == length(population))



##########################
# Dimension reduction: PCA
##########################

n_dims <- 2

# run PCA
# (note: no scaling, since asinh-transformed dimensions are already comparable)
out_PCA <- prcomp(d_sub, center = TRUE, scale. = FALSE)
dims_PCA <- out_PCA$x[, seq_len(n_dims)]
colnames(dims_PCA) <- c("PC_1", "PC_2")
head(dims_PCA)

dim(dims_PCA)
stopifnot(nrow(dims_PCA) == length(population))

colnames(dims_PCA) <- c("dimension_x", "dimension_y")
dims_PCA <- cbind(as.data.frame(dims_PCA), population, type = "PCA")
head(dims_PCA)
str(dims_PCA)



###########################
# Dimension reduction: tSNE
###########################

library(Rtsne)

set.seed(123)
out_Rtsne <- Rtsne(as.matrix(d_sub), dims = n_dims)
dims_Rtsne <- out_Rtsne$Y
colnames(dims_Rtsne) <- c("tSNE_1", "tSNE_2")
head(dims_Rtsne)

dim(dims_Rtsne)
stopifnot(nrow(dims_Rtsne) == length(population))

colnames(dims_Rtsne) <- c("dimension_x", "dimension_y")
dims_Rtsne <- cbind(as.data.frame(dims_Rtsne), population, type = "tSNE")
head(dims_Rtsne)
str(dims_Rtsne)



###########################
# Dimension reduction: UMAP
###########################

library(uwot)

set.seed(123)
out_umap <- umap(d_sub)
dims_umap <- out_umap
colnames(dims_umap) <- c("UMAP_1", "UMAP_2")
head(dims_umap)

dim(dims_umap)
stopifnot(nrow(dims_umap) == length(population))

colnames(dims_umap) <- c("dimension_x", "dimension_y")
dims_umap <- cbind(as.data.frame(dims_umap), population, type = "UMAP")
head(dims_umap)
str(dims_umap)



################
# Generate plots
################

library(ggplot2)

dir_figures <- "plots"


stopifnot(nrow(dims_Rtsne) == nrow(dims_PCA))
stopifnot(nrow(dims_umap) == nrow(dims_PCA))

d_plot <- rbind(dims_PCA, dims_Rtsne, dims_umap)
str(d_plot)

d_plot$population <- as.factor(d_plot$population)
d_plot$type <- factor(d_plot$type, levels = c("PCA", "tSNE", "UMAP"))

colors <- c(rainbow(14), "gray75")


d_plot_PCA <- d_plot[d_plot$type == "PCA", ]
d_plot_PCA_tSNE <- d_plot[d_plot$type %in% c("PCA", "tSNE"), ]


ggplot(d_plot_PCA, aes(x = dimension_x, y = dimension_y, color = population)) + 
  facet_wrap(~ type, scales = "free") + 
  geom_point(size = 0.7, alpha = 0.5) + 
  scale_color_manual(values = colors) + 
  labs(x = "dimension x", y = "dimension y") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        legend.key.height = unit(4, "mm"))

filename <- "Levine_32dim_PCA.png"

ggsave(file.path(dir_figures, filename), height = 3, width = 6)


ggplot(d_plot_PCA_tSNE, aes(x = dimension_x, y = dimension_y, color = population)) + 
  facet_wrap(~ type, scales = "free") + 
  geom_point(size = 0.7, alpha = 0.5) + 
  scale_color_manual(values = colors) + 
  labs(x = "dimension x", y = "dimension y") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        legend.key.height = unit(4, "mm"))

filename <- "Levine_32dim_PCA_tSNE.png"

ggsave(file.path(dir_figures, filename), height = 3, width = 8.5)


ggplot(d_plot, aes(x = dimension_x, y = dimension_y, color = population)) + 
  facet_wrap(~ type, scales = "free") + 
  geom_point(size = 0.7, alpha = 0.5) + 
  scale_color_manual(values = colors) + 
  labs(x = "dimension x", y = "dimension y") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        legend.key.height = unit(4, "mm"))

filename <- "Levine_32dim_PCA_tSNE_UMAP.png"

ggsave(file.path(dir_figures, filename), height = 3, width = 10)



