######################################################################
## EDA for technical test for Junior Data Scientist
## Author: Antonio Ortega
## Date: 28-10-2018
## This script is compatible with a knitr Rnw file
## *save and *load chunks save and load variables
## stored as RData files, useful during compilation of Rnw files
######################################################################

## ---- load_libraries ----
library(ggplot2)
library(viridis)
library(dplyr)
library(magrittr)
library(waffle)
library(tidyr)
library(tibble)
library(ade4)
library(data.table)
library(stringr)
library(FactoMineR)
library(cowplot)
library(pheatmap)
library(kableExtra)
# library(MASS)
# library(scales)
plot_dir <- "plots"
output_data <- "proc_data"
rdata_dir <- "RData"
tables_dir <- "tables"

## ---- load_data ----
datos <- read.table("datos.csv", sep = ",", header=T, stringsAsFactors = F)
datos <- datos[, c(colnames(datos)[1:4], colnames(datos)[5:(ncol(datos)-1)] %>% sort,"Y")]
write(x = kable(x = datos %>% colnames %>% strsplit(., split = "_") %>%
                  lapply(., function(x) x[[1]]) %>%
                  unlist %>% table %>% t,
                format = "latex", digits = 2),
      file = file.path(tables_dir, "data_summary.tex")
      )


datos$edad_integer <- str_match(string = datos$edad,
                                pattern = "\\[\\d{2}-(\\d{2})\\)") %>%
  .[, 2] %>%
  as.integer

train_set <- datos[!is.na(datos["Y"]),]
x_train <- train_set %>% select(-Y)
y_train <- select(train_set, Y) %>% mutate(Y=as.factor(Y))

test_set <- datos[is.na(datos["Y"]),]
x_test <- test_set %>% select(-Y)

## ---- load_data_save ----
datasets <- list(x_train = x_train, y_train = y_train, x_test = x_test)
save("datasets", file = file.path("RData", "datasets.RData"))

## ---- load_data_load ----
load(file = "RData/datasets.RData")
x_train <- datasets$x_train
x_test <- datasets$x_test
y_train <- datasets$y_train


## ---- plots1 ----
race <- as.character(x_train$raza)
race <- table(race)
# Decrease the counts to make them compatible with waffle()
# Counts are not important, just the proportion
race <- 200*race/(max(race)) %>% sort
races <- names(race)
race <- race %>% as.numeric
names(race) <- races
race <- race %>% sort %>% rev
# Customize colors for more clear visualization
# politically incorrect :-)
race_colors <- c(
  "Caucasian"="orange",
  "AfricanAmerican"="black",
  "Hispanic"="olivedrab",
  "Other"="gray",
  "Asian"="yellow"
)[names(race)]

p1 <- waffle(race, rows = 25, flip = T, colors = race_colors)
ggsave(p1, filename = file.path(plot_dir, "race_waffle.png"))

p2 <- ggplot(data=x_train, mapping=aes(x = edad, y=..count../1e3))
p2 <- p2 + geom_bar() + labs(y = "kCounts", x="Age") + ggtitle("Age distribution")
ggsave(p2, filename = file.path(plot_dir, "age_histogram.png"))
p3 <- plot_grid(ncol=2, p1 + theme(legend.position = "top"), p2, labels = "AUTO")
ggsave(p3, filename = file.path(plot_dir, "visualize_categories.png"), width=14)

## ---- preprocess_function ----
preprocess_data <- function(x_train, x_test, etiquettes=FALSE) {
  
  # Every feature is preprocessed in groups according to their type (nominal, ordinal, etc) 
  # Train and test sets have separate pieces of code for "hygiene",
  # even if it implies code cluttering
  # The result for each group is a list with the processed data for the training and test set respectively,
  # in data.frame format
  
  if(etiquettes) {
    ## Preprocess etiquettes
    print("Processing etiquettes")
    train_etiquettes <- x_train[, c("etiqueta_1", "etiqueta_2", "etiqueta_3")]
    test_etiquettes <- x_test[, c("etiqueta_1", "etiqueta_2", "etiqueta_3")]
    rownames(train_etiquettes) <- x_train$identificador
    rownames(test_etiquettes) <- x_test$identificador
    
    # Extract a vector storing all the etiquetas occuring in the TRAINING SET ONLY
    unique_etiquettes <- train_etiquettes %>% unlist %>% unique %>% sort
    
    # Encode the etiquetas in format one-hot
      # Initialize of vector of counts of etiquetas to 0 to all of them
    etiquettes_template <- rep(0, length(unique_etiquettes))
    names(etiquettes_template) <- unique_etiquettes
    
    # For both train and test set
    etiquettes_proc <- lapply(
      list(train = train_etiquettes, test = test_etiquettes),
      function(y) {
        # For every row, apply the same function
        res <- y %>% as.matrix %>% apply(., 1, function(x) {
          # Count how many times each etiqueta appears
          # This returns a table with the counts of the etiquetas appearing in this individual
          # but not the ones not appearing (which obv should be 0)
          local_count <- table(x)
          # Set the counts of the template to 0
          et <- etiquettes_template
          # Drop any etiquette that's not in the training set. Makes sense when analyzing test set
          local_count <- local_count[names(local_count) %in% names(et)]
          # Set the counts of the found etiquetas to the right count
          et[names(local_count)] <- local_count 
          return(et)
          }) %>%
        t %>% unlist
        
        # Format the colnames of the etiquetas so they start with the et_ prefix
        colnames(res) <- paste0("et_", colnames(res))
        # Make the data structure a data.frame of factors
        res <- res %>% apply(.,2,as.factor) %>% as.data.frame
        return(res)
        })
  }
  
  ## Preprocess nominals
  print("Processing nominals")
  # Drop nominal 4
  # Make one-hot encoding using acm.disjonctif
  train_nominals_one_hot <- x_train %>% select(nominal_1:nominal_3) %>%
    acm.disjonctif() %>% apply(.,2,as.factor) %>% as.data.frame
  test_nominals_one_hot <- x_test %>% select(nominal_1:nominal_3) %>%
    acm.disjonctif() %>% apply(.,2,as.factor) %>% as.data.frame
  
  # Drop the nominal categories not present in the training set
  missing_nominal_test <- colnames(train_nominals_one_hot)[!colnames(train_nominals_one_hot) %in% colnames(test_nominals_one_hot)]
  missing_test <- as.data.frame(matrix(0, nrow=nrow(x_test), ncol=length(missing_nominal_test)))
  colnames(missing_test) <- missing_nominal_test
  test_nominals_one_hot <- cbind(test_nominals_one_hot, missing_test)[colnames(train_nominals_one_hot)]
  
  # Make the list
  nominals_proc <- list(train = train_nominals_one_hot, test = test_nominals_one_hot)
  
  ## Preprocess drugs (farmacos)
  print("Processing drugs")
  
  # Select the farmaco features
  train_drugs <- x_train[, grep(pattern = "farmaco", x = colnames(x_train))]
  test_drugs <- x_test[, grep(pattern = "farmaco", x = colnames(x_test))]
  
  # Drop drugs with no variability (non-informative) 
  train_drugs <- train_drugs[,train_drugs %>% apply(.,2,function(x) length(unique(x))) > 1]
  test_drugs <- test_drugs[, colnames(train_drugs)]
  
  # Replace strings with integer and make list
  drugs_proc <- list(train = train_drugs, test = test_drugs) %>% lapply(function(x) {
    x[x=="No"] <- "-1"
    x[x=="Down"] <- "0"
    x[x=="Steady"] <- "1"
    x[x=="Up"] <- "2"
    x <- x %>% apply(.,2, as.integer) %>% as.data.frame
  })
  

  ## Preprocess ordinal
  print("Processing ordinals")
  
  # Select the features
  train_ordinal <- x_train %>% select(ordinal_1:ordinal_2)
  test_ordinal <- x_test %>% select(ordinal_1:ordinal_2)
  
  # Replace strings with integer and make list
  ordinals_proc <- lapply(list(train = train_ordinal, test = test_ordinal), function(x) {
    x[x=="None"] <- "0"
    x[x=="Norm"] <- "1"
    x[x==">7" | x==">200"] <- "2"
    x[x==">8" | x==">300"] <- "3"
    x <- x %>% apply(.,2,as.integer) %>% as.data.frame
    x
  })
  
  ## Preprocess binary
  print("Processing binaries")
  
  # Select the binary features and make them factors of 0 and 1.
  train_binary <- x_train %>% select(binary_1:binary_3) %>%
    apply(.,2, function(x) as.factor(as.integer(as.factor(x)) - 1)) %>% as.data.frame
  test_binary <- x_test %>% select(binary_1:binary_3) %>%
    apply(.,2, function(x) as.factor(as.integer(as.factor(x)) - 1)) %>% as.data.frame
  
  # Make the list
  binary_proc <- list(train = train_binary, test = test_binary)
  
  
  ## Preprocess counter
  print("Processing counters")
  
  # Just make the list (no need to modify them :))
  counter_proc <- list(
    train = x_train %>% select(counter_1:counter_7),
    test = x_test %>% select(counter_1:counter_7)
    )

  ## Preprocess race
  print("Processing race")
  
  # Make the list while one-hot encoding
  # the same way as nominals
  race_proc <- list(
  train = x_train %>% select(raza) %>% acm.disjonctif() %>%
    apply(.,2,as.factor) %>% as.data.frame,
  test = x_test %>% select(raza) %>% acm.disjonctif() %>%
    apply(.,2,as.factor) %>% as.data.frame
  )
  
  ## Preprocess sex
  print("Processing sex")
  
  # Select the sex feature and process them like binaries
  sex_proc <- list(
    train = x_train %>% select(sexo) %>%
      apply(.,2,function(x) as.factor(as.integer(as.factor(x)) - 1)) %>%
      as.data.frame,
    test = x_test %>% select(sexo) %>%
      apply(.,2,function(x) as.factor(as.integer(as.factor(x)) - 1)) %>%
      as.data.frame
  )
  
  ## Preprocess age
  print("Processing age")
  
  # Select the age feature in integer format
  train_age <- x_train$edad_integer
  test_age <- x_test$edad_integer
  
  # Compute the mean age in the TRAINING SET ONLY
  train_age_mean <- mean(train_age, na.rm = T)
  
  # Impute missing ages using this mean on both sets
  train_age[is.na(train_age)] <- train_age_mean
  test_age[is.na(test_age)] <- train_age_mean
  
  # Make the list
  age_proc <- list(train = data.frame(edad_integer = train_age),
                   test = data.frame(edad_integer = test_age))
  
  ## CBIND all the features for each dataset returning a single data.frame
  ## Optionally include etiquetas (huge amount of one-hot columns)
  ## Make a list for each dataset
  datasets <- c("train", "test")
  if(etiquettes) {
  processed_datasets <- lapply(datasets, function(x) {
    cbind(race_proc[[x]], sex_proc[[x]], age_proc[[x]], nominals_proc[[x]], counter_proc[[x]],
          drugs_proc[[x]], ordinals_proc[[x]], binary_proc[[x]], etiquettes_proc[[x]])
    })
  } else {
    processed_datasets <- lapply(datasets, function(x) {
      cbind(race_proc[[x]], sex_proc[[x]], age_proc[[x]], nominals_proc[[x]], counter_proc[[x]],
            drugs_proc[[x]], ordinals_proc[[x]], binary_proc[[x]])
    })
  }
  
  # Name the list and return it
  names(processed_datasets) <- datasets
  return(processed_datasets)
}

## Preprocess using the function above without and with etiquetas
## Store the indes of quantitative and qualitative features
## ---- preprocess ----
processed_datasets <- preprocess_data(x_train, x_test, etiquettes = FALSE)
proc_x_train <- processed_datasets$train
proc_x_test <- processed_datasets$test

quanti <- c("farmaco", "ordinal", "edad_integer", "counter")
quali <- c("raza", "binary", "nominal", "sexo")
quanti_index <- lapply(quanti, function(x) grep(x = colnames(proc_x_train), pattern = x)) %>% unlist
quali_index <- lapply(quali, function(x) grep(x = colnames(proc_x_train), pattern = x)) %>% unlist

## ---- preprocess_etiqueta ----
processed_datasets <- preprocess_data(x_train, x_test, etiquettes = TRUE)
quanti_index_full <- lapply(quanti, function(x) grep(x = colnames(processed_datasets$train), pattern = x)) %>% unlist
quali_index_full <- lapply(c(quali, "et"), function(x) grep(x = colnames(processed_datasets$train), pattern = x)) %>% unlist

## Save the processed datasets with etiquetas to csv files here
write.table(x = processed_datasets$train, file = file.path(output_data, "x_train.csv"),
            sep = ",", row.names = x_train$identificador, quote = F, col.names = T)
write.table(x = processed_datasets$test, file = file.path(output_data, "x_test.csv"),
            sep = ",", row.names = x_test$identificador, quote = F, col.names = T)

## ---- preprocess_save ----
# rm(processed_datasets)
# Save all objects to an RData file
save(list = c("quanti", "quali",
              "quanti_index", "quali_index",
              "quanti_index_full", "quali_index_full",
              "proc_x_train", "proc_x_test", "processed_datasets"),
     file = file.path(rdata_dir, "preprocessed.RData")) 

## ---- preprocess_load ----
# Load them
load(file = file.path(rdata_dir, "preprocessed.RData"))

## ---- PCA ----
# Compute a PCA of the quantitative features in the training set
# with centering and scaling to give equal importance to all features
res.pca <- prcomp(x = proc_x_train[, quanti_index], center = T, scale. = T)
# Store in data.frame data for variance captured by each PC
pca_barplot_data <- data.frame(PC = 1:length(res.pca$sdev),
                               var = res.pca$sdev**2)

# Store in data.frame the PCS, the label and the etiquetas (which are not present in the proc_x_train)
pca_data <- cbind(x_train[, grep(pattern = "et", x = colnames(x_train), invert = T)],
                  as.data.frame(res.pca$x), y_train)

## ---- PCA_visualization ----

# Create barplots showing the variance captured by each PC
# and the cumulative variance
p0 <- plot_grid(ncol = 2, rel_widths = c(1,1),
                ggplot(data = pca_barplot_data, aes(x=PC, y=var)) +
                  geom_bar(stat="identity") + labs(y = "Variance"),
                ggplot(data = pca_barplot_data, aes(x=PC, y=cumsum(var)/sum(var))) +
                  geom_bar(stat="identity") + labs(y = "Fraction cumulative variance")
)

# Create PCA plots facetting and coloring with different categories/label
p1 <- ggplot(
  data = pca_data,
  mapping = aes(x = PC1, y = PC2, col = raza)
  ) + geom_point(size=.1) +
  guides(col=F) + ggtitle("Raza")

p2 <- ggplot(
  data = pca_data,
  mapping = aes(x = PC1, y = PC2, col = sexo)
  ) + geom_point(size=.1) +
  guides(col=F) + ggtitle("Sexo")

p3 <- ggplot(
  data = pca_data,
  mapping = aes(x = PC1, y = PC2, col = edad_integer)
  ) + geom_point(size=.1) +
  guides(col=F) + ggtitle("Edad")

p4 <- ggplot(
  data = pca_data,
  mapping = aes(x = PC1, y = PC2, col = Y)
  ) + geom_point(size=.1) +
  guides(col=F) + ggtitle("Y")

p <- plot_grid(ncol=4, p1,p2,p3,p4, labels="AUTO")

ggsave(
  filename = file.path(plot_dir, "PCA_multicategory.png"),
  plot = p,
  width=20
  )


## ---- MCA_functions ----
# Utility functions to prepare data for MCA and handle its output
prepare_projection_data <- function(x, y, proc_x, quali_index, quanti_index) {
  res.mca <- MCA(
    X = proc_x, ncp = ncp,
    quanti.sup = quanti_index, graph = F
  )
  cats <- apply(proc_x[,quali_index], 2, function(x) nlevels(as.factor(x)))
  mca_df <- extract_mca_data(proc_x, x, cats, res.mca)
  mca_df$obs$Y <- y$Y
  sup_proj_data <- cbind(proc_x[, quanti_index],
                         mca_df$obs %>% select(`Dim 1`:`Dim 5`))
  
  colnames(sup_proj_data)[colnames(sup_proj_data) %>% grep(pattern = "Dim")] <- paste0("MCA_", 1:ncp)
  return(list(mca = res.mca, vars = mca_df$vars, obs = mca_df$obs, proj=sup_proj_data))
}


extract_mca_data <- function(proc_data, x, cats, mca) {
  mca_vars_df <- data.frame(mca$var$coord, Variable = rep(names(cats), cats))
  
  mca_vars_df$Type <- mca_vars_df$Variable %>% as.character %>%
    gsub(pattern = "\\.", replacement = "_", x = .) %>%
    strsplit(x = ., split = "_") %>%
    lapply(function(x) x[[1]]) %>% unlist
  
  mca_obs_df <- cbind(mca$ind$coord,
                      edad_integer = x[,"edad_integer"],
                      x[, lapply(quali, function(y) grep(x = colnames(x), pattern = y)) %>% unlist]
  )
  return(list(vars = mca_vars_df, obs = mca_obs_df))
}

make_mca_plots <- function(train, test, point_size, prefix="") {
  p1 <- ggplot(data=train$vars, 
               aes(x = Dim.1, y = Dim.2, col = Type))
  p1 <- p1 + geom_hline(yintercept = 0, colour = "gray70")
  p1 <- p1 + geom_vline(xintercept = 0, colour = "gray70")
  p1 <- p1 + geom_point(size=2)
  p1 <- p1 + ggtitle("MCA plot of variables using R package FactoMineR")
  p1 <- p1 + scale_color_viridis(discrete = T)
  
  p2 <- ggplot(data = train$obs, aes(x = `Dim 1`, y = `Dim 2`, col = raza))
  p2 <- p2 + geom_point(size=point_size) + facet_wrap(~raza)
  p2 <- p2 + ggtitle("MCA plot facetted by race") 
  
  p3 <- ggplot(data = train$obs, aes(x = `Dim 1`, y = `Dim 2`))
  p3 <- p3 + geom_point(size=point_size) + facet_wrap(~Y)
  p3 <- p3 + ggtitle("MCA plot facetted by Y") 
  
  p4 <- ggplot(data = test$obs, aes(x = `Dim 1`, y = `Dim 2`))
  p4 <- p4 + geom_point(size=point_size) + ggtitle("MCA plot of the test set") 
  
  p5 <- ggplot(data = train$obs, aes(x = `Dim 1`, y = `Dim 2`, col = raza))
  p5 <- p5 + geom_point(size=point_size) + facet_wrap(~Y)
  p5 <- p5 + ggtitle("MCA plot facetted by Y") 
  
  p6 <- ggplot(data = train$obs, aes(x = `Dim 1`, y = `Dim 2`, col = as.factor(edad_integer)))
  p6 <- p6 + geom_point(size=point_size) + facet_wrap(~Y)
  p6 <- p6 + ggtitle("MCA plot facetted by Y")
  p6 <- p6 + scale_color_discrete() + guides(col = guide_legend(title = "Edad"))
  
  p7 <- plot_grid(p3, p4, ncol=2)
  
  # Contribution plot
  # Extract the contributions frm the mca object
  mca_contribution <- train$mca$var$contrib[,1]
  # Compute the variable type for nice coloring in barplot
  type <- strsplit(
    gsub(x = names(mca_contribution), pattern = "\\.", replacement = "_"),
    split = "_") %>%
    lapply(., function(x) x[[1]]) %>% unlist
  
  # Make dataframe
  mca_contribution_df <- data.frame(variable = names(mca_contribution),
                                    type = type,
                                    contrib = mca_contribution,
                                    stringsAsFactors = F) %>%
    arrange(-contrib)
  # Sort features by their contribution
  mca_contribution_df$variable <- factor(mca_contribution_df$variable, levels = mca_contribution_df$variable)
  
  p8 <- ggplot(mca_contribution_df  %>% head(30),
         aes(x = variable, y = contrib, fill = type)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_viridis(discrete = T)
  
  ggsave(width=7, height=7, plot = p1, filename = file.path(plot_dir, paste0(prefix, "mca_variables.png")))
  ggsave(width=7, height=7, plot = p2, filename = file.path(plot_dir, paste0(prefix, "mca_obs_race_facet.png")))
  ggsave(width=7, height=7, plot = p3, filename = file.path(plot_dir, paste0(prefix, "mca_obs_Y_facet.png")))
  ggsave(width=7, height=7, plot = p4, filename = file.path(plot_dir, paste0(prefix, "mca_obs_test.png")))
  ggsave(width=7, height=7, plot = p5, filename = file.path(plot_dir, paste0(prefix, "mca_obs_Y_facet_raza_col.png")))
  ggsave(width=7, height=7, plot = p6, filename = file.path(plot_dir, paste0(prefix, "mca_obs_Y_facet_edad_col.png")))
  ggsave(width=10, height=7, plot = p7, filename = file.path(plot_dir, paste0(prefix, "mca_obs_Y_facet_all.png")))
  ggsave(width=14, height=7, plot = p8, filename = file.path(plot_dir, paste0(prefix, "mca_contrib.png")))
  return(list(p1,p2,p3,p4,p5,p6,p7))
}

## ---- MCA ----

# Perform MCA as implemented in FactoMineR with five dimensions as output
# both with and without etiquetas
# DANGER! THIS CHUNK CAN TAKE SEVERAL MINUTES TO RUN,
# SPECIALLY THE PART WITH ETIQUETAS (many more features)
ncp <- 5

# Without etiquetas
train_sup_proj_data <- prepare_projection_data(x_train, y_train, proc_x_train, quali_index, quanti_index)
unique_cols <- apply(proc_x_test, 2, function(x) length(unique(x)) == 1) %>% which
quanti_index_mca <- lapply(quanti, function(x) grep(x = colnames(proc_x_test[,-unique_cols]), pattern = x)) %>% unlist
quali_index_mca <- lapply(quali, function(x) grep(x = colnames(proc_x_test[,-unique_cols]), pattern = x)) %>% unlist
test_sup_proj_data <- prepare_projection_data(x_test, NULL, proc_x_test[, -unique_cols], quali_index_mca, quanti_index_mca)
  

# With etiquetas
train_sup_proj_data_full <- prepare_projection_data(x_train, y_train, processed_datasets$train, quali_index_full, quanti_index_full)
unique_cols <- apply(processed_datasets$test, 2, function(x) length(unique(x)) == 1) %>% which
quanti_index_full_mca <- lapply(quanti, function(x) grep(x = colnames(processed_datasets$test[,-unique_cols]), pattern = x)) %>% unlist
quali_index_full_mca <- lapply(c("et", quali), function(x) grep(x = colnames(processed_datasets$test[,-unique_cols]), pattern = x)) %>% unlist
test_sup_proj_data_full <- prepare_projection_data(x_test, NULL, processed_datasets$test[, -unique_cols], quali_index_full_mca, quanti_index_full_mca)


## ---- MCA_save ----
save(list = c("train_sup_proj_data", "test_sup_proj_data",
              "train_sup_proj_data_full", "test_sup_proj_data_full"),
     file = file.path(rdata_dir, "mca_dfs.RData"))

## ---- MCA_load ----
load(file.path(rdata_dir, "mca_dfs.RData"))


## ---- MCA_visualization ----

# Visualize the MC results and the contribution of each variable
# to the found dimensions
point_size <- .5

# MCA results plots
plots_without <- make_mca_plots(train_sup_proj_data, test_sup_proj_data, point_size)
plots_with <- make_mca_plots(train_sup_proj_data_full, test_sup_proj_data_full, point_size, prefix = "etiqueta_")
# Combine some plots for easier integration in latex
ggsave(filename = file.path(plot_dir, "mca_variables_combined.png"),
       plot = plot_grid(plots_without[[1]], plots_with[[1]], labels="AUTO"),
       width = 14
       )
ggsave(filename = file.path(plot_dir, "mca_obs_Y_facet_all_combined.png"),
       plot = plot_grid(plots_without[[7]], plots_with[[7]], nrow = 2, labels="AUTO"),
       width = 14, height = 14
)

## ---- supervised_projections ----
# lda_res <- lda(formula = Y ~ ., data = cbind(train_sup_proj_data_full$proj, Y = y_train$Y))
# prop.lda <- lda_res$svd^2/sum(lda_res$svd^2)
# 
# # extra <- train_sup_proj_data_full$proj[, !(colnames(train_sup_proj_data_full$proj) %in% colnames(test_sup_proj_data_full$proj))]
# # extra <- extra %>% apply(.,2,function(x) rep(0, length(x)))
# # 
# # test_data_lda <- cbind(test_sup_proj_data_full$proj, extra)[colnames(train_sup_proj_data_full$proj)]
# 
# plda <- predict(object = lda_res,
#                 newdata = train_sup_proj_data_full$proj)
# 
# dataset <- data.frame(Y = y_train$Y, lda = plda$x)
# ggplot(dataset) + geom_point(aes(lda.LD1, lda.LD2, colour = Y, shape = Y), size = 2.5) + 
#   labs(x = paste("LD1 (", percent(prop.lda[1]), ")", sep=""),
#        y = paste("LD2 (", percent(prop.lda[2]), ")", sep=""))

## lda <- plotRLDF(
##   t(sup_proj_data),
##   labels.y = y_train$Y,
##   trend = TRUE, robust = TRUE
##   )
## str(lda)




## Perform PLSDA
# plsda <- DiscriMiner::plsDA(variables = scale(sup_proj_data),
#                             group = y_train$Y, 
#                             autosel = FALSE, comps = 2)
# 
# # Lots of output, we are interested in the components
# summary(plsda)
# 
# qplot(data=as.data.frame(plsda$components), x=t1, y=t2, geom=c("point"), color=y_train$Y)

## ---- heatmap ----
# Select only some individuals because pheatmap cannot handle to many individuals
idx <- sample(x = 1:nrow(train_sup_proj_data$proj),size = 1e3)
heatmap_data <- train_sup_proj_data$proj[idx, ] %>% as.matrix
heatmap <- pheatmap(heatmap_data,
                    annotation_row  = cbind(select(x_train[idx, ], raza, sexo, edad_integer),
                                            Y = y_train[idx,"Y"]),
                    scale="row", filename = file.path(plot_dir, "heatmap.png"))

# Export counts of categories in counter variables to a table to be integrated in the latex report
write(x = kable(x = train_sup_proj_data$proj[,train_sup_proj_data$proj %>% colnames %>%
                                               grep(pattern = "counter", x = .)] %>%
                  lapply(., function(x) length(table(x))) %>% unlist %>% sort %>% rev),
      file = file.path(tables_dir, "counters_overview.tex"))


## ---- heatmap_save ----
save(list = c("heatmap"),
     file = file.path(rdata_dir, "heatmap.RData"))

## ---- heatmap_load ----
load(file.path(rdata_dir, "heatmap.RData"))

## ---- export_data ----
# Save processed datasets without one-hot encoding data but with mca features
write.table(x = train_sup_proj_data$proj, file = file.path(output_data, "x_train_mca.csv"),
            sep = ",", row.names = x_train$identificador, quote = F, col.names = T)
write.table(x = test_sup_proj_data$proj, file = file.path(output_data, "x_test_mca.csv"),
            sep = ",", row.names = x_test$identificador, quote = F, col.names = T)

# Save processed datasets with one-hot encoding data (full) and the mca features found with them too
write.table(x = train_sup_proj_data_full$proj, file = file.path(output_data, "x_train_mca_full.csv"),
            sep = ",", row.names = x_train$identificador, quote = F, col.names = T)
write.table(x = test_sup_proj_data_full$proj, file = file.path(output_data, "x_test_mca_full.csv"),
            sep = ",", row.names = x_test$identificador, quote = F, col.names = T)

# Save label to a separate file
write.table(x = y_train, file = file.path(output_data, "y_train.csv"),
            sep = ",", row.names = x_train$identificador, quote = F, col.names = T)

## ---- session_info ----
sessionInfo()
