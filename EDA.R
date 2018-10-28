## ---- load_libraries and settings
setwd("/home/antortjim/MEGA/Kaggle/Datrik")
library(ggplot2)
library(dplyr)
library(waffle)
library(tidyr)
library(tibble)
library(ade4)
library(data.table)
library(stringr)
library(FactoMineR)
library(limma)
library(DiscriMiner)
plot_dir <- "plots"
output_data <- "proc_data"

## ---- load_data
datos <- read.table("datos.csv", sep = ",", header=T, stringsAsFactors = F)
datos <- datos[, c(colnames(datos)[1:4], colnames(datos)[5:(ncol(datos)-1)] %>% sort,"Y")]
datos$edad_integer <- str_match(string = datos$edad,
                                pattern = "\\[\\d{2}-(\\d{2})\\)") %>%
  .[, 2] %>%
  as.integer

train_set <- datos[!is.na(datos["Y"]), ]
x_train <- train_set %>% select(-Y)
y_train <- select(train_set, Y) %>% mutate(Y=as.factor(Y))

test_set <- datos[is.na(datos["Y"]), ]
x_test <- test_set %>% select(-Y)
y_test <- test_set %>% select(Y) %>% mutate(Y=as.factor(Y))

## ---- plots1
x_train$raza %>% table
race <- as.character(x_train$raza)
race <- table(race)
race <- 200*race/(max(race)) %>% sort
races <- names(race)
race <- race %>% as.numeric
names(race) <- races
race <- race %>% sort %>% rev
race_colors <- c(
  "Caucasian"="orange",
  "AfricanAmerican"="black",
  "Hispanic"="olivedrab",
  "Other"="gray",
  "Asian"="yellow"
)[names(race)]

p <- waffle(race, rows = 20, colors = race_colors)
ggsave(p, filename = file.path(plot_dir, "race_waffle.png"))

p <- ggplot(data=x_train,
       mapping=aes(x = edad, y=..count../1000)
       ) +
  geom_bar() +
  labs(y = "kCounts", x="Age") +
  ggtitle("Age distribution")
  # scale_x_continuous(breaks = sort(unique(histogram_data$edad)))
p
ggsave(p, filename = file.path(plot_dir, "age_histogram.png"))

## ---- preprocess_function
preprocess_data <- function(x_train, x_test, etiquettes=FALSE) {
  
  if(etiquettes) {
    ## Preprocess etiquettes
    print("Processing etiquettes")
    train_etiquettes <- x_train[, c("etiqueta_1", "etiqueta_2", "etiqueta_3")]
    test_etiquettes <- x_test[, c("etiqueta_1", "etiqueta_2", "etiqueta_3")]
    rownames(train_etiquettes) <- x_train$identificador
    rownames(test_etiquettes) <- x_test$identificador
    unique_etiquettes <- train_etiquettes %>% unlist %>% unique %>% sort
    
    etiquettes_template <- rep(0, length(unique_etiquettes))
    names(etiquettes_template) <- unique_etiquettes
    etiquettes_proc <- lapply(
      list(train = train_etiquettes, test = test_etiquettes),
      function(y) {
        res <- y %>% as.matrix %>% apply(., 1, function(x) {
          local_count <- table(x)
          et <- etiquettes_template
          # drop any etiquette that's not in the training set
          local_count <- local_count[names(local_count) %in% names(et)]
          et[names(local_count)] <- local_count 
          return(et)
          }) %>%
        t %>% unlist
        colnames(res) <- paste0("et_", colnames(res))
        res <- res %>% apply(.,2,as.factor) %>% as.data.frame
        return(res)
        })
    
    # etiquettes_proc$train %>% View
    # OK
  }
  
  ## Preprocess nominals
  print("Processing nominals")
  # Drop nominal 4 for now
  ## Should remove the last column produced with each nominal
  train_nominals_one_hot <- x_train %>% select(nominal_1:nominal_3) %>% acm.disjonctif() %>%
    apply(.,2,as.factor) %>% as.data.frame
  test_nominals_one_hot <- x_test %>% select(nominal_1:nominal_3) %>% acm.disjonctif() %>%
    apply(.,2,as.factor) %>% as.data.frame
  
  missing_nominal_test <- colnames(train_nominals_one_hot)[!colnames(train_nominals_one_hot) %in% colnames(test_nominals_one_hot)]
  missing_test <- as.data.frame(matrix(0, nrow=nrow(x_test), ncol=length(missing_nominal_test)))
  colnames(missing_test) <- missing_nominal_test
  test_nominals_one_hot <- cbind(test_nominals_one_hot, missing_test)[colnames(train_nominals_one_hot)]
  
  nominals_proc <- list(train = train_nominals_one_hot, test = test_nominals_one_hot)
  
  nominals_proc$train %>% class
  # OK
  
  ## Preprocess drugs (farmacos)
  print("Processing drugs")
  
  train_drugs <- x_train[, grep(pattern = "farmaco", x = colnames(x_train))]
  test_drugs <- x_test[, grep(pattern = "farmaco", x = colnames(x_test))]
  
  # Drop drugs with no variability (non-informative) 
  train_drugs <- train_drugs[,train_drugs %>% apply(.,2,function(x) length(unique(x))) > 1]
  test_drugs <- test_drugs[, colnames(train_drugs)]
  
  # Replace strings with integer
  drugs_proc <- list(train = train_drugs, test = test_drugs) %>% lapply(function(x) {
    x[x=="No"] <- "-1"
    x[x=="Down"] <- "0"
    x[x=="Steady"] <- "1"
    x[x=="Up"] <- "2"
    x <- x %>% apply(.,2, as.integer) %>% as.data.frame
  })
  
  drugs_proc$train  %>% class
    
  ## Preprocess ordinal
  print("Processing ordinals")
  
  train_ordinal <- x_train %>% select(ordinal_1:ordinal_2)
  test_ordinal <- x_test %>% select(ordinal_1:ordinal_2)
  ordinals_proc <- lapply(list(train = train_ordinal, test = test_ordinal), function(x) {
    x[x=="None"] <- "0"
    x[x=="Norm"] <- "1"
    x[x==">7" | x==">200"] <- "2"
    x[x==">8" | x==">300"] <- "3"
    x <- x %>% apply(.,2,as.integer) %>% as.data.frame
    x
  })
  
  ordinals_proc$train %>% class
  
  
  ## Preprocess binary
  print("Processing binaries")
  
  train_binary <- x_train %>% select(binary_1:binary_3) %>% apply(.,2, function(x) as.factor(as.integer(as.factor(x)) - 1)) %>% as.data.frame
  test_binary <- x_test %>% select(binary_1:binary_3) %>% apply(.,2, function(x) as.factor(as.integer(as.factor(x)) - 1)) %>% as.data.frame
  binary_proc <- list(train = train_binary, test = test_binary)
  
  binary_proc$train %>% class
  
  ## Preprocess counter
  print("Processing counters")
  
  counter_proc <- list(
    train = x_train %>% select(counter_1:counter_7),
    test = x_test %>% select(counter_1:counter_7)
    )
  counter_proc$train %>% class
  
  ## Preprocess race
  print("Processing race")
  
  race_proc <- list(
  train = x_train %>% select(raza) %>% acm.disjonctif() %>%
    apply(.,2,as.factor) %>% as.data.frame,
  test = x_test %>% select(raza) %>% acm.disjonctif() %>%
    apply(.,2,as.factor) %>% as.data.frame
  )
  
  race_proc$train %>% class
  
  ## Preprocess sex
  print("Processing sex")
  
  sex_proc <- list(
    train = x_train %>% select(sexo) %>% apply(.,2,function(x) as.factor(as.integer(as.factor(x)) - 1)) %>% as.data.frame,
    test = x_test %>% select(sexo) %>% apply(.,2,function(x) as.factor(as.integer(as.factor(x)) - 1)) %>% as.data.frame
  )
  
  ## Preprocess age
  print("Processing age")
  
  train_age <- x_train$edad_integer
  train_age_mean <- mean(train_age, na.rm = T)
  train_age[is.na(train_age)] <- train_age_mean
  test_age <- x_test$edad_integer
  test_age[is.na(test_age)] <- train_age_mean
  
  age_proc <- list(train = data.frame(edad_integer = train_age), test = data.frame(edad_integer = test_age))
  
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
  
  names(processed_datasets) <- datasets
  return(processed_datasets)
}

## ---- preprocess
quantit <- c("farmaco", "ordinal", "edad_integer", "counter")
qualit <- c("raza", "binary", "nominal", "sexo")
# if(etiquettes_in) qualit <- c(qualit, "et")
processed_datasets <- preprocess_data(x_train, x_test, etiquettes = FALSE)
proc_x_train <- processed_datasets$train
proc_x_test <- processed_datasets$test
rm(processed_datasets)
quanti_index <- lapply(quantit, function(x) grep(x = colnames(proc_x_train), pattern = x)) %>% unlist
qualit_index <- lapply(qualit, function(x) grep(x = colnames(proc_x_train), pattern = x)) %>% unlist
processed_datasets <- preprocess_data(x_train, x_test, etiquettes = TRUE)

## ---- PCA
res.pca <- prcomp(x = proc_x_train[, quanti_index], center = T, scale. = T)
barplot(1:length(res.pca$sdev), cumsum(res.pca$sdev))


pca_data <- cbind(x_train[, grep(pattern = "et", x = colnames(x_train), invert = T)],
                  as.data.frame(res.pca$x), y_train)

## ---- PCA_visualization
ggplot(data = pca_data,
       mapping = aes(x = PC1, y = PC2, col = raza)) + geom_point(size=.1)

ggplot(data = pca_data,
       mapping = aes(x = PC1, y = PC2, col = sexo)) + geom_point(size=.1)

ggplot(data = pca_data,
       mapping = aes(x = PC1, y = PC2, col = binary_1)) + geom_point(size=.1)

ggplot(data = pca_data,
       mapping = aes(x = PC1, y = PC2, col = edad_integer)) + geom_point(size=.05)

ggplot(data = pca_data,
       mapping = aes(x = PC1, y = PC2, col = Y)) + geom_point(size=.05)


## ---- MCA_function
extract_mca_data <- function(proc_data, x, cats, mca) {
  mca_vars_df <- data.frame(mca$var$coord, Variable = rep(names(cats), cats))
  
  mca_vars_df$Type <- mca_vars_df$Variable %>% as.character %>%
    gsub(pattern = "\\.", replacement = "_", x = .) %>%
    strsplit(x = ., split = "_") %>%
    lapply(function(x) x[[1]]) %>% unlist
  
  mca_obs_df <- cbind(mca$ind$coord,
                      x[, lapply(qualit, function(y) grep(x = colnames(x), pattern = y)) %>% unlist]
  )
  return(list(vars = mca_vars_df, obs = mca_obs_df))
}

## ---- MCA
ncp <- 5
train_res.mca <- MCA(
  X = proc_x_train, ncp = ncp,
  quanti.sup = quanti_index
)

unique_cols <- apply(proc_x_test, 2, function(x) length(unique(x)) == 1) %>% which
quanti_index_2 <- lapply(quantit, function(x) grep(x = colnames(proc_x_train[,-unique_cols]), pattern = x)) %>% unlist
proc_x_test_mca <- proc_x_test[, -unique_cols]


test_res.mca <- MCA(
  X = proc_x_test_mca, ncp = ncp,
  quanti.sup = quanti_index_2
)


cats <- apply(proc_x_train[,qualit_index], 2, function(x) nlevels(as.factor(x)))
train_mca_df <- extract_mca_data(proc_x_train, x_train, cats, train_res.mca)
train_mca_df$obs$Y <- y_train$Y
cats <- apply(proc_x_test_mca[,qualit_index_2], 2, function(x) nlevels(as.factor(x)))
test_mca_df <- extract_mca_data(proc_x_test_mca, x_test, cats, test_res.mca)

## ---- MCA_visualization
p <- ggplot(data=train_mca_df$vars, 
       aes(x = Dim.1, y = Dim.2, col = Type)) +
  geom_hline(yintercept = 0, colour = "gray70") +
  geom_vline(xintercept = 0, colour = "gray70") +
  geom_point(size=1) +
  # geom_text(aes(label = rownames(mca1_vars_df))) + guides(col=F)
  ggtitle("MCA plot of variables using R package FactoMineR")

p
ggsave(plot = p, filename = file.path(plot_dir, "mca_variables.png"))


p <- ggplot(data = train_mca_df$obs, aes(x = `Dim 1`, y = `Dim 2`, col = raza)) +
  geom_point(size=.1) +
  facet_wrap(~raza)
p
ggsave(plot = p, filename = file.path(plot_dir, "mca_obs_race_facet.png"))

p <- ggplot(data = train_mca_df$obs, aes(x = `Dim 1`, y = `Dim 2`, col = Y)) +
  geom_point(size=.1) +
  facet_wrap(~Y)
p
ggsave(plot = p, filename = file.path(plot_dir, "mca_obs_Y_facet.png"))

p <- ggplot(data = test_mca_df$obs, aes(x = `Dim 1`, y = `Dim 2`)) +
  geom_point(size=.1)
p


## ---- MCA_features
train_sup_proj_data <- cbind(proc_x_train[, quanti_index],
                             train_mca_df$obs %>% select(`Dim 1`:`Dim 5`))
test_sup_proj_data <- cbind(proc_x_test[, quanti_index],
                            test_mca_df$obs %>% select(`Dim 1`:`Dim 5`))

## ---- supervised_projections
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

## ---- Heatmap
idx <- sample(x = 1:nrow(train_sup_proj_data),size = 1e3)
heatmap <- pheatmap::pheatmap(train_sup_proj_data[idx, ] %>% as.matrix,
                   annotation_row  = x_train[idx, ] %>% select(raza, sexo, edad_integer),
                   # annotation_col=y_train$Y[idx],
                   scale="row", filename = "heatmap.png")
save(file = "heatmap.rda", heatmap)

## ---- export_data
# Save processed datasets without one-hot encoding data but with mca features
write.table(x = train_sup_proj_data, file = file.path(output_data, "x_train_mca.csv"),
            sep = ",", row.names = x_train$identificador, quote = F, col.names = T)
write.table(x = test_sup_proj_data, file = file.path(output_data, "x_test_mca.csv"),
            sep = ",", row.names = x_test$identificador, quote = F, col.names = T)
# Save label to a separate file
write.table(x = y_train, file = file.path(output_data, "y_train.csv"),
            sep = ",", row.names = x_train$identificador, quote = F, col.names = T)

# Save processed datasets containing etiquetas too
write.table(x = processed_datasets$train, file = file.path(output_data, "x_train.csv"),
            sep = ",", row.names = x_train$identificador, quote = F, col.names = T)
write.table(x = processed_datasets$test, file = file.path(output_data, "x_test.csv"),
            sep = ",", row.names = x_test$identificador, quote = F, col.names = T)