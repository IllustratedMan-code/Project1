# important -- this makes sure our runs are consistent and reproducible
set.seed(0)

file <- "TCGA_breast_cancer_ERstatus_allGenes.txt"
nfold <- 5
sd_threashold <- 4 # selects only data with sd > 4
# top_num <- 10

# This was relegated to an argument to cross_validation()
# top_num <- 50

header <- scan(file, nlines = 1, sep = "\t", what = character())
data <- read.table(file, skip = 2, header = FALSE, sep = "\t", quote = "", check.names = FALSE)

header[1] <- "gene_id"
names(data) <- header

header2 <- scan(file, skip = 1, nlines = 1, sep = "\t", what = character())

# cleanup - remove genes with sd < 1 (FIXME)
# This looks redundant since the genes above the sd_threshold are removed and the sd_threshold is 4

# compute sd for each gene
data_sd <- sapply(seq(nrow(data)), function(x) {
  as.numeric(sd(data[x, -1]))
})

# add gene names to the sd list
data_sd_names <- cbind(data.frame(data_sd), data[, 1])

# create an "include" list of all those genes where sd > threashold
# this removes
include_list <- data_sd_names[data_sd_names[, 1] > sd_threashold, 2]

Positive <- data[data$gene_id %in% include_list, header2 == "Positive"]
Negative <- data[data$gene_id %in% include_list, header2 == "Negative"]

# define function cross_valid so we can rerun the cross validataion with various parameters
cross_validation <- function(nfold, alg = "centroid", top_num = 50) {
  Positive_groups <- split(sample(colnames(Positive)), 1 + (seq_along(colnames(Positive)) %% nfold))
  Negative_groups <- split(sample(colnames(Negative)), 1 + (seq_along(colnames(Negative)) %% nfold))

  result <- array()

  for (test_group in 1:nfold) {
    testA <- Positive[, colnames(Positive) %in% unlist(Positive_groups[test_group])]
    testB <- Negative[, colnames(Negative) %in% unlist(Negative_groups[test_group])]

    trainingA <- Positive[, !(colnames(Positive) %in% unlist(Positive_groups[test_group]))]
    trainingB <- Negative[, !(colnames(Negative) %in% unlist(Negative_groups[test_group]))]

    # Feature selection --

    # compute t-statistic for each row
    training_t_stat <- data.frame(sapply(seq(nrow(trainingA)), function(x) {
      abs(as.numeric(t.test(trainingA[x, ], trainingB[x, ])$statistic))
    }))

    # add gene id column
    training_t_stat_geneid <- cbind(training_t_stat, rownames(trainingA))
    colnames(training_t_stat_geneid) <- c("t", "id")

    # pick top 50 based on t-statistic
    selected_genes <- head(training_t_stat_geneid[order(-training_t_stat_geneid$t), ], n = top_num)[, 2]

    # narrow down the list of genes based on t-statistic
    testA <- testA[rownames(testA) %in% selected_genes, ]
    testB <- testB[rownames(testB) %in% selected_genes, ]
    trainingA <- trainingA[rownames(trainingA) %in% selected_genes, ]
    trainingB <- trainingB[rownames(trainingB) %in% selected_genes, ]

    centroidA <- rowMeans(trainingA)
    centroidB <- rowMeans(trainingB)

    misclassifiedA <- sum(sapply(testA, function(x) {
      sqrt(sum((x - centroidA)^2)) - sqrt(sum((x - centroidB)^2)) > 0
    }))
    misclassifiedB <- sum(sapply(testB, function(x) {
      sqrt(sum((x - centroidA)^2)) - sqrt(sum((x - centroidB)^2)) < 0
    }))

    if (alg == "GLM") {
      trainingAB <- rbind(
        cbind(data.frame(t(trainingA)), cancer = 0),
        cbind(data.frame(t(trainingB)), cancer = 1)
      )
      testA0 <- data.frame(t(testA))
      testB0 <- data.frame(t(testB))
      model <- glm(cancer ~ ., trainingAB, family = "binomial")
      p <- predict(model, newdata = testA0, type = "response")
      misclassifiedA <<- sum(ifelse(p < 0.5, 0, 1))
      p <- predict(model, newdata = testB0, type = "response")
      misclassifiedB <<- sum(ifelse(p >= 0.5, 0, 1))
    }

    result[test_group] <- (misclassifiedA + misclassifiedB) / (ncol(testA) + ncol(testB))
  }
  paste0("mean=(", round(mean(result), 4), ") ", " sd=(", round(sd(result), 4), ")")
}

kable(data.frame(
  "Top 50 genes (centroid)" = cross_validation(5, top_num = 50),
  "Top 100 genes (centroid)" = cross_validation(5, top_num = 100)
))


print("extra credit")
kable(data.frame(
  "Top 50 genes (GLM)" = cross_validation(5, top_num = 50, alg = "GLM"),
  "Top 100 genes (GLM)" = cross_validation(5, top_num = 100, alg = "GLM")
))
