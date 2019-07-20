#' @title Train and test scClassify model
#'
#' @param exprsMat_train A matrix of expression matrix of training dataset
#' @param cellTypes_train A vector of cell types of training dataset
#' @param exprsMat_test A list or a matrix indicates the expression matrices of the testing datasets
#' @param cellTypes_test A list or a vector indicates cell types of the testing datasets (Optional).
#' @param tree A vector indicates the method to build hierarchical tree, set as "HOPACH" by default.
#' This should be one of "HOPACH" and "HC" (using hclust).
#' @param selectFeatures A vector indicates the method to select features, set as "limma" by default.
#' This should be one or more of "limma", "DV", "DD", "chisq", "BI".
#' @param algorithm A vector indicates the KNN method that are used, set as "WKNN" by default. This
#' should be one or more of "WKNN", "KNN", "DWKNN".
#' @param similarity A vector indicates the similarity measure that are used, set as "pearson" by default.
#' This should be one or more of "pearson",  "spearman", "cosine", "jaccard", "kendall", "binomial", "weighted_rank","manhattan"
#' @param cutoff_method A vector indicates the method to cutoff the correlation distribution. Set as "dynamic" by default.
#' @param weighted_ensemble A logical input indicates in ensemble learning, whether the results is combined by a weighted score for each base classifier.
#' @param k An integer indicates the number of neighbour
#' @param topN An integer indicates the top number of features that are selected
#' @param hopach_kmax An integer between 1 and 9 specifying the maximum number of
#' children at each node in the HOPACH tree.
#' @param pSig A numeric indicates the cutoff of pvalue for features
#' @param prob_threshold A numeric indicates the probability threshold for KNN/WKNN/DWKNN.
#' @param cor_threshold_static A numeric indicates the static correlation threshold.
#' @param cor_threshold_high A numeric indicates the highest correlation threshold
#' @param parallel A logical input indicates whether running in paralllel or not
#' @param ncores An integer indicates the number of cores that are used
#' @param verbose A logical input indicates whether the intermediate steps will be printed
#' @return A \code{list} of results including
#' \itemize{
#' \item{testRes}: A list including the test results for each test dataset and base classifier
#' \item{trainRes}: A list including the training results including the features that are selected and the hierarchical cell type tree
#' \item{ensembleRes}: A list including the ensemble results and scores for the test dataset if ensemble learning is performed.
#' }
#'
#'
#' @author Yingxin Lin
#' @importFrom pbmcapply pbmclapply
#' @export


scClassify <- function(exprsMat_train = NULL,
                       cellTypes_train = NULL,
                       exprsMat_test = NULL,
                       cellTypes_test = NULL,
                       tree = "HOPACH",
                       algorithm = "WKNN",
                       selectFeatures = "limma",
                       similarity = "pearson",
                       cutoff_method = c("dynamic", "static"),
                       weighted_ensemble = FALSE,
                       k = 10,
                       topN = 50,
                       hopach_kmax = 5,
                       pSig = 0.01,
                       prob_threshold = 0.7,
                       cor_threshold_static = 0.5,
                       cor_threshold_high = 0.7,
                       parallel = FALSE,
                       ncores = 1,
                       verbose = FALSE){


  # check input
  if (is.null(exprsMat_train) | is.null(cellTypes_train) | is.null(exprsMat_test)) {
    stop("exprsMat_train or cellTypes_train or exprsMat_test is NULL!")
  }

  tree <- match.arg(tree, c("HOPACH", "HC"), several.ok = FALSE)
  selectFeatures <- match.arg(selectFeatures,
                              c("limma", "DV", "DD", "chisq", "BI"),
                              several.ok = TRUE)

  algorithm <- match.arg(algorithm,
                         c("WKNN", "KNN", "DWKNN"),
                         several.ok = TRUE)

  similarity <- match.arg(similarity, c("pearson",  "spearman",
                                        "cosine", "jaccard", "kendall",
                                        "weighted_rank","manhattan"), several.ok = TRUE)
  cutoff_method <- match.arg(cutoff_method, c("dynamic", "static"))


  # To check if need to run weighted ensemble learning

  if ((length(selectFeatures) > 1 | length(algorithm) > 1 | length(similarity) > 1) ) {
    if (weighted_ensemble) {
      weighted_ensemble <- TRUE
      ensemble <- TRUE

      if (verbose) {
        cat("Performing weighted ensemble learning... \n")
      }

    } else {
      weighted_ensemble <- FALSE
      ensemble <- TRUE

      if (verbose) {
        cat("Performing unweighted ensemble learning... \n")
      }

    }
  } else {
    weighted_ensemble <- FALSE
    ensemble <- FALSE

    if (verbose) {
      cat("Ensemble learning is disabled... \n")
    }

  }



  # QC for the training data set
  zeros <- apply(exprsMat_train, 1, function(x) sum(x == 0)/length(x))
  minPctCell <- min(table(cellTypes_train)/length(cellTypes_train))
  exprsMat_train <- exprsMat_train[zeros <= max(1 - minPctCell, 0.95), ]
  cat("after filtering not expressed genes \n")
  print(dim(exprsMat_train))

  # Add the train dataset in test data list if weighted ensemble learning
  if (weighted_ensemble) {
    if (class(exprsMat_test) == "matrix") {
      exprsMat_test <- list(test = exprsMat_test)
      cellTypes_test <- list(test = cellTypes_test)
    }
    exprsMat_test[[length(exprsMat_test) + 1]] <- exprsMat_train
    cellTypes_test[[length(cellTypes_test) + 1]] <- cellTypes_train
    names(exprsMat_test)[length(exprsMat_test)] <- names(cellTypes_test)[length(cellTypes_test)] <- "train"
  }



  trainRes <- train_scClassify(exprsMat_train,
                               cellTypes_train,
                               tree = tree,
                               selectFeatures = selectFeatures,
                               topN = topN,
                               hopach_kmax = hopach_kmax,
                               pSig = pSig,
                               parallel = parallel,
                               ncores = min(ncores, length(selectFeatures)),
                               verbose = verbose)

  ensemble_methods <- as.matrix(expand.grid(similarity = similarity,
                                            algorithm = algorithm,
                                            features = selectFeatures))

  if (verbose) {
    cat("Predicting using followings parameter combinations: \n")
    print(ensemble_methods)
  }


  ### if there are multiple testing datasets

  if (class(exprsMat_test) == "list") {
    testRes <- list()

    for (testDataset_idx in 1:length(exprsMat_test)) {

      if (verbose) {
        cat("Predicting: ")
        print(names(exprsMat_test)[testDataset_idx])
      }

      if (parallel) {

        predictRes <- pbmcapply::pbmclapply(1:nrow(ensemble_methods), function(em){
          predict_scClassify(exprsMat_test =  exprsMat_test[[testDataset_idx]],
                             trainRes = trainRes,
                             cellTypes_test = cellTypes_test[[testDataset_idx]],
                             k = k,
                             prob_threshold = prob_threshold,
                             cor_threshold_static = cor_threshold_static,
                             cor_threshold_high = cor_threshold_high,
                             features = ensemble_methods[em, 3],
                             algorithm = ensemble_methods[em, 2],
                             similarity = ensemble_methods[em, 1],
                             cutoff_method = cutoff_method,
                             verbose = verbose)
        }, mc.cores = ncores)

        names(predictRes) <- paste(ensemble_methods[, 1],
                                   ensemble_methods[, 2],
                                   ensemble_methods[, 3],
                                   sep = "_")

      }else{
        predictRes <- list()
        for (em in seq(nrow(ensemble_methods))) {

          if (verbose) {
            cat("Using parameters: \n")
            print(ensemble_methods[em, ])
          }

          predictRes[[em]] <- predict_scClassify(exprsMat_test =  exprsMat_test[[testDataset_idx]],
                                                 trainRes = trainRes,
                                                 cellTypes_test = cellTypes_test[[testDataset_idx]],
                                                 k = k,
                                                 prob_threshold = prob_threshold,
                                                 cor_threshold_static = cor_threshold_static,
                                                 cor_threshold_high = cor_threshold_high,
                                                 features = ensemble_methods[em, 3],
                                                 algorithm = ensemble_methods[em, 2],
                                                 similarity = ensemble_methods[em, 1],
                                                 cutoff_method = cutoff_method,
                                                 verbose = verbose)
        }

        names(predictRes) <- paste(ensemble_methods[, 1],
                                   ensemble_methods[, 2],
                                   ensemble_methods[, 3],
                                   sep = "_")
      }
      testRes[[testDataset_idx]] <- predictRes
    }

    names(testRes) <- names(exprsMat_test)


  }else{
    # else only one dataset as a matrix in the test
    if (parallel) {

      predictRes <- pbmcapply::pbmclapply(1:nrow(ensemble_methods), function(em){
        predict_scClassify(exprsMat_test = exprsMat_test,
                           trainRes = trainRes,
                           cellTypes_test = cellTypes_test,
                           k = k,
                           prob_threshold = prob_threshold,
                           cor_threshold_static = cor_threshold_static,
                           cor_threshold_high = cor_threshold_high,
                           features = ensemble_methods[em, 3],
                           algorithm = ensemble_methods[em, 2],
                           similarity = ensemble_methods[em, 1],
                           cutoff_method = cutoff_method,
                           verbose = verbose)
      }, mc.cores = ncores)

      names(predictRes) <- paste(ensemble_methods[, 1],
                                 ensemble_methods[, 2],
                                 ensemble_methods[, 3],
                                 sep = "_")

    }else{
      predictRes <- list()
      for (em in seq(nrow(ensemble_methods))) {
        predictRes[[em]] <- predict_scClassify(exprsMat_test = exprsMat_test,
                                               trainRes = trainRes,
                                               cellTypes_test = cellTypes_test,
                                               k = k,
                                               prob_threshold = prob_threshold,
                                               cor_threshold_static = cor_threshold_static,
                                               cor_threshold_high = cor_threshold_high,
                                               features = ensemble_methods[em, 3],
                                               algorithm = ensemble_methods[em, 2],
                                               similarity = ensemble_methods[em, 1],
                                               cutoff_method = cutoff_method,
                                               verbose = verbose)
      }

      names(predictRes) <- paste(ensemble_methods[, 1],
                                 ensemble_methods[, 2],
                                 ensemble_methods[, 3],
                                 sep = "_")
    }

    testRes <- list(test = predictRes)


  }


  if (ensemble) {

    ensembleRes <- getEnsembleRes(testRes, exclude = NULL, weighted_ensemble = weighted_ensemble)

    return(list(testRes = testRes, trainRes = trainRes, ensembleRes = ensembleRes))

  } else {
    return(list(testRes = testRes, trainRes = trainRes))
  }



}
