#' Testing scClassify model
#'
#' @param exprsMat_test A list or a matrix indicates the expression matrices of the testing datasets
#' @param trainRes A vector of cell types
#' @param cellTypes_test A list or a vector indicates cell types of the testing datasets (Optional).
#' @param k An integer indicates the number of neighbour
#' @param prob_threshold A numeric indicates the probability threshold for KNN/WKNN/DWKNN.
#' @param cor_threshold_static A numeric indicates the static correlation threshold.
#' @param cor_threshold_high A numeric indicates the highest correlation threshold
#' @param features A vector indicates the method to select features, set as "limma" by default.
#' This should be one or more of "limma", "DV", "DD", "chisq", "BI".
#' @param algorithm A vector indicates the KNN method that are used, set as "WKNN" by default. This
#' should be one or more of "WKNN", "KNN", "DWKNN".
#' @param similarity A vector indicates the similarity measure that are used, set as "pearson" by default.
#' This should be one or more of "pearson",  "spearman", "cosine", "jaccard", "kendall", "binomial", "weighted_rank","manhattan"
#' @param cutoff_method A vector indicates the method to cutoff the correlation distribution. Set as "dynamic" by default.
#' @param verbose A logical input indicates whether the intermediate steps will be printed
#' @return list of results
#' @author Yingxin Lin
#'
#' @importFrom proxy as.dist
#' @importFrom proxy dist
#' @import mixtools
#' @export



predict_scClassify <- function(exprsMat_test,
                               trainRes,
                               cellTypes_test,
                               k = 10,
                               prob_threshold = 0.8,
                               cor_threshold_static = 0.5,
                               cor_threshold_high = 0.7,
                               features = "pearson",
                               algorithm = c("WKNN", "KNN", "DWKNN"),
                               similarity = c("pearson",  "spearman",
                                              "cosine", "jaccard", "kendall",
                                              "weighted_rank","manhattan"),
                               cutoff_method = c("dynamic", "static"),
                               verbose = T){

  # check input
  algorithm <- match.arg(algorithm, c("WKNN", "KNN", "DWKNN"))
  similarity <- match.arg(similarity, c("pearson",  "spearman",
                                        "cosine", "jaccard", "kendall",
                                        "weighted_rank","manhattan"))
  cutoff_method <- match.arg(cutoff_method, c("dynamic", "static"))


  if (!features %in% names(trainRes$hierarchyKNNRes)) {
    stop("The selected features are not trained in the provided model!")
  }
  # get the input from train model

  levelModel <- trainRes$hierarchyKNNRes[[features]]$model
  levelHVG <- trainRes$hierarchyKNNRes[[features]]$hvg
  cutree_list <- trainRes$cutree_list

  pred <- list()

  # For each level
  for (i in 1:length(levelModel)) {
    # print(i)
    # If this level is not NULL (Not Level1)
    if (!is.null(levelModel[[i]])) {

      pred_level <- list()

      for (j in 1:length(levelModel[[i]])) {
        # print(paste(i,j))
        # If the model of level i-1, cells labeled as j is NOT "no model"
        if (class(levelModel[[i]][[j]]) != "character") {

          # Select the cells that are going to classified
          # (according to what they are classified in last level)
          predIdx <- which(pred[[i - 1]] == j)

          if (length(predIdx) != 0) {

            # features that are in the test dataset
            common_HVG <- intersect(rownames(exprsMat_test), levelHVG[[i]][[j]])
            exprsMat_toTest <- exprsMat_test[common_HVG, predIdx, drop = FALSE]
            exprsMat_toTest <- exprsMat_toTest[rowSums(exprsMat_toTest) != 0, , drop = FALSE]

            # Calculate the similarity
            corMat <- calculateSimilarity(exprsMat_train = t(levelModel[[i]][[j]]$train[ ,rownames(exprsMat_toTest)]),
                                          exprsMat_test = exprsMat_toTest,
                                          similarity = similarity)

            # Different algorithm
            if (algorithm == "KNN") {
              predRes <- KNNcor(corMat = corMat,
                                subLevelModel = levelModel[[i]][[j]],
                                cutoff_method = cutoff_method,
                                k = k,
                                prob_threshold = prob_threshold,
                                cor_threshold_static = cor_threshold_static,
                                cor_threshold_high = cor_threshold_high,
                                topLevel = F,
                                verbose = verbose)
            }
            if (algorithm == "WKNN") {
              predRes <- WKNNcor(corMat = corMat,
                                 subLevelModel = levelModel[[i]][[j]],
                                 cutoff_method = cutoff_method,
                                 k = k,
                                 prob_threshold = prob_threshold,
                                 cor_threshold_static = cor_threshold_static,
                                 cor_threshold_high = cor_threshold_high,
                                 topLevel = F,
                                 verbose = verbose)
            }

            if (algorithm == "DWKNN") {
              predRes <- DWKNNcor(corMat = corMat,
                                  subLevelModel = levelModel[[i]][[j]],
                                  cutoff_method = cutoff_method,
                                  k = k,
                                  prob_threshold = prob_threshold,
                                  cor_threshold_static = cor_threshold_static,
                                  cor_threshold_high = cor_threshold_high,
                                  topLevel = F,
                                  verbose = verbose)
            }





            pred_level[[j]] <- predRes$predRes


          }

        }

        # Else, the model of (level i-1, cells labeled as j) IS "no model"
        # maintain the same class.
        else {
          predIdx <- which(pred[[i - 1]] == j)
          # check the label of current level based on the label of last level
          pred_level[[j]] <- as.factor(rep(unique(cutree_list[[i]][cutree_list[[i - 1]] == j]), length(predIdx)))
          names(pred_level[[j]]) <- colnames(exprsMat_test)[predIdx]
        }
      }

      # Get the predict results for level i, change it to the numeric

      pred[[i]] <- unlist(lapply(pred_level, function(x){
        cellNames <- names(x)
        x <- as.numeric(as.character(x))
        names(x) <- cellNames
        x
      }))

      # reorder the prediction results to consistent with the exprsMat_test
      pred[[i]] <- pred[[i]][colnames(exprsMat_test)]

      names(pred[[i]]) <- colnames(exprsMat_test)
      pred[[i]] <- as.numeric(as.character(pred[[i]]))
      # there will be NA since they are unassigned from the last level, and therefore are not predicted in this level
      pred[[i]][is.na(pred[[i]])] <- 0
      names(pred[[i]]) <- colnames(exprsMat_test)
    }else{
      # If this level is NULL (Level 1)
      pred[[i]] <- as.factor(rep(1,ncol(exprsMat_test)))
      names(pred[[i]]) <- colnames(exprsMat_test)
    }

  }

  predMat <- do.call(cbind, pred)
  predMat <- sapply(1:ncol(predMat), function(x)predRes(predMat, cutree_list, x))

  predRes <- apply(predMat, 1, function(x){
    unAssIdx <- which(x == "unassigned")
    if (length(unAssIdx) != 0) {
      if (min(unAssIdx) >= 3) {
        x[min(unAssIdx) - 1]
      }else{
        "unassigned"
      }
    }else{
      x[length(x)]
    }
  })

  if (is.null(cellTypes_test)) {
    return(list(predRes = predRes, predLabelMat = predMat))
  }else{

    classify_res <-  ClassifyError(cellTypes_pred = predRes,
                                   cellTypes_test = cellTypes_test,
                                   cellTypes_train = trainRes$cellTypes_train)

    if (verbose) {
      print(table(classify_res)/length(classify_res))
    }

    return(list(predRes = predRes, predLabelMat = predMat, classifyRes = classify_res))

  }

}




# Function to calculate similarity

calculateSimilarity <- function(exprsMat_train,
                                exprsMat_test,
                                similarity = c("pearson",  "spearman",
                                               "cosine", "jaccard", "kendall",
                                               "weighted_rank","manhattan")) {

  similarity <- match.arg(similarity, c("pearson",  "spearman",
                                        "cosine", "jaccard", "kendall",
                                        "weighted_rank","manhattan"))

  if (similarity == "cosine") {

    corMat <- 1 - as.matrix(proxy::dist(t(exprsMat_train), t(exprsMat_test), method = "cosine"))
    corMat[is.na(corMat)] <- min(corMat)

  } else if (similarity == "kendall") {

    corMat <- stats::cor(exprsMat_train, exprsMat_test, method = "kendall")
    corMat[is.na(corMat)] <- -1

  } else if (similarity == "jaccard") {

    corMat <- 1 - as.matrix(proxy::dist(t(exprsMat_train > 0), t(exprsMat_test > 0), method = "Jaccard"))
    corMat[is.na(corMat)] <- min(corMat)

  }else if (similarity == "weighted_rank") {

    corMat <- wtd_rank2(exprsMat_train, exprsMat_test, method = "pearson")
    corMat[is.na(corMat)] <- -1

  }else if (similarity == "manhattan") {

    corMat <- 1 - as.matrix(proxy::dist(t(exprsMat_train), t(exprsMat_test), method = "Manhattan"))
    corMat[is.na(corMat)] <- min(corMat)

  }else if (similarity == "spearman") {

    corMat <- stats::cor(exprsMat_train, exprsMat_test, method = "spearman")
    corMat[is.na(corMat)] <- -1

  }else if (similarity == "pearson") {

    corMat <- stats::cor(exprsMat_train, exprsMat_test, method = "pearson")
    corMat[is.na(corMat)] <- -1

  }else{

    # performing pearson
    corMat <- stats::cor(exprsMat_train, exprsMat_test)
    corMat[is.na(corMat)] <- -1
  }

  return(corMat)
}


# Function to caculate weighted rank correlation
# the codes are modified from the wtd_rank function in dismay package

wtd_rank2 <- function(mat1, mat2 = NULL, method = "pearson") {

  method <- match.arg(method, c("pearson", "spearman"))

  ranks1 <- apply(mat1, 2, function(x) rank(-x, ties.method = "average"))
  # weight the ranks
  # calculate the savage scores
  n1 <- nrow(mat1)
  reciprocals1 <- 1 / seq_len(n1)
  savage1 <- sapply(seq_len(n1), function(i) sum(reciprocals1[i:n1]))
  # replace each rank with the savage score
  savages1 <- ranks1
  savages1[] <- savage1[ranks1]

  if (!is.null(mat2)) {
    ranks2 <- apply(mat2, 2, function(x) rank(-x, ties.method = "average"))
    # weight the ranks
    # calculate the savage scores
    n2 <- nrow(mat2)
    reciprocals2 <- 1 / seq_len(n2)
    savage2 <- sapply(seq_len(n2), function(i) sum(reciprocals2[i:n2]))
    # replace each rank with the savage score
    savages2 <- ranks2
    savages2[] <- savage2[ranks2]

    cor <- stats::cor(savages1, savages2, method = method)
  } else {
    cor <- stats::cor(savages1, method = method)
  }


  # calculate pearson correlation

  return(cor)
}

# Function to get the prediction labels according to the tree
predRes <- function(predMat, cutree_list, level){
  res <- predMat[,level]
  for (i in 1:length(predMat[,level])) {
    if (predMat[i,level] == 0) {
      res[i] <- "unassigned"
    }else{
      res[i] <- paste(names(cutree_list[[level]])[cutree_list[[level]] %in% predMat[i,level]],
                      collapse = "_")
    }

  }
  return(res)
}


# Functions to classify the error for the predicted cell types

ClassifyError <- function(cellTypes_pred, cellTypes_test, cellTypes_train){

  errClass <- c("correct", "correctly unassigned",  "intermediate", "incorrectly unassigned",
                "error assigned", "misclassified")


  if (length(cellTypes_pred) != length(cellTypes_test)) {
    stop("wrong input")
  }
  train_ref <- unique(cellTypes_train)
  # test_ref <- unique(cellTypes_test)
  res <- sapply(1:length(cellTypes_pred), function(i){
    if (cellTypes_test[i] %in% train_ref) {
      if (cellTypes_pred[i] %in% c("unassigned", "Unassigned")) {
        "incorrectly unassigned"
      } else if (cellTypes_pred[i] == "intermediate") {
        "intermediate"
      }else{
        if (cellTypes_test[i] == cellTypes_pred[i]) {
          "correct"
          # }else if(grepl(cellTypes_test[i], cellTypes_pred[i])|grepl("Node", cellTypes_pred[i])){
          #   "intermediate"
          # }
        } else if (length(strsplit(cellTypes_pred[i],"_")[[1]]) >= 2 | grepl("Node", cellTypes_pred[i])) {
          "intermediate"
        }else{
          "misclassified"
        }
      }
    }else{
      if (cellTypes_pred[i] %in% c("unassigned","Unassigned")) {
        "correctly unassigned"
      }else{
        "error assigned"
      }
    }
  })

  res <- factor(res, levels = errClass)
  return(res)
}


