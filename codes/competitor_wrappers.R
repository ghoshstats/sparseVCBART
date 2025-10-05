####################### All competitor wrapper functions !! ###########################
###########===============================================================#############
library(rpart)
library(treeClust)
library(ranger)
library(gbm)
library(np)
library(vcrpart)
library(BART)

# Bart wrapper
##########################################################################################################################
bart_wrapper <- function(Y_train, X_train, Z_train, X_test, Z_test, nd = 1000, burn = 1000, thin = 1){
  
  bart_time <- system.time(
    bart_fit <- BART::wbart(x.train = cbind(X_train, Z_train), 
                            y.train = Y_train, 
                            x.test = cbind(X_test, Z_test),
                            ndpost = nd, nskip = burn, keepevery = thin,
                            printevery=nd*thin+burn+1))["elapsed"]
  sigma_samples <- bart_fit$sigma[-(1:burn)]
  f_samples_train <- bart_fit$yhat.train
  f_samples_test <- bart_fit$yhat.test
  
  # Summarize posterior over regression function
  fit_summary_train <- matrix(nrow = N_train, ncol = 3, dimnames = list(c(), c("MEAN","L95", "U95")))
  fit_summary_test <- matrix(nrow = N_test, ncol = 3, dimnames = list(c(), c("MEAN", "L95", "U95")))
  
  fit_summary_train[,"MEAN"] <- apply(f_samples_train, FUN = mean, MARGIN = 2, na.rm = TRUE)
  fit_summary_train[,"L95"] <- apply(f_samples_train, FUN = quantile, MARGIN = 2, probs = 0.025)
  fit_summary_train[,"U95"] <- apply(f_samples_train, FUN = quantile, MARGIN = 2, probs = 0.975)
  
  fit_summary_test[,"MEAN"] <- apply(f_samples_test, FUN = mean, MARGIN = 2, na.rm = TRUE)
  fit_summary_test[,"L95"] <- apply(f_samples_test, FUN = quantile, MARGIN = 2, probs = 0.025)
  fit_summary_test[,"U95"] <- apply(f_samples_test, FUN = quantile, MARGIN = 2, probs = 0.975)
  
  # Summarize posterior predictive
  ystar_summary_train <- matrix(nrow = N_train, ncol = 3, dimnames = list(c(), c("MEAN", "L95", "U95")))
  ystar_summary_test <- matrix(nrow = N_test, ncol = 3, dimnames = list(c(), c("MEAN", "L95", "U95")))
  
  N_train <- ncol(f_samples_train)
  N_test <- ncol(f_samples_test)
  
  for(i in 1:N_train){
    tmp_ystar_train <- f_samples_train[,i] + sigma_samples * rnorm(nrow(f_samples_train), 0, 1)
    ystar_summary_train[i,"MEAN"] <- mean(tmp_ystar_train)
    ystar_summary_train[i,"L95"] <- quantile(tmp_ystar_train, probs = 0.025)
    ystar_summary_train[i,"U95"] <- quantile(tmp_ystar_train, probs = 0.975)
  }
  
  for(i in 1:N_test){
    tmp_ystar_test <- f_samples_test[,i] + sigma_samples * rnorm(nrow(f_samples_test), 0, 1)
    ystar_summary_test[i,"MEAN"] <- mean(tmp_ystar_test)
    ystar_summary_test[i,"L95"] <- quantile(tmp_ystar_test, probs = 0.025)
    ystar_summary_test[i,"U95"] <- quantile(tmp_ystar_test, probs = 0.975)
  }
  
  return(list(time = bart_time, train_time = bart_time,
              train = list(fit = fit_summary_train, ystar = ystar_summary_train),
              test = list(fit = fit_summary_test, ystar = ystar_summary_test)))
}

##########################################################################################################################
# Wrapper function for BTVCM
honest.rpart.structure <- function(X, Y, method="standard", structY=NULL, leaf.size=3, control=NULL, diameter.test=NULL) {
  n <- nrow(X)
  colnames(X) <- paste("x", 1:ncol(X), sep="")
  pass.diameter.check <- FALSE
  while (!pass.diameter.check) {
    if (method=="random") {
      trainData <- data.frame(Y=rnorm(n), X=X)
    } else if (method == "hybrid") {
      trainData <- data.frame(Y = rnorm(n, 0, sd(structY)) + structY, X = X)
    } else {
      trainData <- data.frame(Y=structY, X=X)
    }
    if (is.null(control)) {
      tree.structure <- rpart(Y ~ . , data=trainData, control=rpart.control(cp = 0, minsplit=leaf.size+1), y = FALSE)
    } else {
      tree.structure <- rpart(Y ~ . , data=trainData, control=control, y = FALSE)
    }
    if (!is.null(diameter.test)) {
      colnames(diameter.test) <- paste("x", 1:ncol(X), sep="")
      where <- rpart.predict.leaves(tree.structure, newdata=data.frame(newdata), type="where")
      if (length(where)==length(unique(where))) return(tree.structure)
    } else {
      return(tree.structure)
    }
    method <- "random"
  }
}

honest.rpart <- function(X, Y, method="standard", structY=NULL, subset=NULL, leaf.size=3, control = NULL, diameter.test=NULL) {
  colnames(X) <- paste("x", 1:ncol(X), sep="")
  trainData <- data.frame(Y, X)
  
  tree <- list()
  tree$train.where <- rep(0, length(Y))
  if (is.null(subset)) {
    newX <- X
    newY <- Y
  } else {
    subset <- sort(subset)
    newX <- X[subset, ]
    newY <- Y[subset]
  }
  if (!is.null(diameter.test)) {
    tree$rpart.tree <- honest.rpart.structure(X, Y, method=method, structY=structY, leaf.size=leaf.size, control=control, diameter.test=diameter.test)
  } else {
    tree$rpart.tree <- honest.rpart.structure(X, Y, method=method, structY=structY, leaf.size=leaf.size, control=control)
  }
  
  tmpWhere <- unique(tree$rpart.tree$where)
  where <- rpart.predict.leaves(tree$rpart.tree, newdata=data.frame(X=newX), type="where")
  tmpPredict <- c()
  for (i in tmpWhere) {
    tmpCount <- sum(where==i)
    if (tmpCount > 0) {
      tmpPredict <- c(tmpPredict, mean(newY[where==i]))
      tree$train.where[subset[where==i]] = i
    } else {
      tmpPredict <- c(tmpPredict, 0)
    }
  }
  tree$predict <- tmpPredict
  names(tree$predict) <- tmpWhere
  return(tree)
}

honest.rpart.predict <- function(tree, newdata) {
  if (is.vector(newdata)) {
    newdata = matrix(newdata, nrow=1)
  }
  colnames(newdata) = paste("x", 1:ncol(newdata), sep="")
  where <- rpart.predict.leaves(tree$rpart.tree, newdata=data.frame(X=newdata), type="where")
  return(as.vector(tree$predict[paste(where)]))
}

honest.rpart.predict.weight <- function(tree, newdata) {
  if (is.vector(newdata)) {
    newdata = matrix(newdata, nrow=1)
  }
  colnames(newdata) = paste("x", 1:ncol(newdata), sep="")
  where <- rpart.predict.leaves(tree$rpart.tree, newdata=data.frame(X=newdata), type="where")
  ans <- t(sapply(where, function(x) return(tree$train.where==x)))*1
  for (i in 1:nrow(ans)){
    w <- sum(ans[i, ])
    if(w>0) {
      ans[i, ] = ans[i, ] / w
    }
  }
  return(ans)
}

predict.boulevard <- function(blv, X) {
  ntree <- length(blv$trees)
  lambda <- blv$lambda
  ans <- rep(0, nrow(X))
  for (b in 1:ntree) {
    ans <- (b-1)/b*ans + lambda/b*honest.rpart.predict(blv$trees[[b]], newdata = X)
  }
  return(ans*(1+lambda)/lambda)
}

predict.boulevard.variance <- function(blv, newdata, narrow = FALSE) {
  ntree <- length(blv$trees)
  lambda <- blv$lambda
  ans <- honest.rpart.predict.weight(blv$trees[[1]], newdata=newdata)
  if (ntree > 1) {
    for (b in 2:ntree) {
      ans <- ans + honest.rpart.predict.weight(blv$trees[[b]], newdata=newdata)
    }
    ans <- ans / ntree
  }
  ans <- apply(ans, 1, function(x) sum(x*x))
  ans <- ans * (1+lambda)^2 * blv$sigma2
  if (narrow) {
    ans <- ans / (1+lambda)^2
  } 
  return(ans)
}

boulevard <- function(X, Y, ntree=1000, lambda = 0.8, subsample=0.8, xtest=NULL, ytest=NULL, leaf.size=10, method="random") {
  n <- nrow(X)
  tree <- list()
  nss <- floor(n * subsample)
  ans <- rep(0, nrow(X))
  if (!is.null(xtest)) {
    predtest <- rep(0, nrow(xtest))
  }
  trainmse <- c()
  testmse <- c()
  for (b in 1:ntree) {
    if (b%%50 == 0) {
      cat("Training Iteration:", b, "\n")
    }
    res <- Y - ans
    if (method=="random") {
      tree[[b]] <- honest.rpart(X, res, method=method, subset=sample(n, nss, replace=FALSE), leaf.size=leaf.size)
    } else if (method=="standard") {
      tree[[b]] <- honest.rpart(X, res, method=method, structY=res, subset=sample(n, nss, replace=FALSE), leaf.size=leaf.size)
    }
    ans <- (b-1)/b*ans + lambda/b*honest.rpart.predict(tree[[b]], newdata = X)
    trainmse <- c(trainmse,mean((ans/lambda*(1+lambda)-Y)^2))
    if (!is.null(xtest)) {
      predtest <- (b-1)/b*predtest + lambda/b*honest.rpart.predict(tree[[b]], newdata=xtest)
      testmse <- c(testmse, mean((predtest/lambda*(1+lambda)-ytest)^2))
    }
  }
  return(list(trees=tree, 
              mse=trainmse, 
              testmse=testmse, 
              lambda=lambda,
              sigma2=trainmse[ntree]))
}


OUT_LOG_FLAG <- TRUE
OUT_LOG_LEVEL <- 0
GLOBAL_MAX <- 1e100

logging <- function(x=NULL, msg="", level=1) {
  if(OUT_LOG_FLAG && level > OUT_LOG_LEVEL) {
    cat(paste("Logging << ", msg, "\n"))
    if(!is.null(x)) {
      print(x)
    }
  }
}

colorize <- function(x) {
  x <- rank(x)
  if (max(x) == min(x)) {
    return(rep(0.5, length(x)))
  } else {
    return((x-min(x)) / (max(x)-min(x)))
  }
}
#### Witchwoods Part ####
crown <- function(z) {
  colnames(z) <- paste("z", 1:ncol(z), sep="")
  return(z)
}

pred.woods <- function(z, model) {
  if (is.vector(z)) {
    z <- matrix(z, nrow=1)
  }
  tmpBeta <- matrix(0, nrow=nrow(z), ncol=model$p)
  for (i in 1:p) {
    tmpBeta[, i] <- predict(model$woods[[1]][[i]], crown(z))
    if (model$woods.len > 1) {
      for (j in 2:model$woods.len) {
        tmpBeta[, i] <- tmpBeta[, i] + model$lambda *  predict(model$woods[[j]][[i]], crown(z))
      }
    }
  }
  return(tmpBeta)
}

#### for ols
ols.diff <- function(y, x, beta) {
  return((y - rowSums(x * beta)) * x)
}

ols.init.beta <- function(y, x) {
  return(matrix(as.vector(lm(y ~ x+0)$coef), 
                nrow=nrow(x), ncol=ncol(x), byrow=TRUE))
}

ols.pred <- function(x, beta) {
  return(rowSums(x * beta))
}

ols.loss <- function(y1, y2) {
  return(mean((y1-y2)^2))
}

#### for logistic
sigmoid <- function(z) {
  return(1/(1+exp(-z)))    
}

lg.pred <- function(x, beta) {
  return(rowSums(x * beta))
}

lg.diff <- function(y, x, beta) {
  return(x*(y-sigmoid(lg.pred(x, beta))))
}

lg.init.beta <- function(y, x) {
  return(matrix(as.vector(glm(y ~ x+0,
                              family = binomial(link = "logit"))$coef), 
                nrow=nrow(x), ncol=ncol(x), byrow=TRUE))
}

lg.loss <- function(y1, y2) {
  return(-mean((1-y1)*(-y2)-log(1+exp(-y2))))
}

grow.woods <- function(u, z, control, is.plinear = FALSE) {
  if (ncol(u) != model$p) {
    logging(msg="Err: dim mismatch")
  }
  if (is.plinear) {
    w_control <- control
    w_control$cp <- 1e10
  }
  tmpWoods <- list()
  if (is.plinear) {
    tmpWoods[[1]] <- rpart(u[, 1] ~ crown(z), control=control, y = FALSE)
    if (model$p > 1) {
      for (i in 2:model$p) {
        tmpWoods[[i]] <- rpart(u[, i] ~ crown(z), control=w_control, y = FALSE)
      }
    }
  } else {
    for (i in 1:model$p) {
      tmpWoods[[i]] <- rpart(u[, i] ~ crown(z), control=control, y = FALSE)
    }
  }
  return(tmpWoods)
}

grow.woods.v2 <- function(u, z, model) {
  if (ncol(u) != model$p) {
    logging(msg="Err: dim mismatch")
  }
  if (model$method == "plinear") {
    w_control <- model$control
    w_control$cp <- 1e10
  }
  tmpSubset <- sample(model$n, model$n * model$subsample)
  tmpDF <- as.data.frame(crown(z))
  tmpWoods <- list()
  if (model$method == "plinear") {
    tmpDF$response = u[, 1]
    tmpWoods[[1]] <- rpart(response ~ ., data=tmpDF, subset=tmpSubset, control=model$control, y = FALSE)
    if (model$p > 1) {
      for (i in 2:model$p) {
        tmpDF$response = u[, i]
        tmpWoods[[i]] <- rpart(response ~ ., data=tmpDF, subset = tmpSubset, control=w_control, y = FALSE)
      }
    }
  } else if (model$method == "boulevard") {
    for (i in 1:model$p) {
      tmpWoods[[i]] <- honest.rpart(Y=u[, i], X=z, structY=u[, i], 
                                    subset=tmpSubset, 
                                    control=model$control)
    }
  } else {
    for (i in 1:model$p) {
      tmpDF$response = u[, i]
      tmpWoods[[i]] <- rpart(response ~., data=tmpDF, subset=tmpSubset, control=model$control, y = FALSE)
    }
  }
  return(tmpWoods)
}

predict.woods <- function(tree, model, newdata = NULL, train=FALSE) {
  if (model$method == "boulevard") {
    return(honest.rpart.predict(tree, newdata = newdata))
    # } else if (train) {
    #     return(predict(tree))
  } else {
    return(predict(tree, newdata = as.data.frame(crown(newdata))))
  }
}

lgd <- function(y, x, z, model) {
  #### set default values
  if (is.null(model$savetrees)) {
    model$savetrees = TRUE
  }
  if (is.null(model$xscale)) {
    model$xscale = FALSE
  }
  if (is.null(model$yscale)) {
    model$yscale = FALSE
  }
  if (is.null(model$dummy)) {
    model$dummy = FALSE
  }
  if (is.null(model$method)) {
    model$method = "ordinary"
  } 
  if (is.null(model$quantile_cut)) {
    model$quantile_cut = 0.05
  }
  if (model$method == "boulevard") {
    if (is.null(model$control$minsplits)) {
      model$control$minsplits = 10
    }
    if (is.null(model$subsample)) {
      model$subsample = 0.5
    }
  }
  
  #### scale predictor / responses
  if (model$xscale) {
    x <- scale(x)
    # print(head(x))
    model$xscale.center = attr(x, "scaled:center")
    model$xscale.scale = attr(x, "scaled:scale")
    if (model$dummy) {
      x <- cbind(1, x)
      model$xscale.center <- c(-1, model$xscale.center)
      model$xscale.scale <- c(1, model$xscale.scale)
    }
    if (model$yscale) {
      y <- scale(y)
      model$yscale.center = attr(y, "scaled:center")
      model$yscale.scale = attr(y, "scaled:scale")
      y <- as.vector(y)
    } else {
      model$yscale.center = 0
      model$yscale.scale = 1
    }
  } else if (model$dummy) {
    x <- cbind(1, x)
  }
  
  if (model$method == "boulevard") {
    x <- x * sqrt(model$n)
  }
  z <- crown(z)
  #### lgb
  model$lc <- c()
  # if (model$method == "boulevard") {
  #     tmpBeta <- model$diff(y, x, 0)
  # } else {
  tmpBeta <- model$init(y, x)
  # }
  model$woods[[1]] <- grow.woods.v2(tmpBeta, z, model)
  model$beta <- matrix(0, nrow=nrow(x), ncol=ncol(x))
  for (i in 1:model$p) {
    model$beta[, i] <- predict.woods(model$woods[[1]][[i]], model, newdata=z)
  }
  model$yhat <- model$pred(x, model$beta)
  model$lc <- c(model$lc, model$loss(y, model$yhat))
  if (model$ntree > 1) {
    for (ntree in 2:model$ntree) {
      gc()
      if (model$savetrees) {
        w_ntree = ntree
      } else {
        w_ntree = 1
      }
      #if (ntree %% 5 == 0) {cat("ITER:", ntree, "\tSave to", w_ntree)} #[modification by SKD to suppress print statements]
      # if (model$method == "boulevard") {
      #     u <- model$diff(y, x, model$beta / ntree * (ntree - 1))
      # } else {
      u <- model$diff(y, x, model$beta)
      # }
      for (i in 1:model$p) {
        if (!is.null(model$gradient_truncate)) {
          quantile_truncate <- c(-model$gradient_truncate, model$gradient_truncate)
        } else {
          quantile_truncate = as.vector(quantile(u[, i], 
                                                 probs=c(model$quantile_cut, 1-model$quantile_cut)))
          u[u[, i] < quantile_truncate[1], i] = quantile_truncate[1]
          u[u[, i] > quantile_truncate[2], i] = quantile_truncate[2]
        }
      }
      model$woods[[w_ntree]] <- grow.woods.v2(u, z, model)
      for (i in 1:model$p) {
        if (model$method == "ordinary") {
          model$beta[, i] <- model$beta[, i] + model$lambda * predict.woods(model$woods[[w_ntree]][[i]], model, newdata=z)
        } else if (model$method == "boulevard") {
          model$beta[, i] <- model$beta[, i] * (ntree - 1) / ntree + 
            model$lambda / ntree * predict.woods(model$woods[[w_ntree]][[i]], model, newdata=z)
        } else {
          model$beta[, i] <- model$beta[, i] + model$lambda * predict.woods(model$woods[[w_ntree]][[i]], model, newdata=z)
        }
      }
      # if (ntree %% 10 == 0) {print(head(model$beta))}
      model$yhat <- model$pred(x, model$beta)
      # plot(y, model$yhat)
      if (model$method=="boulevard") {
        tmpInflate <- lm(y ~ model$pred(x, model$beta) + 0)$coef[1]
        # boxplot(y / model$pred(x, model$beta), breaks=50)
        #if(ntree %% 5 == 0) {cat("\tblvRefLOSS:", model$loss(y, model$pred(x, model$beta * (1 + model$lambda) / model$lambda)))}
        #if(ntree %% 5 == 0) {cat("\tInflate:", tmpInflate)}
        #if(ntree %% 5 == 0) {cat("\tblvInfLOSS:", model$loss(y, model$pred(x, model$beta * tmpInflate)))}
        tmpLoss <- model$loss(y, model$yhat)
      } else {
        tmpLoss <- model$loss(y, model$yhat)
      }
      # tmpLoss <- model$loss(y, model$yhat)
      model$lc <- c(model$lc, tmpLoss)
      #if(ntree %% 5 == 0) {cat("\tLOSS:", tmpLoss, "\n")}
    }
  }
  
  if (model$method == "boulevard") {
    if (!is.null(model$inflate)) {
      model$inflate = lm(y ~ model$pred(x, model$beta) + 0)$coef[1]
    }
    model$beta <- model$beta * sqrt(model$n)
  }
  #### Unscale
  if (model$xscale) {
    for (i in 1:model$p) {
      model$beta[, i] <- model$beta[, i] * model$yscale.scale / model$xscale.scale[i]
    }
    if (model$dummy) model$beta[, 1] <- -model$beta%*%model$xscale.center  + model$yscale.center
  }
  return(model)
}

lgd.predict <- function(x, z, model, flag = FALSE) {
  if (model$dummy) {
    x <- cbind(1, x)
  }
  # print(head(x))
  tmpBeta <- matrix(0, nrow=nrow(x), ncol=ncol(x))
  for (ntree in 1:model$ntree) {
    for (i in 1:model$p) {
      if (model$method == "ordinary") {
        if (ntree == 1) {
          tmpBeta[, i] <- predict.woods(model$woods[[ntree]][[i]], model, newdata = crown(z))
        } else {
          tmpBeta[, i] <- tmpBeta[, i] + model$lambda * predict.woods(model$woods[[ntree]][[i]], model, newdata = crown(z))
        }
      } else if (model$method == "boulevard") {
        tmpBeta[, i] <- tmpBeta[, i] * (ntree-1) / ntree + model$lambda / ntree * predict.woods(model$woods[[ntree]][[i]], model, newdata = crown(z))
      }
    }
  }
  if (model$method == "boulevard") {
    tmpBeta <- tmpBeta * sqrt(model$n)
  }
  #### Unscale
  if (model$xscale) {
    for (i in 1:model$p) {
      tmpBeta[, i] <- tmpBeta[, i] * model$yscale.scale / model$xscale.scale[i]
    }
    if (model$dummy) tmpBeta[, 1] <- -tmpBeta%*%model$xscale.center  + model$yscale.center
  }
  if (model$method == "boulevard") {
    if (flag) {
      tmpBeta = tmpBeta * (1 + model$lambda) / model$lambda 
    }
  }
  return(model$pred(x, tmpBeta))
}

boosted_tvcm_wrapper <- function(Y_train, X_train, Z_train, X_test, Z_test,
                                 intercept = TRUE, lambda = 0.05, ntree = 200, B = 50)
{
  n_obs_train <- nrow(X_train)
  n_obs_test <- nrow(X_test)
  p <- ncol(X_train)
  
  model <- list()
  model$dummy <- intercept
  model$xscale <- TRUE
  model$n <- n_obs_train
  model$p <- p + model$dummy
  model$q <- ncol(Z_train)
  model$diff <- ols.diff
  model$diff <- ols.diff
  model$init <- ols.init.beta
  model$pred <- ols.pred
  model$loss <- ols.loss
  model$lambda <- lambda
  model$subsample <- 0.5
  model$ntree <- ntree
  model$control <- rpart.control(maxdepth = 6, cp=0.0001)
  model$woods <- list()
  model$savetree <- TRUE
  
  time_train <- system.time(  
    tvcm_fit <- tryCatch(lgd(y=Y_train, x=X_train, z=Z_train, model=model),
                         error = function(e){return(NULL)}))["elapsed"]
  
  if(!is.null(tvcm_fit)){
    yhat_train <- tvcm_fit$yhat
    time_test <- system.time(yhat_test <- lgd.predict(X_test, Z_test, tvcm_fit))["elapsed"]
    beta_train <- tvcm_fit$beta
    
    beta_test <- matrix(nrow = nrow(X_test), ncol = 1 + p)
    beta_time <- 0
    tmp_X <- matrix(0, nrow = nrow(X_test), ncol = p)
    beta_time <- system.time(beta_test[,1] <- lgd.predict(tmp_X, Z_test, tvcm_fit))["elapsed"]
    for(j in 1:ncol(X_test)){
      tmp_X <- matrix(0, nrow = nrow(X_test), ncol = p)
      tmp_X[,j] <- 1
      tmp_time <- system.time(beta_test[,j+1] <- lgd.predict(tmp_X, Z_test, tvcm_fit) - beta_test[,1])["elapsed"] # remember to remove the intercept!
      beta_time <- beta_time + tmp_time
    }
    if(B == 0){
      ystar_summary_train <- matrix(nrow = n_obs_train, ncol = 3, dimnames = list(c(), c("MEAN", "L95", "U95")))
      ystar_summary_test <- matrix(nrow = n_obs_test, ncol = 3, dimnames = list(c(), c("MEAN", "L95", "U95")))
      rmse_train <- sqrt(mean( (Y_train - yhat_train)^2))
      ystar_summary_train[,"MEAN"] <- yhat_train
      ystar_summary_test[,"MEAN"] <- yhat_test
      
      ystar_summary_train[,"L95"] <- ystar_summary_train[,"MEAN"] - qnorm(0.975) * rmse_train
      ystar_summary_train[,"U95"] <- ystar_summary_train[,"MEAN"] + qnorm(0.97) * rmse_train
      
      ystar_summary_test[,"L95"] <- ystar_summary_test[,"MEAN"] - qnorm(0.975) * rmse_train
      ystar_summary_test[,"U95"] <- ystar_summary_test[,"MEAN"] + qnorm(0.975 ) * rmse_train
      results <- list(train = list(ystar = ystar_summary_train),
                      test = list(ystar = ystar_summary_test),
                      time = time_train + beta_time,
                      train_time = time_train)      
    } else{
      ######
      # Do the bootstrap now
      ######
      boot_beta <- array(dim = c(n_obs_train + n_obs_test, p + 1, B))
      
      boot_time <- tryCatch(
        {
          system.time(
            {
              for(b in 1:B){
                if(b %% 5 == 0) print(paste("tvcm-boot: b = ", b, Sys.time()))
                boot_train <- sample(1:n_obs_train, size = n_obs_train, replace = TRUE)
                boot_fit <- lgd(y=Y_train[boot_train], x=X_train[boot_train,], z=Z_train[boot_train,], model=model)
                
                boot_beta[1:n_obs_train,,b] <- boot_fit$beta
                
                tmp_X <- matrix(0, nrow = n_obs_train + n_obs_test, ncol = ncol(X_train))
                tmp_int <- lgd.predict(tmp_X, rbind(Z_train, Z_test), boot_fit)
                
                boot_beta[,1,b] <- tmp_int
                for(j in 1:p){
                  tmp_X <- matrix(0, nrow = n_obs_train + n_obs_test, ncol = p)
                  tmp_X[,j] <- 1
                  boot_beta[,j+1,b] <- lgd.predict(tmp_X, rbind(Z_train, Z_test), boot_fit) - tmp_int
                }
              }
            }
          )["elapsed"] 
        }, error = function(e){return(NULL)})
      if(!is.null(boot_time)){
        boot_beta_se <- apply(boot_beta, FUN = sd, MAR = c(1,2))
        
        beta_summary_train <- array(dim = c(n_obs_train, 4, 1 + p), dimnames = list(c(), c("MEAN", "SD", "L95", "U95"), c()))
        beta_summary_test <- array(dim = c(n_obs_test, 4, 1 + p), dimnames = list(c(), c("MEAN", "SD", "L95", "U95"), c()))
        
        ystar_summary_train <- matrix(nrow = n_obs_train, ncol = 3, dimnames = list(c(), c("MEAN", "L95", "U95")))
        ystar_summary_test <- matrix(nrow = n_obs_test, ncol = 3, dimnames = list(c(), c("MEAN", "L95", "U95")))
        
        beta_summary_train[,"MEAN",] <- beta_train
        beta_summary_test[,"MEAN",] <- beta_test
        
        beta_summary_train[,"SD",] <- boot_beta_se[1:n_obs_train,]
        beta_summary_test[,"SD",] <- boot_beta_se[(1 + n_obs_train):(n_obs_train + n_obs_test),]
        
        beta_summary_train[,"L95",] <- beta_summary_train[,"MEAN",] - qnorm(0.975) * beta_summary_train[,"SD",]
        beta_summary_train[,"U95",] <- beta_summary_train[,"MEAN",] + qnorm(0.975) * beta_summary_train[,"SD",]
        
        beta_summary_test[,"L95",] <- beta_summary_test[,"MEAN",] - qnorm(0.975) * beta_summary_test[,"SD",]
        beta_summary_test[,"U95",] <- beta_summary_test[,"MEAN",] + qnorm(0.975) * beta_summary_test[,"SD",]
        
        rmse_train <- sqrt(mean( (Y_train - yhat_train)^2))
        ystar_summary_train[,"MEAN"] <- yhat_train
        ystar_summary_test[,"MEAN"] <- yhat_test
        
        ystar_summary_train[,"L95"] <- ystar_summary_train[,"MEAN"] - qnorm(0.975) * rmse_train
        ystar_summary_train[,"U95"] <- ystar_summary_train[,"MEAN"] + qnorm(0.975) * rmse_train
        
        
        ystar_summary_test[,"L95"] <- ystar_summary_test[,"MEAN"] - qnorm(0.975) * rmse_train
        ystar_summary_test[,"U95"] <- ystar_summary_test[,"MEAN"] + qnorm(0.975 ) * rmse_train
        
        
        results <- list(train = list(beta = beta_summary_train, ystar = ystar_summary_train),
                        test = list(beta = beta_summary_test, ystar = ystar_summary_test),
                        time = time_train + beta_time + boot_time,
                        train_time = time_train,
                        boot_time = boot_time)
      } else results <- NULL  
    }
  } else results <- NULL
  return(results)
  
}


####### Extra trees wrapper ######
##########################################################################################################################

extraTrees_wrapper <- function(Y_train, X_train, Z_train, X_test, Z_test){
  n_obs_train <- nrow(X_train)
  n_obs_test <- nrow(X_test)
  
  p <- ncol(X_train)
  R <- ncol(Z_train)
  
  # Define necessary data frames
  colnames(X_train) <- paste0("X", 1:p)
  colnames(X_test) <- paste0("X",1:p)
  
  colnames(Z_train) <- paste0("Z", 1:R)
  colnames(Z_test) <- paste0("Z", 1:R)
  
  tmp_data_train <- data.frame(Y_train,X_train, Z_train)
  tmp_data_test <- data.frame(Y_test,X_test, Z_test)
  colnames(tmp_data_train)[1] <- "Y"
  colnames(tmp_data_test)[1] <- "Y"
  
  train_time <- system.time(rf_fit <- 
                              tryCatch(ranger(Y~ ., data = tmp_data_train, write.forest = TRUE, splitrule = "extratrees"),
                                       error = function(e){return(NULL)}))["elapsed"]
  if(!is.null(rf_fit)){
    ystar_train <- predict(rf_fit, data = tmp_data_train[,-1])$predictions
    ystar_test <- predict(rf_fit, data = tmp_data_test[,-1])$predictions
    
    rmse_train <- sqrt(mean( (Y_train - ystar_train)^2 ))
    
    ystar_summary_train <- matrix(nrow = n_obs_train, ncol = 3, dimnames = list(c(), c("MEAN", "L95", "U95")))
    ystar_summary_test <- matrix(nrow = n_obs_test, ncol = 3, dimnames = list(c(), c("MEAN", "L95", "U95")))
    ystar_summary_train[,"MEAN"] <- ystar_train
    ystar_summary_test[,"MEAN"] <- ystar_test
    
    ystar_summary_train[,"L95"] <- ystar_summary_train[,"MEAN"] - qnorm(0.975) * rmse_train
    ystar_summary_train[,"U95"] <- ystar_summary_train[,"MEAN"] + qnorm(0.97) * rmse_train
    ystar_summary_test[,"L95"] <- ystar_summary_test[,"MEAN"] - qnorm(0.975) * rmse_train
    ystar_summary_test[,"U95"] <- ystar_summary_test[,"MEAN"] + qnorm(0.975 ) * rmse_train
    
    results <- list(train = list(ystar = ystar_summary_train),
                    test = list(ystar = ystar_summary_test),
                    time = train_time, train_time = train_time) 
    
  } else{
    results <- NULL
  }
  
  return(results)
}


gbm_wrapper <- function(Y_train, X_train, Z_train, X_test, Z_test){
  
  n_obs_train <- nrow(X_train)
  n_obs_test <- nrow(X_test)
  
  p <- ncol(X_train)
  R <- ncol(Z_train)
  
  # Define necessary data frames
  colnames(X_train) <- paste0("X", 1:p)
  colnames(X_test) <- paste0("X",1:p)
  
  colnames(Z_train) <- paste0("Z", 1:R)
  colnames(Z_test) <- paste0("Z", 1:R)
  
  tmp_data_train <- data.frame(Y_train,X_train, Z_train)
  tmp_data_test <- data.frame(Y_test,X_test, Z_test)
  colnames(tmp_data_train)[1] <- "Y"
  colnames(tmp_data_test)[1] <- "Y"
  
  
  train_time <- system.time(
    gbm_fit <- tryCatch(
      gbm(Y ~ ., data = tmp_data_train, distribution = "gaussian",
          interaction.depth = 6, cv.folds = 5, n.trees = 1000)
    )
  )["elapsed"]
  
  
  if(!is.null(gbm_fit)){
    n_tree_opt <- gbm.perf(gbm_fit, method = "cv", plot.it = FALSE)
    ystar_train <- predict(gbm_fit, newdata = tmp_data_train, n.trees = n_tree_opt)
    ystar_test <- predict(gbm_fit, newdata = tmp_data_test, n.trees = n_tree_opt)
    
    rmse_train <- sqrt(mean( (Y_train - ystar_train)^2 ))
    
    ystar_summary_train <- matrix(nrow = n_obs_train, ncol = 3, dimnames = list(c(), c("MEAN", "L95", "U95")))
    ystar_summary_test <- matrix(nrow = n_obs_test, ncol = 3, dimnames = list(c(), c("MEAN", "L95", "U95")))
    ystar_summary_train[,"MEAN"] <- ystar_train
    ystar_summary_test[,"MEAN"] <- ystar_test
    
    ystar_summary_train[,"L95"] <- ystar_summary_train[,"MEAN"] - qnorm(0.975) * rmse_train
    ystar_summary_train[,"U95"] <- ystar_summary_train[,"MEAN"] + qnorm(0.97) * rmse_train
    ystar_summary_test[,"L95"] <- ystar_summary_test[,"MEAN"] - qnorm(0.975) * rmse_train
    ystar_summary_test[,"U95"] <- ystar_summary_test[,"MEAN"] + qnorm(0.975 ) * rmse_train
    
    results <- list(train = list(ystar = ystar_summary_train),
                    test = list(ystar = ystar_summary_test),
                    time = train_time, train_time = train_time) 
  } else{
    results <- NULL
  }
  
  return(results)
}

####### wrapper for kernel smoothing (KS) #########
##########################################################################################################################
# Prepare the outcomes
# beta_summary_ gives the point estimate, bootstrap standard error, and bootstrap CI for each beta
# fit_summary_ gives the point estimate, SE returned by npscoef, and approximate CI for E[y | x,z]
# ypred_summary_ gives the point estimate and upper/lower 95 prediction interval. Here we just point estimate +/- 2 * RMSE.

# Wrapper for kernel smoothing (np::npscoef) in the same style as tvc_wrapper
kernel_smoothing_wrapper <- function(Y_train, X_train, Z_train,
                                     X_test,  Z_test,
                                     B = 50) {
  
  n_obs_train <- nrow(X_train)
  n_obs_test  <- nrow(X_test)
  p <- ncol(X_train)
  
  # Bandwidth selection (train only)
  train_time <- system.time(
    bw <- tryCatch(
      npscoefbw(xdat = X_train, ydat = Y_train, zdat = Z_train),
      error = function(e) NULL
    )
  )["elapsed"]
  
  if (is.null(bw)) return(NULL)
  
  # Fit and predict for train+test (no eydat; we won't ask for residuals/errors)
  ks_fit <- tryCatch(
    npscoef(
      bw,
      betas = TRUE, residuals = FALSE, errors = FALSE,
      txdat = X_train, tydat = Y_train, tzdat = Z_train,
      exdat = rbind(X_train, X_test),
      ezdat = rbind(Z_train, Z_test)
    ),
    error = function(e) NULL
  )
  if (is.null(ks_fit)) return(NULL)
  
  # Predictions
  yhat_train <- ks_fit$mean[seq_len(n_obs_train)]
  yhat_test  <- ks_fit$mean[n_obs_train + seq_len(n_obs_test)]
  
  if (B == 0) {
    # Fast y* intervals from train RMSE
    rmse_train <- sqrt(mean((Y_train - yhat_train)^2))
    q975 <- qnorm(0.975)
    
    ystar_summary_train <- matrix(NA_real_, nrow = n_obs_train, ncol = 3,
                                  dimnames = list(NULL, c("MEAN","L95","U95")))
    ystar_summary_test  <- matrix(NA_real_, nrow = n_obs_test,  ncol = 3,
                                  dimnames = list(NULL, c("MEAN","L95","U95")))
    
    ystar_summary_train[, "MEAN"] <- yhat_train
    ystar_summary_test[,  "MEAN"] <- yhat_test
    
    ystar_summary_train[, "L95"] <- yhat_train - q975 * rmse_train
    ystar_summary_train[, "U95"] <- yhat_train + q975 * rmse_train
    ystar_summary_test[,  "L95"] <- yhat_test  - q975 * rmse_train
    ystar_summary_test[,  "U95"] <- yhat_test  + q975 * rmse_train
    
    return(list(
      train = list(ystar = ystar_summary_train),
      test  = list(ystar = ystar_summary_test),
      time = train_time,
      train_time = train_time
    ))
  }
  
  # -------- Bootstrap branch  --------
  # npscoef$beta has (p + 1) columns: intercept + p slopes
  boot_beta <- array(NA_real_, dim = c(n_obs_train + n_obs_test, p + 1, B))
  
  boot_time <- system.time({
    for (b in 1:B) {
      if (b %% 5 == 0) message(sprintf("np-boot: b = %d at %s", b, Sys.time()))
      boot_index <- sample.int(n_obs_train, size = n_obs_train, replace = TRUE)
      
      tmp_beta <- tryCatch({
        boot_fit <- npscoef(
          bw,
          betas = TRUE, residuals = FALSE, errors = FALSE,
          txdat = X_train[boot_index, , drop = FALSE],
          tydat = Y_train[boot_index],
          tzdat = Z_train[boot_index, , drop = FALSE],
          exdat = rbind(X_train, X_test),
          ezdat = rbind(Z_train, Z_test)
        )
        boot_fit$beta
      }, error = function(e) NULL)
      
      if (!is.null(tmp_beta)) boot_beta[, , b] <- tmp_beta
    }
  })["elapsed"]
  
  boot_beta_se <- apply(boot_beta, c(1, 2), sd, na.rm = TRUE)
  
  beta_summary_train <- array(NA_real_, dim = c(n_obs_train, 4, p + 1),
                              dimnames = list(NULL, c("MEAN","SD","L95","U95"), NULL))
  beta_summary_test  <- array(NA_real_, dim = c(n_obs_test,  4, p + 1),
                              dimnames = list(NULL, c("MEAN","SD","L95","U95"), NULL))
  
  # Means from the full fit
  beta_summary_train[, "MEAN", ] <- ks_fit$beta[seq_len(n_obs_train), ]
  beta_summary_test[,  "MEAN", ] <- ks_fit$beta[n_obs_train + seq_len(n_obs_test), ]
  
  # SDs from bootstrap
  beta_summary_train[, "SD", ] <- boot_beta_se[seq_len(n_obs_train), ]
  beta_summary_test[,  "SD", ] <- boot_beta_se[n_obs_train + seq_len(n_obs_test), ]
  
  # Wald 95% CIs
  q975 <- qnorm(0.975)
  beta_summary_train[, "L95", ] <- beta_summary_train[, "MEAN", ] - q975 * beta_summary_train[, "SD", ]
  beta_summary_train[, "U95", ] <- beta_summary_train[, "MEAN", ] + q975 * beta_summary_train[, "SD", ]
  beta_summary_test[,  "L95", ] <- beta_summary_test[,  "MEAN", ] - q975 * beta_summary_test[,  "SD", ]
  beta_summary_test[,  "U95", ] <- beta_summary_test[,  "MEAN", ] + q975 * beta_summary_test[,  "SD", ]
  
  # y* summaries on train/test using train RMSE
  rmse_train <- sqrt(mean((Y_train - yhat_train)^2))
  ystar_summary_train <- matrix(NA_real_, nrow = n_obs_train, ncol = 3,
                                dimnames = list(NULL, c("MEAN","L95","U95")))
  ystar_summary_test  <- matrix(NA_real_, nrow = n_obs_test,  ncol = 3,
                                dimnames = list(NULL, c("MEAN","L95","U95")))
  
  ystar_summary_train[, "MEAN"] <- yhat_train
  ystar_summary_test[,  "MEAN"] <- yhat_test
  ystar_summary_train[, "L95"]  <- yhat_train - q975 * rmse_train
  ystar_summary_train[, "U95"]  <- yhat_train + q975 * rmse_train
  ystar_summary_test[,  "L95"]  <- yhat_test  - q975 * rmse_train
  ystar_summary_test[,  "U95"]  <- yhat_test  + q975 * rmse_train
  
  list(
    train = list(beta = beta_summary_train, ystar = ystar_summary_train),
    test  = list(beta = beta_summary_test,  ystar = ystar_summary_test),
    time = train_time + boot_time,
    train_time = train_time,
    boot_time = boot_time
  )
}


###### LM wrapper #####
##########################################################################################################################
lm_wrapper <- function(Y_train, X_train, Z_train, X_test, Z_test) {
  n_obs_train <- nrow(X_train)
  n_obs_test  <- nrow(X_test)
  
  p <- ncol(X_train)
  R <- ncol(Z_train)
  
  # Name columns
  colnames(X_train) <- paste0("X", 1:p)
  colnames(X_test)  <- paste0("X", 1:p)
  colnames(Z_train) <- paste0("Z", 1:R)
  colnames(Z_test)  <- paste0("Z", 1:R)
  
  # Training/test frames for predict()
  tmp_data_train <- data.frame(Y = Y_train, X_train, Z_train)
  tmp_data_test  <- data.frame(           X_test,  Z_test)   # <-- no Y_test
  
  # Fit plain linear model: Y ~ X (Z are ignored but harmless in data)
  form <- reformulate(termlabels = colnames(X_train), response = "Y")
  train_time <- system.time(
    lm_fit <- lm(form, data = tmp_data_train)
  )[["elapsed"]]
  
  # Coeff summary
  lm_coef <- as.data.frame(summary(lm_fit)$coef)  # rows: (Intercept), X1..Xp
  
  # Beta summaries: replicate the same coef/SE per row to match your API
  beta_summary_train <- array(NA_real_, dim = c(n_obs_train, 4, 1 + p),
                              dimnames = list(NULL, c("MEAN","SD","L95","U95"), NULL))
  beta_summary_test  <- array(NA_real_, dim = c(n_obs_test,  4, 1 + p),
                              dimnames = list(NULL, c("MEAN","SD","L95","U95"), NULL))
  
  for (k in 1:(p + 1)) {
    est <- lm_coef[k, "Estimate"]
    se  <- lm_coef[k, "Std. Error"]
    beta_summary_train[, "MEAN", k] <- est
    beta_summary_train[, "SD",   k] <- se
    beta_summary_test[ , "MEAN", k] <- est
    beta_summary_test[ , "SD",   k] <- se
  }
  beta_summary_train[, "L95", ] <- beta_summary_train[, "MEAN", ] - qnorm(0.975) * beta_summary_train[, "SD", ]
  beta_summary_train[, "U95", ] <- beta_summary_train[, "MEAN", ] + qnorm(0.975) * beta_summary_train[, "SD", ]
  beta_summary_test[ , "L95", ] <- beta_summary_test[ , "MEAN", ] - qnorm(0.975) * beta_summary_test[ , "SD", ]
  beta_summary_test[ , "U95", ] <- beta_summary_test[ , "MEAN", ] + qnorm(0.975) * beta_summary_test[ , "SD", ]
  
  # Predictive summaries (fit/lwr/upr) on train/test
  ystar_summary_train <- predict(lm_fit, newdata = tmp_data_train, interval = "prediction")
  ystar_summary_test  <- predict(lm_fit, newdata = tmp_data_test,  interval = "prediction")
  colnames(ystar_summary_train) <- c("MEAN","L95","U95")
  colnames(ystar_summary_test)  <- c("MEAN","L95","U95")
  
  list(
    train = list(beta = beta_summary_train, ystar = ystar_summary_train),
    test  = list(beta = beta_summary_test,  ystar = ystar_summary_test),
    time = train_time, train_time = train_time
  )
}

######### Wrapper function for the TVC ##########
##########################################################################################################################

tvc_wrapper <- function(Y_train, X_train, Z_train, X_test, Z_test, B = 50) {
  
  n_obs_train <- nrow(X_train)
  n_obs_test  <- nrow(X_test)
  
  p <- ncol(X_train)
  R <- ncol(Z_train)
  
  # Name columns
  colnames(X_train) <- paste0("X", 1:p)
  colnames(X_test)  <- paste0("X", 1:p)
  colnames(Z_train) <- paste0("Z", 1:R)
  colnames(Z_test)  <- paste0("Z", 1:R)
  
  # Train has Y; test does NOT (we can add a dummy Y = NA for consistent columns)
  tmp_data_train <- data.frame(Y = Y_train, X_train, Z_train)
  tmp_data_test  <- data.frame(X_test, Z_test)
  tmp_data_test$Y <- NA_real_   # keep same column set/order for rbind during bootstrap
  # reorder to match tmp_data_train col order exactly
  tmp_data_test  <- tmp_data_test[colnames(tmp_data_train)]
  
  # Formula: intercept-free varying coefficients
  form_list <- rep(NA_character_, p + 1)
  form_list[1] <- paste0("vc(", paste0(colnames(Z_train), collapse = ", "), ")")
  for (j in 1:p) {
    form_list[j + 1] <- paste0("vc(",
                               paste0(colnames(Z_train), collapse = ", "),
                               ", by = X", j, ")")
  }
  tvc_form <- as.formula(paste0("Y ~ -1 + ", paste0(form_list, collapse = "+")))
  
  # Initial fit (with CV to pick cp)
  train_time <- system.time(
    tvc_fit <- tryCatch(
      tvcglm(formula = tvc_form,
             family = gaussian(),
             data   = tmp_data_train,
             control = tvcglm_control(cv = TRUE)),
      error = function(e) NULL)
  )["elapsed"]
  
  if (is.null(tvc_fit)) {
    return(NULL)
  }
  
  # Cross-validated pruning parameter
  tvc_cp <- cvloss(tvc_fit)$cp.hat
  
  # Coefs and predictions on train/test
  beta_hat_train <- predict(tvc_fit, newdata = tmp_data_train, type = "coef")
  beta_hat_test  <- predict(tvc_fit, newdata = tmp_data_test,  type = "coef")
  
  ystar_train <- predict(tvc_fit, newdata = tmp_data_train, type = "response")
  ystar_test  <- predict(tvc_fit, newdata = tmp_data_test,  type = "response")
  
  if (B == 0) {
    # Fast intervals based on train RMSE
    rmse_train <- sqrt(mean((Y_train - ystar_train)^2))
    
    ystar_summary_train <- matrix(NA_real_, nrow = n_obs_train, ncol = 3,
                                  dimnames = list(NULL, c("MEAN", "L95", "U95")))
    ystar_summary_test  <- matrix(NA_real_, nrow = n_obs_test,  ncol = 3,
                                  dimnames = list(NULL, c("MEAN", "L95", "U95")))
    ystar_summary_train[, "MEAN"] <- ystar_train
    ystar_summary_test[,  "MEAN"] <- ystar_test
    
    ystar_summary_train[, "L95"] <- ystar_summary_train[, "MEAN"] - qnorm(0.975) * rmse_train
    ystar_summary_train[, "U95"] <- ystar_summary_train[, "MEAN"] + qnorm(0.975) * rmse_train
    ystar_summary_test[,  "L95"] <- ystar_summary_test[,  "MEAN"] - qnorm(0.975) * rmse_train
    ystar_summary_test[,  "U95"] <- ystar_summary_test[,  "MEAN"] + qnorm(0.975) * rmse_train
    
    results <- list(
      train = list(ystar = ystar_summary_train),
      test  = list(ystar = ystar_summary_test),
      time = train_time,
      train_time = train_time
    )
    
  } else {
    # Bootstrap for beta intervals
    boot_beta <- array(NA_real_, dim = c(n_obs_train + n_obs_test, p + 1, B))
    
    boot_time <- system.time({
      for (b in 1:B) {
        if (b %% 5 == 0) message(sprintf("tvc-boot: b = %d at %s", b, Sys.time()))
        boot_index <- sample.int(n_obs_train, size = n_obs_train, replace = TRUE)
        tmp_beta <- tryCatch({
          boot_fit <- tvcglm(formula = tvc_form,
                             family  = gaussian(),
                             data    = tmp_data_train[boot_index, , drop = FALSE],
                             control = tvcglm_control(cv = FALSE))
          boot_pruned <- prune(boot_fit, cp = tvc_cp)
          predict(boot_pruned,
                  newdata = rbind(tmp_data_train, tmp_data_test),
                  type = "coef")
        }, error = function(e) NULL)
        if (!is.null(tmp_beta)) boot_beta[, , b] <- tmp_beta
      }
    })["elapsed"]
    
    boot_beta_se <- apply(boot_beta, c(1, 2), sd, na.rm = TRUE)
    
    beta_summary_train <- array(NA_real_, dim = c(n_obs_train, 4, p + 1),
                                dimnames = list(NULL, c("MEAN", "SD", "L95", "U95"), NULL))
    beta_summary_test  <- array(NA_real_, dim = c(n_obs_test,  4, p + 1),
                                dimnames = list(NULL, c("MEAN", "SD", "L95", "U95"), NULL))
    
    ystar_summary_train <- matrix(NA_real_, nrow = n_obs_train, ncol = 3,
                                  dimnames = list(NULL, c("MEAN", "L95", "U95")))
    ystar_summary_test  <- matrix(NA_real_, nrow = n_obs_test,  ncol = 3,
                                  dimnames = list(NULL, c("MEAN", "L95", "U95")))
    
    beta_summary_train[, "MEAN", ] <- beta_hat_train
    beta_summary_test[,  "MEAN", ] <- beta_hat_test
    
    beta_summary_train[, "SD", ] <- boot_beta_se[1:n_obs_train, ]
    beta_summary_test[,  "SD", ] <- boot_beta_se[(n_obs_train + 1):(n_obs_train + n_obs_test), ]
    
    beta_summary_train[, "L95", ] <- beta_summary_train[, "MEAN", ] - qnorm(0.975) * beta_summary_train[, "SD", ]
    beta_summary_train[, "U95", ] <- beta_summary_train[, "MEAN", ] + qnorm(0.975) * beta_summary_train[, "SD", ]
    beta_summary_test[,  "L95", ] <- beta_summary_test[,  "MEAN", ] - qnorm(0.975) * beta_summary_test[,  "SD", ]
    beta_summary_test[,  "U95", ] <- beta_summary_test[,  "MEAN", ] + qnorm(0.975) * beta_summary_test[,  "SD", ]
    
    rmse_train <- sqrt(mean((Y_train - ystar_train)^2))
    ystar_summary_train[, "MEAN"] <- ystar_train
    ystar_summary_test[,  "MEAN"] <- ystar_test
    
    ystar_summary_train[, "L95"] <- ystar_summary_train[, "MEAN"] - qnorm(0.975) * rmse_train
    ystar_summary_train[, "U95"] <- ystar_summary_train[, "MEAN"] + qnorm(0.975) * rmse_train
    ystar_summary_test[,  "L95"] <- ystar_summary_test[,  "MEAN"] - qnorm(0.975) * rmse_train
    ystar_summary_test[,  "U95"] <- ystar_summary_test[,  "MEAN"] + qnorm(0.975) * rmse_train
    
    results <- list(
      train = list(beta = beta_summary_train, ystar = ystar_summary_train),
      test  = list(beta = beta_summary_test,  ystar = ystar_summary_test),
      time = train_time + boot_time,
      train_time = train_time,
      boot_time = boot_time
    )
  }
  
  return(results)
}
