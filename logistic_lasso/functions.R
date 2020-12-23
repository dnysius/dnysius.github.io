library(tidymodels)
library(tidyverse)

fit_logistic_lasso <- function(x, y, lambda, beta0 = NULL, eps = 0.0001, max_iter = 100) {
  ## fit_logistic_lasso(x, y, lambda, beta0 = NULL, eps = 0.0001, max_iter = 100)
  ## fits the logistic LASSO on data 'y' based on the feature matrix 'x'. To be
  ## used as a parsnip model in a tidymodels workflow.
  ##
  ## Input:
  ## - x:         A numeric matrix of predictors not including the intercept  
  ##              column. Columns should be named.
  ## - y:         A factor vector of data. Values span the classification classes.
  ## - lambda:    A number denoting the L1-penalty of the logistic regression.
  ## - beta0:     A numeric vector of the initial estimates of the parameters.
  ## - eps:       A numeric value specifying the tolerance level for detecting 
  ##              convergence.
  ## - max_iter:  A numeric value which limits the number of iterations to
  ##              wait for convergence. 
  ## Output:
  ## - A list of: the numeric parameter vector 'beta',
  ##              logical for estimate convergence 'passed',
  ##              numeric intercept term 'intercept',
  ##              numeric penalty 'lambda',
  ##              logical for the algorithm success 'success',
  ##              number of coordinate descent iterations 'iter',
  ##              list of class names 'fct_levels'.
  ##   
  ## Example:
  ## library(tidymodels)
  ## library(tidyverse)
  ## n = 1000
  ## dat <- tibble(x = seq(-20,20, length.out = n),
  ##               w = cos(seq(-pi,pi, length.out = n)),
  ##               y = rbinom(n,size = 1, prob = 1/(1 + exp(-w+2*x)) )%>% as.numeric %>% factor
  ## )
  ## split <- initial_split(dat)
  ## train <- training(split)
  ## test <- testing(split)
  ## x<-data.matrix(train[, -3])
  ## y<-train$y
  ## ret <- fit_logistic_lasso(x, y, lambda=1e-10)
  n <- dim(x)[1] # number of observations
  predictors <- dim(x)[2]  # number of predictors (doesn't include intercept)
  X <- data.matrix(cbind("intercept"=rep(1, n), x))  # with intercept column
  if (is.null(beta0)) {  # Check if initial guess is not given
    # initial vector, and used to store the previous beta approximation
    beta_all <- rep(0, (predictors+1))
    beta_all0<-beta_all
    beta0<-beta_all[2:(predictors+1)]  # features only, without intercept
    beta<-beta0
    names(beta) <- colnames(x)
    intercept<-c("intercept"=rep(0, 1))
    intercept0<-intercept
  } else {  # Initial guess beta0 given (includes intercept)
    beta_all<-beta0
    beta_all0<-beta_all
    beta0<-beta_all[2:(predictors+1)]  # features only, without intercept
    beta<-beta0
    names(beta) <- colnames(x)
    intercept<-c("intercept"=rep(0, 1))
    intercept[1]<-beta_all[1]
    intercept0<-intercept
  }
  pass<-as.logical(rep(0, predictors+1))  # Store estimate convergence indicator
  names(pass)<-c(names(intercept), names(beta))
  fct_levels <- levels(y)
  y <- as.numeric(y) - 1
  
  # Parameters calculated with the initial beta0
  X_beta <- (X %*% beta_all) %>% as.numeric
  p <- 1/(1 + exp(-X_beta))
  w <- p * (1 - p)
  z <- X_beta + (y - p)/w
  # Coordinate descent: update each parameter using single predictor IRLS solution
  # then compare the beta vector approximations to see if convergence is reached.
  for(iter in 1:max_iter) {
    # Update step derived from IRLS objective functions
    for(j in 0:predictors) {
      if (j==0){
        # Don't penalize intercept: just the regular IRLS so lm() with weights.
        # first column of x is the 1-vector
        r0 <- z - x %*% beta
        # Solving for just the intercept
        intercept[1] <- (lm(r0 ~ 1, weights = w) %>% coef)[1]
        beta_all[1] <- intercept
      } else{
        if (predictors>2){
          # L1-penalty with non-intercept coefficients
          # Distributing 'w' so we get a nicer form
          rj = sqrt(w)*(z - intercept[1]) - ((sqrt(w)*x[,-j]) %*% beta[-j])
          
        } else{
          rj = sqrt(w)*(z - intercept[1]) - ((sqrt(w)*x[,-j]) * beta[-j])
        }
        # Take only non-negative values; if negative we force it to be 0.
        S <- sign(t(sqrt(w)*x[, j])%*%rj)
        M <- max(2*abs((sqrt(w)*t(x[, j]))%*%rj)-n*lambda, 0)
        # Single predictor solution
        beta[j] <-  S * M / sum(2*(sqrt(w)*x[,j])**2)
        beta_all[-1] <- beta
      }
      # Check each estimate convergence
      pass[j+1]<-max(abs(beta_all[j+1] - beta_all0[j+1])) < eps
    }
    # All beta components have been updated once through already. 
    # Let's update the weights now.
    X_beta <- (X %*% beta_all) %>% as.numeric
    p <- 1/(1 + exp(-X_beta))
    w <- p * (1 - p)
    w[w==0]<-1e-16
    z <- X_beta + (y - p)/w
    # It's time to check convergence by comparing beta to the previous beta 
    # approximation and by checking if the L2-norm of the gradient of beta 
    # is small enough.
    
    # 1. Check gradient convergence for the entire beta vector
    converged_grad <- all(pass==TRUE)==TRUE
    # 2. Check change in beta vector for coordinate descent convergence
    converged_cd <- max(abs(c(intercept-intercept0, beta - beta0))) < eps
    if ((converged_grad*converged_cd)==TRUE){
      return(list(success=TRUE, passed=pass, iter=iter, intercept=intercept, 
                  beta=beta, lambda=lambda, fct_levels=fct_levels))
    }
    # 3. If no convergence, record the current beta as the previous beta 'beta0'
    # and update all the components in beta again.
    # Store the updated vector and intercept term.
    intercept0<-intercept
    beta0<-beta
    beta_all0<-beta_all
 }
  # If we get here we have exceeded max_iter iterations for convergence.
  warning(paste("Convergence took more than", max_iter, "iterations"))
  return(list(success=FALSE, passed=pass, iter=max_iter, intercept=intercept, 
              beta=beta, lambda=lambda, fct_levels=fct_levels))
}

predict_logistic_lasso <- function(fit, new_x) {
  ## predict_logistic_lasso(fit, new_x) creates predicted data values based on
  ## a logistic LASSO fit. To be used as a parsnip model in a tidymodels workflow.
  ##
  ## Input:
  ## - fit:   A list of the numeric parameter vector 'beta',
  ##          logical for estimate convergence 'passed',
  ##          numeric intercept term 'intercept',
  ##          numeric penalty 'lambda',
  ##          logical for the algorithm success 'success',
  ##          number of coordinate descent iterations 'iter',
  ##          list of class names 'fct_levels'.
  ## - new_x: A numeric matrix containing the new data to be predicted.
  ##
  ## Output:
  ## - A factor vector of the predicted classification data.
  ## 
  ## Example:
  ## library(tidymodels)
  ## library(tidyverse)
  ## n = 1000
  ## dat <- tibble(x = seq(-20,20, length.out = n),
  ##               w = cos(seq(-pi,pi, length.out = n)),
  ##               y = rbinom(n,size = 1, prob = 1/(1 + exp(-w+2*x)) )%>% as.numeric %>% factor
  ## )
  ## split <- initial_split(dat)
  ## train <- training(split)
  ## test <- testing(split)
  ## x<-data.matrix(train[, -3])
  ## y<-train$y
  ## ret <- fit_logistic_lasso(x, y, lambda=1e-10)
  ## new_x<-test[, -3]
  ## new_y<-predict_logistic_lasso(ret, new_x)
  ## sum((as.numeric(new_y)-as.numeric(test$y))!=0)  # number of incorrect predictions# n = 1000
  n <- dim(new_x)[1] # number of observations
  new_x<-data.matrix(cbind("intercept"=rep(1, n), new_x))
  beta_all<-c(fit$intercept, fit$beta)
  numeric_pred <- (new_x %*% beta_all >= 0) %>% as.numeric
  return( fit$fct_levels[numeric_pred + 1] %>% factor )
}

