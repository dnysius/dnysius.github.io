logistic_lasso <- function(mode = "classification", penalty) {
  args <- list(penalty = rlang::enquo(penalty))
  new_model_spec("logistic_lasso",
                 args = args,
                 mode = mode,
                 eng_args = NULL,
                 method = NULL,
                 engine = NULL)
}
set_new_model("logistic_lasso")
set_model_mode(model = "logistic_lasso", mode = "classification")
set_model_engine("logistic_lasso",
                 mode = "classification",
                 eng = "fit_logistic_lasso"
)
set_dependency("logistic_lasso", eng = "fit_logistic_lasso", pkg = "base")
set_model_arg(
  model = "logistic_lasso",
  eng = "fit_logistic_lasso",
  parsnip = "penalty", ## what parsnip will call it
  original = "lambda", ## what we call it!
  func = list(pkg = "dials", fun = "penalty"), ## Use dials::penalty() to set
  has_submodel = FALSE
)
set_encoding(
  model = "logistic_lasso",
  eng = "fit_logistic_lasso",
  mode = "classification",
  options = list(
    predictor_indicators = "traditional",
    compute_intercept = TRUE,
    remove_intercept = TRUE,
    allow_sparse_x = FALSE
  )
)

set_fit(
  model = "logistic_lasso",
  eng = "fit_logistic_lasso",
  mode = "classification",
  value = list(
    interface = "matrix",
    protect = c("x", "y"),
    func = c(fun = "fit_logistic_lasso"),
    defaults = list()
  )
)
set_pred(
  model = "logistic_lasso",
  eng = "fit_logistic_lasso",
  mode = "classification",
  type = "class",
  value = list(
    pre = NULL,
    post = NULL,
    func = c(fun = "predict_logistic_lasso"),
    args = list(
      fit = expr(object$fit),
      new_x = expr(data.matrix(new_data[, names(object$fit$beta)]))
    )
  )
)

pred_logistic_lasso_prob <- function(fit, new_x) {
  # pred_logistic_lasso_prob(fit, new_x) specifies the probability function
  # for the parsnip model.
  n<-dim(new_x)[1]
  new_x<-data.matrix(cbind("intercept"=rep(1, n), new_x))
  beta_all <- c(fit$intercept, fit$beta)
  logit_p = (new_x %*% beta_all) %>% as.numeric
  return( 1/(1 + exp(-logit_p) ))
}

set_pred(
  model = "logistic_lasso",
  eng = "fit_logistic_lasso", 
  mode = "classification",
  type = "prob",
  value = list(
    pre = NULL,
    post = NULL,
    func = c(fun = "pred_logistic_lasso_prob"),
    args = list(
      fit = expr(object$fit),
      new_x = expr(data.matrix(new_data[, names(object$fit$beta)]))
    )
  )
)

update.logistic_lasso <- function(object, penalty=NULL, ...) {
  # update.logistic_lasso(object, penalty=NULL, ...) creates a new spec
  # that has the final parameter.
  if(!is.null(penalty)) {
    object$args <- list(penalty = enquo(penalty))
  }
  new_model_spec("logistic_lasso", args = object$args, eng_args = NULL,
                 mode = "classification", method = NULL, engine = object$engine)
}
show_model_info("logistic_lasso")
