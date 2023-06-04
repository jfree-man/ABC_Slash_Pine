## Local ABC function
## Copy of abc::abc() with MAD scaling removed
################################################################

abc_local <- function (target, param, sumstat, tol, method, hcorr = TRUE, 
          transf = "none", logit.bounds = c(0, 0), subset = NULL, kernel = "epanechnikov", 
          numnet = 10, sizenet = 5, lambda = c(1e-04, 0.001, 0.01), 
          trace = FALSE, maxit = 500, ...) 
{
  call <- match.call()
  if (missing(target)) 
    stop("'target' is missing")
  if (missing(param)) 
    stop("'param' is missing")
  if (missing(sumstat)) 
    stop("'sumstat' is missing")
  if (!is.matrix(param) && !is.data.frame(param) && !is.vector(param)) 
    stop("'param' has to be a matrix, data.frame or vector.")
  if (!is.matrix(sumstat) && !is.data.frame(sumstat) && !is.vector(sumstat)) 
    stop("'sumstat' has to be a matrix, data.frame or vector.")
  if (missing(tol)) 
    stop("'tol' is missing")
  if (missing(method)) 
    stop("'method' is missing with no default")
  if (!any(method == c("rejection", "loclinear", "neuralnet", 
                       "ridge"))) {
    stop("Method must be rejection, loclinear, or neuralnet or ridge")
  }
  if (method == "rejection") 
    rejmethod <- TRUE
  else rejmethod <- FALSE
  if (!any(kernel == c("gaussian", "epanechnikov", "rectangular", 
                       "triangular", "biweight", "cosine"))) {
    kernel <- "epanechnikov"
    warning("Kernel is incorrectly defined. Setting to default (Epanechnikov)")
  }
  if (is.data.frame(param)) 
    param <- as.matrix(param)
  if (is.data.frame(sumstat)) 
    sumstat <- as.matrix(sumstat)
  if (is.list(target)) 
    target <- unlist(target)
  if (is.vector(sumstat)) 
    sumstat <- matrix(sumstat, ncol = 1)
  if (length(target) != dim(sumstat)[2]) 
    stop("Number of summary statistics in 'target' has to be the same as in 'sumstat'.")
  nss <- length(sumstat[1, ])
  cond1 <- !any(as.logical(apply(sumstat, 2, function(x) length(unique(x)) - 
                                   1)))
  if (cond1) 
    stop("Zero variance in the summary statistics.")
  ltransf <- length(transf)
  if (is.vector(param)) {
    numparam <- 1
    param <- matrix(param, ncol = 1)
  }
  else numparam <- dim(param)[2]
  for (i in 1:ltransf) {
    if (sum(transf[i] == c("none", "log", "logit")) == 0) {
      stop("Transformations must be none, log, or logit.")
    }
    if (transf[i] == "logit") {
      if (logit.bounds[i, 1] >= logit.bounds[i, 2]) {
        stop("Logit bounds are incorrect.")
      }
    }
  }
  if (rejmethod) {
    if (!all(transf == "none")) {
      warning("No transformation is applied when the simple rejection is used.", 
              call. = F)
    }
    transf[1:numparam] <- "none"
  }
  else {
    if (numparam != ltransf) {
      if (length(transf) == 1) {
        transf <- rep(transf[1], numparam)
        warning("All parameters are \"", transf[1], "\" transformed.", 
                sep = "", call. = F)
      }
      else stop("Number of parameters is not the same as number of transformations.", 
                sep = "", call. = F)
    }
  }
  gwt <- rep(TRUE, length(sumstat[, 1]))
  gwt[attributes(na.omit(sumstat))$na.action] <- FALSE
  if (missing(subset)) 
    subset <- rep(TRUE, length(sumstat[, 1]))
  gwt <- as.logical(gwt * subset)
  if (!length(colnames(param))) {
    warning("No parameter names are given, using P1, P2, ...")
    paramnames <- paste("P", 1:numparam, sep = "")
  }
  else paramnames <- colnames(param)
  if (!length(colnames(sumstat))) {
    warning("No summary statistics names are given, using S1, S2, ...")
    statnames <- paste("S", 1:nss, sep = "")
  }
  else statnames <- colnames(sumstat)
  scaled.sumstat <- sumstat
  ####################################################################
  ## Remove scaling 
  # for (j in 1:nss) {
  #   scaled.sumstat[, j] <- normalise(sumstat[, j], sumstat[, 
  #                                                          j][gwt])
  # }
  # for (j in 1:nss) {
  #   target[j] <- normalise(target[j], sumstat[, j][gwt])
  # }
  #####################################################################
  sum1 <- 0
  for (j in 1:nss) {
    sum1 <- sum1 + (scaled.sumstat[, j] - target[j])^2
  }
  dist <- sqrt(sum1)
  dist[!gwt] <- floor(max(dist[gwt]) + 10)
  nacc <- ceiling(length(dist) * tol)
  ds <- sort(dist)[nacc]
  wt1 <- (dist <= ds)
  aux <- cumsum(wt1)
  wt1 <- wt1 & (aux <= nacc)
  if (kernel == "gaussian") {
    wt1 <- rep(TRUE, length(dist))
  }
  for (i in 1:numparam) {
    if (transf[i] == "log") {
      if (min(param[, i]) <= 0) {
        cat("log transform: values out of bounds - correcting...")
        x.tmp <- ifelse(param[, i] <= 0, max(param[, 
                                                   i]), param[, i])
        x.tmp.min <- min(x.tmp)
        param[, i] <- ifelse(param[, i] <= 0, x.tmp.min, 
                             param[, i])
      }
      param[, i] <- log(param[, i])
    }
    else if (transf[i] == "logit") {
      if (min(param[, i]) <= logit.bounds[i, 1]) {
        x.tmp <- ifelse(param[, i] <= logit.bounds[i, 
                                                   1], max(param[, i]), param[, i])
        x.tmp.min <- min(x.tmp)
        param[, i] <- ifelse(param[, i] <= logit.bounds[i, 
                                                        1], x.tmp.min, param[, i])
      }
      if (max(param[, i]) >= logit.bounds[i, 2]) {
        x.tmp <- ifelse(param[, i] >= logit.bounds[i, 
                                                   2], min(param[, i]), param[, i])
        x.tmp.max <- max(x.tmp)
        param[, i] <- ifelse(param[, i] >= logit.bounds[i, 
                                                        2], x.tmp.max, param[, i])
      }
      param[, i] <- (param[, i] - logit.bounds[i, 1])/(logit.bounds[i, 
                                                                    2] - logit.bounds[i, 1])
      param[, i] <- log(param[, i]/(1 - param[, i]))
    }
  }
  ss <- sumstat[wt1, ]
  unadj.values <- param[wt1, ]
  statvar <- as.logical(apply(cbind(sumstat[wt1, ]), 2, function(x) length(unique(x)) - 
                                1))
  cond2 <- !any(statvar)
  if (cond2 && !rejmethod) 
    stop("Zero variance in the summary statistics in the selected region. Try: checking summary statistics, choosing larger tolerance, or rejection method.")
  if (rejmethod) {
    if (cond2) 
      warning("Zero variance in the summary statistics in the selected region. Check summary statistics, consider larger tolerance.")
    weights <- rep(1, length = sum(wt1))
    adj.values <- NULL
    residuals <- NULL
    lambda <- NULL
  }
  else {
    if (cond2) 
      cat("Warning messages:\nStatistic(s)", statnames[!statvar], 
          "has/have zero variance in the selected region.\nConsider using larger tolerance or the rejection method or discard this/these statistics, which might solve the collinearity problem in 'lsfit'.\n", 
          sep = ", ")
    if (kernel == "epanechnikov") 
      weights <- 1 - (dist[wt1]/ds)^2
    if (kernel == "rectangular") 
      weights <- dist[wt1]/ds
    if (kernel == "gaussian") 
      weights <- 1/sqrt(2 * pi) * exp(-0.5 * (dist/(ds/2))^2)
    if (kernel == "triangular") 
      weights <- 1 - abs(dist[wt1]/ds)
    if (kernel == "biweight") 
      weights <- (1 - (dist[wt1]/ds)^2)^2
    if (kernel == "cosine") 
      weights <- cos(pi/2 * dist[wt1]/ds)
    if (method == "loclinear") {
      fit1 <- lsfit(scaled.sumstat[wt1, ], param[wt1, ], 
                    wt = weights)
      pred <- t(structure(cbind(fit1$coefficients)[fit1$qr$pivot, 
      ], names = names(fit1$coefficients))) %*% c(1, 
                                                  target)
      pred <- matrix(pred, ncol = numparam, nrow = sum(wt1), 
                     byrow = TRUE)
      residuals <- param[wt1, ] - t(t(structure(cbind(fit1$coefficients)[fit1$qr$pivot, 
      ], names = names(fit1$coefficients))) %*% t(cbind(1, 
                                                        scaled.sumstat[wt1, ])))
      residuals <- cbind(residuals)
      the_m <- apply(residuals, FUN = mean, 2)
      residuals <- sapply(1:numparam, FUN = function(x) {
        residuals[, x] - the_m[x]
      })
      pred <- sapply(1:numparam, FUN = function(x) {
        pred[, x] + the_m[x]
      })
      sigma2 <- apply(as.matrix(residuals), FUN = function(x) {
        sum((x)^2 * weights)/sum(weights)
      }, MARGIN = 2)
      aic <- sum(wt1) * sum(log(sigma2)) + 2 * (nss + 1) * 
        numparam
      bic <- sum(wt1) * sum(log(sigma2)) + log(sum(wt1)) * 
        (nss + 1) * numparam
      if (hcorr == TRUE) {
        fit2 <- lsfit(scaled.sumstat[wt1, ], log(residuals^2), 
                      wt = weights)
        auxaux <- t(structure(cbind(fit2$coefficients)[fit2$qr$pivot, 
        ], names = names(fit2$coefficients))) %*% c(1, 
                                                    target)
        pred.sd <- sqrt(exp(auxaux))
        pred.sd <- matrix(pred.sd, nrow = sum(wt1), ncol = numparam, 
                          byrow = T)
        pred.si <- t(t(structure(cbind(fit2$coefficients)[fit2$qr$pivot, 
        ], names = names(fit2$coefficients))) %*% t(cbind(1, 
                                                          scaled.sumstat[wt1, ])))
        pred.si <- sqrt(exp(pred.si))
        adj.values <- pred + (pred.sd * residuals)/pred.si
        residuals <- (pred.sd * residuals)/pred.si
      }
      else {
        adj.values <- pred + residuals
      }
      colnames(adj.values) <- colnames(unadj.values)
      lambda <- NULL
    }
    if (method == "neuralnet") {
      linout <- TRUE
      param.mad <- c()
      for (i in 1:numparam) {
        param.mad[i] <- mad(param[, i][gwt])
        param[, i] <- normalise(param[, i], param[, i][gwt])
      }
      lambda <- sample(lambda, numnet, replace = T)
      fv <- array(dim = c(sum(wt1), numparam, numnet))
      pred <- matrix(nrow = numparam, ncol = numnet)
      for (i in 1:numnet) {
        fit1 <- nnet(scaled.sumstat[wt1, ], param[wt1, 
        ], weights = weights, decay = lambda[i], size = sizenet, 
        trace = trace, linout = linout, maxit = maxit, 
        ...)
        cat(i)
        fv[, , i] <- fit1$fitted.values
        pred[, i] <- predict(fit1, data.frame(rbind(target)))
      }
      cat("\n")
      pred.med <- apply(pred, 1, median)
      pred.med <- matrix(pred.med, nrow = sum(wt1), ncol = numparam, 
                         byrow = T)
      fitted.values <- apply(fv, c(1, 2), median)
      residuals <- param[wt1, ] - fitted.values
      if (hcorr == TRUE) {
        pred2 <- matrix(nrow = numparam, ncol = numnet)
        fv2 <- array(dim = c(sum(wt1), numparam, numnet))
        for (i in 1:numnet) {
          fit2 <- nnet(scaled.sumstat[wt1, ], log(residuals^2), 
                       weights = weights, decay = lambda[i], size = sizenet, 
                       trace = trace, linout = linout, ...)
          cat(i)
          fv2[, , i] <- fit2$fitted.values
          pred2[, i] <- predict(fit2, data.frame(rbind(target)))
        }
        cat("\n")
        pred.sd <- sqrt(exp(apply(pred2, 1, median)))
        pred.sd <- matrix(pred.sd, nrow = sum(wt1), ncol = numparam, 
                          byrow = T)
        fv.sd <- sqrt(exp(apply(fv2, c(1, 2), median)))
        adj.values <- pred.med + (pred.sd * residuals)/fv.sd
        residuals <- (pred.sd * residuals)/fv.sd
      }
      else {
        adj.values <- pred.med + residuals
      }
      colnames(adj.values) <- colnames(unadj.values)
      for (i in 1:numparam) {
        adj.values[, i] <- adj.values[, i] * param.mad[i]
      }
    }
    if (method == "ridge") {
      param.mad <- c()
      for (i in 1:numparam) {
        param.mad[i] <- mad(param[, i][gwt])
        param[, i] <- normalise(param[, i], param[, i][gwt])
      }
      numnet <- length(lambda)
      fv <- array(dim = c(sum(wt1), numparam, numnet))
      pred <- matrix(nrow = numparam, ncol = numnet)
      mataux <- sqrt(diag(weights))
      paramaux <- as.matrix(mataux %*% param[wt1, ])
      scaledaux <- mataux %*% scaled.sumstat[wt1, ]
      for (parcur in (1:numparam)) {
        fit1 <- lm.ridge(paramaux[, parcur] ~ scaledaux, 
                         lambda = lambda)
        for (i in 1:numnet) {
          fv[, parcur, i] <- drop(cbind(1, scaled.sumstat[wt1, 
          ]) %*% (rbind(coef(fit1))[i, ]))
          pred[parcur, i] <- drop(c(1, target) %*% (rbind(coef(fit1))[i, 
          ]))
        }
      }
      pred.med <- apply(pred, 1, median)
      pred.med <- matrix(pred.med, nrow = sum(wt1), ncol = numparam, 
                         byrow = T)
      fitted.values <- apply(fv, c(1, 2), median)
      residuals <- param[wt1, ] - fitted.values
      if (hcorr == TRUE) {
        pred2 <- matrix(nrow = numparam, ncol = numnet)
        fv2 <- array(dim = c(sum(wt1), numparam, numnet))
        for (parcur in (1:numparam)) {
          lresidaux <- (mataux %*% (log(residuals[, parcur]^2)))
          fit2 <- lm.ridge(lresidaux ~ scaledaux, lambda = lambda)
          for (i in 1:numnet) {
            fv2[, parcur, i] <- drop(cbind(1, scaled.sumstat[wt1, 
            ]) %*% (rbind(coef(fit2))[i, ]))
            pred2[parcur, i] <- drop(c(1, target) %*% 
                                       (rbind(coef(fit2))[i, ]))
          }
        }
        cat("\n")
        pred.sd <- sqrt(exp(apply(pred2, 1, median)))
        pred.sd <- matrix(pred.sd, nrow = sum(wt1), ncol = numparam, 
                          byrow = T)
        fv.sd <- sqrt(exp(apply(fv2, c(1, 2), median)))
        adj.values <- pred.med + (pred.sd * residuals)/fv.sd
        residuals <- (pred.sd * residuals)/fv.sd
      }
      else {
        adj.values <- pred.med + residuals
      }
      colnames(adj.values) <- colnames(unadj.values)
      for (i in 1:numparam) {
        adj.values[, i] <- adj.values[, i] * param.mad[i]
      }
    }
  }
  if (numparam == 1) {
    unadj.values <- matrix(unadj.values, ncol = 1)
    if (method != "rejection") {
      adj.values <- matrix(adj.values, ncol = 1)
      residuals <- matrix(residuals, ncol = 1)
    }
  }
  for (i in 1:numparam) {
    if (transf[i] == "log") {
      unadj.values[, i] <- exp(unadj.values[, i])
      adj.values[, i] <- exp(adj.values[, i])
    }
    else if (transf[i] == "logit") {
      unadj.values[, i] <- exp(unadj.values[, i])/(1 + 
                                                     exp(unadj.values[, i]))
      unadj.values[, i] <- unadj.values[, i] * (logit.bounds[i, 
                                                             2] - logit.bounds[i, 1]) + logit.bounds[i, 1]
      adj.values[, i] <- exp(adj.values[, i])/(1 + exp(adj.values[, 
                                                                  i]))
      adj.values[, i] <- adj.values[, i] * (logit.bounds[i, 
                                                         2] - logit.bounds[i, 1]) + logit.bounds[i, 1]
    }
  }
  abc:::abc.return(transf, logit.bounds, method, call, numparam, 
             nss, paramnames, statnames, unadj.values, adj.values, 
             ss, weights, residuals, dist, wt1, gwt, lambda, hcorr, 
             aic, bic)
}
