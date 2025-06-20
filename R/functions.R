# library(survival)
# library(dplyr)
# library(boot)
# library(ggplot2)
# library(progress)
#' Multiply robust estimator for calculating the principal survival causal effects among always takers, compliers, and never takers
#'
#' @param times a vector of time when the principal survival causal effects (PSCEs) are of interest
#' @param data the dataset
#' @param Xpi_names names of the covariates for the propensity score model
#' @param Xe_names names of the covariates for the principal score model
#' @param Xc_names names of the covariates for the censoring model
#' @param Xt_names names of the covariates for the failure outcome model
#' @param Z_name name of the treatment assignment status
#' @param S_name name of the true treatment receipt status
#' @param U_name name of the observed failuture time
#' @param delta_name name of the censoring indicator
#' @param B number of the bootstrap replications (default 100)
#'
#' @returns The PSCE estimates and their 95\% confidence intervals
#'
#' @examples
#' # example code
#' attach(sim_data)
#' res = mrPStrata(times=c(1,2,3,4,5,6,7,8),
#'                 data = sim_data,
#'                 Xpi_names = c("X1","X2","X3","X4","X5"),
#'                 Xe_names = c("X1","X2","X3","X4","X5"),
#'                 Xc_names = c("X1","X2","X3","X4","X5"),
#'                 Xt_names = c("X1","X2","X3","X4","X5"),
#'                 Z_name = "z",
#'                 S_name = "s",
#'                 U_name ="U",
#'                 delta_name = "delta",
#'                 B=50)
#' print(res)
mrPStrata = function(times, data, Xpi_names, Xe_names, Xc_names, Xt_names, Z_name, S_name, U_name, delta_name, B=100) {
  # times=c(1,2,3,4)
  # data = sim_data
  # Xpi_names = c("X1","X2","X3","X4","X5")
  # Xe_names = c("X1","X2","X3","X4","X5")
  # Xc_names = c("X1","X2","X3","X4","X5")
  # Xt_names = c("X1","X2","X3","X4","X5")
  # Z_name = "z"
  # S_name = "s"
  # U_name = "U"
  # delta_name = "delta"
  # B=100
  # clean dataset and change the colnames
  names(data)[names(data) == Z_name] <- "z"
  names(data)[names(data) == S_name] <- "s"
  names(data)[names(data) == U_name] <- "U"
  names(data)[names(data) == delta_name] <- "delta"
  # specify the working models
  propensity.model = as.formula(paste("z~",paste(Xpi_names,collapse="+"),sep=""))
  principal.model = as.formula(paste("s~",paste(Xpi_names,collapse="+"),sep=""))
  censor.model = as.formula(paste("delta~",paste(Xpi_names,collapse="+"),sep=""))
  failure.model = as.formula(paste("U~",paste(Xpi_names,collapse="+"),sep=""))
  # number of bootstrap replications
  bootstrap = B
  # multiply robust estimates (point + bootstrap samples)
  mr_estimates = BootEst.PI.SA(times,
                              propensity.model,
                              principal.model,
                              censor.model,
                              failure.model,
                              data,
                              xi0=0, xi1=0, eta0=1, eta1=1,
                              estimand = c("NACE", "CACE", "AACE"),
                              bootstrap)
  bound_vec <- function(x) {x[x < 0] <- 0;x[x > 1] <- 1; return(x)}
  # organize the point estimation and 95\% confidence interval
  # compliers
  S_c_1 = data.frame(time=times,
                     point = mr_estimates$estimate$CACE$z1_mono %>% bound_vec,
                     CI_lower = apply(mr_estimates$bootstrap$CACE$z1_mono.B,1,quantile, 0.025, na.rm = T) %>% bound_vec,
                     CI_upper = apply(mr_estimates$bootstrap$CACE$z1_mono.B,1,quantile, 0.975, na.rm = T) %>% bound_vec)
  S_c_0 = data.frame(time=times,
                     point = mr_estimates$estimate$CACE$z0_mono %>% bound_vec,
                     CI_lower = apply(mr_estimates$bootstrap$CACE$z0_mono.B,1,quantile, 0.025, na.rm = T) %>% bound_vec,
                     CI_upper = apply(mr_estimates$bootstrap$CACE$z0_mono.B,1,quantile, 0.975, na.rm = T) %>% bound_vec)
  PSCE_c = data.frame(time=times,
                      point = mr_estimates$estimate$CACE$diff_mono,
                      CI_lower = apply(mr_estimates$bootstrap$CACE$diff_mono.B,1,quantile, 0.025, na.rm = T),
                      CI_upper = apply(mr_estimates$bootstrap$CACE$diff_mono.B,1,quantile, 0.975, na.rm = T))
  # always takers
  S_a_1 = data.frame(time=times,
                     point = mr_estimates$estimate$AACE$z1_mono %>% bound_vec,
                     CI_lower = apply(mr_estimates$bootstrap$AACE$z1_mono.B,1,quantile, 0.025, na.rm = T) %>% bound_vec,
                     CI_upper = apply(mr_estimates$bootstrap$AACE$z1_mono.B,1,quantile, 0.975, na.rm = T) %>% bound_vec)
  S_a_0 = data.frame(time=times,
                     point = mr_estimates$estimate$AACE$z0_mono %>% bound_vec,
                     CI_lower = apply(mr_estimates$bootstrap$AACE$z0_mono.B,1,quantile, 0.025, na.rm = T) %>% bound_vec,
                     CI_upper = apply(mr_estimates$bootstrap$AACE$z0_mono.B,1,quantile, 0.975, na.rm = T) %>% bound_vec)
  PSCE_a = data.frame(time=times,
                      point = mr_estimates$estimate$AACE$diff_mono,
                      CI_lower = apply(mr_estimates$bootstrap$AACE$diff_mono.B,1,quantile, 0.025, na.rm = T),
                      CI_upper = apply(mr_estimates$bootstrap$AACE$diff_mono.B,1,quantile, 0.975, na.rm = T))
  # never takers
  S_n_1 = data.frame(time=times,
                     point = mr_estimates$estimate$NACE$z1_mono %>% bound_vec,
                     CI_lower = apply(mr_estimates$bootstrap$NACE$z1_mono.B,1,quantile, 0.025, na.rm = T) %>% bound_vec,
                     CI_upper = apply(mr_estimates$bootstrap$NACE$z1_mono.B,1,quantile, 0.975, na.rm = T) %>% bound_vec)
  S_n_0 = data.frame(time=times,
                     point = mr_estimates$estimate$NACE$z0_mono %>% bound_vec,
                     CI_lower = apply(mr_estimates$bootstrap$NACE$z0_mono.B,1,quantile, 0.025, na.rm = T) %>% bound_vec,
                     CI_upper = apply(mr_estimates$bootstrap$NACE$z0_mono.B,1,quantile, 0.975, na.rm = T) %>% bound_vec)
  PSCE_n = data.frame(time=times,
                      point = mr_estimates$estimate$NACE$diff_mono,
                      CI_lower = apply(mr_estimates$bootstrap$NACE$diff_mono.B,1,quantile, 0.025, na.rm = T),
                      CI_upper = apply(mr_estimates$bootstrap$NACE$diff_mono.B,1,quantile, 0.975, na.rm = T))
  out = list(Always_Takers = list(S_a_1=S_a_1,S_a_0=S_a_0,PSCE_a=PSCE_a),
             Compliers = list(S_c_1=S_c_1,S_c_0=S_c_0,PSCE_c=PSCE_c),
             Never_Takers = list(S_n_1=S_n_1,S_n_0=S_n_0,PSCE_n=PSCE_n))
  class(out) = "PSCE.MAIN"
  return(out)
}




#' Point estimation for the multiply robust estimator
#'
#' @param times a vector of time when the principal survival causal effects (PSCEs) are of interest
#' @param propensity.model propensity score model
#' @param principal.model principal score model
#' @param censor.model censoring model
#' @param failure.model outcome model
#' @param data dataset
#' @param xi0 sensitivity parameter xi_0
#' @param xi1 sensitivity parameter xi_1
#' @param eta0 sensitivity parameter eta_0
#' @param eta1 sensitivity parameter eta_1
#' @param estimand the estimands of interest
#'
#' @returns The PSCE point estimates
PointEst.PI.SA <- function(times, propensity.model, principal.model, censor.model, failure.model, data,
                      xi0, xi1, eta0, eta1,
                      estimand = c("NACE", "CACE", "AACE")) {
  tm <- Sys.time()
  lower_bound_vec <- function(vals, seq) {
    # Find the position of the first element that is not less than
    # a given value, or length(inc_seq) + 1 if no such element is found.
    # Time complexity: O(m + n) - often faster than O(n log m) by
    # doing elementwise binary search.
    result <- integer(length(vals))
    p <- 1
    q <- 1
    while (p <= length(vals)) {
      while (q <= length(seq) && seq[q] < vals[p])
        q <- q + 1
      result[p] <- q
      p <- p + 1
    }
    return (result)
  }

  upper_bound_vec <- function(vals, seq) {
    # Find the position of the first element that is greater than
    # a given value, or length(inc_seq) + 1 if no such element is found.
    # Time complexity: O(m + n) - often faster than O(n log m) by
    # applying the univariate version of upper_bound
    result <- integer(length(vals))
    p <- 1
    q <- 1
    while (p <= length(vals)) {
      while (q <= length(seq) && seq[q] <= vals[p])
        q <- q + 1
      result[p] <- q
      p <- p + 1
    }
    return (result)
  }

  mono <- function(seq){
    cur <- Inf
    for (i in 1:length(seq)) {
      seq[i] <- min(cur, seq[i])
      cur <- min(cur, seq[i])
    }
    return (seq)
  }

  Z_name <- propensity.model[[2]]
  S_name <- principal.model[[2]]
  Delta_name <- censor.model[[2]]
  U_name <- failure.model[[2]]

  data <- data %>% arrange({{U_name}})
  n <- nrow(data)
  z <- data %>% pull({{Z_name}})
  s <- data %>% pull({{S_name}})
  U <- data %>% pull({{U_name}})
  delta <- data %>% pull({{Delta_name}})
  XC <- model.matrix(censor.model, data)[, -1]
  XT <- model.matrix(failure.model, data)[, -1]

  fit.tp <- glm(propensity.model, family=binomial(link = "logit"), data = data)
  tp_est <- fit.tp$fitted.values  #estimation of tp
  fit.p1 <- glm(principal.model, family = binomial(link = "logit"), data = data[z == 1, ])
  p1_est <- predict(fit.p1, data, type = "response")
  fit.p0 <- glm(principal.model, family = binomial(link = "logit"), data = data[z == 0, ])
  p0_est <- predict(fit.p0, data, type = "response")

  ps01_est <<- pmax(p1_est - p0_est,0)
  ps00_est <<- 1-p1_est
  ps11_est <<- p0_est
  p1.mean <- mean(z*(s-p1_est)/tp_est+p1_est)
  p0.mean <- mean((1-z)*(s-p0_est)/(1-tp_est)+p0_est)
  w00_est <- (1 - p1_est) / (1 - p1.mean)
  w01_est <- (p1_est - p0_est)/(p1.mean - p0.mean)
  w11_est <- p0_est / p0.mean

  eps1 <- exp(xi1 * (times/max(times))^eta1)
  eps0 <- exp(xi0 * (times/max(times))^eta0)
  w1c <- sapply(eps1, function(x) (x * ps01_est + x * ps11_est) / (x * ps01_est + ps11_est))
  w0c <- sapply(eps0, function(x) (x * ps01_est + x * ps00_est) / (x * ps01_est + ps00_est))
  w0n <- sapply(eps0, function(x) (ps01_est + ps00_est) / (x * ps01_est + ps00_est))
  w1a <- sapply(eps1, function(x) (ps01_est + ps11_est) / (x * ps01_est + ps11_est))

  pos1 <- pmin(lower_bound_vec(times, U), n)
  pos2 <- pmax(upper_bound_vec(times, U) - 1, 1)
  time_id <- 1:length(times)

  subsets <- list(
    "00" = s == 0 & z == 0,
    "01" = s == 1 & z == 0,
    "10" = s == 0 & z == 1,
    "11" = s == 1 & z == 1
  )

  call_censor <- call(
    "~",
    as.call(list(
      quote(Surv), U_name, as.call(list(quote(`-`), 1, Delta_name))
    )),
    censor.model[[3]]
  )

  call_failure <- call(
    "~",
    as.call(list(
      quote(Surv), U_name, Delta_name
    )),
    failure.model[[3]]
  )

  calc_cumhaz <- function(ss) {

    formula_censor <- eval(call_censor)
    formula_failure <- eval(call_failure)

    fit <- list(
      T = coxph(formula_failure, data = data[ss, ]),
      C = coxph(formula_censor, data = data[ss, ])
    )
    tmp_func <- function(fit, X, delta) {
      # basehaz <- basehaz(fit, centered = F)
      # basehaz1 <- basehaz$hazard[pmax(upper_bound_vec(U, basehaz$time) - 1L, 1L)]
      # cumhaz <- outer(c(exp(X %*% fit$coefficients)), basehaz1)
      b_est <- fit$coefficients
      s1 <- cumsum(ss[n:1]*exp(as.matrix(X[n:1,])%*%b_est))
      dls <- ss*delta/s1[n:1]
      dls[(1:length(dls))[is.nan(dls)]] <- 0
      cumhaz <- outer(c(exp(X %*% fit$coefficients)), cumsum(dls))
      return (cumhaz)
    }
    return (list(
      T = tmp_func(fit$T, XT, delta),
      C = tmp_func(fit$C, XC, 1-delta)
    ))
  }

  calc_q <- function(cumhaz) {
    q <- diag(NA, n)
    q[,1] = -cumhaz$C[,1] * exp(cumhaz$C[,1] + cumhaz$T[,1])

    for(j in 2:n) {
      q[j:n,j] = q[j:n,j-1] - (cumhaz$C[j:n,j] - cumhaz$C[j:n,j-1]) * exp(cumhaz$C[j:n,j] + cumhaz$T[j:n,j])
    }
    for(i in 1:(n-1)){
      q[i,(i+1):n] = q[i,i] + (1 - delta[i]) * exp(cumhaz$C[i, i] + cumhaz$T[i, i])
    }
    return (q)
  }

  subset_ind <- c(
    "NACE" %in% estimand || "CACE" %in% estimand, # 00
    "AACE" %in% estimand, # 01
    "NACE" %in% estimand, # 10
    "CACE" %in% estimand || "AACE" %in% estimand # 11
  )

  subsets_to_use <- subsets[subset_ind]

  # calculate q only for needed

  cumhaz <- lapply(subsets_to_use, calc_cumhaz)
  q <- lapply(cumhaz, calc_q)
  gc()

  result <- list()

  gen_result <- function(func1, func0) {
    est1 <- mapply(function(x, y, t) func1(x, y, t) / n, pos1, pos2, time_id)
    est0 <- mapply(function(x, y, t) func0(x, y, t) / n, pos1, pos2, time_id)
    z1_mono = mono(est1)
    z0_mono = mono(est0)

    return (list(
      z1 = est1,
      z0 = est0,
      diff = est1 - est0,
      z1_mono = z1_mono,
      z0_mono = z0_mono,
      diff_mono = z1_mono - z0_mono
    ))
  }

  NACE_1 <- function(x, y, t) {
    sum(
      (1 - s[x:n]) * (z[x:n]) / ((1 - p1.mean) * (tp_est[x:n]) * exp(-cumhaz[["10"]]$C[x:n, y]))
    ) + sum(
      (1 - p1_est) / (1 - p1.mean) * (1 - z / tp_est) * exp(-cumhaz[["10"]]$T[, y])
    ) + sum(
      (1 - s) * z * exp(-cumhaz[["10"]]$T[, y]) * q[["10"]][, x] / ((1 - p1.mean) * tp_est)
    )
  }

  NACE_0 <- function(x, y, t) {
    s1 <- w00_est * (1 - s) * (1 - z)
    s2 <- (1 - p0_est) * (1 - tp_est)
    sum(
      w0n[x:n ,t] * s1[x:n] / (s2[x:n] * exp(-cumhaz[["00"]]$C[x:n, y]))
    ) + sum(
      w0n[, t] * w00_est * (1 - (1 - z) / (1 - tp_est)) * exp(-cumhaz[["00"]]$T[, y])
    ) - sum(
      w0n[, t]^2 / eps0[t] * (z * (s - p1_est) / tp_est + (1 - p1_est) * (1 - z) * (p0_est - s) / s2) * exp(-cumhaz[["00"]]$T[, y])
    ) / (1 - p1.mean) + sum(
      w0n[, t] * s1 * exp(-cumhaz[["00"]]$T[, y]) * q[["00"]][, x] / s2
    )
  }

  CACE_1 <- function(x, y, t) {
    s1 <- w01_est * s * z
    s2 <- p1_est * tp_est
    sum(
      w1c[x:n, t] * s1[x:n] / (s2[x:n] * exp(-cumhaz[["11"]]$C[x:n, y]))
    ) + sum(
      w1c[, t] * w01_est * (1 - z / tp_est) * exp(-cumhaz[["11"]]$T[, y])
    ) - sum(
      w1c[, t]^2 / eps1[t] *((1 - z) * (s - p0_est) / (1 - tp_est) + p0_est * z * (p1_est - s) / s2) * exp(-cumhaz[["11"]]$T[, y])
    ) / (p1.mean - p0.mean) + sum(
      w1c[, t] * s1 * exp(-cumhaz[["11"]]$T[, y]) * q[["11"]][, x] / s2
    )
  }

  CACE_0 <-  function(x, y, t) {
    s1 <- w01_est * (1 - s) * (1 - z)
    s2 <- (1 - p0_est) * (1 - tp_est)
    sum(
      w0c[x:n, t] * s1[x:n] / (s2[x:n] * exp(-cumhaz[["00"]]$C[x:n, y]))
    ) + sum(
      w0c[, t] * w01_est * (1 - (1 - z) / (1 - tp_est)) * exp(-cumhaz[["00"]]$T[, y])
    ) + sum(
      w0c[, t]^2 / eps0[t] * (z * (s - p1_est) / tp_est + (1 - p1_est) * (1 - z) * (p0_est - s) / s2) * exp(-cumhaz[["00"]]$T[, y])
    ) / (p1.mean - p0.mean) + sum(
      w0c[, t] * s1 * exp(-cumhaz[["00"]]$T[, y]) * q[["00"]][, x] / s2
    )
  }

  AACE_1 <-  function(x, y, t) {
    s1 <- w11_est * s * z
    s2 <- p1_est * tp_est
    sum(
      w1a[x:n, t] * s1[x:n] / (s2[x:n] * exp(-cumhaz[["11"]]$C[x:n, y]))
    ) + sum(
      w1a[, t] * w11_est * (1 - z / tp_est) * exp(-cumhaz[["11"]]$T[, y])
    ) + sum(
      w1a[, t]^2 / eps1[t] * ((1 - z) * (s - p0_est) / (1 - tp_est) + p0_est * z * (p1_est - s) / s2) * exp(-cumhaz[["11"]]$T[, y])
    ) / (p0.mean) + sum(
      w1a[, t] * s1 * exp(-cumhaz[["11"]]$T[, y]) * q[["11"]][, x] / s2
    )
  }

  AACE_0 <- function(x, y, t) {
    sum(
      s[x:n] * (1 - z[x:n]) / (p0.mean * (1 - tp_est[x:n]) * exp(-cumhaz[["01"]]$C[x:n, y]))
    ) + sum(
      (p0_est) / (p0.mean) * (1 - (1 - z) / (1 - tp_est)) * exp(-cumhaz[["01"]]$T[, y])
    ) + sum(
      s * (1 - z) * exp(-cumhaz[["01"]]$T[, y]) * q[["01"]][, x] / (p0.mean * (1 - tp_est))
    )
  }

  if ("NACE" %in% estimand)
    result$NACE <- gen_result(NACE_1, NACE_0)
  if ("CACE" %in% estimand)
    result$CACE <- gen_result(CACE_1, CACE_0)
  if ("AACE" %in% estimand)
    result$AACE <- gen_result(AACE_1, AACE_0)
  class(result) <- "est.MR.SA"
  return (result)
}

#' Bootstrap confidence intervals for the multiply robust estimator
#'
#' @param times a vector of time when the principal survival causal effects (PSCEs) are of interest
#' @param propensity.model propensity score model
#' @param principal.model principal score model
#' @param censor.model censoring model
#' @param failure.model outcome model
#' @param data dataset
#' @param xi0 sensitivity parameter xi_0
#' @param xi1 sensitivity parameter xi_1
#' @param eta0 sensitivity parameter eta_0
#' @param eta1 sensitivity parameter eta_1
#' @param estimand the estimands of interest
#' @param bootstrap number of bootstrap replications
#'
#' @returns The bootstrap confidence intervals
BootEst.PI.SA <- function(times, propensity.model, principal.model, censor.model, failure.model, data,
                   xi0, xi1, eta0, eta1,
                   estimand = c("NACE", "CACE", "AACE"), bootstrap) {
  mean_result <- PointEst.PI.SA(times, propensity.model, principal.model, censor.model, failure.model, data,
                           xi0, xi1, eta0, eta1, estimand)
  if (bootstrap == 0)
    return (mean_result)
  bs_result <- c()
  for (name in estimand) {
    z1.B <- z0.B <- diff.B <- z1_mono.B <- z0_mono.B <- diff_mono.B <- matrix(nrow = length(times), ncol = bootstrap)
    bs_result[[name]] <- list(
      z1.B = z1.B,
      z0.B = z0.B,
      diff.B = diff.B,
      z1_mono.B = z1_mono.B,
      z0_mono.B = z0_mono.B,
      diff_mono.B = diff_mono.B
    )
  }
  pb <- progress::progress_bar$new(total = bootstrap, format = "(:spin) [:bar] :percent. :elapsed | :eta")
  for (i in 1:bootstrap) {
    gc()
    pb$tick()
    data.B <- slice_sample(data, n = nrow(data), replace = T)
    result <- PointEst.PI.SA(times, propensity.model, principal.model, censor.model, failure.model, data.B,
                             xi0, xi1, eta0, eta1, estimand)
    for (name in estimand) {
      bs_result[[name]]$z1.B[, i] <- result[[name]]$z1
      bs_result[[name]]$z0.B[, i] <- result[[name]]$z0
      bs_result[[name]]$diff.B[, i] <- result[[name]]$diff
      bs_result[[name]]$z1_mono.B[, i] <- result[[name]]$z1_mono
      bs_result[[name]]$z0_mono.B[, i] <- result[[name]]$z0_mono
      bs_result[[name]]$diff_mono.B[, i] <- result[[name]]$diff_mono
    }
  }
  r_val <- list(
    estimate = mean_result,
    bootstrap = bs_result
  )
  class(r_val) <- "est.MR.SA.bootstrap"
  return (r_val)
}


#' Bias-corrected multiply robust estimator of the PSCE under violation of the principal ignorability assumption
#'
#' @param times a vector of time when the principal survival causal effects (PSCEs) are of interest
#' @param data the dataset
#' @param Xpi_names names of the covariates for the propensity score model
#' @param Xe_names names of the covariates for the principal score model
#' @param Xc_names names of the covariates for the censoring model
#' @param Xt_names names of the covariates for the outcome model
#' @param Z_name name of the treatment assignment status
#' @param S_name name of the true treatment receipt status
#' @param U_name name of the observed failuture time
#' @param delta_name name of the censoring indicator
#' @param xi0 sensitivity parameter xi_0
#' @param xi1 sensitivity parameter xi_1
#' @param eta0 sensitivity parameter eta_0
#' @param eta1 sensitivity parameter eta_1
#' @param B number of the bootstrap replications (default 100)
#'
#' @returns The PSCE estimates and their 95\% confidence intervals
mrPStrata_PI_SA = function(times, data, Xpi_names, Xe_names, Xc_names, Xt_names, Z_name, S_name, U_name, delta_name, xi0=0,xi1=0,eta0=1,eta1=1, B=100) {

  # clean dataset and change the colnames
  names(data)[names(data) == Z_name] <- "z"
  names(data)[names(data) == S_name] <- "s"
  names(data)[names(data) == U_name] <- "U"
  names(data)[names(data) == delta_name] <- "delta"
  # specify the working models
  propensity.model = as.formula(paste("z~",paste(Xpi_names,collapse="+"),sep=""))
  principal.model = as.formula(paste("s~",paste(Xpi_names,collapse="+"),sep=""))
  censor.model = as.formula(paste("delta~",paste(Xpi_names,collapse="+"),sep=""))
  failure.model = as.formula(paste("U~",paste(Xpi_names,collapse="+"),sep=""))
  # number of bootstrap replications
  bootstrap = B
  # multiply robust estimates (point + bootstrap samples)
  mr_estimates = BootEst.PI.SA(times,
                               propensity.model,
                               principal.model,
                               censor.model,
                               failure.model,
                               data,
                               xi0, xi1, eta0, eta1,
                               estimand = c("NACE", "CACE", "AACE"),
                               bootstrap)
  bound_vec <- function(x) {x[x < 0] <- 0;x[x > 1] <- 1; return(x)}
  # organize the point estimation and 95\% confidence interval
  # compliers
  S_c_1 = data.frame(time=times,
                     point = mr_estimates$estimate$CACE$z1_mono %>% bound_vec,
                     CI_lower = apply(mr_estimates$bootstrap$CACE$z1_mono.B,1,quantile, 0.025, na.rm = T) %>% bound_vec,
                     CI_upper = apply(mr_estimates$bootstrap$CACE$z1_mono.B,1,quantile, 0.975, na.rm = T) %>% bound_vec)
  S_c_0 = data.frame(time=times,
                     point = mr_estimates$estimate$CACE$z0_mono %>% bound_vec,
                     CI_lower = apply(mr_estimates$bootstrap$CACE$z0_mono.B,1,quantile, 0.025, na.rm = T) %>% bound_vec,
                     CI_upper = apply(mr_estimates$bootstrap$CACE$z0_mono.B,1,quantile, 0.975, na.rm = T) %>% bound_vec)
  PSCE_c = data.frame(time=times,
                      point = mr_estimates$estimate$CACE$diff_mono,
                      CI_lower = apply(mr_estimates$bootstrap$CACE$diff_mono.B,1,quantile, 0.025, na.rm = T),
                      CI_upper = apply(mr_estimates$bootstrap$CACE$diff_mono.B,1,quantile, 0.975, na.rm = T))
  # always takers
  S_a_1 = data.frame(time=times,
                     point = mr_estimates$estimate$AACE$z1_mono %>% bound_vec,
                     CI_lower = apply(mr_estimates$bootstrap$AACE$z1_mono.B,1,quantile, 0.025, na.rm = T) %>% bound_vec,
                     CI_upper = apply(mr_estimates$bootstrap$AACE$z1_mono.B,1,quantile, 0.975, na.rm = T) %>% bound_vec)
  S_a_0 = data.frame(time=times,
                     point = mr_estimates$estimate$AACE$z0_mono %>% bound_vec,
                     CI_lower = apply(mr_estimates$bootstrap$AACE$z0_mono.B,1,quantile, 0.025, na.rm = T) %>% bound_vec,
                     CI_upper = apply(mr_estimates$bootstrap$AACE$z0_mono.B,1,quantile, 0.975, na.rm = T) %>% bound_vec)
  PSCE_a = data.frame(time=times,
                      point = mr_estimates$estimate$AACE$diff_mono,
                      CI_lower = apply(mr_estimates$bootstrap$AACE$diff_mono.B,1,quantile, 0.025, na.rm = T),
                      CI_upper = apply(mr_estimates$bootstrap$AACE$diff_mono.B,1,quantile, 0.975, na.rm = T))
  # never takers
  S_n_1 = data.frame(time=times,
                     point = mr_estimates$estimate$NACE$z1_mono %>% bound_vec,
                     CI_lower = apply(mr_estimates$bootstrap$NACE$z1_mono.B,1,quantile, 0.025, na.rm = T) %>% bound_vec,
                     CI_upper = apply(mr_estimates$bootstrap$NACE$z1_mono.B,1,quantile, 0.975, na.rm = T) %>% bound_vec)
  S_n_0 = data.frame(time=times,
                     point = mr_estimates$estimate$NACE$z0_mono %>% bound_vec,
                     CI_lower = apply(mr_estimates$bootstrap$NACE$z0_mono.B,1,quantile, 0.025, na.rm = T) %>% bound_vec,
                     CI_upper = apply(mr_estimates$bootstrap$NACE$z0_mono.B,1,quantile, 0.975, na.rm = T) %>% bound_vec)
  PSCE_n = data.frame(time=times,
                      point = mr_estimates$estimate$NACE$diff_mono,
                      CI_lower = apply(mr_estimates$bootstrap$NACE$diff_mono.B,1,quantile, 0.025, na.rm = T),
                      CI_upper = apply(mr_estimates$bootstrap$NACE$diff_mono.B,1,quantile, 0.975, na.rm = T))
  out = list(Always_Takers = list(S_a_1=S_a_1,S_a_0=S_a_0,PSCE_a=PSCE_a),
             Compliers = list(S_c_1=S_c_1,S_c_0=S_c_0,PSCE_c=PSCE_c),
             Never_Takers = list(S_n_1=S_n_1,S_n_0=S_n_0,PSCE_n=PSCE_n))
  class(out) = "PSCE.PI.SA"
  return(out)
}

#' Bias-corrected multiply robust estimator of the PSCE under violation of monotonicity
#'
#' @param times a vector of time when the principal survival causal effects (PSCEs) are of interest
#' @param data the dataset
#' @param Xpi_names names of the covariates for the propensity score model
#' @param Xe_names names of the covariates for the principal score model
#' @param Xc_names names of the covariates for the censoring model
#' @param Xt_names names of the covariates for the outcome model
#' @param Z_name name of the treatment assignment status
#' @param S_name name of the true treatment receipt status
#' @param U_name name of the observed failuture time
#' @param delta_name name of the censoring indicator
#' @param zeta sensitivity parameter zeta
#' @param B number of the bootstrap replications (default 100)
#'
#' @returns The PSCE estimates and their 95\% confidence intervals
mrPStrata_MO_SA = function(times, data, Xpi_names, Xe_names, Xc_names, Xt_names, Z_name, S_name, U_name, delta_name, zeta=0.01, B=100) {
  # clean dataset and change the colnames
  names(data)[names(data) == Z_name] <- "z"
  names(data)[names(data) == S_name] <- "s"
  names(data)[names(data) == U_name] <- "U"
  names(data)[names(data) == delta_name] <- "delta"
  # specify the working models
  propensity.model = as.formula(paste("z~",paste(Xpi_names,collapse="+"),sep=""))
  principal.model = as.formula(paste("s~",paste(Xpi_names,collapse="+"),sep=""))
  censor.model = as.formula(paste("delta~",paste(Xpi_names,collapse="+"),sep=""))
  failure.model = as.formula(paste("U~",paste(Xpi_names,collapse="+"),sep=""))
  # number of bootstrap replications
  bootstrap = B
  # multiply robust estimates (point + bootstrap samples)
  mr_estimates = BootEst.MO.SA(times,
                               propensity.model,
                               principal.model,
                               censor.model,
                               failure.model,
                               data,
                               zeta,
                               estimand = c("NACE", "CACE", "AACE","DACE"),
                               bootstrap)
  bound_vec <- function(x) {x[x < 0] <- 0;x[x > 1] <- 1; return(x)}
  # organize the point estimation and 95\% confidence interval
  # compliers
  S_c_1 = data.frame(time=times,
                     point = mr_estimates$estimate$CACE$z1_mono %>% bound_vec,
                     CI_lower = apply(mr_estimates$bootstrap$CACE$z1_mono.B,1,quantile, 0.025, na.rm = T) %>% bound_vec,
                     CI_upper = apply(mr_estimates$bootstrap$CACE$z1_mono.B,1,quantile, 0.975, na.rm = T) %>% bound_vec)
  S_c_0 = data.frame(time=times,
                     point = mr_estimates$estimate$CACE$z0_mono %>% bound_vec,
                     CI_lower = apply(mr_estimates$bootstrap$CACE$z0_mono.B,1,quantile, 0.025, na.rm = T) %>% bound_vec,
                     CI_upper = apply(mr_estimates$bootstrap$CACE$z0_mono.B,1,quantile, 0.975, na.rm = T) %>% bound_vec)
  PSCE_c = data.frame(time=times,
                      point = mr_estimates$estimate$CACE$diff_mono,
                      CI_lower = apply(mr_estimates$bootstrap$CACE$diff_mono.B,1,quantile, 0.025, na.rm = T),
                      CI_upper = apply(mr_estimates$bootstrap$CACE$diff_mono.B,1,quantile, 0.975, na.rm = T))
  # always takers
  S_a_1 = data.frame(time=times,
                     point = mr_estimates$estimate$AACE$z1_mono %>% bound_vec,
                     CI_lower = apply(mr_estimates$bootstrap$AACE$z1_mono.B,1,quantile, 0.025, na.rm = T) %>% bound_vec,
                     CI_upper = apply(mr_estimates$bootstrap$AACE$z1_mono.B,1,quantile, 0.975, na.rm = T) %>% bound_vec)
  S_a_0 = data.frame(time=times,
                     point = mr_estimates$estimate$AACE$z0_mono %>% bound_vec,
                     CI_lower = apply(mr_estimates$bootstrap$AACE$z0_mono.B,1,quantile, 0.025, na.rm = T) %>% bound_vec,
                     CI_upper = apply(mr_estimates$bootstrap$AACE$z0_mono.B,1,quantile, 0.975, na.rm = T) %>% bound_vec)
  PSCE_a = data.frame(time=times,
                      point = mr_estimates$estimate$AACE$diff_mono,
                      CI_lower = apply(mr_estimates$bootstrap$AACE$diff_mono.B,1,quantile, 0.025, na.rm = T),
                      CI_upper = apply(mr_estimates$bootstrap$AACE$diff_mono.B,1,quantile, 0.975, na.rm = T))
  # never takers
  S_n_1 = data.frame(time=times,
                     point = mr_estimates$estimate$NACE$z1_mono %>% bound_vec,
                     CI_lower = apply(mr_estimates$bootstrap$NACE$z1_mono.B,1,quantile, 0.025, na.rm = T) %>% bound_vec,
                     CI_upper = apply(mr_estimates$bootstrap$NACE$z1_mono.B,1,quantile, 0.975, na.rm = T) %>% bound_vec)
  S_n_0 = data.frame(time=times,
                     point = mr_estimates$estimate$NACE$z0_mono %>% bound_vec,
                     CI_lower = apply(mr_estimates$bootstrap$NACE$z0_mono.B,1,quantile, 0.025, na.rm = T) %>% bound_vec,
                     CI_upper = apply(mr_estimates$bootstrap$NACE$z0_mono.B,1,quantile, 0.975, na.rm = T) %>% bound_vec)
  PSCE_n = data.frame(time=times,
                      point = mr_estimates$estimate$NACE$diff_mono,
                      CI_lower = apply(mr_estimates$bootstrap$NACE$diff_mono.B,1,quantile, 0.025, na.rm = T),
                      CI_upper = apply(mr_estimates$bootstrap$NACE$diff_mono.B,1,quantile, 0.975, na.rm = T))
  # defiers
  S_d_1 = data.frame(time=times,
                     point = mr_estimates$estimate$DACE$z1_mono %>% bound_vec,
                     CI_lower = apply(mr_estimates$bootstrap$DACE$z1_mono.B,1,quantile, 0.025, na.rm = T) %>% bound_vec,
                     CI_upper = apply(mr_estimates$bootstrap$DACE$z1_mono.B,1,quantile, 0.975, na.rm = T) %>% bound_vec)
  S_d_0 = data.frame(time=times,
                     point = mr_estimates$estimate$DACE$z0_mono %>% bound_vec,
                     CI_lower = apply(mr_estimates$bootstrap$DACE$z0_mono.B,1,quantile, 0.025, na.rm = T) %>% bound_vec,
                     CI_upper = apply(mr_estimates$bootstrap$NACE$z0_mono.B,1,quantile, 0.975, na.rm = T) %>% bound_vec)
  PSCE_d = data.frame(time=times,
                      point = mr_estimates$estimate$DACE$diff_mono,
                      CI_lower = apply(mr_estimates$bootstrap$DACE$diff_mono.B,1,quantile, 0.025, na.rm = T),
                      CI_upper = apply(mr_estimates$bootstrap$DACE$diff_mono.B,1,quantile, 0.975, na.rm = T))


  out = list(Always_Takers = list(S_a_1=S_a_1,S_a_0=S_a_0,PSCE_a=PSCE_a),
             Compliers = list(S_c_1=S_c_1,S_c_0=S_c_0,PSCE_c=PSCE_c),
             Never_Takers = list(S_n_1=S_n_1,S_n_0=S_n_0,PSCE_n=PSCE_n),
             Defiers = list(S_d_1=S_d_1,S_d_0=S_d_0,PSCE_d=PSCE_d))
  class(out) = "PSCE.MO.SA"
  return(out)
}


#' Point estimation for the multiply robust estimator under violation of monotonicity assumption
#'
#' @param times a vector of time when the principal survival causal effects (PSCEs) are of interest
#' @param propensity.model propensity score model
#' @param principal.model principal score model
#' @param censor.model censoring model
#' @param failure.model outcome model
#' @param data dataset
#' @param zeta the sensitivity parameter zeta
#' @param estimand the estimands of interest
#'
#' @returns The PSCE point estimates
PointEst.MO.SA <- function(times, propensity.model, principal.model, censor.model, failure.model, data,
                           zeta,
                           estimand = c("NACE", "CACE", "AACE", "DACE")) {
  # propensity: Z ~ Xtp
  # principal: S ~ Xps
  # censor: delta ~ XC
  # failure: U ~ XT
  tm <- Sys.time()
  lower_bound_vec <- function(vals, seq) {
    # Find the position of the first element that is not less than
    # a given value, or length(inc_seq) + 1 if no such element is found.
    # Time complexity: O(m + n) - often faster than O(n log m) by
    # doing elementwise binary search.
    result <- integer(length(vals))
    p <- 1
    q <- 1
    while (p <= length(vals)) {
      while (q <= length(seq) && seq[q] < vals[p])
        q <- q + 1
      result[p] <- q
      p <- p + 1
    }
    return (result)
  }

  upper_bound_vec <- function(vals, seq) {
    # Find the position of the first element that is greater than
    # a given value, or length(inc_seq) + 1 if no such element is found.
    # Time complexity: O(m + n) - often faster than O(n log m) by
    # applying the univariate version of upper_bound
    result <- integer(length(vals))
    p <- 1
    q <- 1
    while (p <= length(vals)) {
      while (q <= length(seq) && seq[q] <= vals[p])
        q <- q + 1
      result[p] <- q
      p <- p + 1
    }
    return (result)
  }

  mono <- function(seq){
    cur <- Inf
    for (i in 1:length(seq)) {
      seq[i] <- min(cur, seq[i])
      cur <- min(cur, seq[i])
    }
    return (seq)
  }

  Z_name <- propensity.model[[2]]
  S_name <- principal.model[[2]]
  Delta_name <- censor.model[[2]]
  U_name <- failure.model[[2]]

  data <- data %>% arrange({{U_name}})
  n <- nrow(data)
  z <- data %>% pull({{Z_name}})
  s <- data %>% pull({{S_name}})
  U <- data %>% pull({{U_name}})
  delta <- data %>% pull({{Delta_name}})
  XC <- model.matrix(censor.model, data)[, -1]
  XT <- model.matrix(failure.model, data)[, -1]

  fit.tp <- glm(propensity.model, family=binomial(link = "logit"), data = data)
  tp_est <- fit.tp$fitted.values  #estimation of tp
  fit.p1 <- glm(principal.model, family = binomial(link = "logit"), data = data[z == 1, ])
  p1_est <- predict(fit.p1, data, type = "response")
  fit.p0 <- glm(principal.model, family = binomial(link = "logit"), data = data[z == 0, ])
  p0_est <- predict(fit.p0, data, type = "response")

  ps00_est = 1-p0_est-(p1_est-p0_est)/(1-zeta)
  ps11_est = p1_est-(p1_est-p0_est)/(1-zeta)
  ps01_est = (p1_est-p0_est)/(1-zeta)
  ps10_est = zeta*(p1_est-p0_est)/(1-zeta)
  p1.mean.dr = mean((z*(s-p1_est)/tp_est+p1_est))
  p0.mean.dr = mean(((1-z)*(s-p0_est)/(1-tp_est)+p0_est))
  # print(min(1-(p1.mean.dr-p0.mean.dr)/pmin(p1.mean.dr,1-p0.mean.dr)))
  # if(zeta>min(1-(p1.mean.dr-p0.mean.dr)/pmin(p1.mean.dr,1-p0.mean.dr))){
  #   print("zeta out of bound")
  # }

  p1.mean.zeta = mean((z*(s-p1_est)/tp_est+p1_est)/(1-zeta))
  p0.mean.zeta = mean(((1-z)*(s-p0_est)/(1-tp_est)+p0_est)/(1-zeta))

  p1.mean = mean(z*(s-p1_est)/tp_est+p1_est)
  p0.mean = mean((1-z)*(s-p0_est)/(1-tp_est)+p0_est)
  e01.zeta = p1.mean.zeta-p0.mean.zeta
  e00.zeta = 1-p0.mean-p1.mean.zeta+p0.mean.zeta
  e11.zeta = p1.mean-p1.mean.zeta+p0.mean.zeta
  e10.zeta = p1.mean.zeta-p0.mean.zeta-p1.mean+p0.mean
  w01_est = ps01_est/(e01.zeta)
  w00_est = ps00_est/(e00.zeta)
  w11_est = ps11_est/(e11.zeta)
  w10_est = ps10_est/(e10.zeta)

  xi1 = xi0 = 0
  eta1 = eta0 = 1
  eps1 <- exp(xi1 * (times/max(times))^eta1)
  eps0 <- exp(xi0 * (times/max(times))^eta0)
  w1c <- sapply(eps1, function(x) (x * ps01_est + x * ps11_est) / (x * ps01_est + ps11_est))
  w0c <- sapply(eps0, function(x) (x * ps01_est + x * ps00_est) / (x * ps01_est + ps00_est))
  w0n <- sapply(eps0, function(x) (ps01_est + ps00_est) / (x * ps01_est + ps00_est))
  w1a <- sapply(eps1, function(x) (ps01_est + ps11_est) / (x * ps01_est + ps11_est))

  pos1 <- pmin(lower_bound_vec(times, U), n)
  pos2 <- pmax(upper_bound_vec(times, U) - 1, 1)
  time_id <- 1:length(times)

  subsets <- list(
    "00" = s == 0 & z == 0,
    "01" = s == 1 & z == 0,
    "10" = s == 0 & z == 1,
    "11" = s == 1 & z == 1
  )

  call_censor <- call(
    "~",
    as.call(list(
      quote(Surv), U_name, as.call(list(quote(`-`), 1, Delta_name))
    )),
    censor.model[[3]]
  )

  call_failure <- call(
    "~",
    as.call(list(
      quote(Surv), U_name, Delta_name
    )),
    failure.model[[3]]
  )

  new_outer <- function(a, b) {
    return (structure(list(
      a = a,
      b = b
    ), class = 'outer'))
  }

  `[.outer` <- function(x, i, j) {
    if (missing(i)) {
      return (x$a * x$b[j])
    }
    if (missing(j)) {
      return (x$a[i] * x$b)
    }
    return (x$a[i] * x$b[j])
  }

  calc_cumhaz <- function(ss) {

    formula_censor <- eval(call_censor)
    formula_failure <- eval(call_failure)

    fit <- list(
      T = coxph(formula_failure, data = data[ss, ]),
      C = coxph(formula_censor, data = data[ss, ])
    )
    tmp_func <- function(fit, X, delta) {
      #basehaz <- basehaz(fit, centered = F)
      #basehaz1 <- basehaz$hazard[pmax(upper_bound_vec(U, basehaz$time) - 1L, 1L)]
      #cumhaz <- outer(c(exp(X %*% fit$coefficients)), basehaz1)
      b_est <- fit$coefficients
      s1 <- cumsum(ss[n:1]*exp(as.matrix(X[n:1,]) %*% b_est))
      dls <- ss * delta / s1[n:1]
      dls[(1:length(dls))[is.nan(dls)]] <- 0
      cumhaz <- new_outer(c(exp(X %*% fit$coefficients)), cumsum(dls))
      #cumhaz <- outer(c(exp(X %*% fit$coefficients)), cumsum(dls))
      return (cumhaz)
    }
    return (list(
      T = tmp_func(fit$T, XT, delta),
      C = tmp_func(fit$C, XC, 1 - delta)
    ))
  }

  calc_q <- function(cumhaz) {
    q <- diag(NA, n)
    q[,1] = -cumhaz$C[,1] * exp(cumhaz$C[,1] + cumhaz$T[,1])

    for(j in 2:n) {
      q[j:n,j] = q[j:n,j-1] - (cumhaz$C[j:n,j] - cumhaz$C[j:n,j-1]) * exp(cumhaz$C[j:n,j] + cumhaz$T[j:n,j])
    }
    for(i in 1:(n-1)){
      q[i,(i+1):n] = q[i,i] + (1 - delta[i]) * exp(cumhaz$C[i, i] + cumhaz$T[i, i])
    }
    return (q)
  }

  subset_ind <- c(
    "NACE" %in% estimand || "CACE" %in% estimand, # 00
    "AACE" %in% estimand || "DACE" %in% estimand, # 01
    "NACE" %in% estimand || "DACE" %in% estimand, # 10
    "CACE" %in% estimand || "AACE" %in% estimand # 11
  )

  subsets_to_use <- subsets[subset_ind]

  # calculate q only for needed

  cumhaz <- lapply(subsets_to_use, calc_cumhaz)
  q <- lapply(cumhaz, calc_q)
  gc()

  result <- list()

  gen_result <- function(func1, func0) {
    est1 <- mapply(function(x, y, t) func1(x, y, t) / n, pos1, pos2, time_id)
    est0 <- mapply(function(x, y, t) func0(x, y, t) / n, pos1, pos2, time_id)
    z1_mono = mono(est1)
    z0_mono = mono(est0)

    return (list(
      z1 = est1,
      z0 = est0,
      diff = est1 - est0,
      z1_mono = z1_mono,
      z0_mono = z0_mono,
      diff_mono = z1_mono - z0_mono
    ))
  }

  NACE_1 <- function(x, y, t) {
    s1 <- w00_est * (1 - s) * z
    s2 <- (1 - p1_est) * tp_est
    sum(
      s1[x:n] / (s2[x:n] * exp(-cumhaz[["10"]]$C[x:n, y]))
    ) + sum(
      w00_est * (1 - z / tp_est) * exp(-cumhaz[["10"]]$T[, y])
    ) + sum(
      zeta * ((1 - z) * (s - p0_est) / (1 - tp_est) + (1 - p0_est) * z * (p1_est - s) / s2) * exp(-cumhaz[["10"]]$T[, y])
    ) / (e00.zeta * (1 - zeta)) + sum(
      s1 * exp(-cumhaz[["10"]]$T[, y]) * q[["10"]][, x] / s2
    )
  }

  NACE_0 <- function(x, y, t) {
    s1 <- w00_est * (1 - s) * (1 - z)
    s2 <- (1 - p0_est) * (1 - tp_est)
    sum(
      s1[x:n] / (s2[x:n] * exp(-cumhaz[["00"]]$C[x:n, y]))
    ) + sum(
      w00_est * (1 - (1 - z) / (1 - tp_est)) * exp(-cumhaz[["00"]]$T[, y])
    ) - sum(
      (z * (s - p1_est) / tp_est + (1 - p1_est) * (1 - z) * (p0_est - s) / s2) * exp(-cumhaz[["00"]]$T[, y])
    ) / (e00.zeta * (1 - zeta)) + sum(
      s1 * exp(-cumhaz[["00"]]$T[, y]) * q[["00"]][, x] / s2
    )
  }

  CACE_1 <- function(x, y, t) {
    s1 <- w01_est * s * z
    s2 <- p1_est * tp_est
    sum(
      s1[x:n] / (s2[x:n] * exp(-cumhaz[["11"]]$C[x:n, y]))
    ) + sum(
      w01_est * (1 - z / tp_est) * exp(-cumhaz[["11"]]$T[, y])
    ) - sum(
      ((1 - z) * (s - p0_est) / (1 - tp_est) + p0_est * z * (p1_est - s) / s2) * exp(-cumhaz[["11"]]$T[, y])
    ) / (e01.zeta * (1 - zeta)) + sum(
      s1 * exp(-cumhaz[["11"]]$T[, y]) * q[["11"]][, x] / s2
    )
  }

  CACE_0 <-  function(x, y, t) {
    s1 <- w01_est * (1 - s) * (1 - z)
    s2 <- (1 - p0_est) * (1 - tp_est)
    sum(
      s1[x:n] / (s2[x:n] * exp(-cumhaz[["00"]]$C[x:n, y]))
    ) + sum(
      w01_est * (1 - (1 - z) / (1 - tp_est)) * exp(-cumhaz[["00"]]$T[, y])
    ) + sum(
      (z * (s - p1_est) / tp_est + (1 - p1_est) * (1 - z) * (p0_est - s) / s2) * exp(-cumhaz[["00"]]$T[, y])
    ) / (e01.zeta * (1 - zeta)) + sum(
      s1 * exp(-cumhaz[["00"]]$T[, y]) * q[["00"]][, x] / s2
    )
  }

  AACE_1 <-  function(x, y, t) {
    s1 <- w11_est * s * z
    s2 <- p1_est * tp_est
    sum(
      s1[x:n] / (s2[x:n] * exp(-cumhaz[["11"]]$C[x:n, y]))
    ) + sum(
      w11_est * (1 - z / tp_est) * exp(-cumhaz[["11"]]$T[, y])
    ) + sum(
      ((1 - z) * (s - p0_est) / (1 - tp_est) + p0_est * z * (p1_est - s) / s2) * exp(-cumhaz[["11"]]$T[, y])
    ) / (e11.zeta * (1 - zeta)) + sum(
      s1 * exp(-cumhaz[["11"]]$T[, y]) * q[["11"]][, x] / s2
    )
  }

  AACE_0 <- function(x, y, t) {
    s1 <- w11_est * s * (1 - z)
    s2 <- p0_est * (1 - tp_est)
    sum(
      s1[x:n] / (s2[x:n] * exp(-cumhaz[["01"]]$C[x:n, y]))
    ) + sum(
      w11_est * (1 - (1 - z) / (1 - tp_est)) * exp(-cumhaz[["01"]]$T[, y])
    ) + sum(
      zeta * (z * (p1_est - s) / tp_est + p1_est * (1 - z) * (s - p0_est) / s2) * exp(-cumhaz[["01"]]$T[, y])
    ) / (e11.zeta * (1 - zeta)) + sum(
      s1 * exp(-cumhaz[["01"]]$T[, y]) * q[["01"]][, x] / s2
    )
  }

  DACE_1 <-  function(x, y, t) {
    s1 <- w10_est * (1 - s) * z
    s2 <- (1 - p1_est) * tp_est
    sum(
      s1[x:n] / (s2[x:n] * exp(-cumhaz[["10"]]$C[x:n, y]))
    ) + sum(
      w10_est * (1 - z / tp_est) * exp(-cumhaz[["10"]]$T[, y])
    ) + sum(
      zeta * ((1 - z) * (p0_est - s) / (1 - tp_est) + (1 - p0_est) * z * (s - p1_est) / s2) * exp(-cumhaz[["10"]]$T[, y])
    ) / (e10.zeta * (1-zeta)) + sum(
      s1 * exp(-cumhaz[["10"]]$T[, y]) * q[["10"]][, x] / s2
    )
  }

  DACE_0 <- function(x, y, t) {
    s1 <- w10_est * s * (1 - z)
    s2 <- p0_est * (1 - tp_est)
    sum(
      s1[x:n] / (s2[x:n] * exp(-cumhaz[["01"]]$C[x:n, y]))
    ) + sum(
      w10_est * (1 - (1 - z) / (1 - tp_est)) * exp(-cumhaz[["01"]]$T[, y])
    ) + sum(
      zeta * (z * (s - p1_est) / tp_est + p1_est * (1 - z) * (p0_est - s) / s2) * exp(-cumhaz[["01"]]$T[, y])
    ) / (e10.zeta * (1 - zeta)) + sum(
      s1 * exp(-cumhaz[["01"]]$T[, y]) * q[["01"]][, x] / s2
    )
  }

  if ("NACE" %in% estimand)
    result$NACE <- gen_result(NACE_1, NACE_0)
  if ("CACE" %in% estimand)
    result$CACE <- gen_result(CACE_1, CACE_0)
  if ("AACE" %in% estimand)
    result$AACE <- gen_result(AACE_1, AACE_0)
  if ("DACE" %in% estimand)
    result$DACE <- gen_result(DACE_1, DACE_0)
  result$bound <- min(1-(p1.mean.dr-p0.mean.dr)/pmin(p1.mean.dr,1-p0.mean.dr))
  class(result) <- "est.MR.SA"
  return (result)
}



#' Bootstrap confidence intervals for the bias-corrected multiply robust estimator under violations of monotonicity
#'
#' @param times a vector of time when the principal survival causal effects (PSCEs) are of interest
#' @param propensity.model propensity score model
#' @param principal.model principal score model
#' @param censor.model censoring model
#' @param failure.model outcome model
#' @param data dataset
#' @param zeta the sensitivity parameter zeta
#' @param estimand the estimands of interest
#' @param bootstrap number of bootstrap replications
#'
#' @returns The bootstrap confidence intervals
BootEst.MO.SA <- function(times, propensity.model, principal.model, censor.model, failure.model, data,
                         zeta,
                         estimand = c("NACE", "CACE", "AACE", "DACE"), bootstrap = 50) {
  gc()
  mean_result <- PointEst.MO.SA(times, propensity.model, principal.model, censor.model, failure.model, data,
                           zeta, estimand)
  if (bootstrap == 0)
    return (mean_result)
  bs_result <- c()
  for (name in estimand) {
    z1.B <- z0.B <- diff.B <- z1_mono.B <- z0_mono.B <- diff_mono.B <- matrix(nrow = length(times), ncol = bootstrap)
    bs_result[[name]] <- list(
      z1.B = z1.B,
      z0.B = z0.B,
      diff.B = diff.B,
      z1_mono.B = z1_mono.B,
      z0_mono.B = z0_mono.B,
      diff_mono.B = diff_mono.B
    )
  }
  gc()
  pb <- progress::progress_bar$new(total = bootstrap, format = "(:spin) [:bar] :percent. :elapsed | :eta")
  cnt <- 0
  for (i in 1:bootstrap) {
    gc()
    data.B <- slice_sample(data, n = nrow(data), replace = T)
    result <- PointEst.MO.SA(times, propensity.model, principal.model, censor.model, failure.model, data.B,
                        zeta, estimand)
    if (result$bound < zeta) {
      cnt <- cnt + 1
      cat("Resample: ", cnt, " with bound ", result$bound, "\n")
      i <- i - 1
      next
    }
    pb$tick()
    for (name in estimand) {
      bs_result[[name]]$z1.B[, i] <- result[[name]]$z1
      bs_result[[name]]$z0.B[, i] <- result[[name]]$z0
      bs_result[[name]]$diff.B[, i] <- result[[name]]$diff
      bs_result[[name]]$z1_mono.B[, i] <- result[[name]]$z1_mono
      bs_result[[name]]$z0_mono.B[, i] <- result[[name]]$z0_mono
      bs_result[[name]]$diff_mono.B[, i] <- result[[name]]$diff_mono
    }
    gc()
  }
  r_val <- list(
    estimate = mean_result,
    bootstrap = bs_result
  )
  class(r_val) <- "est.MR.SA.bootstrap"
  return (r_val)
}






#' Plot of the PSCEs and their associated 95\% pointwise confidence intervals
#'
#' @param res an output from mrPStrata
#'
#' @returns The PSCE point estimates and 95\% pointwise confidence intervals
plot.psce <- function(res) {
  time = res$Compliers$S_c_1$time
  estimands_list = c("Never_Takers", "Compliers", "Always_Takers")
  if (class(res) == "PSCE.MO.SA") {estimands_list = c("Never_Takers", "Compliers", "Always_Takers","Defiers")}
  df_list <- list()
  for (name in estimands_list) {
    df_list <- c(df_list, list(data.frame(
      time = time,
      estimand = name,
      group = "Treated",
      type = "Outcome",
      est = res[[name]][[1]][,"point"],
      lwr = res[[name]][[1]][,"CI_lower"],
      upr = res[[name]][[1]][,"CI_upper"]
    )))
    df_list <- c(df_list, list(data.frame(
      time = time,
      estimand = name,
      group = "Control",
      type = "Outcome",
      est = res[[name]][[2]][,"point"],
      lwr = res[[name]][[2]][,"CI_lower"],
      upr = res[[name]][[2]][,"CI_upper"]
    )))
    df_list <- c(df_list, list(data.frame(
      time = time,
      estimand = name,
      group = "Difference",
      type = "Difference",
      est = res[[name]][[3]][,"point"],
      lwr = res[[name]][[3]][,"CI_lower"],
      upr = res[[name]][[3]][,"CI_upper"]
    )))
  }
  df <- do.call(bind_rows, df_list)
  df$type[which(df$type=="Difference")] = "PSCE"
  df$type[which(df$type=="Outcome")] = "Principal Survival Probability"
  df$type <- factor(df$type, levels = c("PSCE", "Principal Survival Probability"))
  # df$estimand <- ifelse(
  #   df$estimand == 'Never_Takers', 'Low-inclined',
  #   ifelse(
  #     df$estimand == 'Compliers', 'Compliers', 'High-inclined'
  #   )
  # )
  df$estimand[which(df$estimand=="Always_Takers")] = "Always takers"
  df$estimand[which(df$estimand=="Never_Takers")] = "Never takers"
  if (class(res) == "PSCE.MO.SA") {
    df$estimand <- factor(df$estimand, levels = c("Always takers", "Compliers","Never takers","Defiers"))
  } else {
    df$estimand <- factor(df$estimand, levels = c("Always takers", "Compliers","Never takers"))
  }
  df$group <- factor(df$group, levels = c("Treated", "Control", "Difference"))

  cbp2 <- c("#D52000", "#005CE6", "#000000", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  ggplot(df,aes(x=factor(time), y=est,color=group)) +
    geom_errorbar(aes(ymin=lwr, ymax=upr,color=group,width=0.5),position=position_dodge(0.6)) +
    geom_point(position=position_dodge(0.6)) +
    geom_hline(data=df %>% filter(type=="PSCE"),aes(yintercept=0), linetype="dashed", color = "red") +
    facet_grid(type~estimand, scale = "free_y") +
    ylab(NULL)  + xlab("Time")+
    scale_fill_manual(values = cbp2) + scale_colour_manual(values=cbp2)
}








