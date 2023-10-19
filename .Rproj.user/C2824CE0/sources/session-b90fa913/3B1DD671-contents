fatpet <- function(target, se, cluster, weights, boot.reps = 1000, prog.bar = FALSE, seed = NULL) {
  mod <- glm(target ~., data = as.data.frame(se), weights = weights )
  a <- try( .cluster.wild.glm4(mod, dat = as.data.frame( cbind(se, cluster) ), weights = weights,
             cluster = ~cluster, boot.reps = boot.reps, prog.bar = prog.bar, seed = seed ), silent = TRUE )
  ## summary(mod)[[ 12 ]]
  a
}



## The following function is a modification of the function cluster.wild.glm() of the package clusterSEs

.cluster.wild.glm4 <-function(mod, dat, cluster, weights, boot.reps = 1000,
                             prog.bar = TRUE, seed = NULL) {

  # compensate for bizarre R formula updating bug
  # thanks to Jason Thorpe for reporting!
  form.old <- update(mod$formula, 1 ~ 1 )
  while ( form.old != mod$formula ) {
    form.old <- mod$formula
    invisible( mod <- update(mod, formula = .~.) )
  }


  if ( is.null(seed) == FALSE ) {                                               # if user supplies a seed, set it

    tryCatch( set.seed(seed),
             error = function(e) { return("seed must be a valid integer") },
             warning = function(w) { return(NA) } )

  }

  form <- mod$formula                                            # what is the formula of this model?
  variables <- all.vars(form)                                    # what variables are in this model?
  clust.name <- all.vars(cluster)                                # what is the name of the clustering variable?
  used.idx <- which(rownames(dat) %in% rownames( mod$model) )      # what were the actively used observations in the model?
  dat <- dat[used.idx,]                                          # keep only active observations
  clust <- as.vector( unlist(dat[[ clust.name ]]) )                  # store cluster index in convenient vector
  G <- length( unique(clust) )                                       # how many clusters are in this model?
  # ind.variables <- attr(mod$terms, "term.labels")              # what independent variables are in this model? (deprecated)
  "%w/o%" <- function(x, y) x[!x %in% y]                         # create a without function (see ?match)
  dv <- variables %w/o% all.vars( update(form, 1 ~ .) )            # what is the dependent variable?
  ind.variables.data <- all.vars( update(form, 1 ~ .) )            # RHS variables in this model (before variable transforms)
  ind.variables.names.full <- names( coefficients(mod) )           # printed names of coefs (w/ intercept)
  ind.variables.names <- rownames(summary(mod)$coefficients)     # printed names of coefs (w/ intercept), neglecting drops
  ind.variables <- ind.variables.names %w/o% "(Intercept)"       # what independent variables are in this model, neglecting drops and intercept?

  # in case dv is wrapped in a function, need to set it to its functional value
  # so that residuals can be added w/o incident
  dat$dv.new <- mod$y                                            # add transformed DV into data set
  form.new <- update(form, dv.new ~ .)                           # substitute in new dV

  # load in a function to create clustered standard errors
  # by Mahmood Arai: http://thetarzan.wordpress.com/2011/06/11/clustered-standard-errors-in-r/
    cl <- function(dat, fm, cluster, weights) {
      #require(sandwich, quietly = TRUE)
      #require(lmtest, quietly = TRUE)
      M <- length( unique(cluster) )
      N <- length(cluster)
      K <- fm$rank
      dfc <- M / (M - 1)
      uj  <- apply( sqrt(weights) * sandwich::estfun(fm), 2, function(x) tapply( x, cluster, sum) )
      vcovCL <- dfc * sandwich::sandwich(fm, meat. = crossprod(uj)/N)
      lmtest::coeftest(fm, vcovCL)
    }

  se.clust <- cl(dat, mod, clust, weights)[ind.variables.names, 2]               # retrieve the clustered SEs
  beta.mod <- coefficients(mod)[ind.variables.names]                   # retrieve the estimated coefficients
  w <- beta.mod / se.clust                                             # calculate the wald test statistics

  # if the null is to be imposed, execute the following code

    p.store <- c()                                                            # store wild boostrapped p-values
    w.store <- matrix(data = NA, nrow = boot.reps, ncol = length(ind.variables) )    # store bootstrapped test statistics

    if ( attr(mod$terms, "intercept") == 1 )  offset <- 1  else offset <- 0

    # no factors, so remove any variables that were dropped from original model
    # (code from http://www.cookbook-r.com/Formulas/Creating_a_formula_from_a_string/)
    form.new <- as.formula(dv.new ~ se)

    if ( prog.bar )  cat("\n")
    for ( j in 1:length(ind.variables) ) {

      if( prog.bar )  cat("Independent variable being bootstrapped: ", ind.variables[j], "\n")

      # run model imposing the null hypothesis
      mod.null <- glm(dv.new ~ 1, data = dat, weights = weights, family = mod$family)
      null.resid <- residuals(mod.null)

      boot.dat <- dat           # copy the data set into a bootstrap resampling dataset
      wald.store <- c()         # create a container for storing the test statistics

      if ( prog.bar )  pb <- txtProgressBar(min = 0, max = boot.reps, initial = 0, style = 3)
      for ( i in 1:boot.reps ) {

        if ( prog.bar )  setTxtProgressBar(pb, value = i)

        # assign wild bootstrap weights
        weight <- c(1, -1)[rbinom(G, size = 1, prob = 0.5)                   # assign wild bootstrap weights
                           + 1][ match( clust, unique(clust) ) ]
        pseudo.resid <- null.resid * weight                                # create pseudo-residuals using weights
        pseudo.dv <- predict(mod.null) + pseudo.resid                     # create pseudo-observations using pseudo-residuals
        boot.dat[,"dv.new"] <- pseudo.dv                                 # create a bootstrap replicate data set (note dv.new)

        boot.mod <- glm(form.new, data = boot.dat, weights = weights, family = mod$family)  # run a model on the bootstrap replicate data (note form.new)

        se.boot <- cl(boot.dat, boot.mod, clust, weights)[offset + j, 2]             # retrieve the bootstrap clustered SE
        beta.boot <- coefficients(boot.mod)[offset + j]                    # store the bootstrap beta coefficient
        wald.store[i] <- beta.boot / se.boot                             # store the bootstrap test statistic

      }
      if (prog.bar)  close(pb)

      p.store[j] <- 1 - sum( abs(w[offset + j]) > abs(wald.store), na.rm = TRUE ) / boot.reps    # calculate the wild bootstrap p-value
      w.store[, j] <- wald.store

    }

    # calculate t-stat for intercept, if present, w/o imposing the null
    if ( attr(mod$terms, "intercept") == 1 ) {

      if ( prog.bar )  cat("Independent variable being bootstrapped:  Intercept (null not imposed)", "\n")

      # don't impose the null for the constant (but still call it null.resid)
      null.resid <- residuals(mod)

      boot.dat <- dat           # copy the data set into a bootstrap resampling dataset
      wald.store <- c()         # create a container for storing the test statistics

      if ( prog.bar )  pb <- txtProgressBar(min = 0, max = boot.reps, initial = 0, style = 3)
      for ( i in 1:boot.reps) {

        if ( prog.bar )  setTxtProgressBar(pb, value = i)

        weight <- c(1, -1)[rbinom(G, size = 1, prob = 0.5)                   # assign wild bootstrap weights
                    + 1][ match( clust, unique(clust) ) ]
        pseudo.resid <- null.resid * weight                                # create pseudo-residuals using weights
        pseudo.dv <- predict(mod) + pseudo.resid                          # create pseudo-observations using pseudo-residuals
        boot.dat[,"dv.new"] <- pseudo.dv                                 # create a bootstrap replicate data set (note dv.new)

        boot.mod <- glm(form.new, data = boot.dat, weights = weights, family = mod$family)  # run a model on the bootstrap replicate data (note form.new)

        se.boot <- cl(boot.dat, boot.mod, clust, weights)[1, 2]                    # retrieve the bootstrap clustered SE
        beta.boot <- coefficients(boot.mod)[1]                           # store the bootstrap beta coefficient
        wald.store[i] <- (beta.boot - beta.mod[1]) / se.boot             # store the bootstrap test statistic

      }
      if ( prog.bar )  close(pb)

      p.store <- c( 1 - ( sum( abs(w[1]) > abs(wald.store), na.rm = TRUE ) / boot.reps ), p.store)    # calculate the wild bootstrap p-value
      w.store <- cbind(wald.store, w.store)

    }  ##  end  if(attr(mod$terms, "intercept") == 1 ){

    p.store

}

