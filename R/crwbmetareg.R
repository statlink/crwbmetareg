crwbmetareg <- function(target, se, dataset, cluster, weights, boot.reps = 1000, prog.bar = FALSE, seed = NULL) {
  dataset <- model.matrix(target ~., data = as.data.frame(dataset))[, -1]
  mod <- glm(target ~., data = as.data.frame( cbind(se, dataset) ), weights = weights )
  ##summary(mod)[[ 12 ]]
  a <- try( .cluster.wild.glm3( mod, dat = as.data.frame( cbind(se, dataset, cluster) ), weights = weights,
             cluster = ~cluster, boot.reps = boot.reps, prog.bar = prog.bar, seed = seed ), silent = TRUE )
  a$p.values
}



## The following function is a modification of the function cluster.wild.glm() of the package clusterSEs

.cluster.wild.glm3 <-function(mod, dat, cluster, weights, ci.level = 0.95, impose.null = TRUE, boot.reps = 1000,
                             report = FALSE, prog.bar = FALSE, output.replicates = FALSE, seed = NULL) {

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

  if (mod$family[1] != "gaussian" | mod$family[2] != "identity") {
    stop("Use only with gaussian family models with a linear link")
  }
  if ( output.replicates  &  impose.null ) {
    stop("Recovering bootstrap replicates requires setting impose.null = FALSE")
  }

  form <- mod$formula                                            # what is the formula of this model?
  variables <- all.vars(form)                                    # what variables are in this model?
  clust.name <- all.vars(cluster)                                # what is the name of the clustering variable?
  used.idx <- which( rownames(dat) %in% rownames(mod$model) )      # what were the actively used observations in the model?
  dat <- dat[used.idx,]                                          # keep only active observations
  clust <- as.vector( unlist( dat[[ clust.name ]] ) )                  # store cluster index in convenient vector
  G<-length(unique(clust))                                       # how many clusters are in this model?
  # ind.variables <- attr(mod$terms, "term.labels")              # what independent variables are in this model? (deprecated)
  "%w/o%" <- function(x, y) x[!x %in% y]                         # create a without function (see ?match)
  dv <- variables %w/o% all.vars( update(form, 1 ~ .) )            # what is the dependent variable?
  ind.variables.data <- all.vars( update(form, 1 ~ .) )            # RHS variables in this model (before variable transforms)
  ind.variables.names.full <- names(coefficients(mod) )           # printed names of coefs (w/ intercept)
  ind.variables.names <- rownames( summary(mod)$coefficients )     # printed names of coefs (w/ intercept), neglecting drops
  ind.variables <- ind.variables.names %w/o% "(Intercept)"       # what independent variables are in this model, neglecting drops and intercept?
  # in case dv is wrapped in a function, need to set it to its functional value
  # so that residuals can be added w/o incident
  dat$dv.new <- mod$y                                            # add transformed DV into data set
  form.new <- update(form, dv.new ~ .)                           # substitute in new dV
  # check to see whether any IVs are factors
  fac <- max(0, min( length(mod$xlevels), 1) )
  # do not impose the null for factor variables
  if ( fac & impose.null ) {
    #cat("\n","\n", "Note: null not imposed (factor variables are present).", "\n", "\n")
    impose.null <- FALSE
  }
  # check whether there are (automatic) interaction terms
  interaction <- max(attr(mod$terms, "order"),  1)
  # do not impose the null for interaction terms
  if ( interaction > 1 & impose.null ) {
    #cat("\n","\n", "Note: null not imposed (interactions are present).", "\n", "\n")
    impose.null <-FALSE
  }

  # check for polynomial terms
  poly.check <- max(unlist(lapply(mod$model, FUN=class))=="poly")
  # do not impose the null for polynomial terms
  if ( poly.check == 1 & impose.null ) {
    #cat("\n","\n", "Note: null not imposed (polynomial terms are present).", "\n", "\n")
    impose.null <- FALSE
  }
  # load in a function to create clustered standard errors
  # by Mahmood Arai: http://thetarzan.wordpress.com/2011/06/11/clustered-standard-errors-in-r/
    cl <- function(dat, fm, cluster, weights) {
      #require(sandwich, quietly = TRUE)
      #require(lmtest, quietly = TRUE)
      M <- length (unique(cluster) )
      N <- length(cluster)
      K <- fm$rank
      dfc <- (M/(M-1))
      uj  <- apply( sqrt(weights) * sandwich::estfun(fm), 2, function(x) tapply( x, cluster, sum) );
      vcovCL <- dfc * sandwich::sandwich(fm, meat. = crossprod(uj)/N)
      lmtest::coeftest(fm, vcovCL)
    }

  se.clust <- cl(dat, mod, clust, weights)[ind.variables.names, 2]               # retrieve the clustered SEs
  beta.mod <- coefficients(mod)[ind.variables.names]                   # retrieve the estimated coefficients
  w <- beta.mod / se.clust                                             # calculate the wald test statistics
  # if the null is to be imposed, execute the following code
  if ( impose.null ) {
    p.store <- c()                                                            # store wild boostrapped p-values
    w.store <- matrix(data = NA, nrow = boot.reps, ncol = length(ind.variables) )    # store bootstrapped test statistics
    if (attr(mod$terms, "intercept") == 1 ) { offset <- 1 } else offset <- 0
    # no factors, so remove any variables that were dropped from original model
    # (code from http://www.cookbook-r.com/Formulas/Creating_a_formula_from_a_string/)
    form.new <- as.formula( paste("dv.new", paste(ind.variables, collapse = " + "), sep=" ~ ") )

    if ( prog.bar )  cat("\n")
    for ( j in 1:length(ind.variables) ) {

      if ( prog.bar )  cat("Independent variable being bootstrapped: ", ind.variables[j], "\n")
      # run model imposing the null hypothesis
      form.null <- as.formula( paste ("dv.new", "~", paste( ind.variables[1:length(ind.variables) %w/o% j], collapse= " + " ) ) )
      mod.null <- glm(form.null, data = dat, weights = weights, family = mod$family)
      null.resid <- residuals(mod.null)
      boot.dat <- dat           # copy the data set into a bootstrap resampling dataset
      wald.store <- c()         # create a container for storing the test statistics
      if ( prog.bar ) pb <- txtProgressBar(min = 0, max = boot.reps, initial = 0, style = 3)
      for ( i in 1:boot.reps ) {
        if (prog.bar)  setTxtProgressBar(pb, value = i)
        # assign wild bootstrap weights
        weight <- c(1, -1)[rbinom(G, size = 1, prob = 0.5)                   # assign wild bootstrap weights
                           + 1][ match( clust, unique(clust) ) ]
        pseudo.resid <- null.resid * weight                                # create pseudo-residuals using weights
        pseudo.dv <- predict(mod.null) + pseudo.resid                     # create pseudo-observations using pseudo-residuals
        boot.dat[,"dv.new"] <- pseudo.dv                                 # create a bootstrap replicate data set (note dv.new)
        boot.mod <- glm(form.new, data = boot.dat, weights = weights, family = mod$family)  # run a model on the bootstrap replicate data (note form.new)
        se.boot <- cl(boot.dat, boot.mod, clust, weights)[offset + j,2]             # retrieve the bootstrap clustered SE
        beta.boot <- coefficients(boot.mod)[offset + j]                    # store the bootstrap beta coefficient
        wald.store[i] <- beta.boot / se.boot                             # store the bootstrap test statistic
      }
      if ( prog.bar )  close(pb)

      p.store[j] <- 1 - sum( abs(w[offset + j]) > abs(wald.store), na.rm = TRUE ) / boot.reps    # calculate the wild bootstrap p-value
      w.store[,j] <- wald.store
    }

    # calculate t-stat for intercept, if present, w/o imposing the null
    if (attr(mod$terms, "intercept") == 1 ) {

      if ( prog.bar )  cat("Independent variable being bootstrapped:  Intercept (null not imposed)", "\n")
      # don't impose the null for the constant (but still call it null.resid)
      null.resid <- residuals(mod)
      boot.dat <- dat           # copy the data set into a bootstrap resampling dataset
      wald.store <- c()         # create a container for storing the test statistics
      if ( prog.bar )  pb <- txtProgressBar(min = 0, max = boot.reps, initial = 0, style = 3)
      for ( i in 1:boot.reps ) {
        if ( prog.bar )  setTxtProgressBar(pb, value = i)
        weight <- c(1, -1)[rbinom(G, size = 1, prob = 0.5)                   # assign wild bootstrap weights
                    + 1][ match(clust, unique(clust)) ]
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

    ci.lo <- NULL
    ci.hi <- NULL
    print.ci <- NULL
    out.ci <- NULL
  # if the null is NOT to be imposed...
  } else {

    if ( prog.bar )  cat("Wild Cluster bootstrapping w/o imposing null...", "\n")
    boot.dat <- dat                                              # copy the data set into a bootstrap resampling dataset
    w.store <- matrix( data = NA, nrow = boot.reps, ncol = length(ind.variables.names) )    # store bootstrapped test statistics
    # keep track of the beta bootstrap replicates for possible output
    rep.store <- matrix( data = NA, nrow = boot.reps, ncol = length(beta.mod) )
    colnames(rep.store) <- ind.variables.names
    resid <- residuals(mod)                                                         # get the residuals for the model
    if ( prog.bar )  pb <- txtProgressBar(min = 0, max = boot.reps, initial = 0, style = 3)
    for ( i in 1:boot.reps ) {

      if ( prog.bar )  setTxtProgressBar(pb, value = i)
      weight <- c(1, -1)[rbinom(G, size = 1, prob = 0.5)                   # assign wild bootstrap weights
                         + 1][ match( clust, unique(clust) ) ]
      pseudo.resid <- resid * weight                                     # create pseudo-residuals using weights
      pseudo.dv <- predict(mod) + pseudo.resid                          # create pseudo-observations using pseudo-residuals
      boot.dat[,"dv.new"] <- pseudo.dv                                 # create a bootstrap replicate data set (note dv.new)
      boot.mod <- glm(form.new, data = boot.dat, weights = weights, family = mod$family)  # run a model on the bootstrap replicate data (note form.new)
      se.boot <- cl(boot.dat, boot.mod, clust, weights)[, 2]                     # retrieve the bootstrap clustered SE
      beta.boot <- coefficients(boot.mod)[ind.variables.names]         # store the bootstrap beta coefficient
      w.store[i,] <- (beta.boot-beta.mod) / se.boot                    # store the bootstrap test statistic
      rep.store[i,] <- beta.boot                                       # store the bootstrap beta for output
    }

    if ( prog.bar )  close(pb)

    comp.fun <- function(vec2, vec1){ as.numeric(vec1>vec2) }                            # a simple function comparing v1 to v2
    p.store.s <- t( apply( X = abs(w.store), FUN = comp.fun, MARGIN = 1, vec1 = abs(w) ) ) # compare the BS test stats to orig. result
    if ( dim(p.store.s)[1] == 1 ) {
      p.store <- 1 - ( sum(p.store.s) / dim(w.store)[1] )                          # calculate the cluster bootstrap p-value
    } else {
      p.store <- 1 - ( colSums(p.store.s) / dim(w.store)[1] )                          # calculate the cluster bootstrap p-value
    }
    # compute critical t-statistics for CIs
    crit.t <- apply(X = abs(w.store), MARGIN = 2, FUN = quantile, probs = ci.level )
    ci.lo <- beta.mod - crit.t * se.clust
    ci.hi <- beta.mod + crit.t * se.clust
    print.ci <- cbind(ind.variables.names, ci.lo, ci.hi)
    print.ci <- rbind(c("variable name", "CI lower", "CI higher"), print.ci)
    out.ci <- cbind(ci.lo, ci.hi)
    rownames(out.ci) <- ind.variables.names
    colnames(out.ci) <- c("CI lower", "CI higher")
  }

  out <- matrix(p.store, ncol = 1)
  colnames(out) <- c("wild cluster BS p-value")
  rownames(out) <- ind.variables.names
  out.p <- cbind(ind.variables.names, round(out, 3))
  out.p <- rbind(c("variable name", "wild cluster BS p-value"), out.p)
  printmat <- function(m) {
    write.table( format(m, justify = "right"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "   ")
  }

  if ( report ) {
    cat("\n", "\n", "Wild Cluster Bootstrapped p-values: ", "\n", "\n")
    printmat(out.p)
    if ( is.null(print.ci) == FALSE ) {
      cat("\n", "Confidence Intervals (derived from bootstrapped t-statistics): ", "\n", "\n")
      printmat(print.ci)
    }

    if (length(ind.variables.names) < length( ind.variables.names.full) ) {
    cat("\n", "\n", "****", "Note: ", length(ind.variables.names.full) - length(ind.variables.names), " variables were unidentified in the model and are not reported.", "****", "\n", sep = "")
    cat("Variables not reported:", "\n", sep = "")
    cat(ind.variables.names.full[!ind.variables.names.full %in% ind.variables.names], sep = ", ")
    cat("\n", "\n")
    }

   }

  out.list <- list()
  out.list[[ "p.values" ]] <- out
  out.list[[ "ci" ]] <- out.ci
  if ( output.replicates )  out.list[[ "replicates" ]] <- rep.store
  return( invisible(out.list) )
}

