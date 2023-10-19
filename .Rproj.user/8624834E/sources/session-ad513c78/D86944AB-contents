binomCI <- function(x, n, a = 0.05) {
    z <- qnorm(1 - a/2)
    mod <- list()
    prop <- x/n
    pq <- prop * (1 - prop)
    lcor <- exp( log(a/2)/n )
    mod <- matrix(nrow = 12, ncol = 2)
    colnames(mod) <- c(paste(a/2, "%", sep = ""), paste(1 - a/2, "%", sep = ""))
    ## Jeffreys
    com <- n - x + 0.5
    low <- qbeta(a/2, x + 0.5, com)
    up <-  qbeta(1 - a/2, x + 0.5, com)
    mod[1, ] <- c(low, up)
    ## Wald
    com <- z * sqrt( pq/n )
    res <- c(prop - com, prop + com)
    if (prop == 0) res <- c(0, 1 - lcor)
    if (prop == 1) res <- c(lcor, 1)
    mod[2, ] <- res
    ## Wald corrected
    com <- z * sqrt( pq/n ) + 0.5/n
    res <- c(prop - com, prop + com)
    if (prop == 0) res <- c(0, 1 - lcor)
    if (prop == 1) res <- c(lcor, 1)
    mod[3, ] <- res
    ## Wald Blyth and Still
    com <- z * sqrt( pq ) / sqrt(n - z^2 - 2 * z/sqrt(n) - 1/n) + 0.5/n
    res <- c(prop - com, prop + com)
    if (prop == 0) res <- c(0, 1 - lcor)
    if (prop == 1) res <- c(lcor, 1)
    mod[4, ] <- res
    ## Agresti and Coull
    theta <- (x + 2) / (n + 4)
    com <- z * sqrt( theta * (1 - theta)/(n + 4) )
    mod[5, ] <- c(theta - com, theta + com)
    ## Wilson
    xb <- x + z^2/2
    nb <- n + z^2
    pb <- xb / nb
    com <- z * sqrt(n) / nb * sqrt( pq + z^2/4/n )
    mod[6, ] <- c(pb - com, pb + com)
    ## Score
    com <- z * sqrt(x - x^2/n + z^2/4)
    mod[7, ] <- c( x + z^2 - com, x + z^2 + com ) / (n + z^2)
    ## Score corrected
    b1 <- x - 0.5   ;   b2 <- x + 0.5
    l1 <- b1 + 0.5 * z^2 - z * sqrt( b1 - b1^2/n + 0.25 * z^2 )
    l2 <- b2 + 0.5 * z^2 + z * sqrt( b2 - b2^2/n + 0.25 * z^2 )
    mod[8, ] <- c(l1, l2) / (n + z^2)
    ## Wald logit
    b <- log( x / (n - x) )
    com <- z / sqrt(x * (1 - prop) )
    res <- exp( c(b - com, b + com) )
    res <- 1 - 1 / ( 1 + res)
    if (prop == 0) res <- c(0, 1 - lcor)
    if (prop == 1) res <- c(lcor, 1)
    mod[9, ] <- res
    ## Wald logit corrected
    pb <- x + 0.5
    qb <- n - x + 0.5
    b <- log( pb / qb )
    com <- z / sqrt( (n + 1) * pb/(n + 1) * ( 1 - pb/(n + 1) ) )
    res <- exp( c(b - com, b + com) )
    mod[10, ] <- 1 - 1/( 1 + res)
    ## Arcsine
    pb <- asin( sqrt(prop) )
    com <- 0.5 * z / sqrt(n)
    res <- c(pb - com, pb + com)
    res <- sin(res)^2
    if (prop == 0) res <- c(0, 1 - lcor)
    if (prop == 1) res <- c(lcor, 1)
    mod[11, ] <- res
    ## Exact binomial
    a1 <- n - x + 1
    d1 <- x * qf(a/2, 2 * x, 2 * a1)
    a2 <- a1 - 1
    d2 <- (x + 1) * qf(1 - a/2, 2 * (x + 1), 2 * a2)
    res <- c( 1 + a1/d1, 1 + a2/d2 )
    res <- 1 / res
    if (prop == 0) res[1] <- 0
    if (prop == 1) res[2] <- 1
    mod[12, ] <- res
    rownames(mod) <- c("Jeffreys", "Wald", "Wald corrected", "Wald-Blyth-Still", "Agresti-Coull", "Wilson",
    "Score", "Score corrected", "Wald logit", "Wald logit corrected", "Arcsine", "Exact binomial")
    list(prop = prop, ci = mod)
}

