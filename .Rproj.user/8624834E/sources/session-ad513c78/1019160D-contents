binomCIs <- function(x, n, a = 0.05) {
  z <- qnorm(1 - a/2)
  mod <- list()
  prop <- x/n
  pq <- prop * (1 - prop)
  lcor <- exp( log(a/2)/n )
  id0 <- which(prop == 0)
  id1 <- which(prop == 1)
  mod[[ 1 ]] <- prop
  ## Jeffreys
  com <- n - x + 0.5
  low <- qbeta(a/2, x + 0.5, com)
  up <-  qbeta(1 - a/2, x + 0.5, com)
  mod[[ 2 ]] <- cbind(low, up)
  colnames(mod[[ 2 ]]) <- c(paste(a/2, "%", sep = ""), paste(1 - a/2, "%", sep = ""))
  ## Wald
  com <- z * sqrt( pq/n )
  res <- cbind(prop - com, prop + com)
  res[id0, 1] <- 0
  res[id0, 2] <- 1 - lcor[id0]
  res[id1, 1] <- lcor[id1]
  res[id1, 2] <- 1
  mod[[3]] <- res
  colnames(mod[[ 3 ]]) <- c(paste(a/2, "%", sep = ""), paste(1 - a/2, "%", sep = ""))
  ## Wald corrected
  com <- z * sqrt( pq/n ) + 0.5/n
  res <- cbind(prop - com, prop + com)
  res[id0, 1] <- 0
  res[id0, 2] <- 1 - lcor[id0]
  res[id1, 1] <- lcor[id1]
  res[id1, 2] <- 1
  mod[[ 4 ]] <- res
  colnames(mod[[ 4 ]]) <- c(paste(a/2, "%", sep = ""), paste(1 - a/2, "%", sep = ""))
  ## Wald BS
  com <- z * sqrt( pq ) / sqrt(n - z^2 - 2 * z/sqrt(n) - 1/n) + 0.5/n
  res <- cbind(prop - com, prop + com)
  res[id0, 1] <- 0
  res[id0, 2] <- 1 - lcor[id0]
  res[id1, 1] <- lcor[id1]
  res[id1, 2] <- 1
  mod[[ 5 ]] <- res
  colnames(mod[[ 5 ]]) <- c(paste(a/2, "%", sep = ""), paste(1 - a/2, "%", sep = ""))
  ## Agresti and Coull
  theta <- (x + 2) / (n + 4)
  com <- z * sqrt( theta * (1 - theta)/(n + 4) )
  res <- cbind(theta - com, theta + com)
  mod[[ 6 ]] <- res
  colnames(mod[[ 6 ]]) <- c(paste(a/2, "%", sep = ""), paste(1 - a/2, "%", sep = ""))
  ## Wilson
  xb <- x + z^2/2
  nb <- n + z^2
  pb <- xb / nb
  com <- z * sqrt(n) / nb * sqrt( pq + z^2/4/n )
  res <- cbind(pb - com, pb + com)
  mod[[ 7 ]] <- res
  colnames(mod[[ 7 ]]) <- c(paste(a/2, "%", sep = ""), paste(1 - a/2, "%", sep = ""))
  ## Score
  com <- z * sqrt(x - x^2/n + z^2/4)
  res <- cbind( x + z^2 - com, x + z^2 + com ) / (n + z^2)
  mod[[ 8 ]] <- res
  colnames(mod[[ 8 ]]) <- c(paste(a/2, "%", sep = ""), paste(1 - a/2, "%", sep = ""))
  ## Score corrected
  b1 <- x - 0.5   ;   b2 <- x + 0.5
  l1 <- b1 + 0.5 * z^2 - z * sqrt( b1 - b1^2/n + 0.25 * z^2 )
  l2 <- b2 + 0.5 * z^2 + z * sqrt( b2 - b2^2/n + 0.25 * z^2 )
  res <- cbind(l1, l2) / (n + z^2)
  mod[[ 9 ]] <- res
  colnames(mod[[ 9 ]]) <- c(paste(a/2, "%", sep = ""), paste(1 - a/2, "%", sep = ""))
  ## Wald logit
  b <- log( x / (n - x) )
  com <- z / sqrt(x * (1 - prop) )
  res <- exp( cbind(b - com, b + com) )
  res <- 1 - 1 / ( 1 + res)
  res[id0, 1] <- 0
  res[id0, 2] <- 1 - lcor[id0]
  res[id1, 1] <- lcor[id1]
  res[id1, 2] <- 1
  mod[[ 10 ]] <- res
  colnames(mod[[ 10 ]]) <- c(paste(a/2, "%", sep = ""), paste(1 - a/2, "%", sep = ""))
  ## Wald logit corrected
  pb <- x + 0.5
  qb <- n - x + 0.5
  b <- log( pb / qb )
  com <- z / sqrt( (n + 1) * pb/(n + 1) * ( 1 - pb/(n + 1) ) )
  res <- exp( cbind(b - com, b + com) )
  res <- 1 - 1/( 1 + res)
  mod[[ 11 ]] <- res
  colnames(mod[[ 11 ]]) <- c(paste(a/2, "%", sep = ""), paste(1 - a/2, "%", sep = ""))
  ## Arcsine
  pb <- asin( sqrt(prop) )
  com <- 0.5 * z / sqrt(n)
  res <- cbind(pb - com, pb + com)
  res <- sin(res)^2
  res[id0, 1] <- 0
  res[id0, 2] <- 1 - lcor[id0]
  res[id1, 1] <- lcor[id1]
  res[id1, 2] <- 1
  mod[[ 12 ]] <- res
  colnames(mod[[ 11 ]]) <- c(paste(a/2, "%", sep = ""), paste(1 - a/2, "%", sep = ""))
  ## Exact binomial
  a1 <- n - x + 1
  d1 <- x * qf(a/2, 2 * x, 2 * a1)
  a2 <- a1 - 1
  d2 <- (x + 1) * qf(1 - a/2, 2 * (x + 1), 2 * a2)
  res <- cbind( 1 + a1/d1, 1 + a2/d2 )
  res <- 1 / res
  res[id0, 1] <- 0
  res[id1, 2] <- 1
  mod[[ 13 ]] <- res
  colnames(mod[[ 13 ]]) <- c(paste(a/2, "%", sep = ""), paste(1 - a/2, "%", sep = ""))
  names(mod) <- c("Props", "Jeffreys", "Wald", "Wald corrected", "Wald-Blyth-Still", "Agresti-Coull", "Wilson",
                  "Score", "Score corrected", "Wald logit", "Wald logit corrected", "Arcsine", "Exact binomial")
  mod
}

