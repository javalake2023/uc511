# HIP.R

#' @name hipX1split
#'
#' @title A title.
#'
#' @description A description.
#'
#' @details The details.
#'
#' @param x1pts A parm.
#' @param HaltonIndex Halton indices for all points in x1Hpts.
#' @param BoxIndex A parm.
#' @param xlevel A parm.
#' @param x1Hpts A parm.
#'
#' @return A list containing the following variables:
#'         HaltonIndex - Updated Halton indices for all points in x1Hpts.

#' @export
hipX1split <- function(x1pts, HaltonIndex, BoxIndex, xlevel, x1Hpts) {

  # Determine points in current box
  inBox <- which(HaltonIndex == BoxIndex)
  x1pts <- x1pts[inBox]

  # Randomly remove one point (if needed)
  F <- inBox
  nF <- base::length(F)

  # odd number of points?
  if ((nF %% 2) == 1) {
    r <- sample(nF, 1)
    HaltonIndex[F[r]] <- NA
    F <- F[-r]
    nF <- (nF - 1)
    x1pts <- x1pts[-r]
  }

  # Partition and update Halton indices
  x1sort <- order(x1pts)
  F <- F[x1sort]
  Hpt1 <- x1Hpts[(BoxIndex + 1)]
  Hpt2 <- x1Hpts[(BoxIndex + 1 + xlevel)]
  I <- order(c(Hpt1, Hpt2))

  # set halton index - top half.
  HaltonIndex[F[1:(floor(nF/2))]] <- (BoxIndex + (I[1] - 1) * xlevel)
  # set halton index - bottom half.
  HaltonIndex[F[((floor(nF/2))+1):nF]] <- (BoxIndex + (I[2] - 1) * xlevel)
  # return the halton index.
  return(HaltonIndex)
}

#' @name hipX2split
#'
#' @title A title.
#'
#' @description A description.
#'
#' @details The details.
#'
#' @param x2pts A parm.
#' @param HaltonIndex Halton indices for all points in x2Hpts.
#' @param BoxIndex A parm.
#' @param xlevel A parm.
#' @param x2Hpts A parm.
#'
#' @return A list containing the following variables:
#'         HaltonIndex - Updated Halton indices for all points in x2Hpts.

#' @export
hipX2split <- function(x2pts, HaltonIndex, BoxIndex, xlevel, x2Hpts) {

  # Determine points in current box
  inBox <- which(HaltonIndex == BoxIndex)
  x2pts <- x2pts[inBox]

  # Randomly remove one or two points (if needed)
  F <- inBox
  # get number of points.
  nF <- base::length(F)

  # not divisible by 3.
  if ((nF %% 3) == 1) {
    r <- sample(nF, 1)
    HaltonIndex[F[r]] <- NA
    F <- F[-r]
    nF <- (nF - 1)
    x2pts <- x2pts[-r]
  } else if ((nF %% 3) == 2) {
    r <- sample(nF, 1)
    HaltonIndex[F[r]] <- NA
    F <- F[-r]
    nF <- (nF - 1)
    x2pts <- x2pts[-r]

    r <- sample(nF, 1)
    HaltonIndex[F[r]] <- NA
    F <- F[-r]
    nF <- (nF - 1)
    x2pts <- x2pts[-r]
  }

  # Partition and update Halton indices
  x2sort <- order(x2pts)
  F <- F[x2sort]
  Hpt1 <- x2Hpts[(BoxIndex + 1)]
  Hpt2 <- x2Hpts[(BoxIndex + 1 + xlevel)]
  Hpt3 <- x2Hpts[(BoxIndex + 1 + (2 * xlevel))]
  I <- order(c(Hpt1, Hpt2, Hpt3))

  # update halton indices - top third.
  HaltonIndex[F[1:(floor(nF/3))]] <- (BoxIndex + (I[1] - 1) * xlevel)
  # update halton indices - middle third.
  HaltonIndex[F[((floor(nF/3))+1):(2*(floor(nF/3)))]] <- (BoxIndex + (I[2] - 1) * xlevel)
  # update halton indices - bottom third.
  HaltonIndex[F[((2*(floor(nF/3)))+1):nF]] <- (BoxIndex + (I[3] - 1) * xlevel)
  # return the halton index.
  return(HaltonIndex)
}

#' @name hipPartition
#'
#' @title A title.
#'
#' @description A description.
#'
#' @details The details.
#'
#' @param pts A parm.
#' @param its A parm.
#'
#' @return A list containing the following variables:
#'         ptsIndex    - Some output.
#'         HaltonIndex - Updated Halton indices for all points in pts.

#' @export
hipPartition <- function(pts, its) {
  # Initialize
  N <- nrow(pts)
  pts <- cbind(pts, 1:N)
  # Assuming HaltonPts is a function that generates Halton points.
  Hpts <- HaltonPts(N)
  #res <- cppRSHalton_br(n = N)
  #Hpts <- res$pts
  HaltonIndex <- rep(0, N)

  # Partitioning parameters
  xlevel <- c(1, 2, 6, 12, 24, 72, 144, 432, 864, 1728, 5184, 10368, 20736)
  #xlevel <- c(1, 2, 6, 12, 24, 72, 144, 432, 864, 1728, 5184)
  # support upto 13 levels.
  diml <- c(1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1)
  #diml <- c(1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1)

  # Partitioning loop
  for (j in (1:its)) {
    # for every level...
    if (diml[j] == 1) {
      level <- (xlevel[j] - 1)
      for (i in 0:level) {
        # go calculate halton indices.
        HaltonIndex <- hipX1split(pts[,1], HaltonIndex, i, xlevel[j], Hpts[,1])
      }
    } else {
      level <- (xlevel[j] - 1)
      for (i in 0:level) {
        # go calculate halton indices.
        HaltonIndex <- hipX2split(pts[,2], HaltonIndex, i, xlevel[j], Hpts[,2])
      }
    }

    # Remove discarded points
    TF <- which(!is.na(HaltonIndex))
    HaltonIndex <- HaltonIndex[TF]
    pts <- pts[TF,]
  }
  ptsIndex <- pts[,3]
  return(list(ptsIndex = ptsIndex, HaltonIndex = HaltonIndex))
}

#' @name hipIndexRandomPermutation
#'
#' @title A title.
#'
#' @description A description.
#'
#' @details The details.
#'
#' @param its A parm.
#'
#' @return A list containing the following variables:
#'         permHaltonIndex - some output.

#' @export
hipIndexRandomPermutation <- function(its) {

  # Partitioning parameters - support upto 13 levels.
  diml <- c(1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1)
  #diml <- c(1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1)
  diml <- diml[1:its]

  # Halton Indices
  b1 <- 2
  b2 <- 3
  J1 <- sum(diml == 1)
  J2 <- sum(diml == 2)
  B <- (b1^J1) * (b2^J2)
  HIM <- matrix(0, nrow = (b2^J2), ncol = (b1^J1))
  H <- HaltonPts(B)

  Hindex <- floor(rbind((b1^(J1) + 1e-12) * H[,1], (b2^(J2) + 1e-12) * H[,2])) + 1
  Hindex <- t(Hindex)

  # Halton Matrix
  for (i in 0:(B - 1)) {
    HIM[((b2^J2) + 1) - Hindex[(i + 1), 2], Hindex[(i + 1), 1]] <- i
  }

  # Permutated Halton Matrix
  step2 <- c(2, 4, 8, 16, 32, 64, 128, 256)
  step3 <- c(3, 9, 27, 81, 243)
  order2 <- c((sample(2) - 1), rep(NA, ((b1^J1) - b1)))
  order3 <- c((sample(3) - 1), rep(NA, ((b2^J2) - b2)))

  for (i in 1:(J1 - 1)) {
    #if (J1 < 2) {break}
    if (J1 < 2) {
      #print("J1 < 2")
      next
    }
    v <- vector()
    L <- sum(!is.na(order2))

    for (j in (1:L)) {
      k <- order2[j]
      P <- (sample(2) - 1)
      s <- step2[i]
      v <- c(v, (k + s * P[1]), (k + s * P[2]))
    }
    order2[1:base::length(v)] <- v
  }
  for (i in 1:(J2 - 1)) {
    #if (J2 < 2) {break}
    if (J2 < 2) {
      #print("J2 < 2")
      next
    }
    v <- vector()
    L <- sum(!is.na(order3))
    for (j in (1:L)) {
      k <- order3[j]
      P <- (sample(3) - 1)
      s <- step3[i]
      v <- c(v, (k + s * P[1]),  (k + s * P[2]),  (k + s * P[3]))
    }
    order3[1:base::length(v)] <- v
  }

  # Transform the matrix
  b2vals <- (HIM[1,] %% (b1^J1))
  b3vals <- (HIM[,1] %% (b2^J2))

  I2 <- vector()
  for (i in 1:(b1^J1)) {
    F <- which(b2vals == order2[i])
    I2 <- c(I2, F)
  }
  I3 <- vector()
  for (i in 1:(b2^J2)) {
    F <- which(b3vals == order3[i])
    I3 <- c(I3, F)
  }

  permHIM <- HIM[I3, I2]
  permHaltonIndex <- rep(0, B)

  for (i in 0:(B - 1)) {
    permHaltonIndex[i + 1] <- permHIM[which(HIM == i)]
  }

  return(list(permHaltonIndex = permHaltonIndex, B = B))
}

# delete at some stage - use cpp version instead.
HaltonPts <- function(n) {
  # Initialize
  bases <- c(2, 3)
  pts <- matrix(0, nrow = n, ncol = 2)

  for (ii in (1:n)) {
    for (i in (1:2)) {
      k <- (ii - 1)
      j <- 1
      xk <- (k %% bases[i]) * (1 / bases[i])

      while (floor(k / (bases[i] ^ (j))) > 0) {
        xk <- xk + (floor(k / (bases[i] ^ (j))) %% bases[i]) * (1 / (bases[i] ^ (j + 1)))
        j <- (j + 1)
      }

      pts[ii, i] <- xk
    }
  }
  return(pts)
}

#' @name HIP
#'
#' @title A title.
#'
#' @description A description.
#'
#' @details The details.
#'
#' @param population A population of points. (is the pop always pairs of points?).
#' @param n The number of points to draw from the population. Default 20.
#' @param iterations The levels of partitioning required. Default 7.
#'
#' @return A list containing the following variables:
#'         sampleI - indices of the sample points to use.
#'         popIndex - ???
#'         Order - ???
#'         HaltonIndex - Halton indices for all population points.

#' @export
HIP <- function(population, n = 20, iterations = 7) {

  # must be numeric.
  validate_parameters("hipPopulation", c(population))
  # must be numeric, greater than 0.
  validate_parameters("hipN", c(n))
  # must be numeric, greater than 1 and less than or equal to 13.
  validate_parameters("hipIterations", c(iterations))

  # Partitioning parameters
  iteration_level <- c(1, 2, 6, 12, 24, 72, 144, 432, 864, 1728, 5184, 10368, 20736)
  if(((base::length(population)/2)/iteration_level[iterations+1]) < 1.0){
    msg <- "uc511(HIP) Population size %s not compatible with number of iteration levels %s. %s. %s."
    msgs <- sprintf(msg, base::length(population)/2, iterations, base::length(population), iteration_level[iterations+1])
    stop(msgs)
  }

  # Partition the population
  partitionResult <- hipPartition(population, iterations)
  popIndex <- partitionResult$ptsIndex
  HaltonIndex <- partitionResult$HaltonIndex

  # Random permutation of Halton indices
  permResult <- hipIndexRandomPermutation(iterations)
  permHaltonIndex <- permResult$permHaltonIndex
  B <- permResult$B

  Order <- rep(0, length(HaltonIndex))

  for (i in 0:(B - 1)) {
    Order[which(HaltonIndex == i)] <- permHaltonIndex[i + 1]
  }

  # Assign unique indices
  ptsInBox <- (length(HaltonIndex) / B)
  if (ptsInBox > 1) {
    for (i in 0:(B - 1)) {
      F <- which(Order == i)
      Order[F] <- (i + B * (sample(ptsInBox) - 1))
    }
  }

  # Sample Indices
  sampleI <- rep(1, n)
  for (i in (1:n)) {
    sampleI[i] <- popIndex[which(Order == (i - 1))]
  }

  return(list(sampleI     = sampleI,
              popIndex    = popIndex,
              Order       = Order,
              HaltonIndex = HaltonIndex))
}


