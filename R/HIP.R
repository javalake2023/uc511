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
  inBox <- base::which(HaltonIndex == BoxIndex)
  x1pts <- x1pts[inBox]

  # Randomly remove one point (if needed)
  F <- inBox
  nF <- base::length(F)

  # odd number of points?
  if ((nF %% 2) == 1) {
    r <- base::sample(nF, 1)
    HaltonIndex[F[r]] <- NA
    F <- F[-r]
    nF <- (nF - 1)
    x1pts <- x1pts[-r]
  }

  # Partition and update Halton indices
  x1sort <- base::order(x1pts)
  F <- F[x1sort]
  Hpt1 <- x1Hpts[(BoxIndex + 1)]
  Hpt2 <- x1Hpts[(BoxIndex + 1 + xlevel)]
  I <- base::order(c(Hpt1, Hpt2))

  # set halton index - top half.
  HaltonIndex[F[1:(base::floor(nF/2))]] <- (BoxIndex + (I[1] - 1) * xlevel)
  # set halton index - bottom half.
  HaltonIndex[F[((base::floor(nF/2))+1):nF]] <- (BoxIndex + (I[2] - 1) * xlevel)
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
  inBox <- base::which(HaltonIndex == BoxIndex)
  x2pts <- x2pts[inBox]

  # Randomly remove one or two points (if needed)
  F <- inBox
  # get number of points.
  nF <- base::length(F)

  # not divisible by 3.
  if ((nF %% 3) == 1) {
    r <- base::sample(nF, 1)
    HaltonIndex[F[r]] <- NA
    F <- F[-r]
    nF <- (nF - 1)
    x2pts <- x2pts[-r]
  } else if ((nF %% 3) == 2) {
    r <- base::sample(nF, 1)
    HaltonIndex[F[r]] <- NA
    F <- F[-r]
    nF <- (nF - 1)
    x2pts <- x2pts[-r]

    r <- base::sample(nF, 1)
    HaltonIndex[F[r]] <- NA
    F <- F[-r]
    nF <- (nF - 1)
    x2pts <- x2pts[-r]
  }

  # Partition and update Halton indices
  x2sort <- base::order(x2pts)
  F <- F[x2sort]
  Hpt1 <- x2Hpts[(BoxIndex + 1)]
  Hpt2 <- x2Hpts[(BoxIndex + 1 + xlevel)]
  Hpt3 <- x2Hpts[(BoxIndex + 1 + (2 * xlevel))]
  I <- base::order(c(Hpt1, Hpt2, Hpt3))

  # update halton indices - top third.
  HaltonIndex[F[1:(base::floor(nF/3))]] <- (BoxIndex + (I[1] - 1) * xlevel)
  # update halton indices - middle third.
  HaltonIndex[F[((base::floor(nF/3))+1):(2*(base::floor(nF/3)))]] <- (BoxIndex + (I[2] - 1) * xlevel)
  # update halton indices - bottom third.
  HaltonIndex[F[((2*(base::floor(nF/3)))+1):nF]] <- (BoxIndex + (I[3] - 1) * xlevel)
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
  N <- base::nrow(pts)
  pts <- base::cbind(pts, 1:N)
  # Assuming HaltonPts is a function that generates Halton points.
  Hpts <- HaltonPts(N)
  #res <- cppRSHalton_br(n = N)
  #Hpts <- res$pts
  HaltonIndex <- base::rep(0, N)

  # Partitioning parameters
  xlevel <- base::c(1, 2, 6, 12, 24, 72, 144, 432, 864, 1728, 5184, 10368, 20736)
  #xlevel <- c(1, 2, 6, 12, 24, 72, 144, 432, 864, 1728, 5184)
  # support upto 13 levels.
  diml <- base::c(1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1)
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
    TF <- base::which(!base::is.na(HaltonIndex))
    HaltonIndex <- HaltonIndex[TF]
    pts <- pts[TF,]
  }
  ptsIndex <- pts[,3]
  return(base::list(ptsIndex = ptsIndex, HaltonIndex = HaltonIndex))
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
  diml <- base::c(1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1)
  #diml <- c(1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1)
  diml <- diml[1:its]

  # Halton Indices
  b1 <- 2
  b2 <- 3
  J1 <- base::sum(diml == 1)
  J2 <- base::sum(diml == 2)
  B <- (b1^J1) * (b2^J2)
  HIM <- base::matrix(0, nrow = (b2^J2), ncol = (b1^J1))
  H <- HaltonPts(B)

  Hindex <- base::floor(base::rbind((b1^(J1) + 1e-12) * H[,1], (b2^(J2) + 1e-12) * H[,2])) + 1
  Hindex <- base::t(Hindex)

  # Halton Matrix
  for (i in 0:(B - 1)) {
    HIM[((b2^J2) + 1) - Hindex[(i + 1), 2], Hindex[(i + 1), 1]] <- i
  }

  # Permutated Halton Matrix
  step2 <- base::c(2, 4, 8, 16, 32, 64, 128, 256)
  step3 <- base::c(3, 9, 27, 81, 243)
  order2 <- base::c((base::sample(2) - 1), base::rep(NA, ((b1^J1) - b1)))
  order3 <- base::c((base::sample(3) - 1), base::rep(NA, ((b2^J2) - b2)))

  for (i in 1:(J1 - 1)) {
    if (J1 < 2) {
      next
    }
    v <- base::vector()
    L <- base::sum(!base::is.na(order2))

    for (j in (1:L)) {
      k <- order2[j]
      P <- (base::sample(2) - 1)
      s <- step2[i]
      v <- c(v, (k + s * P[1]), (k + s * P[2]))
    }
    order2[1:base::length(v)] <- v
  }
  for (i in 1:(J2 - 1)) {
    if (J2 < 2) {
      next
    }
    v <- base::vector()
    L <- base::sum(!base::is.na(order3))
    for (j in (1:L)) {
      k <- order3[j]
      P <- (base::sample(3) - 1)
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
    F <- base::which(b2vals == order2[i])
    I2 <- c(I2, F)
  }
  I3 <- vector()
  for (i in 1:(b2^J2)) {
    F <- base::which(b3vals == order3[i])
    I3 <- c(I3, F)
  }

  permHIM <- HIM[I3, I2]
  permHaltonIndex <- base::rep(0, B)

  for (i in 0:(B - 1)) {
    permHaltonIndex[i + 1] <- permHIM[base::which(HIM == i)]
  }

  return(base::list(permHaltonIndex = permHaltonIndex,
                    B               = B))
}


# delete at some stage - use cpp version instead.
HaltonPts <- function(n) {
  # Initialize
  bases <- c(2, 3)
  pts <- base::matrix(0, nrow = n, ncol = 2)

  for (ii in (1:n)) {
    for (i in (1:2)) {
      k <- (ii - 1)
      j <- 1
      xk <- (k %% bases[i]) * (1 / bases[i])

      while (base::floor(k / (bases[i] ^ (j))) > 0) {
        xk <- xk + (base::floor(k / (bases[i] ^ (j))) %% bases[i]) * (1 / (bases[i] ^ (j + 1)))
        j <- (j + 1)
      }
      pts[ii, i] <- xk
    }
  }
  return(pts)
}


#' @name is_sf_points
#'
#' @title Check if an object is an sf points object.
#'
#' @description Tests if the object passed to the function is a sf points object or not.
#' An internal only function.
#'
#' @details Detect if an object is a sf points object or not.
#'
#' @param x A probable sf points object.
#'
#' @return Either TRUE or FALSE.

is_sf_points <- function(x) {
  base::inherits(x, "sf") && base::inherits(x$geometry, "sfc_POINT")
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
#'         Population - ???

#' @export
HIP <- function(population, n = 20, iterations = 7) {

  # must be numeric, greater than 0.
  validate_parameters("hipN", c(n))
  # must be numeric, greater than 1 and less than or equal to 13.
  validate_parameters("hipIterations", c(iterations))

  # population can be either a matrix or sf point object
  sf_points_object <- FALSE
  if (is_sf_points(population)) {
    # It's an sf points object!
    sf_points_object <- TRUE
    # save original population
    sf_population <- population
    # just get coordinates from the population
    population <- sf::st_coordinates(population)
  } else if (base::is.matrix(population)) {
    # It's not an sf points object.
    # must be numeric.
    validate_parameters("hipPopulation", c(population))

    # only 2-d matrix is supported
    if(base::dim(population)[2] != 2){
      stop(base::c("uc511(HIP) Parameter population must have two dimensions."))
    }
  } else {
    # unsupported object type.
    stop(base::c("uc511(HIP) Unsupported type, the population must be defined as a sf point object or matrix."))
  }

  # Partitioning parameters
  iteration_level <- base::c(1, 2, 6, 12, 24, 72, 144, 432, 864, 1728, 5184, 10368, 20736)
  if(((base::length(population)/2)/iteration_level[iterations+1]) < 1.0){
    msg <- "uc511(HIP) Pop. size %s not compatible with number of iteration levels %s. Try pop. size of %s+ or a smaller iteration level."
    msgs <- base::sprintf(msg, base::length(population)/2, iterations, iteration_level[iterations+1])
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

  Order <- base::rep(0, base::length(HaltonIndex))

  for (i in 0:(B - 1)) {
    Order[base::which(HaltonIndex == i)] <- permHaltonIndex[i + 1]
  }

  # Assign unique indices
  ptsInBox <- (base::length(HaltonIndex) / B)
  if (ptsInBox > 1) {
    for (i in 0:(B - 1)) {
      F <- base::which(Order == i)
      Order[F] <- (i + B * (base::sample(ptsInBox) - 1))
    }
  }

  # Sample Indices
  sampleI <- base::rep(1, n)
  for (i in (1:n)) {
    sampleI[i] <- popIndex[base::which(Order == (i - 1))]
  }
  # add extra column to population sf point object
  sf_points <- NULL
  if(sf_points_object){
    sf_points <- population
  }
  # create PopulationSample from sampleI and return as a sf point object.
  ###selected_points <- sf_points[sampleI, ]
  # Convert the selected points into a new SF points object
  #PopulationSample <- sf::st_as_sf(data.frame(selected_points), coords = c("X1", "X2"))
  ###PopulationSample <- sf::st_as_sf(data.frame(selected_points), coords = c("X", "Y"))
  PopulationSample <- sf_population[sampleI,]

  # need to add the following to the sf_population sf points object:
  # the Order and HaltonIndex.
  PopulationSample$PopulationIndex <- popIndex[sampleI]
  PopulationSample$Order <- Order[sampleI]
  PopulationSample$HaltonIndex <- HaltonIndex[sampleI]

  return(list(sampleI          = sampleI,
              popIndex         = popIndex,
              Order            = Order,
              HaltonIndex      = HaltonIndex,
              Population       = sf_points,
              PopulationSample = PopulationSample))
}
