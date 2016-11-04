smoothMatrix <-
function(M, n=1, direction = 0, torus = FALSE) { # direction = 1 smoothing rows only, 2 cols, 3 diagonal from left-top to right-bottom, 4 diagonal from left-bottom to right-top; can be combined
  Msmooth <- M
  if (n > 0) {
    for (i in 1:n) {
      
      cat('\r', i)
      flush.console()
      
      Mtmp <- Msmooth
      M1 <- M2 <- M3 <- M4 <- M5 <- M6 <- M7 <- M8 <- M*0
      if (torus) {
        if (any(direction %in% c(0,1))) {
          M1 <- cbind(Msmooth[, ncol(M)], Msmooth[, 1:(ncol(M)-1)])
          M5 <- cbind(Msmooth[, 2:ncol(M)], Msmooth[, 1])
        }
        if (any(direction %in% c(0,2))) {
          M3 <- rbind(Msmooth[nrow(M), ], Msmooth[1:(nrow(M)-1), ])
          M7 <- rbind(Msmooth[2:nrow(M), ], Msmooth[1, ])
        }
        if (any(direction %in% c(0,3))) {
          M2 <- rbind(cbind(Msmooth[, ncol(M)], Msmooth[, 1:(ncol(M)-1)])[nrow(M), ], cbind(Msmooth[, ncol(M)], Msmooth[, 1:(ncol(M)-1)])[1:(nrow(M)-1), ])
          M6 <- rbind(cbind(Msmooth[, 2:ncol(M)], Msmooth[, 1])[2:nrow(M), ], cbind(Msmooth[, 2:ncol(M)], Msmooth[, 1])[1, ])
        }
        if (any(direction %in% c(0,4))) {
          M4 <- rbind(cbind(Msmooth[, 2:ncol(M)], Msmooth[, 1])[nrow(M), ], cbind(Msmooth[, 2:ncol(M)], Msmooth[, 1])[1:(nrow(M)-1), ])
          M8 <- cbind(rbind(Msmooth[2:nrow(M), ], Msmooth[1, ])[, ncol(M)], rbind(Msmooth[2:nrow(M), ], Msmooth[1, ])[, 1:(ncol(M)-1)])
        }
      } else {
        if (any(direction %in% c(0,1))) {
          M1 <- cbind(0, Msmooth[, 1:(ncol(M)-1)])
          M5 <- cbind(Msmooth[, 2:ncol(M)], 0)
        }
        if (any(direction %in% c(0,2))) {
          M3 <- rbind(0, Msmooth[1:(nrow(M)-1), ])
          M7 <- rbind(Msmooth[2:nrow(M), ], 0)
        }
        if (any(direction %in% c(0,3))) {
          M2 <- rbind(0, cbind(0, Msmooth[, 1:(ncol(M)-1)])[1:(nrow(M)-1), ])
          M6 <- rbind(cbind(Msmooth[, 2:ncol(M)], 0)[2:nrow(M), ], 0)
        }
        if (any(direction %in% c(0,4))) {
          M4 <- rbind(0, cbind(Msmooth[, 2:ncol(M)], 0)[1:(nrow(M)-1), ])
          M8 <- cbind(0, rbind(Msmooth[2:nrow(M), ], 0)[, 1:(ncol(M)-1)])
        }
      }
      denom <- matrix(9, nrow(M), ncol(M))
      Msmooth <- Mtmp+M1+M2+M3+M4+M5+M6+M7+M8
      Msmooth <- Msmooth/denom
      if (all(Mtmp == Msmooth)) {
        print("convergence")
        break()
      }
    }
  }
  return(Msmooth)
}
