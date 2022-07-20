e1 <- function(cant) {
  n <- 40
  
  tabla <- matrix(data = runif(n*cant,0,1), nrow = n, ncol = cant)
  
  tabla
  
  esperanza <- 20
  desvio_estandar <- 1.82
  data <- apply(tabla, 2, sum)
  
  histogram(data,
            type = "density",
            panel = function(x, ...) {
              panel.histogram(x, ...)
              panel.mathdensity(dmath = dnorm, col = "black", lwd = 2, args = list(mean=esperanza,sd=desvio_estandar))
            })
  
}

e2 <- function(cant, p, k) {
  vector <- seq(from = 1, to = cant)
  
  n_ensayos <- rep(0, cant)
  cons <- rep(0, cant)
  
  for (i in vector) {
    while (cons[i] < k) {
      ensayo <- rbinom(1, 1, p)
      
      n_ensayos[i] <- n_ensayos[i] + 1
      if (ensayo == 1) {
        cons[i] <- cons[i] + 1
      } else {
        cons[i] <- 0
      }
    }    
  }
  
  n_ensayos
}

e3a <- function(p) {
  for (i in vector) {
    exito <- rbinom(1,1,p)
  }
  exito
}

e3b <- function(n) {
  i <- 0
  n_senales <- 0
    
  while (i < n) {
    exito <- rbinom(1,1,0.5)
    n_senales <- n_senales + exito
    i <- i + 1
  }
  
  n_senales
}

e4a <- function(k, s, p) {
  trayectoria <- c()
  while (k < s & k > 0) {
    exito <- rbinom(1,1,p)
    
    if (exito == 1) {
      k <- k + 1
    } else {
      k <- k - 1
    }
    
    trayectoria <- c(trayectoria, k)
  }
  trayectoria
}

e4b <- function(cant, k, s, p) {
  ruinas <- 0
  i <- 0
  while (i < cant) {
    res <- tail(e4a(k, s, p), n = 1)
    if (res == 0) {
      ruinas <- ruinas + 1
    }
    i <- i + 1
  }
  ruinas / cant
}