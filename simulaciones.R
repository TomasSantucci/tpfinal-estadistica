library(markovchain)
library(expm)
library(lattice)
library(grid)

e1 <- function(cant) {
  n <- 40
  
  tabla <- matrix(data = runif(n*cant,0,1), nrow = n, ncol = cant)
  
  tabla
  
  esperanza <- 20
  desvio_estandar <- 1.82
  data <- apply(tabla, 2, sum)
 
  # se asemeja a una normal de esperanza 20 y desvio estandar 1.82 por el TCDL 
  histogram(data,
            type = "density",
            panel = function(x, ...) {
              panel.histogram(x, ...)
              panel.mathdensity(dmath = dnorm, col = "black", lwd = 2, args = list(mean=esperanza,sd=desvio_estandar))
            })
  
}

# Ejercicio 2

ej2Matrix <- function(p) {
  r1 <- c(1-p,p,0,0)
  r2 <- c(1-p,0,p,0)
  r3 <- c(1-p,0,0,p)
  r4 <- c(0,0,0,1)
  P <- matrix(c(r1,r2,r3,r4), nrow = 4, ncol = 4, byrow = TRUE)
  P
}

ej2 <- function(cant, p, k) {
  P <- ej2Matrix(p)
  mc <- new("markovchain", transitionMatrix = P)
  
  trayectorias <- c()
  for (i in 1:cant) {
    trayectorias <- c(trayectorias, simularTrayectoria(mc))
  }
  trayectorias
}

simularTrayectoria <- function(mc) {
  pasos <- 1
  sim <- rmarkovchain(n = 1, object = mc, t0 = "1")
  while (sim != "4") {
    pasos <- pasos + 1
    sim <- rmarkovchain(n = 1, object = mc, t0 = sim)
  }
  pasos
}

# esperanza estimada de la variable:
sum(ej2(1000, 0.6, 3))/1000

# Para calcularla exactamente es:
# obtenemos la matriz Q de P
p <- 0.6
r1 <- c(1-p,p,0)
r2 <- c(1-p,0,p)
r3 <- c(1-p,0,0)
Q <- matrix(c(r1,r2,r3), nrow = 3, ncol = 3, byrow = TRUE)
I <- diag(3)
# calculamos (I-Q)^(-1)
N <- solve(I-Q)
# la cantidad esperada de pasos antes de llegar al estado absorbente es N*1
EN <- N%*%c(1,1,1)
EN[1]

# Ejercicio 3

ej3a <- function(p) {
  # El proceso es basicamente un ensayo de bernoulli
  rbinom(1,1,p)
}

ej3b <- function(n,p) {
  # El proceso es basicamente una binomial
  rbinom(1,n,p)
}

# ej3c: son las esperanzas y variancias conocidas de bernoulli y binomial

# Ejercicio 4

ej4MC <- function(s, p) {
  n <- c(1, rep(0,s-1))
  for (i in 2:s) {
    if (i < s) {
      r <- rep(0, s)
      r[i-1] = 1-p
      r[i+1] = p
    } else {
      r <- c(rep(0,s-1), 1)
    }
    n <- c(n,r)
  }
  
  mat <- matrix(data = n, nrow = s, ncol = s, byrow = TRUE)
  new("markovchain", transitionMatrix = mat)
}

trayectoria <- function(mc, n, t0) {
  lastState <- tail(states(mc), n=1)
  sim <- rmarkovchain(1, object = mc, t0 = t0)
  while ((sim != "1") & (sim != lastState)) {
    sim <- rmarkovchain(1, object = mc, t0 = sim)
  }
  sim
}

ej4b <- function(cant, mc, t0) {
  ruinas <- 0
  for (i in 1:cant) {
    t <- "1"
    t <- trayectoria(mc, 1000, t0)
    if (t == "1") {
      ruinas <- ruinas + 1
    }
  }
  ruinas / cant
}

mc <- ej4MC(50, 0.6)
ej4b(1000, mc, "5")
absorptionProbabilities(mc)

# falta hacerlo de manera teorica


# Ejercicio 5

comb <- function(n, x) {
  factorial(n) / factorial(n-x) / factorial(x)
}

p <- function(i, j, k) {
  t1 <- comb(k, j)
  t2 <- (i/k)^j
  t3 <- (1-i/k)^(k-j)
  t1*t2*t3
}

row <- function(row, k) {
  r <- c()
  for (i in 0:k) {
    r <- c(r, p(row,i,k))
  }
  r
}

makeP <- function(k) {
  v <- c()
  for (i in 0:k) {
    r <- row(i, k)
    v <- c(v, r)
  }
  mat <- matrix(data = v, nrow = k+1, ncol = k+1, byrow = TRUE)
  mat
}

k <- 20
mc <- new("markovchain", transitionMatrix = makeP(k))

# simulacion del proceso
rmarkovchain(100, object = mc, t0 = "10")

# probabilidad de caer en los estados absorbentes
absorptionProbabilities(mc)

# manualmente es N*R
k <- 20
mc <- new("markovchain", transitionMatrix = makeP(k))
mc <- canonicForm(mc)
Q <- mc[2:k+1, 2:k+1]
R <- mc[2:k+1,1:2]
I <- diag(k-1)
N <- solve(I-Q)
N%*%R

# EJERICIO 6

entries <- c(0,0.5,0.5,0,0, 0.2,0.2,0.2,0.2,0.2, 1/3,1/3,0,1/3,0, 0,0,0,0,1, 0,0,0.5,0.5,0)
transitionMat <- matrix(data = entries, nrow = 5, ncol = 5, byrow = TRUE)
mc <- new("markovchain", transitionMatrix = transitionMat)
plot(mc)
hittingProbabilities(mc)
is.irreducible(mc)
recurrentStates(mc)
# la probabilidad de visitar cualquier pagina arrancando de cualquiera es 1
# porque la mc es irreducible, todos los estados son alcanzables para todos

# EJERCICIO 7

r1 <- c(0,0,0,0,1)
r2 <- c(0,8/13,3/13,1/13,1/13)
r3 <- c(1/16,3/16,3/8,1/4,1/8)
r4 <- c(0,1/11,4/11,5/11,1/11)
r5 <- c(0,1/8,1/2,1/8,1/4)
entries <- c(r1,r2,r3,r4,r5)
mat <- matrix(data = entries, nrow = 5, ncol = 5, byrow = TRUE)
mc <- new("markovchain", transitionMatrix = mat)

# Probabilidades despues de 7 semanas
mc[5,] * mc^(7)

is.irreducible(mc)
period(mc)
# Como la cm es irreducible y aperiodica la distribucion invariante
# es unica y coincide con la distribucion limite
steadyStates(mc)

# EJERCICIO 8

eventTimes <- function(lambda, length) {
  eventsAmount <- rpois(1, lambda*length)
  
  t <- runif(eventsAmount, 0, length)
  sort(t)
}

ladderPlot <- function(t, length) {
  x <- c(0, t, length)
  y <- 0:(length(x)-1)
  y[length(y)] <- length(y) - 2
  plot(x, y, type = "s", xlab = "longitud del tubo en km", ylab = "numero de defectos")
}

distanceVariable <- function(t) {
  d <- c()
  for (i in 2:length(t)) {
    d <- c(d, t[i]-t[i-1])
  }
  d
}

# item a y b
t <- eventTimes(0.1, 10)
ladderPlot(t, 10)

# 8c
t <- eventTimes(0.1, 100)
d <- distanceVariable(t)
# falta hacer un mejor analisis de d
# la esperanza de d es aprox 1/lambda
mean(d)

