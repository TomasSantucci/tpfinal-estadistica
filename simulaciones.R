library(markovchain)
library(expm)
library(lattice)
library(grid)

squarem <- function(values) {
  size <- sqrt(length(values))
  matrix(values, nrow = size, ncol = size, byrow = TRUE)
}

e1 <- function(variables, samples) {
  tabla <- matrix(
    data = runif(variables * samples, 0, 1),
    nrow = variables,
    ncol = samples)

  data <- apply(tabla, 2, sum)
  mean <- variables * 1 / 2
  sd <- sqrt(variables * 1 / 12)
 
  # se asemeja a una normal de esperanza 20 y desvio estandar 1.82 por el TCDL 
  histogram(
    data,
    type = "density",
    panel = function(x, ...) {
      panel.histogram(x, ...)
      panel.mathdensity(
        dmath = dnorm,
        col = "black",
        lwd = 2,
        args = list(mean, sd))
    })
}

e1(40, 10000)

# Ejercicio 2

ej2Matrix <- function(p) {
  squarem(c(
    1-p, p, 0, 0,
    1-p, 0, p, 0,
    1-p, 0, 0, p,
    0  , 0, 0, 1))
}

ej2 <- function(repetitions, p) {
  P <- ej2Matrix(p)
  mc <- new("markovchain", transitionMatrix = P)
  replicate(
    repetitions,
    length(simulateMC(mc, "1", c("4"))) - 1)
}

simulateMC <- function(mc, initial, tfs) {
  step <- initial
  steps <- c(step)
  while (!(step %in% tfs)) {
    step <- rmarkovchain(1, mc, t0 = step)
    steps <- c(steps, step)
  }
  steps
}

p <- 0.6

P <- ej2Matrix(p)
mc <- new("markovchain", transitionMatrix = P)
# esperanza estimada de la variable:
mean(ej2(1000, p))

# Para calcularla exactamente es:
# obtenemos la matriz Q de P

Q <- squarem(c(
  1-p, p, 0,
  1-p, 0, p,
  1-p, 0, 0))

I <- diag(3)
# calculamos (I-Q)^(-1)
N <- solve(I - Q)
N
# la cantidad esperada de pasos antes de llegar al estado absorbente es N*1
EN <- N %*% c(1, 1, 1)
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

s <- 19
k <- (s + 1) / 2
absorbents <- c("1", toString(s))

zeroes <- function(n) {
  rep(0, n)
}

ej4Matrix <-function(p) {
  n <- c(1, zeroes(s - 1))
  
  for (i in 2:(s-1)) {
    r <- zeroes(s)
    r[i-1] = 1-p
    r[i+1] = p
    n <- c(n, r)
  }
  
  n <- c(n, zeroes(s - 1), 1)
  
  squarem(n)
}

transitionToRecurrentProbabilities <- function(mc, size) {
  mc <- canonicForm(mc)
  recurrent <- length(recurrentStates(mc))
  
  Q <- mc[(recurrent + 1):size, (recurrent + 1):size]
  B <- mc[(recurrent + 1):size, 1:recurrent]
  
  I <- diag(size - recurrent)
  S <- solve(I - Q)
  
  S %*% B
}

ej4a <- list()
ej4b <- c()
ej4bextra <- c()
for (p in c(0.4, 0.5, 0.6)) {
  matrix <- ej4Matrix(p)
  mc <- new("markovchain", transitionMatrix = matrix)
  
  # trayectorias de ejemplo
  ej4a <- append(ej4a, list(simulateMC(mc, toString(k), absorbents)))
  
  # forma empírica
  r <- replicate(1000, tail(simulateMC(mc, toString(k), absorbents), 1))
  ej4b <- c(ej4b, length(r[r == "1"]) / 1000)
  
  # forma teórica
  G <- transitionToRecurrentProbabilities(mc, size = s)
  
  ej4bextra <- c(ej4bextra, G[k - 1, 1])
}
ej4a
ej4b
ej4bextra

# Ejercicio 5

P <- function(i, j, k) {
  choose(k, j) * (i/k)^j * (1-i/k)^(k-j)
}

makeP <- function(k) {
  v <- c()
  for (i in 0:k)
    for (j in 0:k)
      v <- c(v, P(i, j, k))
  
  squarem(v)
}

k <- 20
mc <- new("markovchain", transitionMatrix = makeP(k))

# simulacion del proceso
rmarkovchain(100, mc, t0 = "10")

# calculado teóricamente
mc <- new("markovchain", transitionMatrix = makeP(k))
transitionToRecurrentProbabilities(mc, size = k + 1)

# EJERICIO 6

m <- squarem(c(
  0,   0.5, 0.5, 0,   0,
  0.2, 0.2, 0.2, 0.2, 0.2,
  1/3, 1/3, 0,   1/3, 0,
  0,   0,   0,   0,   1,
  0,   0,   0.5, 0.5, 0))

mc <- new("markovchain", transitionMatrix = m)
plot(mc)
hittingProbabilities(mc)
is.irreducible(mc)
period(mc)
steadyStates(mc)
recurrentStates(mc)
# la probabilidad de visitar cualquier pagina arrancando de cualquiera es 1
# porque la mc es irreducible, todos los estados son alcanzables para todos

# EJERCICIO 7

m <- squarem(c(
  0,    0,    0,    0,    1,
  0,    8/13, 3/13, 1/13, 1/13,
  1/16, 3/16, 3/8,  1/4,  1/8,
  0,    1/11, 4/11, 5/11, 1/11,
  0,    1/8,  1/2,  1/8,  1/4))

mc <- new("markovchain", transitionMatrix = m)

# Probabilidades despues de 7 semanas
c(0, 0, 0, 0, 1) * mc^(7)

is.irreducible(mc)
period(mc)
# Como la cm es irreducible y aperiodica la distribucion invariante
# es unica y coincide con la distribucion limite
steadyStates(mc)

# EJERCICIO 8

eventTimes <- function(lambda, length) {
  eventsAmount <- rpois(1, lambda * length)
  sort(runif(eventsAmount, 0, length))
}

ladderPlot <- function(t, length) {
  x <- c(0, t, length)
  y <- 0:(length(x) - 1)
  y[length(y)] <- length(y) - 2
  plot(x, y, type = "s", xlab = "longitud del tubo en km", ylab = "numero de defectos")
}

distanceVariable <- function(t) {
  d <- c()
  for (i in 2:length(t))
    d <- c(d, t[i] - t[i-1])
  d
}

# item a y b
t <- eventTimes(0.1, 10)
ladderPlot(t, 10)

# 8c

lambda <- 0.1
len <- 10000

t <- eventTimes(lambda, len)
d <- distanceVariable(t)

myDist <- function(x) {
  qexp(x, lambda)
}

# falta hacer un mejor analisis de d
# la esperanza de d es aprox 1/lambda
mean(d)
histogram(d)
qqmath(d, distribution = myDist)