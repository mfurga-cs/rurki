#!/usr/bin/env Rscript
# Mateusz Furga <mfurga@student.agh.edu.pl>

N = 8
H = 2 / N

e <- function(i) {
  return (function (x) {
    if (x >= H * (i - 1) && x <= H * i) {
      return (1 + N / 2 * (x - 2 * i / N))
    } else if (x > H * i && x <= H * (i + 1)) {
      return (1 - N / 2 * (x - 2 * i / N))
    }
    return (0)
  })
}

de <- function(i) {
  return (function(x) {
    if (x >= H * (i - 1) && x <= H * i) {
      return (N / 2)
    } else if (x > H * i && x <= H * (i + 1)) {
      return (-N / 2)
    }
    return (0)
  })
}

# Gaussian quadrature with 2 points.
# xi = +/- 1/sqrt(3), wi = 1
int <- function(f, a, b) {
  return ((b - a) / 2 * (f((b - a) / 2 * (1 / sqrt(3)) + (a + b) / 2)
                       + f((b - a) / 2 * (-1 / sqrt(3)) + (a + b) / 2)))
}

du_dv <- function(i, j) {
  return (function(x) {
    (de(i)(x) * de(j)(x))
  })
}

B <- function(i, j) {
  # Skip the integration of the function equal to 0.
  if (abs(i - j) > 1) {
    return (0);
  }
  # Intersection of domains.
  start = max(0, H * (i - 1), H * (j - 1))
  end = min(H * (i + 1), H * (j + 1))
  return (int(du_dv(i, j), start, end) - e(0)(i) * e(0)(j))
}

L <- function(i) {
  return (-20 * e(i)(0))
}

solution <- function() {
  X <- matrix(0, nrow = N, ncol = N)
  for (i in 1 : N) {
    for (j in 1 : N) {
      X[i, j] = B(i - 1, j - 1)
    }
  }

  Y <- vector()
  for (i in 1 : N){
    Y = c(Y, L(i - 1))
  }

  U = solve(X, Y)
  return (function (x) {
    r = 0;
    for (i in 1 : N) {
      r = r + U[i] * e(i - 1)(x)
    }
    return (r)
  })
}

main <- function() {
  x = seq(0, 2, length.out = 200 * N)

  # Plot the base functions.
  title <- paste("Funkcje bazowe e(x) dla N =", as.character(N))
  plot(x, mapply(e(0), x), type = "l", main = title, xlab = "x", ylab = "ei(x)")
  for (i in 1 : N) {
    lines(x, mapply(e(i), x), col=i+10)
  }
  
  # Plot the solution u.
  u = solution()
  title <- "Rozwiązanie równania: u(x)"
  plot(x, mapply(u, x), type = "l", main = title, xlab = "x", ylab = "u(x)")
}

main()


