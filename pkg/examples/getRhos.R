x <- seq(from=-5,
         to=5,
         length=100)
Z <- makeBasis(x=x,
               nKnots=4L)
lambdas <- svd(Z, nu=0L, nv=0L)$d^2
getRhos(lambdas=lambdas,
        degrees=c(1, 2, 3, 3.9))

