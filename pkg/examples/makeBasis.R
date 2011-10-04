x <- seq(from=-5,
         to=5,
         length=1000)

Z <- makeBasis(x=x,
               nKnots=4L,
               type="linear")
str(Z)

plot(x,
     Z %*% c(1, -3, 2, 1),
     type="l")
abline(v=attr(Z, "knots"),
       col="gray",
       lty=2)

Z <- makeBasis(x=x,
               nKnots=10L,
               type="cubic")
dimZ <- ncol(Z)

coefs <- rnorm(dimZ)

plot(x,
     Z %*% coefs,
     type="l")
abline(v=attr(Z, "knots"),
       col="gray",
       lty=2)
