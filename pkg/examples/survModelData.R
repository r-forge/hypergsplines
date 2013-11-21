## get some data
library(survival)
pbc <- na.omit(pbc)
pbc$sex <- as.numeric(pbc$sex == "f")

## try the function
md <- survModelData(times=pbc$time,
                    X=
                    as.matrix(subset(pbc,
                                     select=
                                     c(trt,
                                       age,
                                       sex,
                                       ascites,
                                       hepato,
                                       spiders,
                                       edema,
                                       bili,
                                       chol,
                                       albumin,
                                       copper,
                                       alk.phos,
                                       ast,
                                       trig,
                                       platelet,
                                       protime,
                                       stage))),
                    observed=
                    pbc$status == 2,
                    continuous=
                    c(FALSE,
                      TRUE,
                      FALSE,
                      FALSE,
                      FALSE,
                      FALSE,
                      FALSE,
                      TRUE,
                      TRUE,
                      TRUE,
                      TRUE,
                      TRUE,
                      TRUE,
                      TRUE,
                      TRUE,
                      TRUE,
                      FALSE),
                    nKnots=6L,
                    splineType="cubic",
                    gPrior="hyper-g/n")

## look at the results
str(md)
