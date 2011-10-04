## examples from the Digital Library of Mathematical Functions
## (http://dlmf.nist.gov/15.3)

x <- seq(from=-0.023,
         to=1,
         length=1000)
plot(x,
     hyp2f1(a=5,
            b=-10,
            c=1,
            x=x),
     type="l")


x <- seq(from=-1,
         to=0.022,
         length=1000)
plot(x,
     hyp2f1(a=5,
            b=10,
            c=1,
            x=x),
     type="l")


