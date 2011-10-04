## artificial configuration:
config <- c(a=0L, b=1L, c=3L, d=0L)

## get log prior probs:
getLogModelPrior(config=config,
                 type="flat")

getLogModelPrior(config=config,
                 type="exponential")
