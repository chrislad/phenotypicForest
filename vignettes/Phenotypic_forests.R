## ------------------------------------------------------------------------
set.seed(42)

nPhenotype <- 17
nSNP <- 7
nPhenotypeGroups <- 3

randomWord<-function(n, nLetters = 5)
   replicate(n,paste(sample(letters, nLetters, replace = TRUE),sep = '', collapse=''))

df <- data.frame(
  phenotype      = rep(randomWord(nPhenotype), 1, each = nSNP),
  value          = rep(1:nSNP, nPhenotype) + rnorm(nSNP * nPhenotype, mean = 0, sd = 0.1),
  lowerBound     = runif(nPhenotype * nSNP, min = 0.0, max = 0.1),
  upperBound     = runif(nPhenotype * nSNP, min = 0.0, max = 0.1),
  phenotypeGroup = rep(sample(toupper(randomWord(nPhenotypeGroups)), nPhenotype, replace = TRUE), 1, each = nSNP),
  SNP            = paste('rs',rep(sample(100000, nSNP), nPhenotype), sep = '')
  )

df <- within(df, {
  lowerBound <- value - lowerBound
  upperBound <- value + upperBound}
           )

print(head(df))

## ------------------------------------------------------------------------
#phorest(df, connectingLines = TRUE)

