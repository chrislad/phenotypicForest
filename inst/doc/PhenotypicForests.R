## ------------------------------------------------------------------------
randomWord<-function(n, nLetters = 5)
   replicate(n,paste(sample(letters, nLetters, replace = TRUE),sep = '', collapse=''))

toyData <- function(nPhenotype = 17, nSNP = 7, nPhenotypeGroups = 3 ){
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
  df
}

set.seed(42)
nPhenotype <- 17
nSNP <- 7
nPhenotypeGroups <- 3

df <- toyData(nPhenotype, nSNP , nPhenotypeGroups)
print(head(df))

## ----, fig.width = 6,fig.height = 3--------------------------------------
library(phenotypicForest)
phorest(df, connectingLines = TRUE)

## ----, fig.width = 6,fig.height = 3--------------------------------------
 set.seed(42)
 df <- toyData(10,1,3)
 print(df)
 phorest(df)

## ----, fig.width = 6,fig.height = 3--------------------------------------
 df$phenotypeGroup <- NULL # delete phenotypeGroup column
  
 userDefined <- list(
   "JONOSE" = c("xyhvq", "ntdrs", "lsygm", "yzdmo"),
   "RYBYGU" = c("xdzyc", "nkxlv", "tvkra"),
   "GYKUXU" = c("vafxp", "jlazl", "yxqzq"))

  phorest(df, phenotypeGroups = userDefined)

## ----, fig.width = 8,fig.height = 4--------------------------------------
# creating some random SNP groups
set.seed(42)
nSNP <- 200
nSNPGroup <- 4
df <- toyData(17,nSNP,4)
df$value <- rnorm(nrow(df))
tmp <- data.frame(
  SNP = unique(df$SNP),
  SNPGroup = sample(paste0("SNP Group#", 1:nSNPGroup), nSNP, replace = TRUE))
df <- merge(df, tmp, by = "SNP")

head(df)

phorest(
  df,
  largeSNPSet = TRUE,
  title = 'default plot for large SNP sets')

## ----, fig.width = 8,fig.height = 4--------------------------------------
# specifying the aggregating functions
phorest(
  df,
  largeSNPSet = TRUE,
  aggregatingFunction = function(x) mean(x, na.rm = TRUE),
  aggregatingLowerBound = function(x) mean(x, na.rm = TRUE) - sd(x, na.rm = TRUE),
  aggregatingUpperBound = function(x) mean(x, na.rm = TRUE) + sd(x, na.rm = TRUE),
  title = 'mean and standard deviation')

## ----, fig.width = 7, fig.height = 7-------------------------------------
set.seed(42)
nFamily <- 20
nItemPerFamily <- sample(1:6, nFamily,replace = TRUE)
nValues <- 3

df <- data.frame(
  family = rep( randomWord(nFamily), times = nValues * nItemPerFamily),
  item   = rep( randomWord(sum(nItemPerFamily), 3), each = nValues ),
  score  = rep( paste0("V",1:nValues), times = sum(nItemPerFamily)),
  value  = round(500 * runif( sum(nItemPerFamily * nValues)),2))

print(head(df))
 
polarHistogram(df, familyLabel = FALSE)

## ------------------------------------------------------------------------
df[df$item == 'edm',]

## ------------------------------------------------------------------------
100*df[df$item == 'edm',]$value/sum(df[df$item == 'edm',]$value)

## ----, fig.width=7, fig.height=7-----------------------------------------
polarHistogram(df, normalised = FALSE, guides = c(100,300,500,1000))

## ------------------------------------------------------------------------
  # not run
  # p <- phorest(df, columnNames = c("SNP" = "snpid", "phenotype" = "assays))

## ------------------------------------------------------------------------
#  p <- polarHistogram(df)
#  p <- p + ggtitle("put title here") + xlab("label for x-axis") + ylab("label for y-axis")  
#  print(p)

