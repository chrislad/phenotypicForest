#' plots multiple phenotypic forests plots for cross-comparison.
#'
#' The plot is produced from a data frame containing the phenotypes, phenotype groups, SNPs, some values and possibly some upper and lower bounds.
#' By default, \code{phorest()} expects the following columns: phenotype,phenotypeGroup,SNP,SNPGroup,value,lowerBound and upperBound.
#' But this can be over-ridden with the columnNames, phenotypeGroups and SNPGroups arguments.
#' SNP are assumed to be of the form 'rsXXXXX'; this is used for ordering labels in the legend.ddddddd
#'
#' @author Christophe Ladroue
#'
#' @param associations data frame with at least columns: phenotype, value ,lowerBound,upperBound and SNP, or phenotype, value and SNP if largeSNPSet is TRUE. phenotypeGroup and SNPGroup can also be defined here.
#' @param columnNames user-defined columns names
#' @param phenotypeGroups list of phenotypes to be grouped together.
#' @param unassignedGroup default group name.
#' @param transparency transparency of the error bar. A number between 0 and 1. Can be a vector of length nrow(associations).
#' @param legendTransparency logical. Whether to show the legend of transparency.
#' @param connectingLines logical. Whether to show connecting lines for each SNP.
#' @param connectingLinesColour colour of the connecting lines.
#' @param legendSNP logical. Whether to show the legend for SNPs.
#' @param xLabel label on the x-axis
#' @param yLabel label on the y-axis
#' @param title plot title
#' @param warningOn whether warnings are shown.
#' @param verticalLine x intercepts for vertical lines
#' @param largeSNPSet logical. Whether to show an aggregate plot.
#' @param SNPGroups list of SNPs to be grouped together
#' @param aggregatingFunction function to use for aggregating the values of each group of SNPs
#' @param aggregatingLowerBound function to use for aggregating the values of each group of SNPs and get a lower bound
#' @param aggregatingUpperBound function to use for aggregating the values of each group of SNPs and get an upper bound
#' @param transparencyAggregate transparency of the connecting lines
#' @return a ggplot object
#' @export
#' @import ggplot2
#' @importFrom plyr rename
#' @importFrom plyr desc
#' @examples
#'
#'
#' set.seed(42)
#'
#' nPhenotype <- 17
#' nSNP <- 7
#' nPhenotypeGroups <- 3
#'
#' randomWord<-function(n, nLetters = 5)
#'   replicate(n,paste(sample(letters, nLetters, replace = TRUE),sep = '', collapse=''))
#'
#' df <- data.frame(
#'   phenotype      = rep(randomWord(nPhenotype),1,each=nSNP),
#'   value          = rep(1:nSNP,nPhenotype) + rnorm(nSNP*nPhenotype,mean=0,sd=1),
#'   lowerBound     = runif(nPhenotype*nSNP,min=0.0,max=0.1),
#'   upperBound     = runif(nPhenotype*nSNP,min=0.0,max=0.1),
#'   phenotypeGroup = rep(sample(toupper(randomWord(nPhenotypeGroups)),nPhenotype,replace=TRUE),1,each=nSNP),
#'   SNP            = paste('rs',rep(sample(100000,nSNP), nPhenotype), sep = '')
#' )
#'
#' df<-within(df,{
#'   lowerBound <- value - lowerBound
#'   upperBound <- value + upperBound}
#' )
#'
#' print(head(df))
#' p <- phorest(df, connectingLines = TRUE)
#' print(p)

phorest <- function(associations, columnNames = NULL, phenotypeGroups = NULL,
                    unassignedGroup = "unassigned", transparency = 0.6, legendTransparency = FALSE,
                    connectingLines = FALSE, connectingLinesColour = "#707070", legendSNP = TRUE,
                    xLabel = NULL, yLabel = NULL, title = "", warningOn = TRUE, verticalLine = NULL,
                    largeSNPSet = FALSE, SNPGroups = NULL,
                    aggregatingFunction = function(x) median(x, na.rm = TRUE),
                    aggregatingLowerBound = function(x) quantile(x, 0.25, na.rm = TRUE),
                    aggregatingUpperBound = function(x) quantile(x, 0.75, na.rm = TRUE),
                    transparencyAggregate = 0.5) {


  # Necessary columns
  if (!is.null(columnNames)) {
    namesColumn <- names(columnNames)
    names(namesColumn) <- columnNames
    associations <- rename(associations, namesColumn)
  }
  if (largeSNPSet)
    necessary <- c("phenotype", "value", "SNP")
  else
    necessary <- c("phenotype", "value", "lowerBound", "upperBound","SNP")

  if (length(intersect(necessary, names(associations))) != length(necessary))
    stop(paste("I need at least the following columns: ", paste(necessary, sep = "", collapse = ", "), sep = ""))

    ## build a 'group' column if necessary apply the group definition to a
    ## column
    applyLookup <- function(groups, keys, unassigned = "unassigned") {
      lookup <- rep(names(groups), sapply(groups, length, USE.NAMES = FALSE))
      names(lookup) <- unlist(groups, use.names = FALSE)
      p <- lookup[as.character(keys)]
      p[is.na(p)] <- unassigned
      p
    }
    # Check whether association$group exists no group columns
    if (!"phenotypeGroup" %in% names(associations)) {
      if (is.null(phenotypeGroups)) stop("No phenotype groups defined.")
      associations$phenotypeGroup <- applyLookup(phenotypeGroups, associations$phenotype, unassignedGroup)
    }

    if (!is.factor(associations$phenotypeGroup))
      associations$phenotypeGroup <- factor(associations$phenotypeGroup)

    if (largeSNPSet | (("SNPGroup" %in% names(associations) | !is.null(SNPGroups))))
      if (!"SNPGroup" %in% names(associations)) {
        # no SNPGroup column
        associations$SNPGroup <- applyLookup(SNPGroups, associations$SNP, unassignedGroup)
      }

    # sort by groups/phenotypes
    idx <- order(sapply(as.character(unique(associations$SNP)), function(s) as.numeric(tail(strsplit(s, "s")[[1]], 1)), USE.NAMES = FALSE))
    associations <- within(associations, {
      SNP <- factor(SNP, unique(SNP)[idx])
      phenotype <- as.character(phenotype)  # de-factor
    })
    associations <- arrange(associations, phenotypeGroup, desc(phenotype), SNP)
    associations$phenotype <- factor(associations$phenotype, unique(associations$phenotype))


    # adding the transparency
    if (length(transparency) %in% c(1, nrow(associations)))
      associations$transparency <- transparency else if (warningOn)
        warning("transparency has been recycled.")

    associations$transparency <- as.numeric(associations$transparency)

    # Building the plot
    p <- ggplot(associations, aes(x = value, y = phenotype))

    if (largeSNPSet) {
      # if(!is.null(SNPGroups)){
      # p<-p+geom_path(aes(group=SNP,colour=SNPGroup),alpha=transparency)

      # summaryAssociations<-ddply( associations,
      # .(SNPGroup,phenotype,phenotypeGroup), summarise,
      # aggregatedValue=aggregatingFunction(value)) # known bug in Plyr
      # (http://stackoverflow.com/questions/4289602/plyr-summarize-only-calls-global-functions)
      summaryAssociations <- aggregate(value ~ SNPGroup + phenotype + phenotypeGroup,
                                       data = associations, FUN = aggregatingFunction)

      summaryLowerBound <- aggregate(value ~ SNPGroup + phenotype + phenotypeGroup,
                                     data = associations, FUN = aggregatingLowerBound)

      summaryLowerBound <- rename(summaryLowerBound, c(value = "lowerBound"))

      summaryUpperBound <- aggregate(value ~ SNPGroup + phenotype + phenotypeGroup,
                                     data = associations, FUN = aggregatingUpperBound)

      summaryUpperBound <- rename(summaryUpperBound, c(value = "upperBound"))

      summaryAssociations <- merge(summaryAssociations, summaryLowerBound)
      summaryAssociations <- merge(summaryAssociations, summaryUpperBound)

      summaryAssociations$phenotype <- as.character(summaryAssociations$phenotype)  # de-factor
      summaryAssociations <- arrange(summaryAssociations, phenotypeGroup,
                                     desc(phenotype))
      summaryAssociations$phenotype <- factor(summaryAssociations$phenotype,
                                              unique(summaryAssociations$phenotype))
      idx <- order(as.character(unique(summaryAssociations$SNPGroup)))
      summaryAssociations$SNP <- factor(summaryAssociations$SNPGroup,
                                        unique(summaryAssociations$SNPGroup)[idx])

      p <- p + geom_path(aes(x = value, y = phenotype, group = SNPGroup, colour = SNPGroup),
                         alpha = transparencyAggregate, size = 1.2,
                         data = summaryAssociations)

      # building polygon data
      datapoly <- data.frame()

      for (g in unique(summaryAssociations$SNPGroup)) {
        p1 <- subset(summaryAssociations, SNPGroup == g)
        datapoly <- rbind(datapoly, data.frame(SNPGroup = g, value = c(p1$lowerBound,
                                                                       rev(p1$upperBound)), phenotype = c(as.character(p1$phenotype),
                                                                                                          rev(as.character(p1$phenotype))), vertexOrder = 1:(2 *
                                                                                                                                                               nrow(p1))))
      }

      datapoly$phenotype <- as.character(datapoly$phenotype)  # de-factor
      datapoly <- merge(datapoly, subset(associations, select = c("phenotype", "phenotypeGroup"), by = "phenotype"))
      datapoly <- arrange(datapoly, phenotypeGroup, SNPGroup, vertexOrder)
      datapoly$phenotype <- factor(datapoly$phenotype, unique(summaryAssociations$phenotype))

      idx <- order(as.character(unique(datapoly$SNPGroup)))
      datapoly$SNPGroup <- factor(datapoly$SNPGroup, unique(datapoly$SNPGroup)[idx])

      p <- p + geom_polygon(data = datapoly, aes(x = value, y = phenotype, group = SNPGroup, fill = SNPGroup), alpha = 0.2)

      if (!legendSNP)
        p <- p + scale_color_hue(guide = "none") + scale_fill_hue(guide = "none")

      # } else
      # p<-p+geom_path(aes(group=SNP),colour='white',alpha=transparency)


      p <- p + theme(panel.background = element_rect(fill = "black"), panel.grid.major = element_line(colour = "black"),
                    panel.grid.minor = element_line(colour = "black"))
      p <- p + theme(legend.position = "right", legend.direction = "vertical")

      if (!is.null(verticalLine))
        p <- p + geom_vline(xintercept = verticalLine, colour = "white")
    } else {
      if ("SNPGroup" %in% names(associations))
        p <- p + geom_segment(aes(x = lowerBound, xend = upperBound, y = phenotype, yend = phenotype, colour = SNPGroup, alpha = transparency), size = 3)
      else
        p <- p + geom_segment(aes(x = lowerBound, xend = upperBound, y = phenotype, yend = phenotype, colour = SNP, alpha = transparency), size = 3)

        p <- p + geom_point(aes(x = value, y = phenotype), size = 0.7, colour = "#000000")

        if (connectingLines)
          if ("SNPGroup" %in% names(associations))
            p <- p + geom_path(aes(x = value, y = phenotype, group = SNP, colour = SNPGroup))
          else
            p <- p + geom_path(aes(x = value, y = phenotype, group = SNP, colour = SNP))

        # Legends
        if (!legendSNP)
          p <- p + scale_color_hue(guide = "none")

        if (!legendTransparency)
          p <- p + scale_alpha_identity(guide = "none")

        if (!is.null(verticalLine))
          p <- p + geom_vline(xintercept = verticalLine, colour = "black")

    }

    p <- p + facet_grid(phenotypeGroup ~ ., scales = "free", space = "free") + theme(strip.text.y = element_text())

    # Labels
    if (!is.null(xLabel))
      if (is.na(xLabel))
        p <- p + theme(axis.title.x = element_blank())
      else
        p <- p + xlab(xLabel)

    if (!is.null(yLabel))
        if (is.na(yLabel))
            p <- p + theme(axis.title.y = element_blank())
        else p <- p + ylab(yLabel)

    p <- p + ggtitle(title)
    p
}
# transparencyScale: scale_alpha_identity(legend=FALSE) linetype:
# scale_linetype(legend=FALSE)
