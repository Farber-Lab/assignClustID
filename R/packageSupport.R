# Supporting imports and settings for the package

# You can either use 'import X' or 'importFrom X abcdefg'. importFrom  is
# better practice, but for dplyr and ggplot2 we were simply importing so many functions
# that it makes  sense to just import the whole package
#' @import dplyr
#' @importFrom data.table fread
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#' @importFrom tidyr separate
#' @importFrom tidyr drop_na
#' @importFrom tidyselect all_of
#' @importFrom stringr str_replace
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ComplexHeatmap draw
#' @importFrom tibble column_to_rownames
#' @importFrom grid gpar
#' @importFrom grid grid.text
#' @importFrom klippy klippy
#' @importFrom utils head

# Because of some issues with NOTEs on R CMD check and CRAN submission,
# (see here: http://stackoverflow.com/questions/9439256/)
# I have to register stuff used in data.table as non-standard evaluation,
# in order to pass some R CMD check NOTES.
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c("ALT", "DO_mm10", "FILTER", "FORMAT", "ID", "QUAL",
                           "REF", "a1", "a1sort", "a2", "a2sort", "chr",
                           "cluster", "corrGenotype", "genotype", "info",
                           "marker", "position", "toss", "x", "y",
                           "comparisonTable", "comparisonTableNoRevComp"))
}
