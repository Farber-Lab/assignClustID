#' clustSimilarity: main function to compare the similarity between genotypes
#' reported by souporcell and our known genotypes
#'
#' @param pathToVCF A string specifying a path to a VCF file produced by
#'    souporcell. Provide full path including the file name.
#' @param genotypes A data.frame or string. String: current string option
#'    is only "DO_mm10" (default) - automatically retrieves genotypes of 821
#'    DO mice profiled in Farber lab. / Data.frame: You can provide your
#'    genotype information. The first three columns are: marker (marker name),
#'    chr (chromosome including "chr" prefix), position (marker position on a
#'    chromosome). Following columns are assigned IDs. Numerical IDs must have
#'    "X" prefix.
#' @param IDs Optional vector of strings specifying knwon IDs from genotype
#'    data.frame. If IDs are numerical use "X" prefix in front of each ID.
#' @param compareIncomplete Logical variable specifying if only positions with
#'    complete genotype information from known genotypes and souporcell output
#'    should be used to calculate similarity. TRUE = compare markers for which
#'    genotypes are NOT available for all tested individuals, FALSE =
#'    compare only markers with genotypes available across all individials -
#'    both from souporcell and known genotypes (default - FALSE, which means
#'    only complete cases will be considered)
#' @param minQual An optional numerical variable specifying minimal QUAL
#'    threshold for genotypes from souporcell vcf output to be considered
#' @return data.frame containing similarities (proportional matches between
#'    genotypes) between souporcell clusters (rows) and known genotypes
#'    (columns)
#'
#' @export
#' @examples
#' VCFpath = system.file("extdata", "cluster_genotypes.vcf",
#'                             package = "assignClustID")
#' genotypes = DO_mm10
#' similarities = clustSimilarity(pathToVCF=VCFpath, genotypes=genotypes)
clustSimilarity = function(pathToVCF,
                           genotypes,
                           IDs=NULL,
                           compareIncomplete=FALSE,
                           minQual=NULL){
  # load VCF file (and check if chromosomes start with "chr")
  vcf = loadVCF(pathToVCF)

  # quality filter is set, apply it
  if (!is.null(minQual)){
    suppressWarnings({
    vcf = vcf %>%
      mutate(QUAL = as.numeric(QUAL)) %>%
      filter(!is.na(QUAL) & QUAL >= minQual)
    })
  }

  # check for "genotype" input type
  if (!is.data.frame(genotypes)){
    stop("Input to genotypes variable must be a data.frame.")
  }

  genotype_SNPS = genotypes %>%
    select(marker, chr, position)

  # extract genotypes from vcf file:
  # extract genotype for individual clusters
  vcfEdit = vcf %>%
    # remove backgound based on filter value
    dplyr::filter(FILTER == ".") %>%
    # keep only variants that were previously genotyped
    inner_join(genotype_SNPS, by = c("chr", "position")) %>%
    # keep only necessary columns - extract the genotype information:
    # 0 = ref, 1 = alt
    dplyr::select(-c(chr:ID, QUAL:FORMAT)) %>%
    pivot_longer(-c(marker, REF, ALT),
                 names_to = "cluster", values_to = "info") %>%
    tidyr::separate(info, c("genotype", "toss"), sep = ":", extra = "drop") %>%
    dplyr::select(-toss) %>%
    # filter out unknown genotypes
    dplyr::filter(genotype != "./.") %>%
    tidyr::separate(genotype, c("x", "y"), sep = "/") %>%
    # recreate genotypes based on information provided
    mutate(a1 = ifelse(x == "1", ALT, REF)) %>%
    mutate(a2 = ifelse(y == "1", ALT, REF)) %>%
    # final genotype info is alphabetically sorted
    mutate(a1sort = pmin(a1, a2)) %>%
    mutate(a2sort = pmax(a1, a2)) %>%
    mutate(genotype = paste0(a1sort, a2sort)) %>%
    # select only desired variable and convert back to wide format
    dplyr::select(marker, cluster, genotype) %>%
    # add prefix before cluster number
    mutate(cluster = paste0("clust_", cluster)) %>%
    pivot_wider(names_from = cluster, values_from = genotype) %>%
    mutate(marker = as.character(marker))

  if (!compareIncomplete){
    # keep only rows with genotype information for all clusters
    vcfEdit = vcfEdit %>%
      drop_na()
  }


  #---
  # # get all possible genotypes from DO mice (to see if order matters:
  # # e.g. AT/TA)
  # allGens = genotypes %>%
  #   select(-c(chr, position)) %>%
  #   pivot_longer(!marker, names_to = "mouse", values_to = "gen")
  # levels(factor(allGens$gen))
  # # everything as alphabetically sorted only GC and CG can be
  # # either way -> edit that in the table with selected mice
  # # to only CG

  # if specific IDs were provided - select just those
  if (!is.null(IDs)){
    IDgenotypes = genotypes %>%
      select(all_of(c("marker", IDs)))
  } else {
    IDgenotypes = genotypes %>%
      select(-c(position, chr))
  }


  # correct GC to CG to match the souporcell output
  correctedGenotypes = IDgenotypes %>%
    # correct GC values to CG
    pivot_longer(!marker, names_to = "ID", values_to = "genotype") %>%
    mutate(corrGenotype = str_replace(genotype, pattern = "GC",
                                      replacement = "CG")) %>%
    dplyr::select(-genotype) %>%
    pivot_wider(names_from = ID, values_from = corrGenotype) %>%
    # replace unknown genotypes (--) with NA
    mutate(across(where(is.character), ~na_if(., "--")))

  if (!compareIncomplete){
    # keep only rows with genotype information for all clusters
    correctedGenotypes = correctedGenotypes %>%
      drop_na()
  }


  # ==== calculate similarities ===
  # combine known genotypes with observed mutations from vcf file
  testTable = as.data.frame(vcfEdit) %>%
    inner_join(correctedGenotypes)

  # define which columns are souporcell clusters, and which are know genotypes
  # the first column is always "marke" - remove
  clusters = colnames(vcfEdit)[-1]
  knownGen = colnames(correctedGenotypes)[-1]


  # calculate pairwise similarities
  for (i in clusters){
    for (j in knownGen) {
      combo = paste(i, j, sep = "__")
      similarity = simFunc(testTable[,i],
                           testTable[,j],
                           compareIncomplete)
      results = data.frame(combo = combo, similarity = similarity)
      if (!exists("resTable")){
        resTable = results
      } else {
        resTable = rbind(resTable, results)
      }
    }
  }

  # Get correlation matrix
  finalRes = resTable %>%
    tidyr::separate(combo, c("cluster", "ID"), sep = "__") %>%
    pivot_wider(values_from = similarity, names_from = "ID") %>%
    column_to_rownames(var = "cluster")

  return(finalRes)
}

#' plotClustSimilarity: plot results from \code{clustSimilarity} function in
#'    form of heatmap
#' @param simMatrix Output from \code{clustSimilarity} function
#' @param RCBpalette Optional string specifying an RColorBrewer paltette to use
#'    (default: "Purples")
#' @return ComplexHeatmap object
#' @export
#' @examples
#' VCFpath = system.file("extdata", "cluster_genotypes.vcf",
#'                             package = "assignClustID")
#' genotypes = DO_mm10
#' similarities = clustSimilarity(pathToVCF=VCFpath, genotypes=genotypes)
#' plotClustSimilarity(similarities)
plotClustSimilarity = function(simMatrix,
                               RCBpalette="Purples"){
  # --- plot as heatmap
  mat = as.matrix(simMatrix)
  # https://jokergoo.github.io/2022/03/08/support-hcl-colormaps-in-complexheatmap/
  colFun = colorRamp2(seq(min(mat), max(mat), length = 10),
                      hcl_palette = RCBpalette, reverse = TRUE)

  ht = Heatmap(as.matrix(mat),
               col = colFun,
               heatmap_legend_param = list(
                 title = "similarity"),
               #column_title = "ID", column_title_side = "bottom",
               #row_title = "Cluster ID", row_title_side = "right",
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(sprintf("%.3f", mat[i, j]), x, y,
                           gp = gpar(fontsize = 10))
               })

  draw(ht)
  return(ht)
}


#' simFunc: Internal helper function to calculate similarities
#' @param x string vector containing with genotypes from souporcell
#' @param y string vector containing with known genotypes
#' @param compareIncomplete Logical variable specifying if only positions with
#'    complete genotype information from known genotypes and souporcell output
#'    should be used to calculate similarity. TRUE = compare markers for which
#'    genotypes are NOT available for all tested individuals, FALSE =
#'    compare only markers with genotypes available across all individials -
#'    both from souporcell and known genotypes (default - FALSE, which means
#'    only complete cases will be considered)
#' @noRd
simFunc = function(x, y, compareIncomplete){
  # of only complete cases are considered
  if (!compareIncomplete){
    similarity = round(sum(x==y)/length(x), digits = 3)
  } else {
    # find how many comparisons don't have NAs
    nonNA = sum(!is.na(x==y))
    similarity = sum(x==y, na.rm = TRUE)/nonNA
  }
  return(similarity)
}


#' Internal helper function to load VCF file
#'
#' @param pathToVCF A string specifying a path to a VCF file produced by
#'    souporcell. Provide full path including the file name.
#' @return A data.frame - VCF file without the header lines
#' @noRd
loadVCF = function(pathToVCF){
  if(endsWith(pathToVCF,suffix = "vcf") | endsWith(pathToVCF,suffix = "VCF")){
    vcf = as.data.frame(fread(pathToVCF, skip = "#CHROM"))

    vcfEdit = vcf %>%
      # remove backgound based on filter value
      dplyr::filter(FILTER == ".") %>%
      # rename variables to match genotype_SNPS, and attach "chr" in front of chr number
      dplyr::rename("chr"="#CHROM") %>%
      dplyr::rename("position" = "POS")

    # check if the chromosomes start with "chr" if not, add it
    if (sum(startsWith(vcfEdit$chr, prefix = "chr")) == 0){
      vcfEdit = vcfEdit %>%
        mutate(chr = paste0("chr", chr))
    }
  } else {
    stop("The pathToVCF has to end with vcf or VCF. Did you include a full file
         name?")
  }

  return(vcfEdit)

}


