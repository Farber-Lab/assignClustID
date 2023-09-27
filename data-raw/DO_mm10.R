rm(list = ls())
# pick few IDs from the DO genotypes as examples

# make mm10 genotype file
DO_gen = as.data.frame(fread("/Users/kristynakupova/Farber_lab/annotations/DO_mice/DO_genotypes_alleles.csv", header = TRUE))
colnames(DO_gen) = paste0("X", colnames(DO_gen))
DO_gen = DO_gen %>%
  dplyr::rename(marker = XV1)

# load locus descriptions (object named snps)
# and select marker names and positions
load("/Users/kristynakupova/Farber_lab/annotations/DO_mice/snps.gigamuga_mm10.Rdata")
snps_key = snps %>%
  dplyr::select(marker, chr, pos) %>%
  dplyr::rename("position" = "pos")
# merge genotypes and snps
DO_mm10 = snps_key %>%
  inner_join(DO_gen) %>%
  dplyr::select(marker, chr, position, X1, X110, X111, X120, X121, X150)

# document data
usethis::use_data(DO_mm10, overwrite = TRUE)
