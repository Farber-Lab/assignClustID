library(spgs)

# create table for comparisons

rm(list = ls())

# create all possible combinations
combinations1 = c(rep("A", 4), rep("C", 4), rep("G",4), rep("T",4))
combinations2 = rep(c("A", "C", "G", "T"),4)

# 1) combinations
combinations = paste0(combinations1, combinations2)

# 2) fliped order
combinationsFlip = paste0(combinations2, combinations1)

# 3) complement
combinationsComp = complement(combinations, case = "upper")


# 4) reverse complement
combinationsRevComp = reverseComplement(combinations, case = "upper")


# initiate matrix with all possible pairs
comparisonTable = matrix(0, length(combinations), length(combinations))
colnames(comparisonTable) = combinations
rownames(comparisonTable) = combinations

# now set 1 to all possible matches
for (i in seq(1, length(combinations))){

  # direct matches
  comparisonTable[combinations[i], combinations[i]] = 1

  # flipped matches
  comparisonTable[combinations[i], combinationsFlip[i]] = 1

  # complement matches
  comparisonTable[combinations[i], combinationsComp[i]] = 1

  # reverse complement matches
  comparisonTable[combinations[i], combinationsRevComp[i]] = 1

}

# document data
usethis::use_data(comparisonTable, overwrite = TRUE)

