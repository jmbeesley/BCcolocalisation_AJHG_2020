#!/usr/bin/Rscript

# Run hyprcoloc on sets of SNPs within each BCAC region
# Overall BC and ERneg BC risk perfomed separately

library(dplyr)
library(hyprcoloc)

# start dir
wrkdir <- "./"

# BCAC region list
bcac_fm_regions <- read.delim(paste0(wrkdir, "data/BCAC_FM_regions.txt"))

# Colocalisation function for each gene in list
hyprcoloc_gene <- function(gene, x) {

  message(gene)

  # betas
  x1 <- x[, c(1, 3)]
  x1 <- x1 %>% mutate_if(is.null, as.numeric)
  x1 <- x1 %>% mutate_if(is.character, as.numeric)
  rownames(x1) <- rownames(x)
  x2 <- as.matrix(x1)

  # se
  y1 <- x[, c(2, 4)]
  y1 <- y1 %>% mutate_if(is.null, as.numeric)
  y1 <- y1 %>% mutate_if(is.character, as.numeric)
  rownames(y1) <- rownames(x)
  y2 <- as.matrix(y1)

  # other variables
  traits <- colnames(x2)
  variants <- rownames(x2)
  binary.traits <- as.numeric(!grepl("GTEx", traits))

  # run co-localisation analysis
  if (nrow(x2) > 2) {
    res <- hyprcoloc(x2, y2, trait.names = traits, snp.id = variants, binary.outcomes = binary.traits, snpscores = T)
    return(res$results)
  }
}

# Run for overall breast cancer risk data
hyprcoloc_func1 <- function(region) {
  message(paste("Running hyprcoloc for", region))

  # import data
  input_data_dir <- paste0(wrkdir, "output/formatted_for_coloc_analysis/", region)
  indata <- readRDS(paste0(input_data_dir, "/gtex_bcacoverall_association_data.RDS"))
  regiondir <- paste0(wrkdir, "/output/hyprcoloc_results/", region)
  dir.create(regiondir, recursive = T, showWarnings = F)

  regionlist <- list()

  # format input data
  for (i in names(indata)) {
    indata2 <- indata[[i]][, c(1, 4, 5, 7, 8)]
    indata2 <- indata2 %>% filter(beta.BCAC != "NULL")

    # remove duplicated variants
    if (any(duplicated(indata2$variant_id))) {
      dup_variants <- indata2$variant_id[duplicated(indata2$variant_id)]
      x <- indata2 %>% filter(!variant_id %in% dup_variants)
    } else {
      x <- indata2
    }

    rownames(x) <- x[, 1]
    x <- x[, -1]

    # function for each gene
    regionlist[[i]] <- hyprcoloc_gene(i, x)
  }

  z <- bind_rows(regionlist, .id = "GENE")
  z %>% write.table(paste0(regiondir, "/BCAC_overall_hyprcoloc_results.txt"), quote = F, sep = "\t", row.names = F)

}

# Run for ERnegative breast cancer risk data
hyprcoloc_func2 <- function(region) {

  message(paste("Running hyprcoloc for", region))

  # import data
  input_data_dir <- paste0(wrkdir, "output/formatted_for_coloc_analysis/", region)
  indata <- readRDS(paste0(input_data_dir, "/gtex_bcacerneg_association_data.RDS"))
  regiondir <- paste0(wrkdir, "/output/hyprcoloc_results/", region)
  dir.create(regiondir, recursive = T, showWarnings = F)

  regionlist <- list()

  # format input data
  for (i in names(indata)) {
    indata2 <- indata[[i]][, c(1, 4, 5, 7, 8)]
    indata2 <- indata2 %>% filter(beta.BCAC != "NULL")

    # remove duplicated variants
    if (any(duplicated(indata2$variant_id))) {
      dup_variants <- indata2$variant_id[duplicated(indata2$variant_id)]
      x <- indata2 %>% filter(!variant_id %in% dup_variants)
    } else {
      x <- indata2
    }

    rownames(x) <- x[, 1]
    x <- x[, -1]

    # function for each gene
    regionlist[[i]] <- hyprcoloc_gene(i, x)
  }

  z <- bind_rows(regionlist, .id = "GENE")
  z %>% write.table(paste0(regiondir, "/BCAC_erneg_hyprcoloc_results.txt"), quote = F, sep = "\t", row.names = F)
}

# loop through regions
lapply(bcac_fm_regions$FM_region, hyprcoloc_func1)
lapply(bcac_fm_regions$FM_region, hyprcoloc_func2)
