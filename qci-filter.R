#!/usr/bin/env Rscript
# load pkgs----
pkg <- lapply(c(
  "readxl",
  "tidyverse",
  "argparser",
  "cli",
  "openxlsx",
  "pbapply"
), function(x) {
  suppressWarnings(suppressMessages(library(x, character.only = TRUE)))
})
#suppress dplyr messages
options(dplyr.summarise.inform = FALSE)
# ----
# options----
get_options <- function() {
  p <- arg_parser("Filter variants according to the heritance")
  p <- add_argument(
    parser = p,
    arg = c("--files"),
    nargs = Inf,
    help = "Files of annotated variant tables.")
  p <- add_argument(
    parser = p,
    arg = c("--output"),
    nargs = Inf,
    help = "Name of output files.")
  # p <- add_argument(
  #   parser = p,
  #   arg = c("--qci"),
  #   help = "If table(s) exported from QCI.",
  #   flag = TRUE,)
  p <- add_argument(
    parser = p,
    arg = c("--threshold"),
    type = "numeric",
    help = "Threshold for lab frequency filter.")
  p <- add_argument(
    parser = p,
    arg = c("--mode"),
    default = "mono",
    help = "How to filter table. Takes duo, trio or quadro")
  parse_args(p)
}
#----
# read files----
read_multi <- function(file) {
  if (tools::file_ext(file) %in% c("xlsx", "xls")) {
    cli_alert_info("Reading excel for sample: {file}")
    readxl::read_xlsx(file)
  } else if (tools::file_ext(file) %in% c("tsv", "gz")) {
    cli_alert_info("Reading tsv for sample: {file}")
    suppressWarnings(
      readr::read_delim(file, "\t",
        show_col_types = FALSE,
        col_types = cols(.default = "c"))
    )
  } else {
    cli_alert_warning(
      "File extension(s) \"{tools::file_ext(file)}\" can not be handled!")
    cli_abort("Aborting")
  }
}
#----
# get density based elbow----
get_delbow_points <- function(x, iter, probs = 0.05) {
  counts <- hist(x, breaks = iter, plot = FALSE)$counts
  breaks <- hist(x, breaks = iter, plot = FALSE)$breaks
  names(counts) <- breaks[-length(breaks)]
  peak_indx <- c(F,diff(sign(c(diff(counts))))==-2,F) %>% which()
  topcounts <- counts[peak_indx]
  return(quantile(as.numeric(names(topcounts)), probs = probs))
}
#----
# threshold----
calc_threshold <- function(table, threshold = NA, iter = nrow(table)/10) {
  if (is.na(threshold)) {
    cli_alert_info("Calculating threshold.")
    thres <- get_delbow_points(as.double(table$Lab_frequency), iter)
    # thres <- pbsapply(b, function(x) {quantile(as.numeric(x$Lab_frequency), na.rm = T, probs = 0.035)})
  } else if (is.double(threshold)) {
    thres <- threshold
  } else {
    thres <- 10
    cli_alert_warning("Threshold could not calculated. Default is 10")
  }
  cat_bullet(thres)
  return(thres)
}
#----
# default filter----
filter_default <- function(table) {
  # this function apply default filters.
  # selects only relevant columns.
  # selects only low frequency (lower than average)
  # selects full OMIM phenotype
  # base_https <- "https://varsome.com/variant/hg38/"
  # anno_mode <- "?annotation-mode=germline"
  if (sum(is.na(table$Lab_frequency)) > 0) {
    table <- table %>%
      mutate(Lab_frequency = replace_na(Lab_frequency, 0))
  }
  table <- table %>%
    filter(!is.na(`OMIM Phenotype`)) %>%
    arrange(Lab_frequency)
  # mutate(query = paste0(base_https, paste(Chromosome, Position, `Reference Allele`, `Sample Allele`, sep = "-"), anno_mode))
  # empty positions are '-', not empty. this is a temporary workaround.
  # table$query <- sub("---", "--", table$query)
  return(table)
}
#----
# strict filter----
filter_strict <- function(table, threshold, classification_terms = c("Pathogenic")) {
  # terms only pathogenic
  t <- table %>%
    filter_default() %>%
    filter(Lab_frequency < threshold) %>%
    filter((Classification %in% classification_terms))
  cli::cat_bullet("Applying strict filter.")
  cli::cat_bullet("# of final strict list is ", nrow(t))
  return(t)
}
#----
# loose filter----
filter_loose <- function(table, threshold, classification_terms = c("Uncertain Significance", "Likely Pathogenic")) {
  # terms other than benign
  genotype_col <- paste0(table$`Case Samples`[1], " - Genotype")
  t <- table %>% 
    filter_default() %>%
    filter(Lab_frequency < threshold) %>%
    filter((Classification %in% classification_terms)) %>%
    filter(str_detect(`OMIM Phenotype`, "recessive")) %>%
    filter(.data[[genotype_col[[1]]]] %in% c("Hom")) %>%
    arrange(Lab_frequency)
  cli::cat_bullet("Applying loose filter.")
  cli::cat_bullet("# of final loose list is ", nrow(t))
  return(t)
}
#----
# compound filter----
filter_compound <- function(table, threshold, classification_terms = c("Uncertain Significance", "Likely Pathogenic")) {
  # terms other than benign
  genotype_col <- paste0(table$`Case Samples`[1], " - Genotype")
  t <- table %>% 
    filter_default() %>%
    filter(Lab_frequency < threshold) %>%
    filter((Classification %in% classification_terms)) %>%
    filter(str_detect(`OMIM Phenotype`, "recessive")) %>%
    filter(.data[[genotype_col[[1]]]] %in% c("Het")) %>% 
    filter(duplicated(`Gene Symbol`) | duplicated(`Gene Symbol`, fromLast = TRUE)) %>%
    arrange(`Gene Symbol`)
  cli::cat_bullet("Applying compound filter.")
  cli::cat_bullet("# of final compound list is ", nrow(t))
  return(t)
}
#----
# de novo filter----
filter_denovo <- function(table, threshold, classification_terms = c("Uncertain Significance", "Likely Pathogenic")) {
  # terms other than benign
  genotype_col <- paste0(table$`Case Samples`[1], " - Genotype")
  t <- table %>% 
    filter_default() %>% 
    filter(Lab_frequency < threshold) %>%
    filter((Classification %in% classification_terms)) %>%
    filter(.data[[genotype_col[[1]]]] %in% c("Het")) %>%
    filter(str_detect(`OMIM Phenotype`, "dominant")) %>%
    arrange(Lab_frequency)
  cli::cat_bullet("Applying de novo filter.")
  cli::cat_bullet("# of final de novo list is ", nrow(t))
  return(t)
}
#----
# x-linked filter----
filter_xlinked <- function(table, threshold, classification_terms = c("Uncertain Significance", "Likely Pathogenic")) {
  # terms other than benign
  t <- table %>% 
    filter_default() %>%
    filter(Lab_frequency < threshold) %>%
    filter((Classification %in% classification_terms)) %>% 
    filter(Chromosome == "X") %>% 
    filter(str_detect(`OMIM Phenotype`, "X-linked")) %>%
    arrange(Lab_frequency)
  cli::cat_bullet("Applying X linked filter.")
  cli::cat_bullet("# of final X linked list is ", nrow(t))
  return(t)
}
#----
# export----
export_table <- function(
    default = default,
    strict = strict, 
    loose = loose, 
    compound = compound, 
    denovo = denovo, 
    x_linked = x_linked, 
    names) {
  for(i in 1:length(names)){
    list_cat <- list(
      default = default[[i]],
      strict = strict[[i]], 
      loose = loose[[i]], 
      compound = compound[[i]], 
      denovo = denovo[[i]], 
      x_linked = x_linked[[i]])
    suffix <- '_filtered.xlsx'
    filename <- paste0(names[i], suffix)
    cli_alert_info("Exporting: {filename}")
    write.xlsx(list_cat, file = filename)
  }
}
#----
# parse options----
argv <- get_options()
files <- argv$files
output <- argv$output
if (!is.na(output) & length(files) != length(output)) {
  stop("Number of files and outputs should be same. Aborting.")
}
threshold <- argv$threshold
# qci <- argv$qci
mode <- argv$mode
#----

input_list <- pblapply(files, read_multi)
names(input_list) <- str_extract(files, "MG[0-9]+(.[0-9]+)+")
if (any(is.na(names(input_list)))) {
  cli_abort("There are NA values in sample names. Aborting.")
} else {
  cli_alert_info("Sample name(s) for provided table(s): {names(input_list)}")
}
# apply filters to all tables
default <- pblapply(input_list, function(x) {
  filter_default(x)
  })
strict <- pblapply(input_list, function(x) {
  thres <- calc_threshold(x)
  filter_strict(x, thres)
  })
loose <- pblapply(input_list, function(x) {
  thres <- calc_threshold(x)
  filter_loose(x, thres)
  })
compound <- pblapply(input_list, function(x) {
  thres <- calc_threshold(x)
  filter_compound(x, thres)
  })
denovo <- pblapply(input_list, function(x) {
  thres <- calc_threshold(x)
  filter_denovo(x, thres)
  })
x_linked <- pblapply(input_list, function(x) {
  thres <- calc_threshold(x)
  filter_xlinked(x, thres)
  })

export_table(
 default = default,
 strict = strict,
 loose = loose,
 compound = compound,
 denovo = denovo,
 x_linked = x_linked,
 names = names(input_list))
