#Relative fit analysis
library(ape)
library(geiger)
library(tidyverse)
library(arbutus)
library(flipR)
library(parallel)

#Load data, look at tree
tree <- read.tree("data/crayfish.nodelabels.tre") %>% drop.tip(c("CSETO", "CCRYP", "PFALL", "OAUST"))
plot(tree)

#Standard Error function
standard_error <- function(x) sd(x) / sqrt(length(x))

log_normalize <- function(x){
  res <- x
  if (!is.na(x)) {
    res <- log(x)
  } 
  res
}

tpm_matrix <- read.delim("TPM_data.tsv") %>% as.data.frame() %>% pivot_longer(!c(Orthogroup)) %>% mutate(value = log(value)) %>% pivot_wider(names_from = "name", values_from = "value")
tpm_matrix[tpm_matrix == -Inf] <- 0
tpm_mean <- tpm_matrix %>%
  group_by(Orthogroup) %>%
  transmute(
    PPALL = PPALL_A,
    CRUST = tryCatch(mean(CRUST_A:CRUST_B, na.rm = TRUE), error = function(x)696969),
    OINCO = tryCatch(mean(OINCO_A:OINCO_B, na.rm = TRUE), error = function(x)696969),
    PHORS = tryCatch(mean(PHORS_A:PHORS_B, na.rm = TRUE), error = function(x)696969),
    PLUCI = tryCatch(mean(PLUCI_A:PLUCI_B, na.rm = TRUE), error = function(x)696969),
    CDUBI = tryCatch(mean(CDUBI_A:CDUBI_C, na.rm = TRUE), error = function(x)696969),
    CGRAY = tryCatch(mean(CGRAY_A:CGRAY_C, na.rm = TRUE), error = function(x)696969),
    CHAMU = tryCatch(mean(CHAMU_A:CHAMU_C, na.rm = TRUE), error = function(x)696969),
    CNERT = tryCatch(mean(CNERT_A:CNERT_C, na.rm = TRUE), error = function(x)696969),
    CTENE = tryCatch(mean(CTENE_A:CTENE_C, na.rm = TRUE), error = function(x)696969)
  ) %>% ungroup() %>% filter(!Orthogroup %in% c("OG0009623", "OG0009630","OG0009637","OG0009651","OG0009672","OG0009674","OG0009678","OG0009725","OG0009744","OG0009770","OG0009791","OG0009940", "OG0009867"))
tpm_mean[tpm_mean == 696969] <- NA

rpkm_matrix <- read.delim("FPKM_data.tsv") %>% as.data.frame() %>% pivot_longer(!c(Orthogroup)) %>% mutate(value = log(value)) %>% pivot_wider(names_from = "name", values_from = "value")
rpkm_matrix[rpkm_matrix == -Inf] <- 0

rpkm_mean <- rpkm_matrix %>%
  group_by(Orthogroup) %>%
  transmute(
    PPALL = PPALL_A,
    CRUST = tryCatch(mean(CRUST_A:CRUST_B, na.rm = TRUE), error = function(x)696969),
    OINCO = tryCatch(mean(OINCO_A:OINCO_B, na.rm = TRUE), error = function(x)696969),
    PHORS = tryCatch(mean(PHORS_A:PHORS_B, na.rm = TRUE), error = function(x)696969),
    PLUCI = tryCatch(mean(PLUCI_A:PLUCI_B, na.rm = TRUE), error = function(x)696969),
    CDUBI = tryCatch(mean(CDUBI_A:CDUBI_C, na.rm = TRUE), error = function(x)696969),
    CGRAY = tryCatch(mean(CGRAY_A:CGRAY_C, na.rm = TRUE), error = function(x)696969),
    CHAMU = tryCatch(mean(CHAMU_A:CHAMU_C, na.rm = TRUE), error = function(x)696969),
    CNERT = tryCatch(mean(CNERT_A:CNERT_C, na.rm = TRUE), error = function(x)696969),
    CTENE = tryCatch(mean(CTENE_A:CTENE_C, na.rm = TRUE), error = function(x)696969)
  ) %>% ungroup() %>% filter(!Orthogroup %in% c("OG0009623", "OG0009630","OG0009637","OG0009651","OG0009672","OG0009674","OG0009678","OG0009725","OG0009744","OG0009770","OG0009791","OG0009940", "OG0009867"))

rpkm_mean[rpkm_mean == 696969] <- NA

#Average expression values across individuals, generate SE table as well
tpm_avg_dat <- tpm_mean %>% group_by(Orthogroup)
rpkm_avg_dat <- rpkm_mean %>% group_by(Orthogroup)

tpm_SE <- tpm_matrix %>% group_by(Orthogroup) %>% summarise(Orthogroup, SE = standard_error(across(CDUBI_A:PPALL_A))) %>% filter(!Orthogroup %in% c("OG0009623", "OG0009630","OG0009637","OG0009651","OG0009672","OG0009674","OG0009678","OG0009725","OG0009744","OG0009770","OG0009791","OG0009940", "OG0009867"))
tpm_SE[is.na(tpm_SE)] <- 0

rpkm_SE <- rpkm_matrix %>% group_by(Orthogroup) %>% summarise(Orthogroup, SE = standard_error(across(CDUBI_A:PPALL_A))) %>% filter(!Orthogroup %in% c("OG0009623", "OG0009630","OG0009637","OG0009651","OG0009672","OG0009674","OG0009678","OG0009725","OG0009744","OG0009770","OG0009791","OG0009940", "OG0009867"))
rpkm_SE[is.na(rpkm_SE)] <- 0

#Now need to flip tables and properly format
format_expr_data <- function (avgdat) {
  temp <- avgdat %>% pull(Orthogroup)
  avgdat <- avgdat %>% ungroup() %>% select(!Orthogroup)
  dat <- flip(avgdat)
  colnames(dat) <- temp
  dat
}

#Use treedata() to trim any extra values
tpm_avg_tree_data <- format_expr_data(tpm_avg_dat) %>% as.matrix()
tpm_input <- list(tpm_avg_tree_data, tpm_SE)

rpkm_avg_tree_data <- format_expr_data(rpkm_avg_dat) %>% as.matrix()
rpkm_input <- list(rpkm_avg_tree_data, rpkm_SE)

#Run relative fit
runFC <- function ( tree_dat, SE ){
  fitResults <- vector(mode = "list", length = 179)
  for(j in 1:179){
    temp <- tree_dat[,j] 
    temp <- temp[!is.na(temp)]
    tdf <- treedata(tree, temp, sort = TRUE)
    phy <- tdf$phy
    data <- tdf$data
    fitBM <- fitContinuous(phy, data, SE[[2]][[j]], model = "BM")
    fitOU <- fitContinuous(phy, data, SE[[2]][[j]], model = "OU")
    fitEB <- fitContinuous(phy, data, SE[[2]][[j]], model = "EB")
    aic <- c(fitBM$opt[["aic"]], fitOU$opt[["aic"]], fitEB$opt[["aic"]])
    fit <- ifelse(min(aic) == aic[1], list(c(fitBM, model = "BM")), 
                  ifelse(min(aic) == aic[2], list(c(fitOU, model = "OU")), 
                         list(c(fitEB, model = "EB"))))
    fitResults[j] <- fit
    print(j)
  }
  fitResults
}


model_count <- function (fit) {
  ou = 0
  bm = 0
  eb = 0
  for(f in fit){
    vec <- f
    ifelse(vec$model == "OU", ou <- ou + 1, ifelse(vec$model == "BM", bm <- bm + 1, eb <- eb + 1))
  }
  df <- data.frame(OU = ou, BM = bm, EB = eb)
  b <- df %>% pivot_longer(c(OU, BM, EB), names_to = "model")
  b
}

class_fix <- function(fitobj){
  class(fitobj) <- "gfit"
  fitobj
}

run_arb <- function (fits){
  arby <- mclapply(fits, function(i) try(arbutus(i), TRUE), mc.cores = 16)
  arby <- arby[sapply(arby, function(x) !inherits(x, "try-error"))]
  arby_df <- map_df(arby, function(i) try(pvalue_arbutus(i), TRUE))
  arby_df
}

#Run best fit
total_process_best <- function (tree_list, nor){
  fit <- runFC(tree_list[[1]], tree_list[[2]])
  fit_name <- paste0("arbutus/best/best_fit_", nor)
  saveRDS(fit, file = fit_name)
  df <- model_count(fit)
  aic_name <- paste0("arbutus/best/AIC_", nor, ".png")
  df %>% ggplot(aes(model, value)) + geom_col() + theme_classic()
  ggsave(aic_name)
  result <- lapply(fit, class_fix) %>% run_arb()
  rds_name <- paste0("arbutus/best/best_pvals_", nor)
  saveRDS(result, file = rds_name)
  result %>% select(!m.sig) %>% pivot_longer(cols = everything(), names_to = "tstat") %>%
    ggplot(aes(value)) + geom_histogram(aes(y = ..density..)) + facet_wrap(~tstat, nrow = 1) + theme_bw()
  pval_name <- paste0("arbutus/best/best_arbutus_", nor, ".png")
  ggsave(pval_name)
}

total_process_best(tpm_input, "tpm_log")
total_process_best(rpkm_input, "rpkm_log")
