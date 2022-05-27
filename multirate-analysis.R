#Relative fit analysis
library(ape)
library(geiger)
library(tidyverse)
library(arbutus)
library(flipR)
library(parallel)
library(OUwie)

#Load data, look at tree
tree <- read.tree("data/crayfish.nodelabels.tre")
nodes <- tree$node.label
plot.phylo(tree, show.node.label = TRUE)
matrix <- read.delim("data/orthogroups.TMM.EXPR.matrix") %>% as.data.frame() %>% rename(Gene = X)

#Make OUwie df
reg <- c(2,1,1,1,2,2,1,1,2,1,1,2,2,2)
tips <- tree$tip.label
OUwie_df <- data.frame(Genus_species = tips, Reg = reg)

#Standard Error function
standard_error <- function(x) sd(x) / sqrt(length(x))

#Average expression values across individuals, generate SE table as well
avg_dat <- matrix %>% group_by(Gene) %>% transmute("CCRYP" = mean(CCRYP_KC8782),
                                                   "CDUBI" = mean(CDUBI_KC8779:CDUBI_KC8781),
                                                   "CGRAY" = mean(CGRAY_KC8769:CGRAY_KC8771),
                                                   "CHAMU" = mean(CHAMU_KC8582:CHAMU_KC8602),
                                                   "CNERT" = mean(CNERT_KC8482:CNERT_KC8484),
                                                   "CRUST" = mean(CRUST_KC8774:CRUST_KC8775),
                                                   "CSETO" = mean(CSETO_KC8659:CSETO_KC8661),
                                                   "CTENE" = mean(CTENE_KC8583:CTENE_KC8587),
                                                   "OAUST" = mean(OAUST_KC8577:OAUST_KC8580),
                                                   "OINCO" = mean(OINCO_KC8597:OINCO_KC8598),
                                                   "PFALL" = mean(PFALL_KC8497:PFALL_KC8499),
                                                   "PHORS" = mean(PHORS_KC8784:PHORS_KC8785),
                                                   "PLUCI" = mean(PLUCI_KC8540:PLUCI_KC8783),
                                                   "PPALL" = mean(PPALL_KC8786))

dat_SE <- matrix %>% group_by(Gene) %>% summarise(Gene, SE = standard_error(across(CCRYP_KC8782:PPALL_KC8786)))

#Now need to flip tables and properly format
format_expr_data <- function (avgdat) {
  temp <- avgdat %>% pull(Gene)
  avgdat <- avgdat %>% ungroup() %>% select(!Gene)
  dat <- flip(avgdat)
  colnames(dat) <- temp
  dat
}

#Use treedata() to trim any extra values
avg_tree_data <- format_expr_data(avg_dat) %>% as.matrix()
input <- list(avg_tree_data, dat_SE)

#Run relative fit
runFC <- function ( tree_dat, SE ){
  fitResults <- vector(mode = "list", length = 3560)
  tdf <- treedata(tree, tree_dat, sort = TRUE)
  phy <- tdf$phy
  phy$node.label <- nodes
  data <- tdf$data
  for(j in 1:3560){
    fitBM <- fitContinuous(phy, data[,j], SE[[2]][[j]], model = "BM")
    fitOU <- fitContinuous(phy, data[,j], SE[[2]][[j]], model = "OU")
    fitEB <- fitContinuous(phy, data[,j], SE[[2]][[j]], model = "EB")
    OUwie_df2 <- OUwie_df %>% mutate(X = (data[,j]))
    fitBMS <- tryCatch(OUwie(phy, OUwie_df2, model = "BMS"), error = function(x)list(AIC = Inf))
    aic <- c(fitBM$opt[["aic"]], fitOU$opt[["aic"]], fitEB$opt[["aic"]], fitBMS$AIC)
    fit <- ifelse(min(aic) == aic[1], list(c(fitBM, model = "BM")), 
                  ifelse(min(aic) == aic[2], list(c(fitOU, model = "OU")), 
                         ifelse(min(aic) == aic[3], list(c(fitEB, model = "EB")),
                                list(c(fitBMS, model = "BMS")))))
    fitResults[j] <- fit
  }
  fitResults
}

model_count <- function (fit) {
  ou = 0
  bm = 0
  eb = 0
  bms = 0
  for(f in fit){
    vec <- f
    ifelse(vec$model == "OU", ou <- ou + 1, ifelse(vec$model == "BM", bm <- bm + 1, ifelse(vec$model == "EB", eb <- eb + 1, bms <- bms + 1)))
  }
  df <- data.frame(OU = ou, BM = bm, EB = eb, BMS = bms)
  b <- df %>% pivot_longer(c(OU, BM, EB, BMS), names_to = "model")
  b
}

class_fix <- function(fitobj){
  if(length(fitobj) == 5) class(fitobj) <- "gfit"
  if(length(fitobj) == 27) class(fitobj) <- "OUwie"
  fitobj
}

run_arb <- function (fits){
  arby <- mclapply(fits, function(i) try(arbutus(i), TRUE), mc.cores = 16)
  arby <- arby[sapply(arby, function(x) !inherits(x, "try-error"))]
  arby_df <- map_df(arby, function(i) try(pvalue_arbutus(i), TRUE))
  arby_df
}

#Run best fit
total_process_best <- function (tree_list){
  fit <- runFC(tree_list[[1]], tree_list[[2]])
  fit_name <- "arbutus/best/best_fit_bms"
  saveRDS(fit, file = fit_name)
  df <- model_count(fit)
  aic_name <- "arbutus/best/AIC_bms.png"
  df %>% ggplot(aes(model, value)) + geom_col() + theme_classic()
  ggsave(aic_name)
  result <- lapply(fit, class_fix) %>% run_arb()
  rds_name <- "arbutus/best/best_pvals_bms"
  saveRDS(result, file = rds_name)
  result %>% select(!m.sig) %>% pivot_longer(cols = everything(), names_to = "tstat") %>%
    ggplot(aes(value)) + geom_histogram(aes(y = ..density..)) + facet_wrap(~tstat, nrow = 1) + theme_bw()
  pval_name <- "arbutus/best/best_arbutus_bms.png"
  ggsave(pval_name)
}

total_process_best(input)
