## Just OU Analysis

#Relative fit analysis
library(ape)
library(geiger)
library(tidyverse)
library(arbutus)
library(flipR)
library(parallel)

#Load data, look at tree
tree <- read.tree("data/crayfish.nodelabels.tre")
plot.phylo(tree, show.node.label = TRUE)
matrix <- read.delim("data/orthogroups.TMM.EXPR.matrix") %>% as.data.frame() %>% rename(Gene = X)

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

run_arb <- function (fits){
  arby <- vector("list", length = length(fits))
  count = 1
  for(f in fits){
    class(f) <- "gfit"
    arby[[count]] <- arbutus(f)
    count = count + 1
  }
  arby_df <- map_df(arby, pvalue_arbutus)
  arby_df
}

runOU <- function ( tree_dat, SE ){
  fitResults <- vector(mode = "list", length = 3560)
  tdf <- treedata(tree, tree_dat, sort = TRUE)
  phy <- tdf$phy
  data <- tdf$data
  for(j in 1:3560){
    fitOU <- fitContinuous(phy, data[,j], SE[[2]][[j]], model = "OU")
    fitResults[[j]] <- fitOU
  }
  fitResults
}

total_process_OU <- function (tree_list){
  fit <- runOU(tree_list[[1]], tree_list[[2]])
  fit_name <- "arbutus/OU/OU_fit"
  saveRDS(fit, file = fit_name)
  result <- run_arb(fit)
  rds_name <- "arbutus/OU/OU_pvals"
  saveRDS(result, file = rds_name)
  result %>% pivot_longer(cols = everything(), names_to = "tstat") %>%
    ggplot(aes(value)) + geom_histogram(aes(y = ..density..)) + facet_wrap(~tstat, nrow = 1) + theme_bw()
  pval_name <- "arbutus/OU/OU_arbutus.png"
  ggsave(pval_name)
}

total_process_OU(input)