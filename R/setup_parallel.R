library(Rcpp)
sourceCpp("src/bcf_clean_overpar_linear.cpp")
library(bcf)
library(ggplot2)
library(latex2exp)
library(rpart)
library(rpart.plot)
library(partykit)

.get_chain_tree_files = function(tree_path, chain_id, no_output = FALSE){
  if (is.null(tree_path) | no_output){
    out <- list(
      "con_trees" = toString(character(0)), 
      "mod_trees" = toString(character(0))
    )
  } else{
    out <- list("con_trees" = paste0(tree_path,'/',"con_trees.", chain_id, ".txt"), 
                "mod_trees" = paste0(tree_path,'/',"mod_trees.", chain_id, ".txt"))
  }
  return(out)
}