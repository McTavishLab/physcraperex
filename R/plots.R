# trenex2pdf
#
# Write nexus tree files to pdf.
#
#' @param treefile A character vector with the name of the nexus or newick file
#' @param studyid TreeBASE study id

trenex2pdf <- function(treefile,
                       treefilesdir = "../data-raw/trees/treebase/S1091/",
                       treepdfsdir = "../data-raw/trees/treebase/S1091-pdf/",
                       height = 7, ...){
  # all treebase trees are in nexus format
  tree <- ape::read.nexus(file= paste0(treefilesdir, "/", treefile))

  outputfile <- paste0(treepdfsdir, gsub(".nex", ".pdf", treefile))
  pdf(file = outputfile, height)
  ape::plot.phylo(ape::ladderize(tree), cex = 0.5, ...)
  graphics::mtext(paste("TreeBASE tree:", gsub(".nex", "", treefile)))
  dev.off()
}

# get_file_names
#
# Write nexus or newick files to pdf or png.
#
#' @param studyid TreeBASE study id
#' @param file_dir Character vector of a directory with one or more tree files

get_file_names <- function(studyid, file_dir = "data-raw/trees/treebase"){
  # dir_original <- getwd()
  print(getwd())
  treefiles <- paste0(studyid, "-treefiles.csv")
  system(paste0("ls ", studyid, " > ", treefiles))
  return(paste0(file_dir, "/", treefiles))
}

# Cummulative sum of tip values for each node/branch down to the root
#
#
#' @param tree a phylo object
#' @param values named numeric character vector with names corresponding to tip labels
#' @author Luna L. Sanchez Reyes
#' @return a vector
#'

sum_tips <- function(tree, values) {
    # tree <- ape::reorder.phylo(tree, "postorder")  # no need to reorder the whole tree
    res <- numeric(max(tree$edge))
    res[1:ape::Ntip(tree)] <- values[match(names(values), tree$tip.label)]
    for (i in ape::postorder(tree))  { # ape postorder doesn't include root
         tmp <- tree$edge[i,1]
         # print(i)
         res[tmp] <- res[tmp] + res[tree$edge[i, 2]]
         # print(res)
   }
   res
}


# Plot a tree with branches colored according to molecular data.
#
#' @param tree a phylo object
#' @param values named numeric character vector with names corresponding to tip labels
#' @author Emily Jane McTavish
#' @return a vector
#' @example
#' treefile = 'data/pg_2827_tree6577/run_pg_2827tree6577_run4/RAxML_bestTree.2020-07-31'
#' otu_info_file = 'data/pg_2827_tree6577/outputs_pg_2827tree6577/otu_info_pg_2827tree6577.csv'

plot_branches <- function(treefile, otu_info_file){

  phytree = ape::read.tree(treefile)

  tip_info = utils::read.csv(otu_info_file, sep = '\t', row.names = 1, stringsAsFactors = FALSE)

  inout <- as.vector(tip_info$X.physcraper.ingroup)
  names(inout) <- rownames(tip_info)

  outgroups = c()
  for (lab in phytree[["tip.label"]]) {
    status = inout[names(inout) == lab]
    if (status == "False")
    {outgroups = c(outgroups, lab)}
  }

  pruned <- ape::drop.tip(phytree, phytree$tip.label[match(outgroups, phytree$tip.label)])

  st <- as.vector(tip_info$X.physcraper.status)
  names(st) <- rownames(tip_info)

  spp <- as.vector(tip_info$X.ot.ottTaxonName)
  names(spp) <- rownames(tip_info)


  new = c()
  spp_labels = c()
  for (lab in pruned[["tip.label"]]) {
    status = st[names(st) == lab]
    val = grepl("new", status, fixed = TRUE)
    val = 100*as.integer(val)
    names(val) <- lab
    new = c(new, val)
    tax = spp[names(spp) == lab]
    spp_labels = c(spp_labels, tax)
  }


  phytools::plotBranchbyTrait(pruned, new, mode=c("tips"),legend=FALSE, show.tip.label = FALSE, palette=colorRampPalette(c("black", "red")), type='fan')

  pruned[["tip.label"]] <- spp_labels
  names(new) <- pruned[["tip.label"]]
  phytools::plotBranchbyTrait(pruned, new, mode=c("tips"), legend=FALSE, show.tip.label = FALSE, palette=colorRampPalette(c("black", "red")), cex=0.5, edge.width=3, type='unrooted')

}


## set working directory to path to https://github.com/McTavishLab/physcraper_example

#treefile = 'data/pg_2827_tree6577/run_pg_2827tree6577_run4/RAxML_bestTree.2020-07-31'
#otu_info_file = 'data/pg_2827_tree6577/outputs_pg_2827tree6577/otu_info_pg_2827tree6577.csv'

#treefile = 'pg_55_local_new/outputs_pg_55tree5864_ndhf/physcraper_pg_55tree5864_ndhf.tre'
#otu_info_file = 'pg_55_local_new/outputs_pg_55tree5864_ndhf/otu_info_pg_55tree5864_ndhf.csv'


#treefile = 'ot_350/outputs_ot_350Tr53297/physcraper_ot_350Tr53297.tre'
#otu_info_file = 'ot_350/outputs_ot_350Tr53297/otu_info_ot_350Tr53297.csv'
