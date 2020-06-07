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
