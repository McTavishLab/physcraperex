% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_branches_method1}
\alias{plot_branches_method1}
\title{Plot a tree with branches colored according to molecular data, method 1}
\usage{
plot_branches_method1(
  x,
  tip_label = "otu",
  drop_outgroup = TRUE,
  ladderize_tree = TRUE,
  color = "red",
  ...
)
}
\arguments{
\item{x}{A list from get_tip_values}

\item{tip_label}{A character vector. Can be one of "otu" or "taxon"}

\item{drop_outgroup}{Boolean}

\item{ladderize_tree}{Boolean}
}
\value{
a plot
}
\description{
Plot a tree with branches colored according to molecular data, method 1
}
\examples{
treefile = 'data/pg_2827_tree6577/run_pg_2827tree6577_run4/RAxML_bestTree.2020-07-31'
otufile = 'data/pg_2827_tree6577/outputs_pg_2827tree6577/otu_info_pg_2827tree6577.csv'
}
\author{
Emily Jane McTavish
}
