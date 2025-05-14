## example of evolutionary accumulation modelling (EvAM) for multidrug resistance (MDR)
# uses HyperTraPS-CT to explore tuberculosis MDR evolution

##########
##### Import libraries

library(readxl)       # to read Excel datafile
library(dplyr)        # for data wrangling
library(phytools)     # to read external tree data
library(ggplot2)      # }
library(ggpubr)       # } for visualisation
library(ggraph)       # }

# download and load HyperTraPS-CT package
remotes::install_github("StochasticBiology/hypertraps-ct@bioconductor")
library(hypertrapsct)

##########
##### Import data

# get profiles of drug resistance/susceptibility
# Supplementary Table 4 of Casali et al., https://www.nature.com/articles/ng.2878.s3
system("wget https://static-content.springer.com/esm/art%3A10.1038%2Fng.2878/MediaObjects/41588_2014_BFng2878_MOESM35_ESM.xls")
o.df = read_excel("41588_2014_BFng2878_MOESM35_ESM.xls")

# get phylogeny linking isolates
# Supplementary Data Set 1 of Casali et al., https://www.nature.com/articles/ng.2878.s3
system("wget https://static-content.springer.com/esm/art%3A10.1038%2Fng.2878/MediaObjects/41588_2014_BFng2878_MOESM34_ESM.txt")
tree = read.tree("41588_2014_BFng2878_MOESM34_ESM.txt")

##########
##### Curate data

# extract isolate ID and resistance profiles to our ten drugs (discarding mutation info)
# remove any incomplete profiles and recast as binary strings
col.interest = c("Isolate", "INH", "RIF", "PZA", "EMB", "STR", "AMI", "CAP", "MOX", "OFL", "PRO")
df = o.df[,which(colnames(o.df) %in% col.interest)]
missing.rows = unique(which(df == ".", arr.ind=TRUE)[,1])
final.df = df[-missing.rows,]
final.df = final.df %>%
  mutate(across(-Isolate, ~ ifelse(. == "R", 1, 0)))
final.df = as.data.frame(final.df)

##########
##### Evolutionary accumulation modelling

# reconstruct ancestral states using parsimony picture and extract before-after transitions
src.data = curate.tree(tree, final.df)

# fit model using HyperTraPS-CT (will take an hour or so with this parameterisation)
# scaling factor of 1000 just adjusts "time" units for easier inference and interpretation
fitted.model = HyperTraPS(src.data$dests, 
                          initialstates = src.data$srcs, 
                          starttimes = src.data$times*1000, endtimes = src.data$times*1000, 
                          penalty = 1, seed = 1, length = 5, kernel = 3)

# example visualisations
fitted.model$featurenames = col.interest[2:11]
ggarrange(
  ggarrange(plotHypercube.curated.tree(src.data),
            plotHypercube.influencegraph(fitted.model, cv.thresh = 0.3) + 
              scale_edge_width_continuous(range = c(1,2)),
            ncol=1),
  plotHypercube.sampledgraph2(fitted.model, edge.label.size=3, edge.label.angle = "none", node.labels=FALSE,
                            no.times=TRUE, small.times=TRUE) + theme(legend.position = "none"),
  ncol=2, widths=c(1,2))

# output results to image file for publication
sf = 2
png("results-grid.png", width=600*sf, height=600*sf, res=72*sf)
ggarrange(
  ggarrange(plotHypercube.curated.tree(src.data),
            plotHypercube.sampledgraph2(fitted.model, edge.label.size=3, edge.label.angle = "none", node.labels=FALSE,
                                        no.times=TRUE, small.times=TRUE) + theme(legend.position = "none") +   theme(plot.margin = unit(c(0, 0., 0, 0.), "cm")),
   nrow=1, widths=c(1,2)),
  ggarrange(plotHypercube.influencegraph(fitted.model, cv.thresh = 0.3) + 
              scale_edge_width_continuous(range = c(1,2)),
  ggplot() + geom_blank() + theme_void(),
  nrow=1, widths=c(1,2)),
  
  nrow=2, heights=c(2.5,1.3))
dev.off()

# some examples of predictions given trained model
predictNextStep(fitted.model, c(0,0,0,0,0,0,0,0,0,0))
predictNextStep(fitted.model, c(1,1,0,0,1,0,0,0,0,0))
plotHypercube.prediction(predictNextStep(fitted.model, c(1,1,0,0,1,0,0,0,0,0)), max.size = 10)
