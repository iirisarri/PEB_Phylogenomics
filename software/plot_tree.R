library("ggtree")
library("ggimage")
library("phytools")

# set working directory
setwd("~/path-to-treefiles/")

# read tree (newick format)
tree <- read.tree("partitioned.treefile.tre")

# reroot on the branch leading to "Callorhinchus_milii"
tree <- reroot(tree, interactive = T)

# get phylopic images for the tips
d <- ggimage::phylopic_uid(tree$tip.label)

# plot rectangular tree with phylopics
ggtree(tree) %<+% d + 
  geom_tiplab(aes(image=uid), geom="phylopic", offset=0.3, cex = 0.07) +
  geom_tiplab(aes(label=label), offset = 0.02) + xlim(NA, 1) +
  geom_nodelab(aes(label=label), cex=3, nudge_x = -0.04, nudge_y = 0.4)

# circular tree with phylopics as labels
ggtree(tree,layout = "circular") %<+% d + 
  geom_tiplab(aes(image=uid), geom="phylopic", offset=0.1) +
  geom_tiplab(aes(label=NA), offset = .2) + xlim(NA, 0.4)
