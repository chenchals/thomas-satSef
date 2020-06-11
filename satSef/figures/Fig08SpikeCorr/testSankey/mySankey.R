# From https://www.r-graph-gallery.com/323-sankey-diagram-with-the-networkd3-library.html
# Load package
library(networkD3)

# Load energy projection data
# URL <- "https://cdn.rawgit.com/christophergandrud/networkD3/master/JSONdata/energy.json"
# Energy <- jsonlite::fromJSON(URL)


# Now we have 2 data frames: a 'links' data frame with 3 columns (from, to, value), and a 'nodes' data frame that gives the name of each node.
# head( Energy$links )
# head( Energy$nodes )

# Thus we can plot it
# p <- sankeyNetwork(Links = Energy$links, Nodes = Energy$nodes, Source = "source",
#                   Target = "target", Value = "value", NodeID = "name",
#                   units = "TWh", fontSize = 12, nodeWidth = 30)
# p


conns <- jsonlite::fromJSON("sankErrorTimingAccuErrorOtherFast.json")

conns$links$value <- conns$links$counts/sum(conns$links$counts)
conns$links$linkColor <- c("Minus","Minus","Plus","Plus","Minus","Minus","Plus","Plus")

conns$links$absMean <- abs(conns$links$rhoValMean)


p1 <- sankeyNetwork(Links = conns$links, Nodes = conns$nodes, Source = "source",
                   Target = "target", Value = "counts", NodeID = "name",
                   fontSize = 12, nodeWidth = 10,LinkGroup = "linkColor")
p1


p2


# save the widget
# library(htmlwidgets)
# saveWidget(p1, file=paste0( getwd(), "sankeyConn.html"))
