aa_path_network
================
Kimberly Ledger
3/11/2022

this code creates a network plot of pathogen co-infections found in
Amblyomma americanum

libraries

``` r
library(tidyverse)
library(igraph)
library(intergraph)
library(ggnetwork)
```

read in the data - this is an edge list where each row represents a
connection between a microbe in the column (micro1) and another column
(micro2)

``` r
edge_data <- read.csv("AA_pathogen_edgelist.csv")
```

make object with just co-infections

``` r
links <- edge_data[,c(4,5)]
```

create network object

``` r
network <- graph_from_data_frame(d=links, directed = F)
```

plot it

``` r
plot(network)
```

![](aa_coinfection_network_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

plot in a circle

``` r
plot(network, layout = layout_in_circle(network))
```

![](aa_coinfection_network_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

make an adjacency matrix

``` r
adj_mat <- get.adjacency(network)
adj_mat
```

    ## 10 x 10 sparse Matrix of class "dgCMatrix"
    ##                                                   
    ## Borrelia_lonestari          . .  1 . . . .  . .  5
    ## Babesia_sp                  . .  . . . . .  . .  8
    ## Theileria_cervi             1 .  . 1 1 . 2  . . 31
    ## Ehrlichia_ewingii           . .  1 . . . .  . 2  7
    ## Rickettsia_parkeri          . .  1 . . . .  . .  3
    ## Babesia_sp_Coco             . .  . . . . .  . .  7
    ## Ehrlichia_chaffeensis       . .  2 . . . .  . 1  4
    ## Hepatozoon_sp_A             . .  . . . . .  . . 12
    ## Ehrlichia_sp_PanolaMountain . .  . 2 . . 1  . .  4
    ## Rickettsia_amblyommatis     5 8 31 7 3 7 4 12 4  .

try making a prettier plot - weight the connections by coinfections

``` r
plot(network, vertex.label.color="black", 
     vertex.label.cex=.5,vertex.label.dist=0.2, 
     edge.curved=0.5, edge.width = edge.betweenness(network),
     layout=layout_in_circle(network))
```

![](aa_coinfection_network_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
#plot(network, vertex.label.color="black", 
#     vertex.label.cex=.5,vertex.label.dist=0.2, 
#     edge.curved=0.5, edge.width = edge.betweenness(network), edge.color = edge.betweenness(network),
#     layout=layout_in_circle(network))
```
