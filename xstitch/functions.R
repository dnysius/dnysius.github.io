library(imager)
library(tidyverse)
library(tidymodels)
library(sp)
library(cowplot)
library(dmc)

change_resolution <- function(image_df, x_size)
{
  ## change_resolution(image_df, x_size) subsamples an image to produce
  ## a lower resolution image. Any non-coordinate columns in the data
  ## frame are summarized with their most common value in the larger
  ## grid cell.
  ##
  ## Input:
  ## - image_df: A data frame in wide format. The x-coordinate column MUST
  ##             be named 'x' and the y-coordinate column MUST be named 'y'.
  ##             Further columns have no naming restrictions.
  ## - x_size:   The number of cells in the x-direction. The number of cells
  ##             in the vertical direction will be computed to maintain the
  ##             perspective. There is no guarantee that the exact number
  ##             of cells in the x-direction is x_size
  ##
  ## Output:
  ## - A data frame with the same column names as image_df, but with fewer
  ##   entries that corresponds to the reduced resolution image.
  ##
  ## Example:
  ##   library(imager)
  ##   library(dplyr)
  ##   fpath <- system.file('extdata/Leonardo_Birds.jpg',package='imager')
  ##   im <- load.image(fpath)
  ##   im_dat<- as.data.frame(im,wide = "c") %>% rename(R = c.1, G = c.2, B = c.3) %>%
  ##            select(x,y,R,G,B)
  ##   agg_image <- change_resolution(im_dat, 50)

  if(!require(sp)) {
    stop("The sp packages must be installed. Run install.packages(\"sp\") and then try again.")
  }
  if(!require(dplyr)) {
    stop("The dplyr packages must be installed. Run install.packages(\"dplyr\") and then try again.")
  }

  sp_dat <- image_df
  gridded(sp_dat) = ~x+y

  persp = (gridparameters(sp_dat)$cells.dim[2]/gridparameters(sp_dat)$cells.dim[1])
  y_size = floor(x_size*persp)
  orig_x_size = gridparameters(sp_dat)$cells.dim[1]
  orig_y_size = gridparameters(sp_dat)$cells.dim[2]

  x_res = ceiling(orig_x_size/x_size)
  y_res = ceiling(orig_y_size/y_size)

  gt = GridTopology(c(0.5,0.5), c(x_res, y_res),
                    c(floor(orig_x_size/x_res), floor(orig_y_size/y_res)))
  SG = SpatialGrid(gt)
  agg = aggregate(sp_dat, SG, function(x) names(which.max(table(x)))[1] )
  agg@grid@cellsize <- c(1,1)
  df <- agg %>% as.data.frame %>% rename(x = s1, y = s2)  %>% select(colnames(image_df))

  return(df)

}

process_image <- function(image_file_name, k_list){
  ## process_image(image_file_name, k_list) performs k-means clustering on an
  ## image for different numbers of cluster centers.
  ##
  ## Input:
  ## - image_file_name: A string containing the PNG or JPEG image file name.
  ## - k_list:          A list of the number of centers in each clustering.
  ##                    Must be greater than or equal to 1.
  ##
  ## Output:
  ## - A list or tibble of tibbles containing summaries of the clusterings
  ##   on various levels.
  ## Example:
  ##   library(imager)
  ##   library(tidyverse)
  ##   library(tidymodels)
  ##   library(dmc)
  ##   fpath <- system.file('extdata/parrots.png', package='imager')
  ##   cluster_info <- process_image(fpath, c(2:10))
  im <- imager::load.image(image_file_name)
  tidy_dat <- as.data.frame(im, wide ="c") %>%
    rename(R = c.1, G = c.2, B = c.3) %>%
    select(x, y, R, G, B)
  dat <- select(tidy_dat, c(-x, -y))
  kclusts <-
    tibble(k_list) %>%
    mutate(
      kclust = map(k_list, ~kmeans(x=dat, centers=.x, nstart=2)),
      augmented=map(kclust, augment, tidy_dat),
      glanced = map(kclust, glance),
      tidied = map(kclust, tidy)
    )
  clusters <-
    kclusts %>%
    unnest(cols = c(tidied))
  # Specifying colors
  R<-as.data.frame(clusters["R"])[[1]]
  G<-as.data.frame(clusters["G"])[[1]]
  B<-as.data.frame(clusters["B"])[[1]]
  n<-length(R)
  hex<-tibble(as.character(c(1:n)))
  dmc_id<-tibble(as.character(c(1:n)))
  dmc_name<-tibble(as.character(c(1:n)))
  for (i in 1:n){
    c<-rgb(R[i], G[i], B[i])
    dmc_colour<-dmc(c, visualize=FALSE)
    hex[i, 1] <- dmc_colour[[3]]
    dmc_id[i, 1]<-dmc_colour[[1]]
    dmc_name[i, 1]<-paste(dmc_colour[[2]]," (", dmc_colour[[1]], ")",sep="")
  }
  # Clusters with colours
  clusters <-
    clusters %>%
    mutate("dmc_name"=dmc_name[[1]],
           "hex"=hex[[1]],
           "dmc"=dmc_id[[1]]) %>%
    select(-kclust,-augmented,-glanced,-size)
  # Clustering summaries per level
  cluster_info<-list(kclusts, clusters)
  cluster_info
}

scree_plot <- function(cluster_info){
  ## scree_plot(cluster_info) produces a scree plot for evaluating the number
  ## of clusters to use for the cross-stitch pattern.
  ##
  ## Input:
  ## - cluster_info: A list containing tibbles of clustering summaries.
  ##
  ## Output:
  ## - Produces a scree plot.
  ##
  ## Example:
  ##   library(imager)
  ##   library(tidyverse)
  ##   library(tidymodels)
  ##   library(dmc)
  ##   fpath <- system.file('extdata/parrots.png', package='imager')
  ##   cluster_info <- process_image(fpath, c(2:10))
  ##   scree_plot(cluster_info)
  kclusts <- cluster_info[[1]]
  clusterings <- kclusts %>%
    unnest(cols = c(glanced))
  ggplot(clusterings, aes(k_list, tot.withinss)) +
    ggtitle("Scree Plot") +
    xlab("k (number of clusters)") + ylab("total within sum of squares") +
    geom_line() +
    geom_point() +
    theme(legend.position="bottom",
          legend.box="vertical")
}

colour_strips <- function(cluster_info){
  ## colour_strips(cluster_info) produces the labelled DMC color strips of the
  ## cross-stitch pattern.
  ## Input:
  ## - cluster_info: A list containing tibbles of clustering summaries.
  ##
  ## Output:
  ## - Produces a grid of labelled color strips.
  ##
  ## Example:
  ##   library(imager)
  ##   library(tidyverse)
  ##   library(tidymodels)
  ##   library(dmc)
  ##   fpath <- system.file('extdata/parrots.png', package='imager')
  ##   cluster_info <- process_image(fpath, c(2:10))
  ##   colour_strips(cluster_info)
  clusters<-cluster_info[[2]]
  square <- function(hex, label_size, k){
    dat<-tibble(x=c(5, 5), y=c(8, 2), l=c(hex, paste("k=",k,sep="")))
    ggplot(dat, aes(x, y, label=l))  +
      geom_text() + geom_label(label.size=.1, alpha=.2) +
      coord_fixed(xlim=c(0,10), ylim = c(0,10)) + theme_void() +
      theme(plot.background = element_rect(fill = hex))
  }
  t <- tibble(colours = clusters$hex,
              k = clusters$k_list,
              squares = purrr::map2(colours, k, ~ square(.x, 24/length(colours), .y)))
  plot_grid(plotlist = t$squares)
  
}

make_pattern <- function(cluster_info, k, x_size,
                         black_white = FALSE, background_colour = NULL){
  ## make_pattern(cluster_info, k, x_size, black_white, background_colour)
  ## produces the plot of the cross-stitch pattern.
  ## Input:
  ## - cluster_info:      A list containing tibbles of clustering summaries.
  ## - k:                 The number of cluster centers chosen for the pattern.
  ## - x_size:            The approximate total number of possible stitches in
  ##                      the horizontal direction.
  ## - black_white:       A logical variable which sets the pattern in black and
  ##                      white if TRUE or colour if FALSE (default).
  ## - background_colour: A character which is the DMC color of the background
  ##                      to be excluded from the pattern.
  ## Output:
  ## - Produces a plot of the cross-stitch pattern.
  ##
  ## Example:
  ##   library(imager)
  ##   library(tidyverse)
  ##   library(tidymodels)
  ##   library(sp)
  ##   library(cowplot)
  ##   library(dmc)
  ##   fpath <- system.file('extdata/parrots.png', package='imager')
  ##   cluster_info <- process_image(fpath, c(2:10))
  ##   make_pattern(cluster_info, 3)
  kclusts <- cluster_info[[1]]
  clusters<-cluster_info[[2]]
  assignments <-
    kclusts %>%
    unnest(cols = c(augmented)) %>%
    filter(k_list==k) %>%
    select(-kclust, -glanced, -tidied)
  
  clusters_colours <-
    clusters %>%
    filter(k_list==k) %>%
    select(cluster, hex, dmc_name, dmc)
  
  coloured_img_df<-assignments %>%
    filter(k_list==k) %>%
    inner_join(clusters_colours, by=c(".cluster"="cluster")) %>%
    mutate(cluster=as.numeric(as.character(.cluster))) %>%
    select(x,y, cluster, hex, dmc_name, dmc)
  
  # Subsample image
  sampled_im<-change_resolution(coloured_img_df, x_size)
  sampled_im<-sampled_im %>%
    mutate(cluster=as.numeric(cluster))
  
  # Plotting
  if (!is.null(background_colour)){
    sampled_im<-sampled_im %>%
      filter(dmc!=background_colour)
  } else {
    
  }
  
  cluster_frame <- clusters_colours %>%
    select(cluster, hex, dmc_name)
  
  if (black_white == TRUE){
    sampled_im %>% 
      ggplot(aes(x, y, shape=factor(cluster))) +
      guides(shape = guide_legend(ncol = 1)) +
      geom_point() + 
      scale_y_reverse() +
      scale_shape_manual(name="DMC",
                         values = cluster_frame %>% select(cluster, cluster) %>% deframe,
                         label =  cluster_frame %>% select(cluster, dmc_name) %>% deframe)+
      theme_half_open() +
      background_grid(major="xy",
                      minor="xy",
                      size.major = 1,
                      size.minor = 1,
                      color.major = "black",
                      color.minor = "grey85") +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.box="vertical",
            legend.title=element_text(size=8),
            legend.text=element_text(size=7)) +
      labs(x="", y="") +
      coord_fixed(ratio=max(sampled_im$y)/max(sampled_im$x)) +
      ggtitle("Cross-stitch Pattern")
    
  } else {
    sampled_im %>% 
      ggplot(aes(x, y, col=factor(cluster), shape=factor(cluster))) +
      geom_point() + 
      scale_y_reverse() +
      scale_colour_manual(name="DMC",
                          values = cluster_frame %>% select(cluster, hex) %>% deframe,
                          label =  cluster_frame %>% select(cluster, dmc_name) %>% deframe) +
      scale_shape_manual(name="DMC",
                         values = cluster_frame %>% select(cluster, cluster) %>% deframe,
                         label =  cluster_frame %>% select(cluster, dmc_name) %>% deframe)+
      theme_half_open() +
      background_grid(major="xy",
                      minor="xy",
                      size.major = 1,
                      size.minor = 1,
                      color.major = "black",
                      color.minor = "grey85") +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.box="vertical",
            legend.title=element_text(size=8),
            legend.text=element_text(size=7)) +
      labs(x="", y="") +
      coord_fixed(ratio=max(sampled_im$y)/max(sampled_im$x)) +
      guides(col = guide_legend(ncol = 1)) +
      ggtitle("Cross-stitch Pattern")
  }
}