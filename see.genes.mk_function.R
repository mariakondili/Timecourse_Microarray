
see.genes.mk <- function (data, edesign=data$edesign,
                          group.cols = c(3:ncol(edesign)), 
                          cluster.data = 1, cluster.method = "hclust",
                          k = clust_num, dis = dis, 
                          show.fit = FALSE, summary.mode = "median",
                          k.mclust = FALSE, distance = "cor", 
                          agglo.method = "ward.D", step.method = "backward", min.obs = 3, 
                          alfa = 0.05,nvar.correction = FALSE, show.lines = TRUE, iter.max = 500, 
                          color.mode = "rainbow", cexlab = 1, legend = TRUE, newX11 = TRUE, 
                          ylim = NULL,  main = NULL,item = "genes", ...) 
{
  
  
  # see.genes.mk(sig_data[[comparison]],
  #              names.groups = colnames(edesign)[3:4],
  #              show.fit = F,
  #              dis = mydesign.mat$dis,
  #              cluster.method="hclust",
  #              cluster.data = 1, k = clust_num ) 
  
  source("plotGroups.mk_function.R")
  library(gridExtra)
  
  edesign = data$edesign
  time_var = edesign$Time
  repvect = edesign$Replicate
  group.cols = c(3:ncol(edesign))
  groups_var = edesign[,group.cols]
  narrays <- length(time_var)
  
  if (!is.null(dim(data))) {
      dat <- as.data.frame(data)
      clusterdata <- data
  } else {
      clusterdata <- data[[ cluster.data ]]
      dat <- as.data.frame(data$sig.profiles)
  }
  
  
  if (nrow(dat) > 1) {
    dat <- as.data.frame(dat[,(ncol(dat) - length(time_var) + 1):ncol(dat)]) 
    # only 1 column,when length(time_var)
    # all columns if nrow(time_var) is used 
    #--------------------------------- NAs REMOVAL ---------------------------------
    count.noNa <- function(x) (length(x) - length(x[is.na(x)]))
    # count_na <- function(x) (length(x[is.na(x)]))
    dat <- dat[which(apply(as.matrix(dat), 1, count.noNa) >= length(unique(repvect))),]
    clusterdata <- dat
  
    #--------------------------------- NAs TREATMENT ---------------------------------
    # If cluster.data is beta or t.values, they are changed for 0s.
    # With profiles, for mean value of the same type of replicate.
  
    if (any(is.na(clusterdata))) {
    
        if (cluster.method == "kmeans" || cluster.method == "Mclust") {
          if (all(cluster.data != 1, cluster.data != "sig.profiles"))   {
              clusterdata[is.na(clusterdata)] <- 0
          } else {
            mean.replic <- function(x) { tapply(as.numeric(x), repvect, mean, na.rm = TRUE) }
            MR <- t(apply(clusterdata, 1, mean.replic))
            if (any(is.na(MR))) {
                row.mean <- t(apply(MR, 1, mean, na.rm = TRUE))
                MRR <- matrix(row.mean, nrow(MR), ncol(MR))
                MR[is.na(MR)] <- MRR[is.na(MR)]
            }
        
            data.noNA <- matrix(NA, nrow(clusterdata),ncol(clusterdata))
            u.repvect <- unique(repvect)
        
            for (i in 1:nrow(clusterdata)) {
                for (j in 1:length(u.repvect)) {
                  data.noNA[i, repvect == u.repvect[j]] = MR[i,u.repvect[j]]
                }
              }
            clusterdata <- data.noNA
          } # close "else"
        }
      } # close if any is.na
  
    
    #--------------------------------------------------------------------------
    # CLUSTERING
    #-------------------------------------------------------------------------
  
    if (!is.null(clusterdata)) {
        k <- min(k, nrow(dat), na.rm = TRUE)
        if (cluster.method == "hclust") {
            if (distance == "cor") {
              dcorrel <- matrix(rep(1, (nrow(clusterdata)^2)), 
                                nrow(clusterdata), 
                                nrow(clusterdata)) - cor(t(clusterdata), use = "pairwise.complete.obs")
              clust <- hclust(as.dist(dcorrel), method = agglo.method)
              c.algo.used = paste(cluster.method, "cor",  agglo.method, sep = "_")
            } else {
              clust <- hclust(dist(clusterdata, method = distance), 
                              method = agglo.method)
              c.algo.used = paste(cluster.method, distance, 
                                  agglo.method, sep = "_")
            }
            cut <- cutree(clust, k = k)
        }
        else if (cluster.method == "kmeans") {
          cut <- kmeans(clusterdata, k, iter.max)$cluster
          c.algo.used = paste("kmeans", k, iter.max, sep = "_")
        }
        else if (cluster.method == "Mclust") {
          library(mclust)
          if (k.mclust) {
            my.mclust <- Mclust(clusterdata)
            k = my.mclust$G
          } else {
            my.mclust <- Mclust(clusterdata, k)
          }
          cut <- my.mclust$class
          c.algo.used = paste("Mclust", k, sep = "_")
        }
        else stop("Invalid cluster algorithm")
        # if (newX11) 
        #  X11()
        groups_var <- as.matrix(groups_var)
        colnames(groups_var) <- colnames(edesign)[group.cols]
        # if (k <= 4) 
        #   par(mfrow = c(2, 2))
        # else if (k <= 6) 
        #   par(mfrow = c(3, 2))
        # else if (k > 6 & k <= 9) 
        #   par(mfrow = c(3, 3))
        # if (k > 9)
        #   par(mfrow = c(5,4))
        # for (i in 1:(k)) {
        #   X11(); PlotProfiles(data = dat[cut == i, ], repvect = repvect, 
        #                        main = i, ylim = ylim, color.mode = color.mode, 
        #                        cond = rownames(edesign), item = item, ...)
        # }
        # if (newX11) 
        #   X11()
        # if (k <= 4) {
        #   par(mfrow = c(2, 2))
        #   cexlab = 0.6
        # }
        # else if (k <= 6) {
        #   par(mfrow = c(3, 2))
        #   cexlab = 0.6
        # }
        # else if (k > 6 & k <= 9 ) {
        #   par(mfrow = c(3, 3))
        #   cexlab = 0.35
        # }
        # else if (k > 9) {
        #   par(mfrow = c(5, 4))
        #   cexlab = 0.30
        # }
       
        plots_list <- list()
        for ( j in 1:k ) {
           #X11() #--> Each graph will appear in a Window
           cl_g <- which(cut == j)
           p <- plotGroups.mk(data = dat[cl_g, ], edesign =edesign, show.fit = show.fit,
                       dis = dis, step.method = step.method, min.obs = min.obs,
                       alfa = alfa, nvar.correction = nvar.correction,
                       show.lines = show.lines, time = time_var,  
                       groups.vector = groups_var, repvect = repvect, 
                       summary.mode = summary.mode, xlab = "Time", 
                       ylab = "Expression LPS/NS", 
                       main = paste("Cluster", j, sep = " "), 
                       ylim = ylim, cexlab = cexlab,legend = legend,item = item)
           plots_list[[j]] <- p
          }
      } else {
        print("warning: impossible to compute hierarchical clustering")
        c.algo.used <- NULL
        cut <- 1
       }
  
  # close if : line 41 
  } else {
    print("warning: NULL data. No visualization possible")
    c.algo.used <- NULL
    cut <- NULL
  }
  
  # else if (nrow(dat) == 1) {
  #     if (newX11) 
  #       X11()
  #       PlotProfiles(data = dat, repvect = repvect, main = NULL, 
  #                     ylim = ylim, color.mode = color.mode, cond = rownames(edesign), 
  #                   ...)
  #     if (newX11) 
  #         #X11()
  #         plots_list <- PlotGroups.mk(data = dat, show.fit = show.fit, dis = dis, 
  #                    step.method = step.method, min.obs = min.obs, alfa = alfa, 
  #                    nvar.correction = nvar.correction, show.lines = show.lines, 
  #                    time = time, groups = groups, repvect = repvect, 
  #                    summary.mode = summary.mode, xlab = "time", 
  #                    main = main, ylim = ylim, cexlab = cexlab, legend = legend, 
  #                   ...)
  #         c.algo.used <- NULL
  #         cut <- 1
  # 
  # } 
  
  
  OUTPUT <- list(cut, c.algo.used, groups, plots_list)
  names(OUTPUT) <- c("cut", "cluster.algorithm.used", "groups","plots_list")

  OUTPUT
}