plotGroups.mk <- function (data, edesign = edesign, show.fit = FALSE, dis = NULL, 
                           step.method="backward", min.obs = 2, 
                           alfa = 0.05, nvar.correction = FALSE,
                           show.lines = TRUE, time=edesign$Time,
                           group.cols = c(3:ncol(edesign)), 
                           groups.vector = groups.vector, 
                           repvect = edesign$Replicate,  
                           summary.mode = "median", xlab="Time", 
                           ylab = NULL, cex.xaxis = 1, ylim = NULL, 
                           main = NULL, cexlab = 0.8, legend = TRUE, 
                           sub = NULL, item = NULL) 
{
  ##> test with:
  # data <- dat[cut == j, grep(names.groups[1],colnames(dat),value=T)] #>> only counts of LPS, j=1
  # time_var = edesign$Time[which(edesign$LPS == 1)]

  groups_var <- as.data.frame(groups.vector)
  
  ## lps_samples <- rownames(groups_var)[which(groups_var$LPS == 1)]
  ## ns_samples  <- rownames(groups_var)[which(groups_var$NS == 1)]
  
  if (!is.vector(data)) { #--> yes !
      if (summary.mode == "representative") {
        distances <- apply(as.matrix(dist(data, diag = TRUE, 
                                          upper = TRUE)), 1, sum)
        representative <- names(distances)[distances == min(distances)]
        yy <- as.numeric(data[rownames(data) == representative, ])
        sub <- paste("Representative:", representative)
      }
      else if (summary.mode == "median") {
        yy <- apply(as.matrix(data), 2, median, na.rm = TRUE)
        if (is.null(sub)) {
          # subtitle of plot below  axis
          sub <- paste("Median profile of", nrow(data), item, sep = " ")
        }
      }
      else stop("not valid summary.mode")
      if (dim(data)[1] == 1) {
        sub <- rownames(data)
      }
  
  } else if (length(data) != 0) {
    yy <- as.numeric(data)
    sub <- rownames(data)
  } else stop("empty data")
  
  
  if (is.null(ncol(groups_var))) {
    ncol = 1
    legend = FALSE
    codeg = "group"
  } else {
    ncol = ncol(groups_var)
    codeg <- as.character(colnames(groups_var))
  }
  
  reps <- i.rank(repvect) ## for transforming column in vector to be accessed with [1] in i.rank 
  y <- vector(mode = "numeric", length = length(unique(reps)))
  x <- vector(mode = "numeric", length = length(unique(reps)))
  g <- matrix(nrow = length(unique(reps)), ncol = ncol)
  
  for (i in 1:length(y)) {
    y[i] <- mean(yy[reps == i], na.rm = TRUE)
    x[i] <- mean(time_var[reps == i])
    for (j in 1:ncol) {
      g[i, j] <- mean(groups_var[reps == i, j])
    }
  }
  if (is.null(ylim)) {
      ylim = c(min(as.numeric(yy), na.rm = TRUE), 
               max(as.numeric(yy), na.rm = TRUE))
  }
  abcissa <- x
  xlim = c(min(abcissa, na.rm = TRUE), max(abcissa, na.rm = TRUE) * 1.3)
  color1 <- as.numeric(sort(factor(colnames(groups_var)))) + 1
  color2 <- groups_var
  
  for (j in 1:ncol) {
      color2[, j] <- color2[, j] * j
  }
  color2 <- as.vector(apply(color2, 1, sum) + 1)
  
  # plot(x = time_var, y = yy, pch = 19, xlab = xlab, ylab = ylab,
  #      xaxt = "n", main = main, sub = sub, ylim = ylim,
  #      xlim = xlim, cex = cexlab, col = color2 )
  # axis(1, at = unique(abcissa), labels = unique(abcissa), cex.axis = cex.xaxis)

  
  library(ggplot2)
  yy.mlt        <- reshape2::melt(yy,value.name = "expression_value")
  yy.mlt$time   <- as.character(time_var)
  yy.mlt$group  <- ifelse(edesign$Mort==1,"Mort","Viv")
  
  # library(ggpubr)
  # p <- ggscatter(data=yy.mlt, x="time", y="expression_value", 
  #                 color = "group",shape=19, 
  #                 size=2,title = main , 
  #                 xlab=xlab, ylab=ylab,
  #                 rug = TRUE , palette = "jco" ) 
  #                 # rug = bins on y axis for each value ,shows density of value range.
  # or : 
  
  p <- ggplot(data=yy.mlt, aes(x=time, y=expression_value,colour=group)) + 
       geom_point( size=3) + 
       ggtitle(main) + 
       theme(axis.text=element_text(size=10)) +
       geom_rug() +
       stat_summary( aes(group=group), geom="line",fun="mean") + 
       theme(plot.margin=unit(c(1, 0, 1, 0), "points")) 
       ##+ xlab(xlab) + ylab(ylab) + 
  
  if (show.fit) {
    rm <- matrix(yy, nrow = 1, ncol = length(yy))
    rownames(rm) <- c("ratio medio")
    colnames(rm) <- rownames(dis)
    fit.y <- T.fit(rm, design = dis, step.method = step.method, 
                   min.obs = min.obs, alfa = alfa, 
                   nvar.correction = nvar.correction)
    betas <- fit.y$coefficients
  }
  
  
  # for (i in 1:ncol(groups_var)) {
  #   group <- g[, i]
  #   if ((show.fit) && !is.null(betas)) {
  #     li <- c(2:6)
  #     a <- reg.coeffs(coefficients = betas,
  #                     groups.vector = groups.vector,
  #                     group = colnames(groups)[i] )
  #     a <- c(a, rep(0, (7 - length(a))))
  #     curve(a[1] + a[2] * x + a[3] * (x^2) + a[4] * (x^3) +
  #           a[5] * (x^4) + a[6] * (x^5) + a[7] * (x^5), from = min(time),
  #           to = max(time), col = color1[i], add = TRUE,
  #           lty = li[i] )
  #   }
  #   if (show.lines) {
  #     lx <- abcissa[group != 0]
  #     ly <- y[group != 0]
  #     ord <- order(lx)
  #     lxo <- lx[ord]
  #     lyo <- ly[ord]
  #     lines( lxo, lyo, col = color1[i] )
  #   }
  # }
  
  ## Create a df for Lines Graph --replace for_loop:
  # if (show.lines) {
      # group_timeline <- list()
      # for (i in 1:ncol(groups_var)) {
      #   group <- g[, i]
      #   lx <- abcissa[group != 0]
      #   ly <- y[group != 0]
      #   ord <- order(lx)
      #   lxo <- lx[ord]
      #   lyo <- ly[ord]
      #   t   <- names(groups_var[i])
      #   # Line for one group:
      #   group_timeline[[i]] <- data.frame("time"=as.character(lxo), "expression_value"=lyo, "group"=t )
      # }
      # group_timeline <- Reduce(f=rbind, x = group_timeline)
      # p + geom_line(data = group_timeline ,
      #               aes(x=factor(time),y=expression_value,col=group))
    
   # }
   
  # op <- par(bg = "white")
  # if (legend) 
  #   legend(max(abcissa, na.rm = TRUE) * 1.02, ylim[1], legend = codeg, 
  #          text.col = color1, col = color1, cex = cexlab, lty = 1, 
  #          yjust = 0)
  # par(op)
  
  return( p )
}