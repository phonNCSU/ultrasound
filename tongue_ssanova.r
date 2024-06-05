#######################################################################################
# tongue_ssanova.r                                revised May 7, 2019
# Jeff Mielke
# functions for SSANOVA comparisons of tongue traces in polar coordinates using gss 
#######################################################################################
#
# BASIC COMMAND TO GENERATE AN SSANOVA PLOT (IF 'phone' IS THE NAME OF YOUR FACTOR)
#  ss <- polar.ssanova(data, 'phone')
#
# BASIC COMMAND TO PLOT THE RAW DATA
#  show.traces(data)
#
# TO PLOT TO FILE, SEPARATING BY TWO DIFFERENT FACTORS (COLUMNS IN YOUR DATA FRAME):
#  cairo_pdf('my_ssanova_pdf.pdf', h=4.5, w=5, onefile=T)
#    ss.by.C <- polar.ssanova(data, 'consonant')
#    ss.by.V <- polar.ssanova(data, 'vowel')
#  dev.off()
#
# TO HIGHLIGHT RAW DATA FOR THE LEVEL ('I'):
#  show.traces(data, c('I'))
#
# DATA FILE SHOULD BE ORGANIZED LIKE THIS (MULTIPLE COLUMNS CAN BE USED INSTEAD OF word):
#
# word,token,X,Y
# dog,1,307,262
# dog,1,311,249
# dog,1,315,240
# dog,2,308,261
# dog,2,311,250
# dog,2,314,249
# cat,1,307,240
# dog,2,311,250
# dog,2,314,259
# ...
#
#######################################################################################
#
# polar.ssanova() ARGUMENTS (ALL OPTIONAL EXCEPT data):
#
#           data: your tongue tracings (minimally including columns X and Y and a
#                 column with a factor)
#       data.cat: the factor to use to categorize the data (defaults to 'word')
#          scale: how much to scale the axis values (e.g. to convert from pixels to 
#                 centimeters)
#  origin.method: how to choose the origin for calculating polar coordinates
#          debug: whether to generate the cartesian and non-transformed polar plots too
#       plotting: whether to plot anything (or just return the result of the test)
#           main: the main title for the plot
#        CI.fill: whether to indicate confidence intervals with shading (like ggplot) 
#                 or with dotted lines (like the earlier SSANOVA code).  
#                 Defaults to FALSE (dotted lines)
#       printing: if TRUE, different splines use different line types, so that the
#                 figure can be printed in black and white.
#           flip: whether to flip the Y values (useful for plotting data from images
#                 in cartesian coordinates, but ignored if using polar coordinates)
# cartesian.only: used by cart.ssanova()
#       is.polar: if TRUE, the data is already in polar coordinates
#
#######################################################################################
#
#  cart.ssanova() SAME AS polar.ssanova() BUT DOESN'T USE POLAR COORDINATES
#
#######################################################################################
#
#   show.traces() ARGUMENTS (ALL OPTIONAL EXCEPT data): 
#
#           data: your tongue tracings (minimally including columns X and Y and a
#                 column with a factor)
#       data.cat: the factor to use to categorize the tongues (defaults to 'word')
#   to.highlight: a list of factor levels to plot while muting the other levels
#        to.plot: a list of factor levels to plot, excluding the rest (defaults to all)
#    token.label: the factor to use to identify individual tokens (defaults to 'token')
#           flip: whether to flip the Y values (useful for plotting data from images)
#           main: the main title for the plot
#       overplot: whether to add the traces to an existing plot
#       is.polar: if TRUE, the data is already in polar coordinates
#         origin: used if the data is in polar coordinates already
#
#######################################################################################

library(gss)
library(plyr)

# CONVERT POLAR COORDINATES TO CARTESIAN COORDINATES
# make.cartesian <- function(tr, origin=c(0,0), flip=TRUE){    
#     X <- apply(tr, 1, function(x,y) origin[1]-x[2]*cos(x[1]))
#     Y <- apply(tr, 1, function(x,y) x[2]*sin(x[1])-origin[2])
#     # this was added 12/13/2021 to prevent positive unflipped values from coming back with negative Y values
#     if (flip){
#         xy <- cbind(X, Y)
#     }else{
#         xy <- cbind(X, -Y)
#         xy <- cbind(X, Y)
#     }
#     return(xy)
# }


#CONVERT POLAR COORDINATES TO CARTESIAN COORDINATES
make.cartesian <- function(tr, origin=c(0,0), flip=TRUE){    
  # X <- apply(tr, 1, function(x,y) origin[1]-x[2]*cos(x[1]))
  # Y <- apply(tr, 1, function(x,y) x[2]*sin(x[1])-origin[2])
  X <- origin[1]-tr[,2]*cos(tr[,1])
  Y <- tr[,2]*sin(tr[,1])-origin[2]
  # this was added 12/13/2021 to prevent positive unflipped values from coming back with negative Y values
  if (flip){
    xy <- cbind(X, Y)
  }else{
    xy <- cbind(X, -Y)
  }
  return(xy)
}


# make.cartesian <- function(tr, origin=c(0,0), flip=TRUE){    
#     # this was added 12/13/2021 to prevent positive unflipped values from coming back with negative Y values
#     # origin[1] = -origin[1]
#     # xy[,1] = -xy[,1]

#     if (flip){
#         # X <- apply(tr, 1, function(x,y) origin[1]-x[2]*cos(x[1]))
#         # Y <- apply(tr, 1, function(x,y) -x[2]*sin(x[1])-origin[2])
#         X <- origin[1]-tr[,2]*cos(tr[,1])
#         Y <- -tr[,2]*sin(tr[,1])-origin[2]
#         xy <- cbind(X, Y)
#     }else{
#         # X <- apply(tr, 1, function(x,y) origin[1]-x[2]*cos(x[1]))
#         # Y <- apply(tr, 1, function(x,y) -x[2]*sin(x[1])-origin[2])
#         X <- origin[1]-tr[2]*cos(tr[1])
#         Y <- -tr[,2]*sin(tr[,1])-origin[2]
#         xy <- cbind(X, -Y)
#     }
#     return(xy)
# }

#CONVERT CARTESIAN COORDINATES TO POLAR COORDINATES
make.polar <- function(data.xy, origin=c(0,0), flip=TRUE){
    # xy <- cbind(data.xy$X, data.xy$Y)
    xy <- data.xy[,c('X','Y')]
    # all_r <- apply(xy, 1, function(x) sqrt((x[1]-origin[1])^2 + (x[2]-origin[2])^2))
    all_r <- sqrt((xy[,1]-origin[1])^2 + (xy[,2]-origin[2])^2)

    # to handle origin higher than tongue root

    origin[1] = -origin[1]
    xy[,1] = -xy[,1]

    if (flip){
        # all_theta <- apply(xy, 1, function(x,y) atan2(x[2]-origin[2], x[1]-origin[1]))
        all_theta <- atan2(xy[,2]-origin[2], xy[,1]-origin[1])
    }else{
        all_theta <- atan2(origin[2]-xy[,2], xy[,1]-origin[1])
    #     all_theta <- -atan2(xy[,2]-origin[2], -xy[,1]-origin[1])
    #     # all_theta <-   -apply(xy, 1, function(x,y) atan2(x[2]-origin[2], x[1]+origin[1]))
    #     # all_theta <-   -apply(xy, 1, function(x,y) atan2(x[2]-origin[2], -x[1]-origin[1]))
    }

    data.tr <- data.xy
    data.tr$X <- all_theta
    data.tr$Y <- all_r
    return(data.tr)
}
# through 1/23/22 adding pi after calculating theta
# from 1/24/22 if not flipping, changing sign of x and then changing sign of theta and not adding pi

#RESCALE DATA FROM PIXELS TO CENTIMETERS
us.rescale<-function(data, usscale, X='X', Y='Y'){
    data[,c(X)] <- data[,c(X)]*usscale
    data[,c(Y)] <- data[,c(Y)]*usscale
    data
}

#SELECT AN APPROPRIATE ORIGIN FOR THE DATA

# if my_ssanova is an ssanova result based on my_data (with flip=TRUE), then these two commands should return the same thing:
#   my_ssanova$origin*c(1,-1)
#   select.origin(my_data$X, my_data$Y, my_data$token, flip=TRUE)
   


select.origin <- function(Xs, Ys, tokens, method='xmean_ymin', flip){

    require(plyr)
    
    if (method=='xmean_ymin'){
        # if (mean(Ys)>0){
        if (flip){
            return(c(mean(Xs), max(Ys)*1.01))
        }else{
            return(c(mean(Xs), min(Ys)-0.01*diff(range(Ys))))
        }
    }else if (method=='yextremes'){
        if (mean(Ys)>0){
            return(c(mean(Xs[which(Ys==min(Ys))]), max(Ys)*1.01))
        }else{
            return(c(mean(Xs[which(Ys==max(Ys))]), min(Ys)-0.01*diff(range(Ys))))
        }
    }else if (method=='xmean_ymean'){
        origin = c(mean(Xs), mean(Ys))
        origin[2] = min(origin[2], Ys[Xs>origin[1]])
        origin[1] = max(origin[1], Xs[Ys<origin[2]])
        return(origin)
    }else if (method=='xmid_ymid'){
        origin = c(mean(range(Xs)), mean(range(Ys)))
        origin[2] = min(origin[2], Ys[Xs>origin[1]])
        origin[1] = max(origin[1], Xs[Ys<origin[2]])
        return(origin)
    }else if (method=='x75_y25'){
        origin = c(min(Xs)+0.75*diff(range(Xs)), min(Ys)+0.25*diff(range(Ys)))
        origin[2] = min(origin[2], Ys[Xs>origin[1]])
        origin[1] = max(origin[1], Xs[Ys<origin[2]])
        return(origin)
    }else if (method=='x50_y25'){
        origin = c(min(Xs)+0.50*diff(range(Xs)), min(Ys)+0.25*diff(range(Ys)))
        origin[2] = min(origin[2], Ys[Xs>origin[1]])
        origin[1] = max(origin[1], Xs[Ys<origin[2]])
        return(origin)
    }else if (method=='x75_y50'){
        origin = c(min(Xs)+0.75*diff(range(Xs)), min(Ys)+0.5*diff(range(Ys)))
        origin[2] = min(origin[2], Ys[Xs>origin[1]])
        origin[1] = max(origin[1], Xs[Ys<origin[2]])
        return(origin)
    }else if (method=='mid_arc'){
        ss_data <- data.frame(token=tokens, X=Xs, Y=Ys)
        ss_ends <- ddply(ss_data, .(token), summarize, firstX=X[which(Y==max(Y[X<mean(X)]))], firstY=max(Y[X<mean(X)]), lastX=max(X), lastY=Y[which(X==max(X))])
        ss_first <- c(median(ss_ends$firstX), median(ss_ends$firstY))
        ss_last <- c(median(ss_ends$lastX), median(ss_ends$lastY))
        ss_mid <- ss_first/2 + ss_last/2 
        return(ss_mid)   
    }else if (method=='mid_arc_mod'){
        ss_data <- data.frame(token=tokens, X=Xs, Y=Ys)
        ss_ends <- ddply(ss_data, .(token), summarize, firstX=X[which(Y==max(Y[X<mean(X)]))], firstY=max(Y[X<mean(X)]), lastX=max(X), lastY=Y[which(X==max(X))])
        ss_first <- c(median(ss_ends$firstX), median(ss_ends$firstY))
        ss_last <- c(median(ss_ends$lastX), median(ss_ends$lastY))
        #ss_mid <- ss_first/3 + 2*ss_last/3 
        ss_mid <- ss_first/2 + ss_last/2 
        return(c(ss_mid[1], max(Ys)*1.01))
    }else if (method=='xmax_ymin'){
        print(paste('using origin method xmax_ymin'))
        # if (mean(Ys)>0){
        if (flip){
            return(c(max(Xs), max(Ys)*1.01))
        }else{
            return(c(max(Xs), min(Ys)-0.01*diff(range(Ys))))
        }
    }else if (method=='xminmax_ymin'){
        print(paste('using origin method xminmax_ymin'))
        # if (mean(Ys)>0){
        if (flip){
            return(c(min(ddply(data.frame(X=Xs,Y=Ys,token=tokens), .(token), summarize, maxX=max(X))$maxX)*1.01, max(Ys)*1.01))
        }else{
            return(c(min(ddply(data.frame(X=Xs,Y=Ys,token=tokens), .(token), summarize, maxX=max(X))$maxX)*1.01, min(Ys)*1.01))
        }
    }else if (method=='xmaxmin_ymin'){
        print(paste('using origin method xmaxmin_ymin'))
        # if (mean(Ys)>0){
        if (flip){
            return(c(max(ddply(data.frame(X=Xs,Y=Ys,token=tokens), .(token), summarize, minX=min(X))$minX)*0.99, max(Ys)*1.01))
        }else{
            return(c(max(ddply(data.frame(X=Xs,Y=Ys,token=tokens), .(token), summarize, minX=min(X))$minX)*0.99, min(Ys)*1.01))
        }
    }

    return(c(mean(Xs), max(Ys)*1.01))
}

crop.result <- function(data, ss.result, data.cat){
    new.result <- c()
    for (data.level in unique(ss.result[,data.cat])){
        sub.data <- data[data[,data.cat]==data.level,]
        sub.result <- ss.result[ss.result[,data.cat]==data.level,]
        new.result <- rbind(new.result, subset(sub.result, sub.result$X>=min(sub.data$X) & sub.result$X<=max(sub.data$X)))
    }
    new.result
}

#PERFORM THE SSANOVA AND RETURN THE RESULTING SPLINES AND CONFIDENCE INTERVALS
#expand.grid + predict scheme based on http://www.ling.upenn.edu/~joseff/papers/fruehwald_ssanova.pdf
tongue.ss <- function(data, data.cat='word', flip=FALSE, length.out=1000, SDs=1.96, alpha=1.4, crop=FALSE, drop_levels=FALSE, verbose=FALSE){    

    require(gss)

    if (flip==TRUE){
        data$Y <- -data$Y
    }
    data$tempword <- data[,data.cat]
    if (drop_levels){
        data$tempword <- droplevels(data$tempword)
    }
    #print(summary(lm(Y ~ tempword * X, data=data)))
    # plot(data$X, data$Y)
    ss.model <- ssanova(Y ~ tempword + X + tempword:X, data=data, alpha=alpha)
    ss.result <- expand.grid(X=seq(min(data$X), max(data$X), length.out=length.out), tempword=levels(data$tempword))
    ss.result$ss.Fit <- predict(ss.model, newdata=ss.result, se=T)$fit
    ss.result$ss.cart.SE  <- predict(ss.model, newdata=ss.result, se=T)$se.fit
    #print(names(ss.result))
    #print(aggregate(ss.Fit ~ tempword, FUN=mean, data=ss.result))
    #print(aggregate(ss.cart.SE ~ tempword, FUN=mean, data=ss.result))
    ss.result$ss.upper.CI.X <- ss.result$X
    ss.result$ss.upper.CI.Y <- ss.result$ss.Fit + SDs*ss.result$ss.cart.SE
    ss.result$ss.lower.CI.X <- ss.result$X
    ss.result$ss.lower.CI.Y <- ss.result$ss.Fit - SDs*ss.result$ss.cart.SE
    names(ss.result)[which(names(ss.result)=='tempword')] <- data.cat

    if(verbose){
        print(summary(ss.result))
    }

    if (crop){
        ss.result <- crop.result(data, ss.result, data.cat)
    }

    #ss.result
    list(model=ss.model, result=ss.result)
}

replot.tongue.ss <- function(what_to_plot, lwd=3, main=NULL, CI.fill=FALSE, printing=FALSE, show.legend=T, xlab='X', ylab='Y', xlim=NULL, ylim=NULL, 
                             overplot=FALSE, Fit.palette=NULL, CI.palette=NULL, color.alpha=0.25, Fit.v=0.75, CI.v=0.75, ltys=1:100, axes=T, palate=NULL, 
                             x.axis.step=NULL, y.axis.step=NULL,  legend.pos='bottomright'){

    if (is.null(main)){
        main <- what_to_plot$main
    }
    if (is.null(xlim)){
        xlim <- what_to_plot$xlim
    }
    if (is.null(ylim)){
        ylim <- what_to_plot$ylim
    }
    if (is.null(x.axis.step)&is.null(y.axis.step)){
        axes=TRUE
    }else{
        axes=FALSE
    }
    
    plot.tongue.ss(ss.result=what_to_plot$pol.cart, data.cat=names(what_to_plot$pol.cart)[2], lwd=lwd, main=main, CI.fill=CI.fill, printing=printing, show.legend=show.legend, 
                   plot.labels=c(main,xlab,ylab), overplot=overplot, xlim=xlim, ylim=ylim, Fit.palette=Fit.palette, CI.palette=CI.palette, 
                   color.alpha=color.alpha, Fit.v=Fit.v, CI.v=CI.v, ltys=ltys, axes=axes, origin=what_to_plot$origin, palate=palate, legend.pos=legend.pos)
    if(!axes){
        if (is.null(x.axis.step)){
            axis(1)
        }else{
            axis(1, seq(ceil(xlim[1]),floor(xlim[2]), by=x.axis.step))
        }
        if (is.null(y.axis.step)){
            axis(2)
        }else{
            axis(2, seq(ceil(ylim[1]),floor(ylim[2]), by=y.axis.step))
        }  
        box()
    }
}

#PLOT THE SSANOVA RESULTS
plot.tongue.ss <- function(ss.result, data.cat='word', lwd=3, main='', CI.fill=FALSE, printing=FALSE, show.legend=T, 
                           plot.labels=c(main,'X','Y'), overplot=FALSE, xlim=NULL, ylim=NULL, Fit.palette=NULL, CI.palette=NULL, 
                           color.alpha=0.25, Fit.v=0.75, CI.v=0.75, ltys=1:100, axes=T, origin=NULL, palate=NULL, 
                           legend.pos='bottomright'){  

    ltys=1:100
    n_categories <- length(levels(ss.result[,data.cat]))
    # print (paste('n_categories', n_categories))
    if (is.null(Fit.palette)){
        Fit.palette <- rainbow(n_categories, v=Fit.v)
    }
    if (is.null(CI.palette)){
        CI.palette <- rainbow(n_categories, alpha=color.alpha, v=CI.v)
    }
    # print (Fit.palette)
    xrange = range(c(ss.result$X, ss.result$ss.lower.CI.X, ss.result$ss.upper.CI.X, origin[1], palate$X))
    yrange = range(c(ss.result$ss.Fit, ss.result$ss.lower.CI.Y, ss.result$ss.upper.CI.Y, origin[2], palate$Y))

    # xrange = range(c(ss.result$X, ss.result$ss.lower.CI.X, ss.result$ss.upper.CI.X))
    # yrange = range(c(ss.result$ss.Fit, ss.result$ss.lower.CI.Y, ss.result$ss.upper.CI.Y))

    # if (!is.null(origin)){
    #     xrange <- range(c(xrange, origin[1]))
    #     yrange <- range(c(yrange, origin[2]))
    # }

    if (is.null(xlim)){
        xlim <- xrange
    }
    if (is.null(ylim)){
        ylim <- yrange
    }
    if (!overplot){
        if (axes){
            main <- plot.labels[1]
        # }else{
        #     main <- ''
        }
        plot(0, 0, xlim=xlim, ylim=ylim,xlab=plot.labels[2],ylab=plot.labels[3], main=main, type='n', axes=axes)
    }

    if (!is.null(palate)){
        points(palate$X,palate$Y, type='l')
    }

    if (printing){
        for (i in 1:n_categories){
            w=levels(ss.result[,data.cat])[i]
            subdata <- ss.result[ss.result[,data.cat]==w,]
            if (CI.fill==TRUE){
                polygon(c(subdata$ss.upper.CI.X, rev(subdata$ss.lower.CI.X)),
                        c(subdata$ss.upper.CI.Y, rev(subdata$ss.lower.CI.Y)),
                        col=CI.palette[i], border=F)
            }else{
                lines(subdata$ss.upper.CI.X, subdata$ss.upper.CI.Y, type='l', col=Fit.palette[i], lty=3)
                lines(subdata$ss.lower.CI.X, subdata$ss.lower.CI.Y, type='l', col=Fit.palette[i], lty=3)
            }
            lines(subdata$X, subdata$ss.Fit, type='l', col=Fit.palette[i], lwd=lwd, lty=ltys[i])
            }
        if (show.legend){
            # if (CI.fill){
            #     lege<-legend(xlim[1]+0.8*diff(ylim), ylim[1]+0.3*diff(ylim), c(levels(ss.result[,data.cat])), lwd=lwd, col='white', lty=ltys[1:n_categories])
            #     for (i in 1:length(Fit.palette)){
            #         lines(lege$rect$left + lege$rect$w * c(0.17,0.6), rep(lege$text$y[i], 2), col=CI.palette[i], lend='butt', lwd=12)
            #     }
            #     for (i in 1:length(Fit.palette)){
            #         lines(lege$rect$left + lege$rect$w * c(0.17,0.6), rep(lege$text$y[i], 2), col=Fit.palette[i], lwd=lwd, lty=ltys[i])
            #     }
            # }else{
                # legend(xlim[1]+0.8*diff(ylim), ylim[1]+0.3*diff(ylim), c(levels(ss.result[,data.cat])), lwd=lwd, col=Fit.palette, lty=ltys[1:n_categories])
                legend(legend.pos, c(levels(ss.result[,data.cat])[levels(ss.result[,data.cat])%in%ss.result[,data.cat]]), 
                    lwd=lwd, col=Fit.palette[levels(ss.result[,data.cat])%in%ss.result[,data.cat]], 
                    lty=ltys[1:n_categories][levels(ss.result[,data.cat])%in%ss.result[,data.cat]])
            # }


        }
    }else{
        # print ('PLOTTING')
        # print (Fit.palette)
        for (i in 1:n_categories){
            w=levels(ss.result[,data.cat])[i]
            subdata <- ss.result[ss.result[,data.cat]==w,]
            if (CI.fill==TRUE){
                polygon(c(subdata$ss.upper.CI.X, rev(subdata$ss.lower.CI.X)),
                        c(subdata$ss.upper.CI.Y, rev(subdata$ss.lower.CI.Y)),
                        col=CI.palette[i], border=F)
                }else{
                lines(subdata$ss.upper.CI.X, subdata$ss.upper.CI.Y, type='l', col=Fit.palette[i], lty=3)
                lines(subdata$ss.lower.CI.X, subdata$ss.lower.CI.Y, type='l', col=Fit.palette[i], lty=3)
                }
            lines(subdata$X, subdata$ss.Fit, type='l', col=Fit.palette[i], lwd=lwd, lty=1)
            }
        if (show.legend){
            legend(legend.pos, c(levels(ss.result[,data.cat])[levels(ss.result[,data.cat])%in%ss.result[,data.cat]]), 
                lwd=lwd, lty=1, col=Fit.palette[levels(ss.result[,data.cat])%in%ss.result[,data.cat]])
        }
    }
}

guess.data.cat <- function(data, data.cat){
    
}


reshow.traces <- function(what_to_plot, to.highlight=c(''), to.plot=c(''), token.label='token', flip=TRUE, main=NULL, 
                            overplot=FALSE, is.polar=FALSE, origin=c(0,0), show.legend=TRUE, xlim=NULL, ylim=NULL,
                            trace.palette=NULL, ghost.palette=NULL, color.alpha=0.1, color.v=0.75, palate=NULL, xlab='X', ylab='Y', 
                            x.axis.step=NULL, y.axis.step=NULL, legend.pos='bottomright', drop_levels=FALSE){

    if (is.null(main)){
        main <- what_to_plot$main
    }
    if (is.null(xlim)){
        xlim <- what_to_plot$xlim
    }
    if (is.null(ylim)){
        ylim <- what_to_plot$ylim
    }
    if (is.null(x.axis.step)&is.null(y.axis.step)){
        axes=TRUE
    }else{
        axes=FALSE
    }

    show.traces(data=what_to_plot$data, data.cat=names(what_to_plot$pol.cart)[2], to.highlight=to.highlight, to.plot=to.plot, token.label=token.label, flip=flip,
                main=main, overplot=overplot, is.polar=is.polar, origin=what_to_plot$origin, show.legend=show.legend, 
                xlim=xlim, ylim=ylim, trace.palette=trace.palette, ghost.palette=ghost.palette, color.alpha=color.alpha, color.v=color.v, palate=palate,
                xlab=xlab, ylab=ylab, axes=axes, legend.pos=legend.pos, drop_levels=drop_levels)
    if(!axes){
        if (is.null(x.axis.step)){
            axis(1)
        }else{
            axis(1, seq(ceil(xlim[1]),floor(xlim[2]), by=x.axis.step))
        }
        if (is.null(y.axis.step)){
            axis(2)
        }else{
            axis(2, seq(ceil(ylim[1]),floor(ylim[2]), by=y.axis.step))
        }  
        box()
    }
}


interpolate_polar_data = function(data.polar, interpolate){
    data.polar.interpolated = c()        
    for (tid in unique(data.polar$token_id)){
        subdata = subset(data.polar, token_id==tid)
        rownames(subdata) = NULL
        subdata_interp_XY <- spline(subdata$X, subdata$Y, xout=seq(subdata$X[1], subdata$X[nrow(subdata)], length.out=interpolate))
        interp_token = subset(subdata, select=-c(X,Y))[1,]
        interp_token = suppressWarnings(cbind(interp_token, data.frame(X=subdata_interp_XY$x, Y=subdata_interp_XY$y)))
        data.polar.interpolated = rbind(data.polar.interpolated, interp_token)
    }
    data.polar.interpolated
}

#PLOT THE ORIGINAL DATA
show.traces <- function(data, data.cat='word', to.highlight=c(''), to.plot=c(''), token.label='token', flip=TRUE, main='', 
                        overplot=FALSE, is.polar=FALSE, origin=c(0,0), show.legend=TRUE, xlim=NULL, ylim=NULL,
                        trace.palette=NULL, ghost.palette=NULL, color.alpha=0.1, color.v=0.75, palate=NULL, xlab='X', ylab='Y', 
                        axes=TRUE, legend.pos='bottomright', drop_levels=FALSE, interpolate=NULL){ 
    if (sum(!names(data)%in%c('token','X','Y'))==1 & !data.cat%in%names(data)){
        data.cat <- names(data)[!names(data)%in%c('token','X','Y')]
        warning(paste('Using column \"',data.cat,'" to group the data.\nTo avoid this warning, use "show.traces(data, \'',data.cat,'\')"',sep=''))
    }
    #print(data.cat)
    show.cat <- function(data, data.cat, w, col){
        subdata <- data[data[,data.cat]==w,]
        subdata[,token.label] <- factor(subdata[,token.label])
        tokens <- levels(subdata[,token.label])
        for (t in tokens){
            token <- subdata[subdata[,token.label]==t,]
            lines(token$X,token$Y,col=col)
        }
    }
    if (flip){
        data$Y <- -data$Y
        if(!is.null(palate)){
            palate$Y <- -palate$Y
        }
    }
    if (is.polar){
        data[,c('X','Y')] <- make.cartesian(data[,c('X','Y')], origin=origin, flip=flip)
    }


    if (!is.null(interpolate)){
        data.polar = data
        data.polar[,c('X','Y')] <- make.polar(data.polar[,c('X','Y')], origin=origin, flip=flip)
        data = interpolate_polar_data(data.polar, interpolate)
        data[,c('X','Y')] <- make.cartesian(data[,c('X','Y')], origin=origin, flip=flip)
        remove(data.polar)
    }




    if (drop_levels){
        data[,data.cat] <- droplevels(data[,data.cat])
    }

    categories <- levels(data[,data.cat])
    n_categories <- length(categories)
    if (is.null(trace.palette)){
        trace.palette <- rainbow(n_categories, v=color.v)
    }
    if (is.null(ghost.palette)){
        ghost.palette <- rainbow(n_categories, alpha=color.alpha, v=color.v)
    }

    if (overplot==FALSE){
        if (is.null(xlim)){
            xlim <- range(rbind(data$X, palate))
        }
        if (is.null(ylim)){
            ylim <- range(rbind(data$Y, palate))
        }
        plot(0,0, type='n', xlim=xlim, ylim=ylim,xlab=xlab,ylab=ylab, main=main, axes=axes)
    }

    if (!is.null(palate)){
        points(palate$X,palate$Y, type='l')
    }

    for (i in 1:n_categories){
        w=levels(data[,data.cat])[i]
        if (w%in%to.plot >= mean(categories%in%to.plot)){
            if (w%in%to.highlight >= mean(categories%in%to.highlight)){
                show.cat(data, data.cat, w, col=trace.palette[i])
            }else{
                show.cat(data, data.cat, w, col=ghost.palette[i])
            }
        }
    }
    if(show.legend){
        legend(legend.pos, categories[categories%in%data[,data.cat]], lwd=1, col=trace.palette[categories%in%data[,data.cat]])
    }
}

#CALCULATE AN SSANOVA IN POLAR COORDINATES AND THEN PLOT IT BACK IN CARTESIAN COORDINATES
polar.ssanova <- function(data, data.cat='word', scale=1, origin.method='xmean_ymin', debug=FALSE, plotting=TRUE, main='', 
                          CI.fill=FALSE, printing=FALSE, flip=TRUE, cartesian.only=FALSE, is.polar=FALSE, show.legend=TRUE, 
                          plot.labels=c(main,'X','Y'), overplot=FALSE, xlim=NULL, ylim=NULL, lwd=3, alpha=1.4, crop=FALSE,
                          Fit.palette=NULL, CI.palette=NULL, color.alpha=0.25, Fit.v=0.75, CI.v=0.75, ltys=1:100, axes=T,
                          origin=NULL, palate=NULL, token.label='token', length.out=1000, SDs=1.96, legend.pos='bottomright',
                          drop_levels=FALSE, interpolate=NULL){
    #origin <- c(NULL, NULL)
    if (sum(!names(data)%in%c(token.label,'X','Y'))==1 & !data.cat%in%names(data)){
        data.cat <- names(data)[!names(data)%in%c(token.label,'X','Y')]
        warning(paste('Using column \"',data.cat,'" to group the data.\nTo avoid this warning, use "polar.ssanova(data, \'',data.cat,'\')"',sep=''))
    }
    if (flip==TRUE){
    #    data$Y <- -data$Y
        if(!is.null(palate)){
            palate$Y <- -palate$Y
        }
        if (is.null(ylim)){
            ylim = range(-data$Y, na.rm=T)
        }

    }else if (is.null(ylim)){
        ylim = range(data$Y, na.rm=T)
    }

    if (is.null(xlim)){
        xlim = range(data$X, na.rm=T)
    }

    data.scaled <- us.rescale(data, scale)
    if (cartesian.only){
        tongue.ss.return <- tongue.ss(data.scaled, data.cat=data.cat, flip=flip, alpha=alpha, crop=crop, length.out=length.out, SDs=SDs, drop_levels=drop_levels)
        ss.pol.cart <- tongue.ss.return$result
        ss.cart <- ss.pol.cart
        ss.polar <- ss.pol.cart
    }else{
        if (is.polar){
            #origin <- select.origin(data.scaled$X, data.scaled$Y, method=origin.method)
            origin <- c(0,0)
            print (origin)
            data.polar <- data.scaled
        }else{
            if (is.null(origin)){
                origin <- select.origin(data.scaled$X, data.scaled$Y, data.scaled[,token.label], method=origin.method, flip=flip)
            }
            print(paste('origin is',paste(origin, collapse=', ')))
            # print(summary(data.scaled$Y))
            data.polar <- make.polar(data.scaled, origin, flip)
            # print('DATA POLAR')
            # print(head(data.polar,15))
        }

        if (!is.null(interpolate)){
            data.polar = interpolate_polar_data(data.polar, interpolate)
        }

        tongue.ss.return <- tongue.ss(data.polar, data.cat=data.cat, alpha=alpha, crop=crop, length.out=length.out, SDs=SDs, drop_levels=drop_levels)

        ss.polar <- tongue.ss.return$result
        ss.cart <- c()

        ss.pol.cart <- ss.polar
        ss.pol.cart[,c('X','ss.Fit')] <- make.cartesian(ss.polar[,c('X','ss.Fit')], origin=origin, flip=flip)
        ss.pol.cart[,c('ss.cart.SE')] <- NA
        ss.pol.cart[,c('ss.upper.CI.X','ss.upper.CI.Y')] <- make.cartesian(ss.polar[,c('ss.upper.CI.X','ss.upper.CI.Y')], origin=origin, flip=flip)
        ss.pol.cart[,c('ss.lower.CI.X','ss.lower.CI.Y')] <- make.cartesian(ss.polar[,c('ss.lower.CI.X','ss.lower.CI.Y')], origin=origin, flip=flip)
    }

    if (plotting){
        if (!is.null(origin)){
            origin <- c(origin[1], ifelse(flip, -origin[2], origin[2]))
        }
        if (debug){
            tongue.ss.return <- tongue.ss(data.scaled, data.cat=data.cat, flip=T, crop=crop, length.out=length.out, SDs=SDs, drop_levels=drop_levels)
            ss.cart <- tongue.ss.return$result
            plot.tongue.ss(ss.cart, data.cat, main=main, CI.fill=CI.fill, printing=printing, show.legend=show.legend, 
                           plot.labels=plot.labels, overplot=overplot, xlim=xlim, ylim=ylim, lwd=lwd,
                           Fit.palette=Fit.palette, CI.palette=CI.palette, color.alpha=color.alpha, Fit.v=Fit.v, CI.v=CI.v, ltys=ltys,
                           axes=axes, origin=origin, palate=palate, legend.pos=legend.pos)
            plot.tongue.ss(ss.polar, data.cat, main=main, CI.fill=CI.fill, printing=printing, show.legend=show.legend, 
                           plot.labels=plot.labels, overplot=overplot, xlim=xlim, ylim=ylim, lwd=lwd,
                           Fit.palette=Fit.palette, CI.palette=CI.palette, color.alpha=color.alpha, Fit.v=Fit.v, CI.v=CI.v, ltys=ltys,
                           axes=axes, origin=origin, palate=palate, legend.pos=legend.pos)
        }
        plot.tongue.ss(ss.pol.cart, data.cat, main=main, CI.fill=CI.fill, printing=printing, show.legend=show.legend, 
                       plot.labels=plot.labels, overplot=overplot, xlim=xlim, ylim=ylim, lwd=lwd,
                       Fit.palette=Fit.palette, CI.palette=CI.palette, color.alpha=color.alpha, Fit.v=Fit.v, CI.v=CI.v, ltys=ltys,
                           axes=axes, origin=origin, palate=palate, legend.pos=legend.pos)
    }
    #return ss.pol.cart
    list(polar=ss.polar, ss.cart=ss.cart, pol.cart=ss.pol.cart, model=tongue.ss.return$model, 
         data=data, origin=origin, xlim=xlim, ylim=ylim, main=main)
}

#CALCULATE AN SSANOVA IN CARTESIAN COORDINATES (NOT ADVISED FOR ULTRASOUND DATA)
cart.ssanova <- function(data, data.cat='word', scale=1, origin.method='xmean_ymin', debug=FALSE, plotting=TRUE, main='', 
                         CI.fill=FALSE, printing=FALSE, flip=TRUE, show.legend=TRUE, plot.labels=c(main,'X','Y'), 
                         overplot=FALSE, xlim=NULL, ylim=NULL, lwd=3, alpha=1.4, crop=FALSE, 
                         Fit.palette=NULL, CI.palette=NULL, color.alpha=0.25, Fit.v=0.75, CI.v=0.75, ltys=1:100, axes=T, 
                         length.out=1000, SDs=1.96, legend.pos='bottomright', drop_levels=FALSE){
    polar.ssanova(data=data, data.cat=data.cat, scale=scale, origin.method=origin.method, debug=debug, plotting=plotting, main=main,
                  CI.fill=CI.fill, printing=printing, flip=flip, cartesian.only=TRUE, show.legend=show.legend, plot.labels=plot.labels, 
                  overplot=overplot, xlim=xlim, ylim=ylim, lwd=lwd, alpha=alpha, crop=crop,
                  Fit.palette=Fit.palette, CI.palette=CI.palette, color.alpha=color.alpha, Fit.v=Fit.v, CI.v=CI.v, ltys=ltys,
                  axes=axes, length.out=length.out, SDs=SDs, legend.pos=legend.pos, drop_levels=drop_levels)
}


rotateXY <- function(data, angle, center=NULL){
    #rotate the data so that the occlusal plane is horizontal
    M <- cbind(data$X,data$Y)
    alpha <- -angle*pi/180
    rotm <- matrix(c(cos(alpha),sin(alpha),-sin(alpha),cos(alpha)),ncol=2)
    if (is.null(center)){
        M2 <- t(rotm %*% ( t(M)-c(M[1,1],M[1,2]))+c(M[1,1],M[1,2]))
    }else{
        M2 <- t(rotm %*% ( t(M)-center)+center)
    }

    data[,c('X','Y')] <- M2
    data
}

separate_token <- function(data, split_into=c('repetition', 'phone_number', 'token_number', 'phone', 'time'), sep='_', remove=FALSE){
    require(tidyr)
    data <- separate(data, 'token', split_into, sep=sep, remove=remove)
    for (column_name in split_into){
        if (column_name=='phone'){
            data$phone <- as.factor(data$phone)
        }else{
            data[,column_name] <- as.numeric(data[,column_name])
        }
    }
    data
}

word_from_filename <- function(filename){
    number_index <- regexpr('[0-9]', filename)[[1]]
    substr(filename, 1, number_index-1)
}

token_from_filename <- function(filename){
    number_index <- regexpr('[0-9]', filename)[[1]]
    gsub('.jpg', '', substr(filename, number_index, nchar(filename)))
}

separate_filename <- function(data){
    words <- as.factor(unlist(lapply(paste(data$Filename), word_from_filename)))
    tokens <- as.factor(unlist(lapply(paste(data$Filename), token_from_filename)))
    separate_token(data.frame(word=words, token=tokens, data))
}

rotateCounterclockwise <- function(data, angle, origin=c(0,0)){
    M <- cbind(data$X, data$Y)
    alpha <- -angle*pi/180
    rotm <- matrix(c(cos(alpha), sin(alpha),-sin(alpha), cos(alpha)), ncol=2)
    M2 <- t(rotm %*% ( t(M) - origin ) + origin )
    data.frame(X=M2[,1], Y=M2[,2])
}

contiguousequal <- function(data, value, position) {
    if(data[position] != value)
        return(rep(FALSE, length(data)))
    id <- cumsum(c(1, as.numeric(diff(data) != 0)))
    id == id[position]
}



difference_plot <- function(ss, level1, level2, main=NULL, flip=TRUE){

    require(pracma)
    require(geosphere)

    if (is.null(main)){
        main <- paste0('/', level1, '/ vs. /', level2, '/')
    }
    first_fit <- ss$polar[ss$polar[,2]==level1,]
    second_fit <- ss$polar[ss$polar[,2]==level2,]

    first_fit_cart <- ss$pol.cart[ss$pol.cart[,2]==level1,]
    second_fit_cart <- ss$pol.cart[ss$pol.cart[,2]==level2,]

    fitXrange <- c(max(min(first_fit$X),min(second_fit$X)), min(max(first_fit$X),max(second_fit$X)))

    first_fit_subset <- first_fit$X>=fitXrange[1] & first_fit$X<=fitXrange[2]
    second_fit_subset <- second_fit$X>=fitXrange[1] & second_fit$X<=fitXrange[2]

    xout <- seq(fitXrange[1], fitXrange[2], length.out=100)

    first_interp <- approx(first_fit[first_fit_subset,'X'],
                           first_fit[first_fit_subset,'ss.Fit'], xout=xout)
    second_interp <- approx(second_fit[second_fit_subset,'X'],
                            second_fit[second_fit_subset,'ss.Fit'], xout=xout)

    if (flip==TRUE){
        first_interp_cart <- make.cartesian(data.frame(X=first_interp$x, Y=first_interp$y), origin=c(ss$origin[1],-ss$origin[2]), flip=flip)  
        second_interp_cart <- make.cartesian(data.frame(X=second_interp$x, Y=second_interp$y), origin=c(ss$origin[1],-ss$origin[2]), flip=flip)
    }else{
        first_interp_cart <- make.cartesian(data.frame(X=first_interp$x, Y=first_interp$y), origin=c(ss$origin[1],ss$origin[2]), flip=flip)  
        second_interp_cart <- make.cartesian(data.frame(X=second_interp$x, Y=second_interp$y), origin=c(ss$origin[1],ss$origin[2]), flip=flip)
    }

    first_interp_cart <- list(x=first_interp_cart[,1], y=first_interp_cart[,2])
    second_interp_cart <- list(x=second_interp_cart[,1], y=second_interp_cart[,2])

    interp_diff = second_interp$y - first_interp$y
    second_higher <- interp_diff > 0

    polygon_rle = rle(second_higher)
    polygon_starts = c(1,1+polygon_rle$lengths[1:length(polygon_rle$lengths)-1])
    polygon_ends = cumsum(polygon_rle$lengths)

    which_span = rep(NA,100)
    span_mass = c()
    for(i in 1:length(polygon_starts)){
        which_span[polygon_starts[i]:polygon_ends[i]] <- i
        span_mass = c(span_mass, sum(interp_diff[polygon_starts[i]:polygon_ends[i]]))
    }

    # polygon_boundaries <- which(diff(second_higher)!=0)
    # polygon_starts <- c(1, polygon_boundaries+1)
    # polygon_ends <- c(polygon_boundaries, length(second_higher))

    largest_positive_displacement = which(span_mass==max(span_mass))
    first_fit_polygon_subset <- first_fit$X>=first_interp$x[polygon_starts[largest_positive_displacement]] & first_fit$X<=first_interp$x[polygon_ends[largest_positive_displacement]]
    second_fit_polygon_subset <- second_fit$X>=first_interp$x[polygon_starts[largest_positive_displacement]] & second_fit$X<=first_interp$x[polygon_ends[largest_positive_displacement]]
    fit_diff = second_fit[second_fit_polygon_subset,'ss.Fit'] - first_fit[first_fit_polygon_subset,'ss.Fit']
    fit_angle = first_fit[first_fit_polygon_subset,'X']
    mean_diff_angle = sum(fit_diff/sum(fit_diff)*fit_angle)
    mean_diff_dist = approx(first_interp$x, interp_diff,mean_diff_angle)$y
    arrow_start_angle = approx(first_interp$x, xout=mean_diff_angle)$x
    arrow_start_radius = approx(first_interp$x, first_interp$y, xout=mean_diff_angle)$y

    deltaX = - ( mean(second_fit_cart$X) - mean(first_fit_cart$X) )
    deltaY = mean(second_fit_cart$ss.Fit) - mean(first_fit_cart$ss.Fit)
     
    deltaX = - ( mean(second_interp_cart$x) - mean(first_interp_cart$x) )
    deltaY = mean(second_interp_cart$y) - mean(first_interp_cart$y)
     


    first_fit_polygon <- rbind(first_fit_cart[,c('X','ss.Fit')], 
                               c(max(ss$pol.cart$X), min(ss$pol.cart$ss.Fit)), 
                               c(min(ss$pol.cart$X), min(ss$pol.cart$ss.Fit)))

    second_fit_polygon <- rbind(second_fit_cart[,c('X','ss.Fit')], 
                               c(max(ss$pol.cart$X), min(ss$pol.cart$ss.Fit)), 
                               c(min(ss$pol.cart$X), min(ss$pol.cart$ss.Fit)))

    first_fit_polygon[,2] <- first_fit_polygon[,2]-mean(ss$pol.cart$ss.Fit)
    second_fit_polygon[,2] <- second_fit_polygon[,2]-mean(ss$pol.cart$ss.Fit)

    centroid_deltaXY = centroid(second_fit_polygon) - centroid(first_fit_polygon)
    centroid_deltaXY[1] = -centroid_deltaXY[1]

    contour_displacement_angle = (atan2(deltaY,deltaX))
    polygon_displacement_angle = (atan2(centroid_deltaXY[2],centroid_deltaXY[1]))

    # polygon(first_fit_polygon[,1], first_fit_polygon[,2], border='blue', col='gray')
    # print(first_fit_polygon)

    if (flip==TRUE){
        arrow_start <- make.cartesian(data.frame(X=arrow_start_angle, Y=arrow_start_radius), origin=c(ss$origin[1],-ss$origin[2]), flip=flip)  
    }else{
        arrow_start <- make.cartesian(data.frame(X=arrow_start_angle, Y=arrow_start_radius), origin=c(ss$origin[1],ss$origin[2]), flip=flip)  
    }


    main = paste(main,'\n',round(mean_diff_dist*cos(mean_diff_angle), 1), round(mean_diff_dist*sin(mean_diff_angle),1), round(mean_diff_angle,2),
                           round(deltaX,1), round(deltaY,1), round(contour_displacement_angle,2),
                           round(centroid_deltaXY[1],1), round(centroid_deltaXY[2],1), round(polygon_displacement_angle,2))

    plot(0,0,type='n', xlim=ss$xlim, ylim=ss$ylim, main=main, xlab='X', ylab='Y')

    total_area <- 0

    # print (polygon_starts)
    for (i in 1:length(polygon_starts)){
        polygon_x <- c(first_interp_cart$x[polygon_starts[i]:polygon_ends[i]], rev(second_interp_cart$x[polygon_starts[i]:polygon_ends[i]]))
        polygon_y <- c(first_interp_cart$y[polygon_starts[i]:polygon_ends[i]], rev(second_interp_cart$y[polygon_starts[i]:polygon_ends[i]]))
        polygon(polygon_x, polygon_y, border=NA, col=ifelse(second_higher[polygon_starts[i]],rgb(1,0,0,a=0.25),rgb(0,0,1,a=0.25)))
        total_area <- total_area + abs(polyarea(polygon_x, polygon_y))
    }
    

    # print(c(arrow_start_angle, arrow_start_radius, arrow_start))
    # angle_to_second_interp <- atan2(second_interp_cart$y-first_interp_cart$y, second_interp_cart$x-first_interp_cart$x)
    # dist_to_second_interp <- sqrt((second_interp_cart$y-first_interp_cart$y)^2+(second_interp_cart$x-first_interp_cart$x)^2)

    # angle_to_second_interp_mean <- atan2(mean(second_interp_cart$y-first_interp_cart$y), mean(second_interp_cart$x-first_interp_cart$x))
    # dist_to_second_interp_mean <- sqrt((mean(second_interp_cart$y-first_interp_cart$y))^2+(mean(second_interp_cart$x-first_interp_cart$x))^2)


    points(first_fit_cart[,c('X','ss.Fit')],type='l',col='black')
    points(second_fit_cart[,c('X','ss.Fit')],type='l',col=rgb(1,0,0))

    arrows(arrow_start[1], arrow_start[2], 
           arrow_start[1]-mean_diff_dist*cos(mean_diff_angle), 
           arrow_start[2]+mean_diff_dist*sin(mean_diff_angle), length=0.1, lwd=3, col='red')

    arrows(ss$origin[1], ss$origin[2], 
           ss$origin[1]-3*deltaX, 
           ss$origin[2]+3*deltaY, length=0.1, lwd=3, col='black')

    arrows(ss$origin[1], ss$origin[2], 
           ss$origin[1]-3*centroid_deltaXY[1], 
           ss$origin[2]+3*centroid_deltaXY[2], length=0.1, lwd=3, col='blue')


    points(first_interp_cart,type='l',col='black',lwd=3)
    points(second_interp_cart,type='l',col='red',lwd=3)

    legend('bottomright', c(level1, level2), lwd=3, col=c('black','red'))
    list(total_area=total_area, displacement_angle=mean_diff_angle, displacement_distance=mean_diff_dist,
        contour_displacement_angle=contour_displacement_angle,
        contour_displacement_magnitude=sqrt(deltaY^2+deltaX^2),
        polygon_displacement_angle=polygon_displacement_angle,
        polygon_displacement_magnitude=sqrt(centroid_deltaXY[2]^2+centroid_deltaXY[1]^2))

}



# require(fdasrvf)

# x = curve_pair_align(t(data.frame(first_interp_cart$x, first_interp_cart$y)),t(data.frame(second_interp_cart$x, second_interp_cart$y)))

# t(data.frame(first_interp_cart$x, first_interp_cart$y))

# library(elasdics)

# library(elasdics)

# xx = align_curves(data.frame(first_interp_cart), data.frame(second_interp_cart))


# xx$data_curve2_aligned[,c('x','y')] - data.frame(first_interp_cart)

