#' @importFrom grDevices colorRampPalette
#' @importFrom graphics image points
#' @importFrom utils head
NULL

#' Plots a 2D Density 
#'
#' This function tries to replicate 'smart' 2D density plots 
#' as they are found in flow cytometry data analysis software.
#' These plots commonly work with millions of points. In dense
#' regions, they will plot these as a heatmap; in sparse regions,
#' they will show the actual points.
#'
#' @param x vector containing the X coordinates of the points.
#' @param y vector containing the Y coordinates of the points.
#'  which are passed on to the plotting function so the plot can be
#'  customized.
#' @param nbins how many bins to use for histogram.
#' @param max.points when to start plotting heatmap instead of 
#'  individual points.
#' @param colfunc color mapping for heatmap.
#' @param tf.fun function to apply to histogram values before 
#'  plotting them, \code{log} can be useful.
#' @param pch point character to use for individual points.
#' @param cex point character expansion to use for individual points.
#' @param color.levels how many different colors to use for heatmap.
#' @param ... further arguments to be passed to the plotting function,
#'  such as \code{xlim} or \code{xlab}. 
#' @export
densityplot2d <- function( x, y, nbins=50, max.points=5, 
	colfunc=colorRampPalette(c("gray","darkblue","lightblue", "green", "yellow","orange", "red")), tf.fun=identity, pch=19, cex=0.3, color.levels=20,
	... ){
	x <- as.numeric(x)
	y <- as.numeric(y)
	if( any( is.na(x) ) ){
                stop("NA values in X!")
        }
        if( any( is.na(y) ) ){
		stop("NA values in Y!")
	}
        x.bin <- seq(floor(min(x)-.001), ceiling(max(x)+.001), length=nbins)
	y.bin <- seq(floor(min(y)-.001), ceiling(max(y)+.001), length=nbins)

	x.int <- as.integer(findInterval( x, x.bin ))
        y.int <- as.integer(findInterval( y, y.bin ))

        pix.ry <- length(y.bin) 

        pix <- x.int*pix.ry + y.int
        ptbl <- tabulate( pix )

        my.freq.x <- floor((1:length(ptbl))/pix.ry)
        my.freq.y <- (1:length(ptbl)) %% pix.ry

        my.freq <- cbind(my.freq.x, my.freq.y, ptbl )
        my.freq <- my.freq[my.freq.x > 0 & my.freq.y > 0,]

        freq2D <- diag(nbins-1)*0
        freq2D[cbind(my.freq[,1], my.freq[,2])] <- my.freq[,3]
        f.int <- freq2D[cbind(x.int,y.int)]	
	freq2D[freq2D <= max.points] <- NA

	freq2D <- tf.fun(freq2D)

        image( head(x.bin,-1)+diff(x.bin)/2, head(y.bin,-1)+diff(y.bin)/2, freq2D,
		useRaster=FALSE,
		col=colfunc(color.levels), ... )

        points( x[f.int<=max.points], y[f.int<=max.points],
		pch=pch, cex=cex, col=colfunc(color.levels)[1] )

        #x.points <- c(0,1,2,5,10,20,50,100,200,500,1000)
        #axis( 1, log( 1+x.points, base=2 ), as.character(x.points) )
        #axis( 2, log( 1+x.points, base=2 ), as.character(x.points) )

}

#' @export
confregion2d <- function(M, max.points=1000, alpha=0.05){
	requireNamespace("depth",quietly=TRUE)
	if( nrow(M) > max.points ){
		M <- M[sample(1:nrow(M),max.points,replace=TRUE),]
	}
	d <- apply( M, 1, function(x) depth::depth(x,M) )
	dt <- quantile( d, alpha )
	M <- M[d>dt,]
	h <- chull( M )
	M[c(h,h[1]),]
}

