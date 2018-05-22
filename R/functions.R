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
#' @examples
#' densityplot2d( facsdata[,1], facsdata[,2], tf.fun=sqrt, nbin=100 )
#' pc <- princomp(facsdata)
#' densityplot2d( pc$scores[,1], pc$scores[,2], tf.fun=sqrt, nbin=100 )
#' densityplot2d( pc$scores[,1], pc$scores[,2], z=facsdata[,2],tf.fun=sqrt, nbin=100, colfunc=colorRampPalette(c("black","red")) )
#' @export
densityplot2d <- function( x, y, z=NULL, nbins=50, max.points=5, 
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

	x.args <- list(...)

	if( "xlim" %in% names(x.args) ){
		x.bin <- seq(x.args$xlim[1],x.args$xlim[2],length=nbins+1)
		x.args[["xlim"]] <- NULL
	} else {
	        x.bin <- seq(floor(min(x)-.001), ceiling(max(x)+.001), length=(nbins+1))
	}

	if( "ylim" %in% names(x.args) ){
		y.bin <- seq(x.args$ylim[1],x.args$ylim[2],length=nbins+1)
		x.args[["ylim"]] <- NULL
	} else {
		y.bin <- seq(floor(min(y)-.001), ceiling(max(y)+.001), length=(nbins+1))
	}


	x.int <- as.integer(findInterval( x, x.bin ))
        y.int <- as.integer(findInterval( y, y.bin ))

        pix <- (x.int-1)*(nbins) + (y.int)


	density.matrix <- tabulate(pix)


	density.matrix <- c(density.matrix,rep(0,nbins^2-length(density.matrix)))
	density.matrix <- matrix( density.matrix, nrow=nbins, byrow=TRUE )
       	points.not.in.image <- density.matrix[cbind(x.int,y.int)] <= max.points

	#stopifnot( all( density.matrix[cbind(x.int,y.int)] > 0 ) )

	density.matrix[density.matrix <= max.points] <- NA

	if( !is.null( z ) ){
		z.matrix <- tapply( z, pix, mean )[as.character(1:(nbins^2))]
		z.matrix <- matrix( z.matrix, nrow=nbins, byrow=TRUE )
	}

	do.call(plot,c(list(NA,
		xlim=range(x.bin), 
		ylim=range(y.bin)),x.args))

	M <- tf.fun( density.matrix )
	zlim <- range(M[is.finite(M)])

 	M <- (M - zlim[1L])/diff(zlim)
	col <- colfunc( color.levels )

	if( !is.null(z) ){
		zp <- z[points.not.in.image]
		zlim <- range( z.matrix[is.finite(z.matrix)] )
		z.matrix <- 1-(z.matrix-zlim[1L])/diff(zlim)
       		zi <- floor((color.levels - 1e-05) * z.matrix + 1e-07)
		zi[zi < 0 | zi >= color.levels] <- NA
		zp <- 1-(zp-zlim[1L])/diff(zlim)
		zp <- floor((color.levels - 1e-05) * zp + 1e-07)
		zp[zp < 0 | zp >= color.levels] <- NA
	}

       	Mi <- floor((color.levels - 1e-05) * M + 1e-07)
	Mi[Mi < 0 | Mi >= color.levels] <- NA

	if( is.null(z) ){
		#1D lookup
		Mc <- col[Mi+1L]
	} else {
		#2D lookup
		col <- sapply(rev(col),function(x)
			colorRampPalette(c("white",x))(color.levels+1) )
		col <- col[2:(color.levels+1),,drop=FALSE]
		Mc <- col[cbind(c(Mi+1L),c(zi+1L))]
	}

	dim(Mc) <- dim(M)

	rasterImage( t(Mc)[nrow(Mc):1, ,drop=FALSE],
		min(x.bin), min(y.bin), max(x.bin), max(y.bin), interpolate=FALSE )
	if( is.null( z ) ){
        	points( x[points.not.in.image], y[points.not.in.image],
			pch=pch, cex=cex, col=colfunc(color.levels)[1] )
	} else {
        	points( x[points.not.in.image], y[points.not.in.image],
			pch=pch, cex=cex, col=col[1,zp+1] )
	}


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

