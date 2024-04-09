#By convention, the variable contribution plot has a circle around the variables that has a radius of 1. Here’s some code to make one.
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
}


##Combining correlogram with the significance test
## Computing the p-value of correlations
## To compute the matrix of p-value, a custom R function is used 
# ... : further arguments to pass to the native R cor.test function
#tutorial: http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram
# mat : is a matrix of data

cor.mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(mat[, i], mat[, j], ...)
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
}


# create a function that is the opposite of intersect
outersect <- function(x, y) {
    sort(c(setdiff(x, y),
           setdiff(y, x)))
}