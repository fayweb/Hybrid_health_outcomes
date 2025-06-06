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


# testing distributions
# Define function to be used to test, get the log lik and aic
tryDistrib <- function(x, distrib){
    # deals with fitdistr error:
    fit <- 
        tryCatch(MASS::fitdistr(x, distrib), error=function(err) "fit failed")
    return(list(fit = fit,
                loglik = tryCatch(fit$loglik, error=function(err) "no loglik computed"), 
                AIC = tryCatch(fit$aic, error=function(err) "no aic computed")))
}


findGoodDist <- function(x, distribs, distribs2){
    l =lapply(distribs, function(i) tryDistrib(x, i))
    names(l) <- distribs
    print(l)
    listDistr <- lapply(distribs2, function(i){
        if (i %in% "t"){
            fitdistrplus::fitdist(x, i, start = list(df =2))
        } else {
            fitdistrplus::fitdist(x,i)
        }}
    ) 
    par(mfrow=c(2,2))
    denscomp(listDistr, legendtext=distribs2)
    cdfcomp(listDistr, legendtext=distribs2)
    qqcomp(listDistr, legendtext=distribs2)
    ppcomp(listDistr, legendtext=distribs2)
    par(mfrow=c(1,1))
}

# Function to italicize y-axis and legend labels
italics_y <- function(ggplot_object, labels) {
    
    ggplot_object +
        scale_y_discrete(labels = labels) +
        scale_fill_manual(values = color_mapping, labels = labels) +
        theme(
            axis.text.y = element_markdown(),  # Apply markdown to y-axis text
            legend.text = element_markdown()   # Apply markdown to legend text
        )
}


# Function to italicize x-axis and legend labels
italics_x <- function(ggplot_object, labels) {
    
    ggplot_object +
        scale_x_discrete(labels = labels) +
        scale_fill_manual(values = color_mapping, labels = labels) +
        theme(
            axis.text.x = element_markdown(),  # Apply markdown to y-axis text
            legend.text = element_markdown()   # Apply markdown to legend text
        )
}

# Modified function to save each table in its own folder
save_table_all_formats <- function(table_object, table_name, output_dir = paste0(tables)) {
    
    # Create table-specific folder
    table_folder <- file.path(output_dir, table_name)
    dir.create(table_folder, recursive = TRUE, showWarnings = FALSE)
    
    # Create base filename (now inside the table folder)
    base_path <- file.path(table_folder, table_name)
    
    # Save in all formats
    tryCatch({
        # HTML (best for Google Docs)
        gtsave(table_object, paste0(base_path, ".html"))
        cat("✓ Saved", table_name, "as HTML\n")
        
        # Word document
        gtsave(table_object, paste0(base_path, ".docx"))
        cat("✓ Saved", table_name, "as DOCX\n")
        
        # PNG (high quality image)
        gtsave(table_object, paste0(base_path, ".png"), 
               vwidth = 1200, vheight = 800)
        cat("✓ Saved", table_name, "as PNG\n")
        
        # PDF (best for LaTeX)
        gtsave(table_object, paste0(base_path, ".pdf"))
        cat("✓ Saved", table_name, "as PDF\n")
        
        # LaTeX code
        table_object %>%
            as_latex() %>%
            writeLines(paste0(base_path, ".tex"))
        cat("✓ Saved", table_name, "as TEX\n")
        
        cat("✅ All formats saved in folder:", table_folder, "\n\n")
        
    }, error = function(e) {
        cat("❌ Error saving", table_name, ":", e$message, "\n")
    })
}


