#================
#
# This function fixed the error when running CV using distribution "coxph" in the original 
# gbm package by Greg Ridgeway.
#
# For manual or reference, please check http://cran.r-project.org/web/packages/gbm/index.html
#
# @Author Wei-Yi Cheng
#
#================



gbm.cvrun = function (formula = formula(data), distribution = "bernoulli", 
    data = list(), weights, var.monotone = NULL, n.trees = 100, 
    interaction.depth = 1, n.minobsinnode = 10, shrinkage = 0.001, 
    bag.fraction = 0.5, train.fraction = 1, cv.folds = 0, seed=92, keep.data = TRUE, 
    verbose = TRUE) 
{
    set.seed(seed)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "weights", "offset"), names(mf), 
        0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf$na.action <- na.pass
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    Terms <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    w <- model.weights(mf)
    offset <- model.offset(mf)
    var.names <- attributes(Terms)$term.labels
    x <- model.frame(terms(reformulate(var.names)), data, na.action = na.pass)
    response.name <- as.character(formula[[2]])
    if (is.character(distribution)) 
        distribution <- list(name = distribution)
    cv.error <- NULL
    if (cv.folds > 1) {
        if (distribution$name == "coxph") 
            i.train <- 1:floor(train.fraction * nrow(y))
        else i.train <- 1:floor(train.fraction * length(y))
        cv.group <- sample(rep(1:cv.folds, length = length(i.train)))
        cv.error <- rep(0, n.trees)
        for (i.cv in 1:cv.folds) {
            if (verbose) 
                cat("CV:", i.cv, "\n")
            i <- order(cv.group == i.cv)
	    if (distribution=="coxph"){
		 gbm.obj <- gbm.fit(x[i.train, , drop = FALSE][i,
                , drop = FALSE], y[i.train,][i,], offset = offset[i.train][i],
                distribution = distribution, w = ifelse(w ==
                  NULL, NULL, w[i.train][i]), var.monotone = var.monotone,
                n.trees = n.trees, interaction.depth = interaction.depth,
                n.minobsinnode = n.minobsinnode, shrinkage = shrinkage,
                bag.fraction = bag.fraction, train.fraction = mean(cv.group !=
                  i.cv), keep.data = FALSE, verbose = verbose,
                var.names = var.names, response.name = response.name)
	    } else {
            	gbm.obj <- gbm.fit(x[i.train, , drop = FALSE][i, 
                , drop = FALSE], y[i.train][i], offset = offset[i.train][i], 
                distribution = distribution, w = ifelse(w == 
                  NULL, NULL, w[i.train][i]), var.monotone = var.monotone, 
                n.trees = n.trees, interaction.depth = interaction.depth, 
                n.minobsinnode = n.minobsinnode, shrinkage = shrinkage, 
                bag.fraction = bag.fraction, train.fraction = mean(cv.group != 
                  i.cv), keep.data = FALSE, verbose = verbose, 
                var.names = var.names, response.name = response.name)
	    }
            cv.error <- cv.error + gbm.obj$valid.error * sum(cv.group == 
                i.cv)
        }
        cv.error <- cv.error/length(i.train)
    }
    gbm.obj <- gbm.fit(x, y, offset = offset, distribution = distribution, 
        w = w, var.monotone = var.monotone, n.trees = n.trees, 
        interaction.depth = interaction.depth, n.minobsinnode = n.minobsinnode, 
        shrinkage = shrinkage, bag.fraction = bag.fraction, train.fraction = train.fraction, 
        keep.data = keep.data, verbose = verbose, var.names = var.names, 
        response.name = response.name)
    gbm.obj$Terms <- Terms
    gbm.obj$cv.error <- cv.error
    gbm.obj$cv.folds <- cv.folds
    return(gbm.obj)
}
