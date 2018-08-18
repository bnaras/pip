#' Compute the prediction statistics
#'
#' Compute three-day predicted usage, waste, remaining and shortage
#' using the actual usage and predicted usage.
#'
#' @param y is the number of units used at the current day i
#' @param t_pred is the sum of the predicted number of units for days
#'     i+1, i+2, i+3
#' @param initial_expiry_data the number of units that can be used the
#'     current day and the next (2-length vector)
#' @param initial_collection_data is the number of units to collect
#'     for days i, i+1, i+2
#' @param start is when we start the clock. So for example, in this
#'     default invocation the shelf storage (initial_collection_data)
#'     was 60, 60, 60, on days 10, 11, and 12. So the evaluation of
#'     the model begins at day 11 but the decision to order collection
#'     starts on day 13.
#' @return a list with four components, \code{x} is the number of
#'     units to collect, \code{r} is a matrix of two columns
#'     indicating remaining units usable on day i+1 and day i+2,
#'     \code{w} is waste, and \code{s} is shortage.
#' @export
compute_prediction_statistics <- function(y, #actual usage arranged according to date
                                          t_pred,
                                          initial_expiry_data = c(0, 0),
                                          initial_collection_data = c(60, 60, 60),
                                          start = 10) {
    N <- length(y)## remove day column
    r <- matrix(0, nrow = N , ncol = 2)
    x <- numeric(N + 3)
    w <- s <- numeric(N)
    x[seq.int(start, start + 2)] <- initial_collection_data
    for (i in seq.int(start, N)) {
        if(i == start){
            w[i] <- pos(initial_expiry_data[1] - y[i])
            r[i, 1] <- pos(initial_expiry_data[1] + initial_expiry_data[2] - y[i] - w[i])
            s[i] <- pos(y[i] - initial_expiry_data[1] - initial_expiry_data[ 2] - x[i])
            r[i, 2] <- pos(x[i] - pos(y[i] - initial_expiry_data[1] - initial_expiry_data[2]))
        }else{
            w[i] <- pos(r[i - 1 , 1] - y[i])
            r[i, 1] <- pos(r[i - 1, 1] + r[i - 1, 2] - y[i] - w[i])
            s[i] <- pos(y[i] - r[i - 1, 1] - r[i - 1, 2] - x[i])
            r[i, 2] <- pos(x[i] - pos(y[i] - r[i - 1, 1] - r[i - 1, 2]))
        }
        x[i + 3] <- floor(pos(t_pred[i] - x[i + 1] - x[i + 2] - r[i, 1] - r[i, 2] + 1))
    }
    return(list(x = x, r = r, w = w, s = s))
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#' Run a fit on a given dataset
#'
#' The dataset should have one more column than the predictors and
#' response. This additional column, which is typically a date or day
#' number, has to be named "date" or "day"
#' @param c0 the c0 value
#' @param history_window Number of days to look back
#' @param penalty_factor penalty for shortage specified by doctors
#' @param initial_expiry_data the number of units expiring in an day
#'     and the day after, a 2-length vector
#' @param initial_collection_data is the number of units that will be
#'     collected for the first three days when the prediction begins
#' @param start the day the model is evaluated??. Default 10
#' @param data the dataset
#' @param date_column the name of the date or day number column as a
#'     regex, default is "day|date" i.e. day or date
#' @param response_column the name of the response column, default is
#'     "plt_used"
#' @param show_progress a TRUE/FALSE flag for showing progress,
#'     default TRUE
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @examples
#' build_model(c0 = 30,
#'            history_window = 200,
#'            penalty_factor = 15,
#'            start = 10,
#'            response_column = "plt.used",
#'            data = day3_data,
#'            show_progress = FALSE)
#'
#' @export
#'
build_model <- function(c0, ## the minimum number of units to keep on shelf (c)
                        history_window, ## Number of days to look back
                        penalty_factor = 15, ## the penalty factor for shortage, specified by doctor
                        start = 10,   ## the day you start evaluating??
                        initial_expiry_data = c(0, 0),
                        initial_collection_data = c(60, 60, 60),
                        data, ## The data set
                        date_column = "day|date",  ## we assume column is named date or day by default
                        response_column = "plt_used",
                        show_progress = TRUE) {

    ## Check that we have enough history!
    n <- nrow(data)
    if (history_window > n) {
        stop(sprintf("build_model: history_window(%d) must be less than data rows(%d)",
                     history_window, n))
    }
    d <- data[seq.int(n - history_window + 1, n), ]
    data_column_names <- names(data)
    number_of_columns <- ncol(data)
    resp_var_index <- grep(response_column, data_column_names)
    date_var_index <- grep(date_column, data_column_names)
    if (length(resp_var_index) < 1 || length(date_var_index) < 1) {
        stop("Response and/or date columns not found!")
    }
    pred_var_indices <- setdiff(seq.int(2L, number_of_columns), c(resp_var_index, date_var_index))
    predictor_names <- data_column_names[pred_var_indices]

    nfolds <- 8
    ## Is this always 200? or does it go as far as the history window?
    lambds <- 200 - seq(0, 100) * 2

    N <- history_window
    p <- length(predictor_names)

    foldid <- create_folds(N, nfolds)
    predictions_cv <- w_cv <- r_cv <- matrix(NA, nrow = N, ncol = length(lambds))

    XX <- as.matrix(d[, pred_var_indices])
    y <- d[, resp_var_index]
    #notice that the repsonse in the proceccessed data have been shifted backward by 1 day
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
    #y_today has date matches the date in the data frame
    xMat <- cbind(1, XX)
    #the optimization always assume that we make prediction by the end of the day, say 23:59:59
    pb <- if (show_progress) utils::txtProgressBar(min = 0, max = nfolds, style = 3) else NULL
    for (k in seq_len(nfolds)) {
        coeffs <- single_lpSolve(d, lamb = lambds, ind = which(foldid != k), start = start, c = c0)
        pred1 <- xMat %*% coeffs
        for (l in seq_along(lambds)){
            rr <- compute_prediction_statistics(y,
                                                t_pred = pred1[, l],
                                                initial_expiry_data = c(0, 0),
                                                initial_collection_data = y[seq.int(start + 1L, start + 3L)],
                                                start = start)
            ## The 30 above is the parameter passed, i.e. min number of units to keep on shelf.
            ## So I changed it to c0
            w_cv[foldid == k, l] <- rr$w[foldid == k]
            r_cv[foldid == k, l] <- rr$r[foldid == k, 1] + rr$r[foldid == k, 2]
        }
        if (show_progress) utils::setTxtProgressBar(pb, k)
    }
    if (show_progress) close(pb)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
    ## Since we skip for start days and then we see the results only three days later,
    ## we see waste is start + 5
    first_day_waste_seen <- start + 5 ## 15 is this number

    cv_loss <- apply(w_cv, 2, function(x) sum( x[-(seq_len(first_day_waste_seen))])) +
        apply(r_cv, 2, function(x) sum(((pos(penalty_factor - x))^2) [-(seq_len(first_day_waste_seen))]) )
    cv_loss <- cv_loss + 2 * (length(cv_loss):1)
    index <- which.min(cv_loss)
    coefs <- as.numeric(single_lpSolve(d, lamb = lambds[index], start = start, c = c0))
    names(coefs) <- c("intercept", predictor_names)
    list(
        lambda = lambds[index],
        w = w_cv[ , index],
        r = r_cv[ , index],
        first_day_waste_seen = first_day_waste_seen,
        coefs = coefs
    )
}

#'
#' Predict the next three day sum and units to order if available
#' @param model the trained model
#' @param new_data a new data frame with the same predictor names as in
#'     the training data, even if one row
#' @return the predicted three day total of how many units to collect
#' @export
#'
predict_three_day_sum <- function(model, ## the trained model
                                  new_data) { ## the new dataset with specified number of features
    predictor_names <- names(model$coefs)[-1]
    common_names <- intersect(predictor_names, names(new_data))
    if (length(common_names) != length(predictor_names)) {
        stop("New data set column names don't match training data set!")
    }
    new_data <- as.matrix(new_data[, predictor_names])
    ceiling(sum(new_data %*% model$coefs[-1] ) + model$coefs[1])
}




