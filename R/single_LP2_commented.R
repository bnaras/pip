#This is the function tranining a model with fixed lambda(1)
#1:p: betas
#p+1:p+N:w
#p+N+1:p+2N:r1
#p+2N+1:p+3N:r2
#p+3N+1:p+4N:x
#p+4N+1:p+5N: t
#p+5N+1:intercept
#p+5N+2:p+5N+1+p:l1 bounds
#########


#' Solve an LP problem
#' @param d the data frame
#' @param lamb the lambda value
#' @param start the start date
#' @param c the value for c_0
#' @param buffer a value to be added to account for training in the initial stages, default 10
#' @param ind the indices for cross validation used only if non-null. (Should be sorted!!)
#' @importFrom lpSolveAPI make.lp add.constraint set.bounds lp.control delete.constraint
#' @importFrom lpSolveAPI set.objfn get.variables
#' @export
#'
single_lpSolve <- function(d, lamb, start = 10, c = 30, buffer = 10, ind = NULL) {
    doing_cv <- !is.null(ind)
    ## ind, if specified, is sorted!!! We make use of that below

    XX <- d[, 2:44]
    y <- d[, 45]

    N <- nrow(XX)
    p <- ncol(XX)

    if (!doing_cv) ind <- 1:N

    pind <- p + ind
    ## Remember y is usage, so y[i] is known usage on day i
    ## p betas, N ws, N r1s, N r2s, N xs, N ts, 1 intercept, p bounds
    N_var <- p + 5*N + 1 + p

    my.lp <- lpSolveAPI::make.lp(0,N_var)

    obj_coefficients <- rep(0,N_var)
    ## Set coefs of w to be 1's (eq 9)
    obj_coefficients[pind] <- 1
    lpSolveAPI::set.objfn(my.lp,obj_coefficients)

    ## Constrain first week sum of coefs of beta to be zero??  Why?
    ## Reply: There is non-identifiability issue with the day of week
    ## information. We can add the same value its all coefficients of
    ## these seven features, and shift the intercept correspondingly,
    ## the final model will be equivalent. Wihout loss of generality,
    ## we constrain the week sum of coefs to be 0, it leads to better
    ## intepretability and solves the identifiability issue. That's
    ## why I do not need to remove the first feature. If you have
    ## removed the first feature, then this constraint is no longer
    ## needed. However, later, when we apply L1 penalty, I have left
    ## out the day of week features(7)+average PLT consumption
    ## feature(1)(because we know that then are important and a linear
    ## model with them will not lead to overfitting), we need to
    ## change that accordingly as a result. I have removed the first
    ## two lines.
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[1:7] <- 1
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, type = "=", 0)

    ##w >=0
    for (i in 1:N){
        ## Setting each coeff of w to be non-negative: (eq 13)
        constraint_coefficients <- rep(0,N_var)
        constraint_coefficients[i+p] <- 1
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, type = ">=", 0)
    }

    ##start
    ## model evaluation begins on day start + 1 but collection
    ## begins start + 3
    ##
    ## What is this 10 in the next below??
    ## reply: our model started making prediction at day start+1,
    ## thus we can only decide on the new units
    ## come in at day start+4, the dates 1, 2, 3's new units are given by
    ## x_{start+1} = y_{start+1}+c+buffer; x_{start+2} = y_{start+2};x_{start+3} = y_{start+3}
    rhs_offset <- c(c + buffer, 0, 0)
    for (i in seq_len(3L)) {
        constraint_coefficients <- rep(0,N_var)
        constraint_coefficients[start+i+p+3*N] <- 1 ## x_{start + 1}
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", y[start+i] + rhs_offset[i])
    }

    ## Set
    ## Setting r1 to zero for start day?? Reply: yes.
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[start+p+N] <- 1 ## r_{start}(1)
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=",0)

    ## Setting r2 to zero for start day??  Reply: yes. The initial
    ## values for x_{start+1}, x_{start+2}, x_{start+3} have made sure
    ## the feasibility of this problem. The model traininig part here
    ## does not involve any initial information we input as a
    ## practioner, including remining bags, new bags will be arriving
    ## for the next three days.
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[start+p+N*2] <- 1 ## r_{start}(2)
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=",0)

    ##prediction
    ## Eqn 11
    ## Total need t_i. Z is the covariate matrix (p columns)
    for(i in (start+1):N) {
        constraint_coefficients <- rep(0,N_var)
        constraint_coefficients[i+p+4*N] <- 1 ## t_i
        constraint_coefficients[1:p] <- -as.numeric(XX[i,]) ## p features on day i
        constraint_coefficients[p+1+5*N] <- -1 ## negate intercept
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", 0)
    }

    for (i in (start+1):N){

        ##remaining
        ##r1_{i} - r2_{i-1}- r_{i-1} + w_{i} >= -y_{i}
        ##w_{i} + r2_{i} + r1_{i}- r1_{i-1}-r2_{i-1} - x_{i} == -y_{i}
        ## r2_{i} - x_{i} <= 0

        ## Eqn 14
        constraint_coefficients <- rep(0,N_var)
        constraint_coefficients[i+p+N] <- 1 #r_i(1)
        constraint_coefficients[i+p] <- 1 #w_i
        constraint_coefficients[i+p+N-1] <- -1 #r_{i-1}(1)
        constraint_coefficients[i+p+2*N-1] <- -1 # r_{i-1}(2)
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", -y[i])

        ## Eqn 16
        constraint_coefficients <- rep(0,N_var)
        constraint_coefficients[i+p+N] <- 1 # r_i(1) >= 0
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)

        ## Eqn 15
        constraint_coefficients <- rep(0,N_var)
        constraint_coefficients[i+p+N] <- 1 ## r_i(1)
        constraint_coefficients[i+p+N*2] <- 1 ## r_i(2)
        constraint_coefficients[i+p] = 1 ## w_i
        constraint_coefficients[i+p+3*N] <- -1 ## x_i
        constraint_coefficients[i+p+N-1] <- -1 ## r_{i-1}(1)
        constraint_coefficients[i+p+N*2-1] <- -1 ## r_{i-1}(2)
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", -y[i])

        ## Eqn 17. Why is it missing c_0? Maybe because of start day??
        ## Reply: Yes! Later from day start+4->last day, we have added
        ## the c there. You have noticed it there. Here we only need
        ## to initialize start+1->start+3, but I was being lazy and
        ## had not written them in a seperate small loop.

        constraint_coefficients <- rep(0,N_var)
        constraint_coefficients[i+p+N*2] <- 1 ## r_i(2)
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)

        ##wasted
        ##w_{i} - r1_{i-1} >=  - y_{i}
        ## Eqn 12
        constraint_coefficients <- rep(0,N_var)
        constraint_coefficients[i+p] <- 1 ## w_i
        constraint_coefficients[i+p+N-1] <- -1 ## r_{i-1}(1)
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", -y[i])
        ##constraint_coefficients[i+p+N-1] = -1 ## repeated from 2 lines above
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "<=", 0)
    }

    ## Only for cross-validation
##    if (doing_cv) {  ## Fix per Lexi's email
    for (i in ind) {
        if(i >= start+4){
            constraint_coefficients <- rep(0,N_var)
            constraint_coefficients[i+p+N*2] = 1 ## r_i(2)
            lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", c)
        }
    }
##    }

    for (i in seq.int(start+4, N)) {

        ##collecting
        ##x_{i+3} +  x_{i+2}+x_{i+1}  + r2_{i} + r1_{i} -t_i== 0:start+N

        ## Eqn 18??  Error in transcribing equation 8 to 18.
        ## Implements equation 8
        constraint_coefficients <- rep(0,N_var)
        constraint_coefficients[i+3*N+p] <- 1 ## x_{i}
        if (i %in% ind) {
            constraint_coefficients[i+3*N-1+p] <- 1 ## x_{i-1}
            constraint_coefficients[i+3*N-2+p] <- 1 ## x_{i-2}
            constraint_coefficients[i+p+N-3] <- 1 ## r_{i-3}(1)
            constraint_coefficients[i+p+N*2-3] <- 1 ## r_{i-3}(2)
            constraint_coefficients[i+p+N*4-3] <- -1 ## t_{i-3}
            lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", 0)
        } else {
            lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", y[i])

        }
    }
    ## End of Only for cross-validation

    ##write.lp(my.lp,'model.lp',type='lp')
    for (j in 1:p) {
        ## penalty term constraints
        constraint_coefficients <- rep(0,N_var)
        constraint_coefficients[j] <- 1 ## beta_j
        constraint_coefficients[j+N*5+p+1] <- 1 ## l1 constraint
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)

        ## penalty term constraints
        constraint_coefficients[j] <- -1  ## beta_j
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)
    }

    lpSolveAPI::set.bounds(my.lp, lower = rep(-1000,ncol(my.lp)), upper = rep(1000,ncol(my.lp)))
    lpSolveAPI::lp.control(my.lp,sense='min')

    ##write.lp(my.lp,'model.lp',type='lp')
    coeffients_matrix = matrix(0, ncol = length(lamb), nrow = p+1)

    ##nConstraints <- nrow(my.lp) + 1  ## we add a constraint in the loop following..
    for (i in seq_along(lamb)) {
        ## What is this 9 in the second line below?
        ##Reply: here, we did not penalize day of week and the
        ## plg avg features, so there are 8 features left out, that is why we start
        ## from the 9th predictor
        constraint_coefficients <- rep(0,N_var)
        constraint_coefficients[(p + 5*N + 1 + 9):(p + 5*N + 1 + p)] <- 1
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "<=", lamb[i])

        ##print(my.lp)
        ##assignInMyNamespace("lplist", c(lplist, list(my.lp)))

        status <- solve(my.lp)
        if (status > 0) {
            stop(sprintf("lpSolveStatus is %d", status))
        }

        coeffs <- c(lpSolveAPI::get.variables(my.lp)[5*N+p+1],
                    lpSolveAPI::get.variables(my.lp)[1:p])
        coeffients_matrix[ , i] = coeffs  ## all the beta_j's

        ## Remove last constraint
        ##lpSolveAPI::delete.constraint(my.lp, nConstraints)

        ## aa = as.numeric(lpSolveAPI::get.variables(my.lp)[1:p])
        ## ## Store the coefficients for each lambda
        ## coeffients_matrix[-1,i] = aa  ## all the beta_j's
        ## coeffients_matrix[1, i] = as.numeric(lpSolveAPI::get.variables(my.lp)[5*N+p+1])  ## the intercept
        ## ## Compute the predictions t_i
        ## predictions[,i] = apply(XX,1, function(x) sum(x*coeffients_matrix[-1,i])+coeffients_matrix[1,i])
    }


    ## aa <- as.numeric(lpSolveAPI::get.variables(my.lp)[1:p])  ## all the beta_j's
    ## coeffs <- c(as.numeric(lpSolveAPI::get.variables(my.lp)[5*N+p+1]),aa)  ## the intercept
    ## coeffs
    coeffients_matrix
}

