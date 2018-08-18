##lplist <- list()

#' Solve an LP problem
#' @param d the data frame
#' @param ind representing the fold
#' @param lambds the lambdas to use in penalization
#' @param start the start date
#' @param c the value for c_0
#' @importFrom lpSolveAPI make.lp add.constraint set.bounds lp.control
#' @importFrom lpSolveAPI set.objfn get.variables
#'
#' @export
#'
cv_lpSolve <- function(d, ind, lambds, start = 10, c = 30){
    XX <- d[, 2:44]
    y <- d[, 45]
    N =  nrow(XX)
    p =  ncol(XX)
    predictions = matrix(NA, N,length(lambds))

    ## Remember y is usage, so y[i] is known usage on day i
    ## p betas, N ws, N r1s, N r2s, N xs, N ts, 1 intercept, p bounds
    N_var = p + 5*N + 1 + p
    my.lp <- lpSolveAPI::make.lp(0,N_var)

    obj_coefficients = rep(0,N_var)
    ## Only use coefficients for observations in the fold! (eq 9.)
    obj_coefficients[(p+ind)] = 1
    lpSolveAPI::set.objfn(my.lp,obj_coefficients)

    ## Constrain  first week sum of coefs of beta to be zero??
    # Why?
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[1:7] = 1
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, type = "=", 0)

    ##w >=0
    for(i in 1:N){
        ## Since w_i >= 0, using more than the ind values should not matter here
        ## especially because we are minimizing sum(w_i)!
        ## Setting each coeff of w to be non-negative: (eq 13)
        constraint_coefficients <- rep(0,N_var)
        constraint_coefficients[i+p] = 1
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, type = ">=", 0)
    }

    ##start
    ## model evaluation begins on day start + 1 but collection begins start + 3
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[start+1+p+3*N] = 1 ## x_{start + 1}
    ## What is this 10 in the next below?
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", y[start+1]+c+10)

    ## This seems to be matching the next two days.??
    for(i in 2:3){
        constraint_coefficients <- rep(0,N_var)
        constraint_coefficients[start+i+p+3*N] = 1 ## x_{start + 2}, x_{start + 3}
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", y[start+i])
    }

    ## Setting r1 to zero for start day??
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[start+p+N] = 1 ## r_{start}(1)
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=",0)

    ## Setting r2 to zero for start day??
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[start+p+N*2] = 1 ## r_{start}(2)
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=",0)

    ##prediction
    ## Eqn 11
    ## Total need t_i. Z is the covariate matrix (p columns)
    for(i in (start+1):N){
        constraint_coefficients <- rep(0,N_var)
        constraint_coefficients[i+p+4*N] = 1 ## t_i
        constraint_coefficients[1:p] = -as.numeric(XX[i,]) ## p features on day i
        constraint_coefficients[p+1+5*N] = -1 ## negate intercept
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", 0)
    }

    for(i in (start+1):N){
        ##remaining
        ##r1_{i} - r2_{i-1}- r_{i-1} + w_{i} >= -y_{i}
        ##w_{i} + r2_{i} + r1_{i}- r1_{i-1}-r2_{i-1} - x_{i} == -y_{i}
        ## r2_{i} - x_{i} <= 0

        ## Eqn 14
        constraint_coefficients <- rep(0,N_var)
        constraint_coefficients[i+p+N] = 1 #r_i(1)
        constraint_coefficients[i+p] = 1 #w_i
        constraint_coefficients[i+p+N-1] = -1 #r_{i-1}(1)
        constraint_coefficients[i+p+2*N-1] = -1 # r_{i-1}(2)
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", -y[i])

        ## Eqn 16
        constraint_coefficients <- rep(0,N_var)
        constraint_coefficients[i+p+N] = 1 # r_i(1) >= 0
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)

        ## Eqn 15
        constraint_coefficients <- rep(0,N_var)
        constraint_coefficients[i+p+N] = 1 ## r_i(1)
        constraint_coefficients[i+p+N*2] = 1 ## r_i(2)
        constraint_coefficients[i+p] = 1 ## w_i
        constraint_coefficients[i+p+3*N] = -1 ## x_i
        constraint_coefficients[i+p+N-1] = -1 ## r_{i-1}(1)
        constraint_coefficients[i+p+N*2-1] = -1 ## r_{i-1}(2)
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", -y[i])

        ## Eqn 17. Why is it missing c_0? Maybe because of start day??
        constraint_coefficients <- rep(0,N_var)
        constraint_coefficients[i+p+N*2] = 1 ## r_i(2)
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)

        ##wasted
        ##w_{i} - r1_{i-1} >=  - y_{i}
        ## Eqn 12
        constraint_coefficients <- rep(0,N_var)
        constraint_coefficients[i+p] = 1 ## w_i
        constraint_coefficients[i+p+N-1] = -1 ## r_{i-1}(1)
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", -y[i])
        ## constraint_coefficients[i+p+N-1] = -1 ## repeated from 2 lines above
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "<=", 0)
    }

    ## Only for cross-validation
    for(i in ind){
        if(i >= start+4){
            constraint_coefficients <- rep(0,N_var)
            constraint_coefficients[i+p+N*2] = 1 ## r_i(2)
            lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", c)
        }
    }
    ## End of Only for cross-validation

    for(i in (start+4):N){
        ##collecting
        ##x_{i+3} +  x_{i+2}+x_{i+1}  + r2_{i} + r1_{i} -t_i== 0:start+N

        ## Eqn 18??  Error in transcribing equation 8 to 18.
        ## Implements equation 8

        constraint_coefficients <- rep(0,N_var)
        if(i%in%ind){
            constraint_coefficients[i+3*N+p] = 1
            constraint_coefficients[i+3*N+p] = 1 ## x_{i}
            constraint_coefficients[i+3*N-1+p] = 1 ## x_{i-1}
            constraint_coefficients[i+3*N-2+p] = 1 ## x_{i-2}
            constraint_coefficients[i+p+N-3] = 1 ## r_{i-3}(1)
            constraint_coefficients[i+p+N*2-3] = 1 ## r_{i-3}(2)
            constraint_coefficients[i+p+N*4-3] = -1 ## t_{i-3}
            lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", 0)
        } else{
            constraint_coefficients[i+3*N+p] = 1 ## x_i
            lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", y[i])
        }
    }

    ##write.lp(my.lp,'model.lp',type='lp')
    coeffients_matrix = matrix(0,ncol = length(lambds), nrow = p+1)
    for(j in 1:p){
        ## penalty term constraints
        constraint_coefficients <- rep(0,N_var)
        constraint_coefficients[j] = 1 ## beta_j
        constraint_coefficients[j+N*5+p+1] = 1 ## l1 constraint
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)

        ## penalty term constraints
        constraint_coefficients[j] = -1 ## beta_j
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)
    }

    lpSolveAPI::set.bounds(my.lp, lower = rep(-1000,ncol(my.lp)), upper = rep(1000,ncol(my.lp)))
    lpSolveAPI::lp.control(my.lp,sense='min')

    for(i in 1:length(lambds)){
        constraint_coefficients <- rep(0,N_var)

        ## What is this 9 below?
        constraint_coefficients[(p+5*N+1+9):(p+5*N+1+p)] = 1
        ## This assumes that the lambda's are ordered, else results will be wrong!!!
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "<=", lambds[i])
        ##assignInMyNamespace("lplist", c(lplist, list(my.lp)))
        ##print(my.lp)
        status <- solve(my.lp)
        if (status > 0) {
            stop(sprintf("lpSolveStatus is %d", status))
        }

        aa = as.numeric(lpSolveAPI::get.variables(my.lp)[1:p])
        ## Store the coefficients for each lambda
        coeffients_matrix[-1,i] = aa  ## all the beta_j's
        coeffients_matrix[1,i] = as.numeric(lpSolveAPI::get.variables(my.lp)[5*N+p+1])  ## the intercept
        ## Compute the predictions t_i
        predictions[,i] = apply(XX,1, function(x) sum(x*coeffients_matrix[-1,i])+coeffients_matrix[1,i])
    }
    return(predictions = predictions)
}
