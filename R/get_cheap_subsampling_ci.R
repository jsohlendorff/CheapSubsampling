### get_cheap_subsampling_ci.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: May 28 2024 (08:38) 
## Version: 
## Last-Updated: May 28 2024 (08:51) 
##           By: Thomas Alexander Gerds
##     Update #: 2
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
get_cheap_subsampling_ci <- function(est, boot_est, size, n_val, alpha) {
    b_val <- length(boot_est)
    s_val <- sqrt(mean((est - boot_est)^2))
    tq <- stats::qt(1 - alpha / 2, df = b_val)
    list(
        estimate = est,
        cheap_lower = est - tq * sqrt((size) / (n_val - size)) * s_val,
        cheap_upper = est + tq * sqrt((size) / (n_val - size)) * s_val
    )
}




######################################################################
### get_cheap_subsampling_ci.R ends here
