# Variable, global to package's namespace. 
# This function is not exported to user space and does not need to be documented.
MYPKGOPTIONS <- settings::options_manager(cache_dir = '.Rcache')

# User function that gets exported:

#' Set or get options for fcutils package
#' 
#' @param ... Option names to retrieve option values or \code{[key]=[value]} pairs to set options.
#'
#' @section Supported options:
#' The following options are supported
#' \itemize{
#'  \item{\code{cache_dir}}{(\code{character};.cache_dir) Path used as local file cache}
#' }
#' @import settings
#' @export
fcutils_options <- function(...){
  # protect against the use of reserved words.
  stop_if_reserved(...)
  MYPKGOPTIONS(...)
}

#' Reset global options for fcutils
#'
#' @export
fcutils_options_reset <- reset(MYPKGOPTIONS)