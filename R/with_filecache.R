CACHE_DIR <- '.Rcache'
#' Execute an expression, while saving result to a local file
#' @param expr R expression whose result should be cached
#' @param filename A unique name for cached item
#' @param cache_dir The directory in which to cache result. Defaults to '.Rcache' (to change default, use \code{fcutils_options})
#' @param parse (optional) Whether to parse expression before executing (default TRUE for character exprs, otherwise FALSE)
#' @return The contents of the cache-file, if that file exists. Otherwise, result of expr
#' @examples
#' res <- with_filecache(some_long_running_function(),
#'                      filename = 'cache_this_result.Rds',
#'                      cache_dir = 'name_of_this_document')
#' @export
with_filecache <- function(expr, filename, cache_dir=fcutils_options()$cache_dir, parse=NULL) {
  ## prepare expr for evaluation
  expr <- substitute(expr)
  if (is.null(parse)) {
    if ("character" %in% class(expr)) {
      parse=TRUE;
      warning("Detected a character expr; consider wrapping in {} instead of quotes.");
    } else {
      parse=FALSE;
    }
  }
  
  cache_file <- file.path(cache_dir, filename)
  if (!dir.exists(cache_dir)) {
    warning('Cache directory does not exist. Creating one.')
    dir.create(cache_dir, recursive=TRUE)
  }
  if (file.exists(cache_file)) {
    try({obj <- readRDS(cache_file)})
    if (!inherits(obj, 'try-error')) {
      return(obj)
    } else {
      warning('Error reading RDS file -- re-executing expr')
    }
  }
  
  ## evaluate expr
  if (parse) {
    obj = eval(parse(text=expr));	
  } else {
    obj = eval(expr)
  }
  
  saveRDS(obj, file = cache_file, compress = TRUE)
  return(obj)
}