
#' modify coloring in DT::datatable json
#' based on work by `NicE`, posted here on SO: https://stackoverflow.com/a/33524422/3457743
#' @importFrom htmlwidgets JS
#' @export
color_from_middle <- function (data, color1, color2) 
{
  max_val <- max(abs(data))
  htmlwidgets::JS(sprintf("isNaN(parseFloat(value)) || value < 0 ? 'linear-gradient(90deg, transparent, transparent ' + (50 + value/%s * 50) + '%%, %s ' + (50 + value/%s * 50) + '%%,%s  50%%,transparent 50%%)': 'linear-gradient(90deg, transparent, transparent 50%%, %s 50%%, %s ' + (50 + value/%s * 50) + '%%, transparent ' + (50 + value/%s * 50) + '%%)'",
                          max_val,color1,max_val,color1,color2,color2,max_val,max_val))
}

