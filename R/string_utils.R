
#' Find all common substrings of two strings
#' Based on https://stackoverflow.com/questions/28261825/longest-common-substring-in-r-finding-non-contiguous-matches-between-the-two-str
#' @param a first string
#' @param b second string
#' @import stringi
common_substrings <- function(a, b) { 
  sbstr_locations <- stringi::stri_locate_all_regex(drop(attr(adist(a, b, counts=TRUE), "trafos")), "M+")[[1]]
  cmn_sbstr <- stringi::stri_sub(longest_string(c(a, b)), sbstr_locations)
}

#' Return the longest string of a vector/list of strings s
#' @param s character or list of strings
longest_string <- function(s) {
  s[which.max(nchar(s))]
}

#' given a list of strings a, find the longest substring common to all strings
#' @param a list or character vector of strings
longest_common_substring <- function(a) {
  a <- unlist(na.omit(a))

  # cross a with itself
  x <- purrr::cross2(a, a) %>%
    # exclude self-matches
    purrr::keep(~ .[[1]] != .[[2]]) %>%
    unique() %>%
    # for each pair, find common substrings
    purrr::map(~ common_substrings(.[[1]], .[[2]])) %>%
    # summarise unique sets of substrings
    unique() %>%
    na.omit() %>%
    purrr::map(~ .[nchar(.) > 1])
  
  # identify all unique substrings shared by at least 2 records
  all_substrings <- x %>% unlist() %>% unique()
  
  # keep the longest of all substrings present in all original strings (a+/b)
  all_substrings %>%
    purrr::keep(~ all(stringr::str_detect(unique(c(a, b)), pattern = .))) %>%
    longest_string(.)
}
