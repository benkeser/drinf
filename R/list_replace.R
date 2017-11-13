#' Helper function to replace values in a list
#' @param x a list
#' @param list indices of list
#' @param values what values to replace with
#' @return x with list values replaced
#' @export
list_replace <- function(x, list, values){
    x[[list]] <- values; x
}