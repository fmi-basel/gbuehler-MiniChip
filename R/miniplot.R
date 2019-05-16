#' @title miniplot
#'
#' @description this makes a silly plot
#'
#' @param number1 the first number
#' @param number2 the second number
#'
#' @return a plot
#'
#' @examples
#' miniplot(1,2)
#'
#' @importFrom graphics plot
#'
#' @export
miniplot <- function(number1,number2){
  graphics::plot(number1,number2)
}
