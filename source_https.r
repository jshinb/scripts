## 2015-10-29
## taken from: http://www.r-bloggers.com/source_https-sourcing-an-r-script-from-github-over-https/
## to source R-files from my github repository: 'https://github.com/jshinb/scripts'

source_https <- function(url, ...) {
  # load package
  require(RCurl)
 
  # parse and evaluate each .R script
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
  })
}
 
# Example
#source_https("https://raw.github.com/tonybreyal/Blog-Reference-Functions/master/R/bingSearchXScraper/bingSearchXScraper.R",
#             "https://raw.github.com/tonybreyal/Blog-Reference-Functions/master/R/htmlToText/htmlToText.R")
