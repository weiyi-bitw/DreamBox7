# sourcing R files using https from github
# reference: http://tonybreyal.wordpress.com/2011/11/24/source_https-sourcing-an-r-script-from-github/
#
source_https <- function(url,...){
	#load package
	require(RCurl)

	# parse and evaluate each .R script
	sapply(c(url, ...), function(u) {
		eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
	})
}
