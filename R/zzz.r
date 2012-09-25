.onAttach <- function(libname, pkgname){
	DBVer <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), fields="Version")
	packageStartupMessage(
	"\n",
	"===============================\n",
	"\n",
	paste(pkgname, DBVer), "\n",
	"\n",
	"===============================\n",
	"\n",
	"Welcome to DreamBox7 package!\n\n",
	" -- Brought to you by Wei-Yi Cheng and Team Attractor Metagenes.\n\n",
	"===============================\n\n",
	"Aren't Attractor Metagenes awesome? If you agree, please recommend our paper at  \n",
	"http://arxiv.org/abs/1204.6538 to everybody and help us further spread the idea!\n",
	"We hope one day cancer patients can be profiled by several metagenes indicating \n",
	"important biomolecular events of the disease such as chromosomal instability, \n",
	"invasiveness, or lymphocyte responses, instead of mutually exclusive subtypes.\n",
	"This may lead to a more comprehensive understanding of cancer, a better prognosis,\n",
	"and even new druggable targets! Let's get it started! \n",
	"                                                           Wei-Yi Cheng 09/24/2012\n",
	"\n"
	)
}


