
.onAttach <- function(libname, pkgname) {

  packageStartupMessage("")
	packageStartupMessage("Welcome to CHiCAGO - version ", utils::packageDescription("Chicago", fields="Version"))
	packageStartupMessage('If you are new to CHiCAGO, please consider reading the vignette through the command: vignette("Chicago").')
	packageStartupMessage('NOTE: Default values of tlb.minProxOEPerBin and tlb.minProxB2BPerBin changed as of Version 1.1.5. No action is required unless you specified non-default values, or wish to re-run the pipeline on old chicagoData objects. See news(package="Chicago")')

}