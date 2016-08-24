
.onAttach <- function(libname, pkgname) {

	packageStartupMessage("Welcome to CHiCAGO - version ", utils::packageDescription("Chicago", fields="Version"))
	packageStartupMessage('If you are new to CHiCAGO, please consider reading the vignette through the command: vignette("Chicago").')
	packageStartupMessage('NOTE: tlb.minProxOEPerBin and tlb.minProxB2BPerBin defaults changed as of Version 1.1.5. See news(package="Chicago")')

}