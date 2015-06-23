
.onAttach <- function(libname, pkgname) {

	packageStartupMessage("Welcome to CHiCAGO - version ", utils::packageDescription("Chicago", fields="Version"))
	packageStartupMessage('If you are new to CHiCAGO, please consider reading the vignette through the command: vignette("Chicago").')
}

