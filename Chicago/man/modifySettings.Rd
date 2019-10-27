\name{modifySettings}
\alias{modifySettings}
\title{
Modify Settings
}
\description{
Modify the settings in a \code{chicagoData} object.
}
\usage{
modifySettings(cd, designDir=NULL, settings=list(), settingsFile=NULL)
}
\arguments{
  \item{cd}{A \code{chicagoData} object.}
  \item{designDir}{The new location of the design directory, e.g "~/resources/path" or NULL if not modified.}
  \item{settings}{A named list containing settings to modify.}
  \item{settingsFile}{The location of a file containing settings or NULL if not provided. Each row should contain the name of a setting, followed by whitespace, followed by the value of that setting.}
}
\details{
\code{cd}'s settings are updated. For a list of available settings, see \code{\link{defaultSettings}}.  
}
\value{
An object of class \code{chicagoData}.
}
\author{
Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett
}

\seealso{
\code{\link{setExperiment}}, \code{\link{defaultSettings}}
}
\examples{
designDir <- file.path(system.file("extdata", package="Chicago"), "unitTestDesign")

cd <- setExperiment(designDir)

##Suppose I want to change the number of samples drawn for dispersion estimation 
print(settings(cd)$brownianNoise.subset)

cd <- modifySettings(cd, settings=list(brownianNoise.subset = 10))
print(settings(cd)$brownianNoise.subset)
}
