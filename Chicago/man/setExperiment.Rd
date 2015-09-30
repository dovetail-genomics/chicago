\name{setExperiment}
\alias{setExperiment}
\title{
Set Experiment
}
\description{
Creates a template CHiCAGO experiment object. This should be the first function called.
}
\usage{
setExperiment(designDir = "", settings = list(),
      settingsFile = NULL, def.settings=defaultSettings())
}
\arguments{
  \item{designDir}{
  The location of the design directory, e.g "~/resources/path". (Should not end with a slash.)
}
  \item{settings}{
  A named list containing settings to apply. Setting \code{names(settings)[1]} is set to (settings)[[1]], and so on. This overrides anything specified in \code{settingsFile}, or in \code{def.settings}.
}
  \item{settingsFile}{
  The location of a file containing settings. Each row should contain the name of a setting, followed by whitespace, followed by the value of that setting. Overrides anything specified in \code{def.settings}.
}
  \item{def.settings}{
  These are the default settings.
}
}
\details{
For a list of settings, see \code{\link{defaultSettings}}.
}
\value{
An object of class \code{chicagoData}.
}
\author{
Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett
}

\seealso{
\code{\link{defaultSettings}}
}
\examples{
  designDir <- file.path(system.file("extdata", package="Chicago"), "unitTestDesign")

  cd <- setExperiment(designDir)
}
