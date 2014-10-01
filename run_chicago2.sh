## EDIT ME
chicagopath="/bi/scratch/wingetts/chicago_wrapper/chicagov2/"
R --slave --no-restore-data --args "$@" < ${chicagopath}/production_line_CHiCAGO2.R
