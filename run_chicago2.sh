SCRIPT=$(readlink -f $0)
chicagopath=`dirname $SCRIPT`

R --slave --no-restore-data --args "$@" < ${chicagopath}/production_line_CHiCAGO2.R
