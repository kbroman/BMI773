LEC=steps2rr

FIGS=Figs/data_dict.pdf

include ../Makefile_lectures.mk

Figs/data_dict.pdf: R/spreadsheets.R
	cd $(<D);R $(R_OPTS) -e "source('$(<F)')"
