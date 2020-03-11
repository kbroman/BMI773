LEC=eda

FIGS=Figs/il3.pdf \
	 Figs/adipose_weight.pdf \
	 Figs/body_weight.pdf

include ../Makefile_lectures.mk

Figs/%.pdf: R/%.R
	cd $(<D);R $(R_OPTS) -e "source('$(<F)')"
