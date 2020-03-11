LEC=eda

FIGS=Figs/il3.pdf \
	 Figs/adipose_weight.pdf \
	 Figs/body_weight.pdf \
	 Figs/gough_weights1.pdf

include ../Makefile_lectures.mk

Figs/%.pdf: R/%.R
	cd $(<D);R $(R_OPTS) -e "source('$(<F)')"

Figs/gough_weights1.pdf: R/gough_weights.R
	cd $(<D);R $(R_OPTS) -e "source('$(<F)')"
