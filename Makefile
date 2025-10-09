R := "/c/Program Files/R/R-4.3.3/bin/x64/Rscript.exe" --vanilla

all: normalize de go goi qc gonet

normalize:
	$(R) scripts/01_normalize.R

de:
	$(R) scripts/02_de_limma.R

go:
	$(R) scripts/03_go_enrichment.R

goi:
	$(R) scripts/04_plots_and_goi.R

qc:
	$(R) scripts/05_qc_and_volcano.R

gonet:
	$(R) scripts/06_go_networks.R
