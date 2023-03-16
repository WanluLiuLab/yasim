main.pdf: main.tex coefficient.tex plot.timestamp
	pdflatex main

plot.timestamp: plot.R data.csv
	Rscript plot.R

coefficient.tex: gen_table.py data.csv
	python gen_table.py
