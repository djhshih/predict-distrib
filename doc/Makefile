all: predict-distrib.pdf
	

%.pdf: %.qmd
	quarto render $<

clean:
	rm -rf *.pdf

