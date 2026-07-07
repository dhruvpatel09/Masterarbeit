PDF=main.pdf

all:
	latexmk -pdf -interaction=nonstopmode main.tex

clean:
	latexmk -C
	rm -f *.run.xml *.bcf *.synctex.gz
