
all: latex/gradus.pdf

latex/gradus.pdf: latex/gradus.tex
	(cd latex && pdflatex gradus && bibtex gradus && pdflatex gradus)
