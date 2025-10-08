
all: latex/gradus.pdf

latex/gradus.pdf: latex/gradus.tex
	(cd latex && pdflatex gradus && bibtex gradus && pdflatex gradus)

diff:
	rm -rf diffs
	cp -r latex diffs
	git show submitted:latex/gradus.tex > /tmp/.submitted-gradus.tex
	latexdiff --type=BOLD /tmp/.submitted-gradus.tex latex/gradus.tex | \
		sed --posix \
			's/\\providecommand{\\DIFadd}\[1\]{{\\bf \#1}}/\\providecommand{\\DIFadd}\[1\]{{\\color{red} \\bf \#1}}/' \
			> diffs/gradus.tex
	(cd diffs && pdflatex gradus && bibtex gradus && pdflatex gradus)
	cp diffs/gradus.pdf diff.pdf

rediff:
	(cd diffs && pdflatex gradus && bibtex gradus && pdflatex gradus)
	cp diffs/gradus.pdf diff.pdf
