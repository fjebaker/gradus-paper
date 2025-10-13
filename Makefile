
BUILD_DIR := build
LATEX_COMPILE := pdflatex -interaction=nonstopmode
PROJECT := gradus

BUNDLE := $(PROJECT)-bundle

all: $(BUILD_DIR)/$(PROJECT).pdf

bundle: $(PROJECT)-bundle.tar.gz

$(BUNDLE).tar.gz: $(BUILD_DIR)/$(PROJECT).pdf
	rm -rf $(BUNDLE)
	mkdir -p $(BUNDLE)
	cp -r latex/* $(BUNDLE)/
	cp $(BUILD_DIR)/$(PROJECT).bbl $(BUNDLE)/
	cp $(BUILD_DIR)/$(PROJECT).pdf $(BUNDLE)/reference.pdf
	cp bundle-Makefile $(BUNDLE)/Makefile
	tar czf $(BUNDLE).tar.gz $(BUNDLE)

$(BUILD_DIR)/$(PROJECT).pdf: latex/$(PROJECT).tex
	mkdir -p $(BUILD_DIR)
	cp -r latex/* $(BUILD_DIR)
	(cd $(BUILD_DIR) && \
		$(LATEX_COMPILE) $(PROJECT) && \
		bibtex $(PROJECT) && \
		$(LATEX_COMPILE) $(PROJECT) && \
		$(LATEX_COMPILE) $(PROJECT) \
	)

diff:
	rm -rf diffs
	cp -r latex diffs
	git show submitted:latex/$(PROJECT).tex > /tmp/.submitted-$(PROJECT).tex
	latexdiff --type=BOLD /tmp/.submitted-$(PROJECT).tex latex/$(PROJECT).tex | \
		sed --posix \
			's/\\providecommand{\\DIFadd}\[1\]{{\\bf \#1}}/\\providecommand{\\DIFadd}\[1\]{{\\color{red} \\bf \#1}}/' \
			> diffs/$(PROJECT).tex
	(cd diffs && $(LATEX_COMPILE) $(PROJECT) && bibtex $(PROJECT) && $(LATEX_COMPILE) $(PROJECT))
	cp diffs/$(PROJECT).pdf diff.pdf

rediff:
	(cd diffs && $(LATEX_COMPILE) $(PROJECT) && bibtex $(PROJECT) && $(LATEX_COMPILE) $(PROJECT))
	cp diffs/$(PROJECT).pdf diff.pdf
