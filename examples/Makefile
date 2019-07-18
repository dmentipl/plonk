#
# Make notebook, rst, html, slides, etc. from markdown

.PHONY:
help:
	@echo "Usage:"
	@echo "  make <target>"
	@echo
	@echo "  all              Build all"
	@echo "  clean            Clean up"
	@echo "  docs             Make .rst for docs"
	@echo "  help             Display this help"
	@echo "  html             Make .html notebook for viewing"
	@echo "  notebook         Make .ipynb notebook for running"
	@echo "  slides           Make .html reveal.js slides"

tutorial.ipynb: tutorial.md
	jupytext --to ipynb tutorial.md

tutorial.html: tutorial.ipynb
	jupyter nbconvert --to html --execute tutorial.ipynb  --ExecutePreprocessor.timeout=180

tutorial.slides.html: tutorial.ipynb
	jupyter nbconvert --to slides --execute tutorial.ipynb --ExecutePreprocessor.timeout=180

tutorial.rst: tutorial.ipynb
	jupyter nbconvert --to rst --execute tutorial.ipynb --ExecutePreprocessor.timeout=180

.PHONY: all clean docs html notebook slides

all: notebook docs html slides

docs: tutorial.rst

html: tutorial.html

notebook: tutorial.ipynb

slides: tutorial.slides.html

clean:
	rm -rf tutorial.ipynb tutorial.html tutorial.slides.html tutorial.rst tutorial_files
