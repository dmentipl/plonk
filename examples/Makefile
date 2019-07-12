#
# Make notebook, html, slides from markdown

.PHONY:
help:
	@echo "Usage:"
	@echo "  make <target>"
	@echo
	@echo "  notebook         Make .ipynb notebook"
	@echo "  html             Make .html notebook"
	@echo "  slides           Make .html reveal.js slides"
	@echo "  help             Display this help"

tutorial.ipynb: tutorial.md
	jupytext --to ipynb tutorial.md

tutorial.html: tutorial.ipynb
	jupyter nbconvert --to html --execute tutorial.ipynb --output tutorial.html --ExecutePreprocessor.timeout=180

tutorial.slides.html: tutorial.ipynb
	jupyter nbconvert --to slides --execute tutorial.ipynb --ExecutePreprocessor.timeout=180

.PHONY:
notebook: tutorial.ipynb

.PHONY:
slides: tutorial.slides.html

.PHONY:
html: tutorial.html

.PHONY:
clean:
	rm -f tutorial.ipynb tutorial.html tutorial.slides.html
