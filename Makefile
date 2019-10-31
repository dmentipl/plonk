#
# Plonk Makefile
#
# Daniel Mentiplay, 2019
#

.PHONY: help
help: ## Display this help
	@echo "Makefile targets:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  %-20s %s\n", $$1, $$2}'

.PHONY: development
development: ## Setup Plonk for development
	@echo ">>> Install Plonk in 'conda development' mode"
	@echo
	@conda env create --file environment.yml
	@conda develop --name plonk-dev .
	@echo "To use, type: 'conda activate plonk-dev'"

.PHONY: conda
conda: ## Build Conda package
	@echo
	@echo ">>> Build Conda package"
	@echo
	@conda build conda

.PHONY: docs
docs: ## Build documentation
	@echo
	@echo ">>> Build documentation"
	@echo
	@make -C docs html

.PHONY: test
test: ## Run tests
	@echo
	@echo ">>> Run tests"
	@echo
	@python -m coverage run -m unittest discover && coverage html
	@echo
	@echo ">>> Check formatting with isort and black"
	@echo
	@isort --check-only -rc
	@black --check --skip-string-normalization .

.PHONY: clean
clean: ## Clean temporary build files
	@\rm -rf .coverage htmlcov
