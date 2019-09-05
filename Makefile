#
# Plonk Makefile
#
# Daniel Mentiplay, 2019
#

.PHONY: clean conda development docs help

development: ## Setup Plonk for development
	@make -C splash install
	@echo
	@echo ">>> Build Splash Cython extension"
	@echo
	@python setup.py build_ext --inplace
	@echo
	@echo ">>> Install Plonk in 'conda development' mode"
	@echo
	@conda develop .

help: ## Display this help
	@echo "Makefile targets:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  %-20s %s\n", $$1, $$2}'

conda: ## Build Conda package
	@echo
	@echo ">>> Build Conda package"
	@echo
	@make -C splash conda

docs: ## Build documentation
	@echo
	@echo ">>> Build documentation"
	@echo
	@make -C docs html

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

clean: ## Clean up
	@make -C splash clean
	@\rm -rf .coverage htmlcov
	@\rm -rf splash/splash.c splash/splash.*.so
