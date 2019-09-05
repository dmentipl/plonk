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
	@make -C splash install
	@echo
	@echo ">>> Build Splash Cython extension"
	@echo
	@python setup.py build_ext --inplace
	@echo
	@echo ">>> Install Plonk in 'conda development' mode"
	@echo
	@conda develop .

.PHONY: conda
conda: ## Build Conda package
	@echo
	@echo ">>> Build Conda package"
	@echo
	@make -C splash conda

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
	@make -C splash clean
	@\rm -rf .coverage htmlcov
	@\rm -rf splash/splash.c
	@\rm -rf splash/libsplash.so

.PHONY: distclean
distclean: ## Clean all generated files
	@make -C splash distclean
	@\rm -rf splash/splash.c
	@\rm -rf splash/*.so
