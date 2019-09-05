#
# Compile Splash Fortran source to produce libsplash.so.
#
# The main targets are:
#   - build:   to build the shared object library
#   - install: to build and install the library to INSTALL_DIR
#   - conda:   to build the library when building the Plonk Conda package
#
# Daniel Mentiplay, 2019
#

################################################################################

# The Splash shared object library is installed to INSTALL_DIR.
# By default this is `~/anaconda/lib`. Set INSTALL_DIR to install elsewhere.

ifeq (X$(INSTALL_DIR), X)
    ifeq (X$(CONDA_PREFIX), X)
        CONDA_PREFIX = $(HOME)/anaconda
    endif
    INSTALL_DIR = $(CONDA_PREFIX)/lib
endif

# Splash compile options

ifndef $(SYSTEM)
    SYSTEM = gfortran
    FC = $(SYSTEM)
endif
ifndef $(OPENMP)
    OPENMP = yes
endif
ifndef $(DOUBLEPRECISION)
	DOUBLEPRECISION = no
endif
ifndef $(DEBUG)
    DEBUG = no
endif

################################################################################

SPLASH_GIT_COMMIT_HASH = 7cdbdd8b52eca6d2a78570f8647223af609742b3
LIBSPLASH = libsplash.so
SPLASH_DIR = splash

SHELL := /bin/bash
UNAME = ${shell uname}

FC_BIN = $(shell command -v $(FC))
FC_VERSION = $(shell $(FC) --version | head -1 | awk '{print $$NF}')

.PHONY: help
help: ## Display this help
	@echo "Makefile targets:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  %-20s %s\n", $$1, $$2}'

.PHONY: print
print:
	@echo
	@echo ">>> Compiling Splash modules"
	@echo ">>> Fortran compiler binary: $(FC_BIN)"
	@echo ">>> Fortran compiler version: $(FC_VERSION)"
	@echo

$(LIBSPLASH): $(SPLASH_DIR) print
	@echo
	@echo ">>> Linking Splash into a shared object library"
	@echo
	make -C $(SPLASH_DIR) \
		SYSTEM=$(SYSTEM) \
		DOUBLEPRECISION=$(DOUBLEPRECISION) \
		OPENMP=$(OPENMP) \
		DEBUG=$(DEBUG) \
		libsplash \
		&& cp $(SPLASH_DIR)/build/$(LIBSPLASH) .

.PHONY: build
build: print $(LIBSPLASH) ## Build libsplash

.PHONY: install
install: $(INSTALL_DIR) $(LIBSPLASH) ## Install libsplash to INSTALL_DIR
	@echo
	@echo ">>> Copying libsplash.so to INSTALL_DIR=$(INSTALL_DIR)"
	@echo
	cp -f $(LIBSPLASH) $(INSTALL_DIR)
ifeq ($(UNAME), Darwin)
	@echo
	@echo ">>> Fixing macOS library path issues with install_name_tool"
	@echo
	install_name_tool -id @rpath/$(LIBSPLASH) $(INSTALL_DIR)/$(LIBSPLASH)
endif

$(INSTALL_DIR):
	@echo
	@echo ">>> INSTALL_DIR=$(INSTALL_DIR) does not exist"
	@echo
	@exit 1

.PHONY: conda
conda: $(LIBSPLASH) ## Build libsplash for Plonk Conda package
	@echo
	@echo ">>> Copying libsplash.so to PREFIX/lib=$(PREFIX)/lib"
	@echo ">>> Only use this Makefile target for building the Conda package"
	@echo
	cp -f $(LIBSPLASH) $(PREFIX)/lib

$(SPLASH_DIR):
	@echo
	@echo ">>> Fetching Splash source code"
	@echo
	@git clone https://github.com/danieljprice/splash.git
	@cd splash && git checkout $(SPLASH_GIT_COMMIT_HASH)

.PHONY: clean
clean: ## Clean temporary build files
	@echo
	@echo ">>> Cleaning up temporary build files"
	@echo
	@\rm -rf $(SPLASH_DIR)

.PHONY: distclean
distclean: ## Clean all generated files
	@echo
	@echo ">>> Remove all generated files"
	@echo
	@\rm -rf $(SPLASH_DIR)
	@\rm -f $(LIBSPLASH)
