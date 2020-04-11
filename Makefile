## Makefile for generating R packages.
## original by 2011 Andrew Redd
## https://gist.github.com/halpo/1374344

## Modified by kaz-yos for a different directory structure.
## To be put in the package directory.
## The build and check output files are put in the package directory.
## These files as well as Makefile have to be .Rbuildignore.
## writing rules
## https://www.gnu.org/software/make/manual/html_node/Rules.html#Rules

## Extract the package name and version from the DESCRIPTION file.
PKG_NAME=$(shell    grep -i ^package: DESCRIPTION | cut -d : -d \  -f 2)
PKG_VERSION=$(shell grep -i ^version: DESCRIPTION | cut -d : -d \  -f 2)

## Define files to check for updates
R_FILES   := $(wildcard R/*.R)
TST_FILES := $(wildcard tests/testthat/*.R)
SAS_FILES := $(wildcard tests/reference_results/sas-*.sas)
SRC_FILES := $(wildcard src/*) $(addprefix src/, $(COPY_SRC))
VIG_FILES := $(wildcard vignettes/*)
PKG_FILES := DESCRIPTION NAMESPACE NEWS.md $(R_FILES) $(TST_FILES) $(SRC_FILES) $(VIG_FILES)


## .PHONY to allow non-file targets (file targets should not be here)
## https://www.gnu.org/software/make/manual/html_node/Phony-Targets.html
.PHONY: test lint winbuild vignettes readme pkgdown build check check_devtools revdep install sas clean


### Define targets

## test just runs testthat scripts. No dependencies.
test: NAMESPACE
	Rscript -e 'options(width = 120); options("testthat.progress.max_fails" = Inf); devtools::test()' | tee test-all.txt

## lintr
# https://github.com/jimhester/lintr#continuous
lint: NAMESPACE
	Rscript -e 'options(width = 120); lintr::lint_package()' | tee lint_package.txt

## winbuild always build regardless of file update status
## Links to results e-mailed (no useful output locally)
winbuild:
	Rscript -e "devtools::check_win_devel()"
	Rscript -e "devtools::check_win_release()"

## Build vignettes in inst/doc
vignettes: ${VIG_FILES}
	Rscript -e "devtools::build_vignettes()"

## Build README.md
readme: README.Rmd
	Rscript -e "rmarkdown::render('README.Rmd', output_format = 'md_document', output_file = 'README.md')"

## Build website
## https://pkgdown.r-lib.org
pkgdown: vignettes readme
	Rscript -e "pkgdown::build_site()"

## build depends on the *.tar.gz file, i.e., its own product.
## *.tar.gz file is defined seprately to prevent build execution on every invocation.
build: $(PKG_NAME)_$(PKG_VERSION).tar.gz

## (file target) The *.tar.gz file depends on package files including NAMESPACE,
## and build *.tar.gz file from these.
$(PKG_NAME)_$(PKG_VERSION).tar.gz: $(PKG_FILES)
	Rscript -e "devtools::build(pkg = '.', path = '.', manual = TRUE)"

## (file target) NAMESPACE depends on *.R files, and excecute roxygen2 on these.
## methods::is() is not automatically loaded by roxygen2 version 4
NAMESPACE: $(R_FILES)
	Rscript -e "devtools::document('.')"

## check requires the *.tar.gz file, and execute strict tests on it.
check: $(PKG_NAME)_$(PKG_VERSION).tar.gz
	R CMD check --as-cran ./$(PKG_NAME)_$(PKG_VERSION).tar.gz | tee cran-check.txt

check_devtools: $(PKG_NAME)_$(PKG_VERSION).tar.gz
	Rscript -e "options(width = 120); devtools::check(pkg = '.', check_dir = '.', manual = TRUE)" | tee cran-check.txt

## revdep requires the *.tar.gz file, and execute strict tests on it.
## https://github.com/r-lib/revdepcheck
revdep: $(PKG_NAME)_$(PKG_VERSION).tar.gz
	Rscript -e "devtools::revdep()"
	Rscript -e "options(width = 120); revdepcheck::revdep_check(num_workers = 6); revdepcheck::revdep_summary()" | tee revdep_check.txt

revdep_clean:
	Rscript -e "revdepcheck::revdep_reset()"

## install requires the *.tar.gz file, and execute installation using it.
install: $(PKG_NAME)_$(PKG_VERSION).tar.gz
	Rscript -e "devtools::install('.')"


## run sas analyses
# https://stackoverflow.com/questions/1789594/how-do-i-write-the-cd-command-in-a-makefile
# It should only run when 01_generate_sas_data.R is updated.
sas_data: tests/reference_results/01_generate_sas_data.R
	cd tests/reference_results/ ; \
	Rscript 01_generate_sas_data.R 2>&1 | tee 01_generate_sas_data.R.txt

# It should only run when either one of these are updated.
sas_scripts: tests/reference_results/02_generate_sas_macro_calls.R tests/reference_results/02_generate_sas_macro_calls_helpers.R
	cd tests/reference_results/ ; \
	Rscript 02_generate_sas_macro_calls.R 2>&1 | tee 02_generate_sas_macro_calls.R.txt

sas_scripts_clean:
	-for f in $(SAS_FILES); do \
	rm $$f ; \
	done;

# This one can only run on a linux server with sas.
# This depends on up-to-date sas_data and sas_scripts.
# These two must depend on files so that they do not run every time.
# This fails without harm on a system without the sas command.
sas: sas_data sas_scripts
	-cd tests/reference_results/ ; \
	for f in $(subst tests/reference_results/,,$(SAS_FILES)) ; do \
	sas $$f ; \
	done;

# It depends on sas which should be ok on an environment without sas.
sas_extract: sas
	for f in $(subst .sas,.lst,$(SAS_FILES)); do \
	cat $${f} | grep " cde \| nde \| cde=nde \| [pt]nde \| nie \| [pt]nie \| total effect " > $${f%lst}txt ; \
	done;

sas_clean:
	-for f in $(subst .sas,.log,$(SAS_FILES)); do \
	rm $$f ; \
	done;
	-for f in $(subst .sas,.lst,$(SAS_FILES)); do \
	rm $$f ; \
	done;


## clean has no dependency, and execute removal of make output files.
clean: revdep_clean
	-rm    -f $(PKG_NAME)_$(PKG_VERSION).tar.gz
	-rm -r -f $(PKG_NAME).Rcheck
	-rm -r -f man/*.Rd
	-rm -r -f NAMESPACE


## Define a target "list" that just prints the names of files.
.PHONY: list
list:
	@echo "R files:"
	@echo $(R_FILES)
	@echo
	@echo "SAS files:"
	@echo $(SAS_FILES)
	@echo
	@echo "Test files:"
	@echo $(TST_FILES)
	@echo
	@echo "Source files:"
	@echo $(SRC_FILES)
	@echo
	@echo "Vignettes:"
	@echo $(VIG_FILES)
