ROOT_DIR := $(shell dirname $(realpath $(lastword ${MAKEFILE_LIST})))

report/%.html: code/%.r
	mkdir -p "$(dir $@)"
	Rscript -e "pdf.options(encoding = 'CP1250'); ezknitr::ezspin('$<', out_dir = '${@D}', fig_dir = 'figures/$*/', wd = '${ROOT_DIR}', doc = '^## ', move_intermediate_file = FALSE)"
