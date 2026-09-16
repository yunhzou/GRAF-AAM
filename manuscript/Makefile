PYTHON ?= python3
.PHONY: all figures verify bundle
all:
	MANUSCRIPT_PYTHON="$(PYTHON)" bash scripts/build.sh
figures:
	$(PYTHON) scripts/build_figures.py
verify:
	$(PYTHON) scripts/validate_outputs.py
bundle:
	$(PYTHON) scripts/package_bundle.py
