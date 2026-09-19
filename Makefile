.DEFAULT_GOAL := help
PKG = riptide-ffa

check-sdist: ## Build and smoke-test the source distribution in a temporary virtualenv
	@set -eu; \
	tmpdir="$$(mktemp -d)"; \
	on_exit() { \
		status=$$?; \
		if [ "$$status" -eq 0 ]; then \
			printf '\033[1;32m✅ Source distribution smoke test passed\033[0m\n'; \
		else \
			printf '\033[1;31m❌ Source distribution smoke test failed\033[0m\n' >&2; \
		fi; \
		rm -rf "$$tmpdir"; \
		exit "$$status"; \
	}; \
	trap on_exit EXIT; \
	python -m build --sdist --outdir "$$tmpdir/dist"; \
	python -m venv "$$tmpdir/venv"; \
	"$$tmpdir/venv/bin/python" -m pip install "$$tmpdir"/dist/*.tar.gz; \
	cd "$$tmpdir"; \
	"$$tmpdir/venv/bin/python" -c "import riptide; print(riptide.__version__)"; \
	"$$tmpdir/venv/bin/rffa" --help >/dev/null; \
	"$$tmpdir/venv/bin/rseek" --help >/dev/null

install: ## Install the package in editable mode with dev dependencies
	pip install -e .[dev]

# GLORIOUS hack to autogenerate Makefile help
# This simply parses the double hashtags that follow each Makefile command
# https://marmelab.com/blog/2016/02/29/auto-documented-makefile.html
help: ## Print this help message
	@echo "Makefile help for ${PKG}"
	@echo "===================================================================="
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

test: ## Run the unit tests and print a coverage report
	pytest --cov=src/ --cov-report=term-missing

.PHONY: check-sdist install help test
