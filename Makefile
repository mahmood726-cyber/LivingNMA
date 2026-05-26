.PHONY: help install test cov lint fmt clean pre-commit-install pre-commit-run
PY ?= python

help:	## Show this help message
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "  \033[36m%-22s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

install:	## Install runtime + test deps + pre-commit
	$(PY) -m pip install --upgrade pip
	if [ -f requirements.txt ]; then $(PY) -m pip install -r requirements.txt; fi
	$(PY) -m pip install ruff pytest pytest-cov pre-commit
	pre-commit install

test:	## Run pytest (terse)
	$(PY) -m pytest -q --tb=short

cov:	## Run pytest with coverage
	$(PY) -m pytest -q --cov=. --cov-report=term-missing --cov-report=html

lint:	## Run ruff check
	$(PY) -m ruff check .

fmt:	## Run ruff format + auto-fix
	$(PY) -m ruff format .
	$(PY) -m ruff check --fix .

clean:	## Remove caches and coverage artifacts
	rm -rf .pytest_cache .ruff_cache htmlcov coverage.xml .coverage
	find . -type d -name __pycache__ -prune -exec rm -rf {} +

pre-commit-install:	## Install pre-commit hooks
	pre-commit install

pre-commit-run:	## Run all pre-commit hooks
	pre-commit run --all-files
