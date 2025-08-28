PY=python3
PIP=pip

.PHONY: setup test fmt lint build clean

setup:
	$(PIP) install -U pip
	$(PIP) install -e .[dev]

test:
	$(PY) -m pytest -q

fmt:
	$(PY) -m black src tests

lint:
	$(PY) -m ruff check src tests

build:
	$(PY) -m build

clean:
	rm -rf build dist *.egg-info .pytest_cache __pycache__
