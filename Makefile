PYTHON ?= .venv/bin/python

.PHONY: init analisis pdfs test check

init:
	python3 -m venv .venv
	$(PYTHON) -m pip install --upgrade pip
	$(PYTHON) -m pip install -r requirements.txt

analisis:
	$(PYTHON) scripts/run_analysis.py

pdfs:
	$(PYTHON) scripts/make_pdfs.py

test:
	$(PYTHON) -m unittest discover -s tests

check:
	$(PYTHON) -m py_compile scripts/run_analysis.py scripts/make_pdfs.py
	$(PYTHON) -m unittest discover -s tests
