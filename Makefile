.PHONY: default black flake8 test test-v test-vv

default: black flake8

black:
	black -l 79 .

black-check:
	black -l 79 --check .

flake8:
	flake8 --ignore=E126,W504,E704,W503,E123,E24,E121,E226,E203,W605  .

test:
	pytest

test-v:
	pytest -v

test-vv:
	pytest -vv
