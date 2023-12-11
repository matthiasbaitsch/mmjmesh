all: test book apidoc

book:
	cd demo && quarto render --to html

apidoc:
	cd docs && julia make.jl

.PHONY: test
test:
	julia test/runtests.jl