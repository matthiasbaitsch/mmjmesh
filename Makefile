all: test book apidoc

book:
	cd demo && quarto render --to html

book-publish:
	cd demo && quarto publish gh-pages --no-prompt

apidoc:
	cd docs && julia --project make.jl

.PHONY: test
test:
	julia --project test/runtests.jl

clean:
	find . -name \*.cov | xargs rm

