all: test book apidoc

guide:
	cd doc/guide && quarto render --to html

guide-publish:
	cd doc/guide && quarto publish gh-pages --no-prompt

apidoc:
	cd doc/apidoc && julia --project make.jl

.PHONY: test
test:
	julia --project test/runtests.jl

clean:
	find . -name \*.cov | xargs rm
	rm -rf demo/_book demo/.quarto demo/_freeze

