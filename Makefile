all: test guide apidoc

guide:
	cd doc/guide && quarto render --to html

guide-publish:
	cd doc/guide && quarto publish gh-pages --no-prompt

apidoc:
	cd doc/apidoc && julia --project make.jl

.PHONY: test
test:
	julia --color=yes --project -e 'using Pkg; Pkg.test()'

clean:
	find . -name \*.cov | xargs rm
	rm -rf doc/guide/_book doc/guide/.quarto doc/guide/_freeze
	rm doc/*/Manifest.toml

