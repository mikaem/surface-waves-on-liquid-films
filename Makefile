all:
	jupyter nbconvert surfacewaves.ipynb --to html
	python postprocess.py
	cp surfacewaves.html docs/index.html
