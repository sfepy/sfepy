
all:
	python setup.py build_ext --inplace
clean:
	python setup.py clean
htmldocs:
	python setup.py htmldocs
pdfdocs:
	python setup.py pdfdocs
doxygendocs:
	python setup.py doxygendocs