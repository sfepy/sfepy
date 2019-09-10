
PYTHON=python$(shell grep python_version site_cfg.py | cut -c 19)

all clean: site_cfg.py
site_cfg.py:
	cp site_cfg_template.py site_cfg.py
	$(info Set correct Python version in site_cfg.py before running `make`.)
	$(info ----------------------------------------------------------------)

all:
	$(PYTHON) setup.py build_ext --inplace
clean:
	$(PYTHON) setup.py clean
htmldocs:
	$(PYTHON) setup.py htmldocs
pdfdocs:
	$(PYTHON) setup.py pdfdocs
doxygendocs:
	$(PYTHON) setup.py doxygendocs
