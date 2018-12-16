Release Tasks
=============

A few notes on what to do during a release.

Things to check before a release
--------------------------------

#. synchronize module documentation (dry run)::

     $ ./script/sync_module_docs.py doc/src/ . -n

#. regenerate gallery page and examples::

    $ script/gen_gallery.py -l ../doc-devel
    $ rm -rf doc/examples/
    $ cp -a gallery/examples/ doc/

#. create temporary/testing tarball::

     $ python setup.py sdist

#. check in-place build::

     $ # unpack the tarball
     $ # cd into

     $ python setup.py build_ext --inplace
     $ ./test_install.py

#. check that documentation can be built::

     $ # copy site_cfg.py
     $ python setup.py htmldocs
     $ firefox doc/_build/html/index.html

   or use::

     $ cd doc/
     $ make html
     $ firefox _build/html/index.html

   try also::

     $ # copy gallery/images
     $ python setup.py pdfdocs

#. check installed build::

     $ pip install . --user
     $ cd
     $ sfepy-run run_tests
     $ rm -r output/

   then remove the installed files so that they do not interfere with
   the local build

#. create final tarball

   * update doc/release_notes.rst, with the help of::

     $ python script/gen_release_notes.py 2018.2

   * update doc/news.rst, doc/archived_news.rst
   * change version number (sfepy/version.py) so that previous release
     tarball is not overwritten!
   * set ``is_release = True`` in site_cfg.py
   * update pdfdocs::

     $ python setup.py pdfdocs

   * create tarball::

     $ python setup.py sdist

#. tag the release using::

     $ git tag release_XXXX.X

Useful Git commands
-------------------

* log ::

    git log --pretty=format:"%s%n%b%n" --topo-order --reverse release_2016.4..HEAD

* who has contributed since <date>::

    git log --after=<date> | grep Author | sort | uniq
    git log release_2012.1..HEAD | grep Author | sort -k3 | uniq
    git shortlog -s -n release_2012.3..HEAD

    git rev-list --committer="Name Surname" --since=6.months.ago HEAD | wc
    git rev-list --author="Name Surname" --since=6.months.ago HEAD | wc
    # ?no-merges

* misc::

    git archive --format=tar HEAD | gzip > name.tar.gz

Web update and file uploading
-----------------------------

* upload the tarball to https://sfepy.org/doc-devel/downloads.html

* make a pull request with the updated version in
  ``sfepy-feedstock/recipe/meta.yaml`` from a fork
  (e.g. https://github.com/rc/sfepy-feedstock) of
  https://github.com/conda-forge/sfepy-feedstock.

* publish development docs also as new release docs

* send announcement to

  * sfepy@python.org, scipy-dev@python.org,
    scipy-user@python.org, numpy-discussion@python.org,
    python-announce-list@python.org
