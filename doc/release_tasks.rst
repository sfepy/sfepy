Release Tasks
=============

A few notes on what to do during a release.

Things to check before a release
--------------------------------

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

     $ python setup.py pdfdocs

#. check installed build::

     $ python setup.py install --root=<some path>
     $ cd
     $ runTests.py
     $ rm -r output/

   then remove the installed files so that they do not interfere with
   the local build

#. create final tarball

   * change version number (sfepy/version.py) so that previous release
     tarball is not overwritten!
   * set ``is_release = True`` in site_cfg.py
   * update pdfdocs::

     $ python setup.py pdfdocs

   * create tarball::

     $ python setup.py sdist

Useful Git commands
-------------------

* log ::

    git log --pretty=format:"%s%n%b%n" release_2009.1..HEAD
    git shortlog -s -n release_2010.1..HEAD

* who has contributed since <date>::

    git log --after=<date> | grep Author | sort | uniq
    git log release_2009.1..HEAD | grep Author | sort | uniq

    git rev-list --committer="Name Surname" --since=6.months.ago HEAD | wc
    git rev-list --author="Name Surname" --since=6.months.ago HEAD | wc
    # ?no-merges

* misc::

    git archive --format=tar HEAD | gzip > name.tar.gz

Web update and file uploading
-----------------------------

* publish development docs as new release docs:

* add RELEASE_NOTES.txt as <version>_RELEASE_NOTES.txt via SVN to the
  Google site

* upload the tarball to http://code.google.com/p/sfepy/downloads/list

  * make it featured, un-feature the previous release
  * update download link at http://code.google.com/p/sfepy/wiki/Downloads

* update news at http://code.google.com/p/sfepy/

  * update wiki page http://code.google.com/p/sfepy/wiki/FrontPage

    * archive last news entry to
      http://code.google.com/p/sfepy/wiki/ArchivedNews

  * administer...

* send announcement to

  * sfepy-devel@googlegroups.com, scipy-dev@scipy.org,
    scipy-user@scipy.org, python-announce-list@python.org

  * ntc@list.zcu.cz, kme@list.zcu.cz
