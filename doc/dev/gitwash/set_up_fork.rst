.. _set-up-fork:

==================
 Set up your fork
==================

First you follow the instructions for :ref:`forking`.

Overview
========

::

   git clone git@github.com:your-user-name/sfepy.git
   cd sfepy
   git remote add upstream git://github.com/sfepy/sfepy.git

In detail
=========

Clone your fork
---------------

#. Clone your fork to the local computer with ``git clone
   git@github.com:your-user-name/sfepy.git``
#. Investigate.  Change directory to your new repo: ``cd sfepy``. Then
   ``git branch -a`` to show you all branches.  You'll get something
   like::

      * master
      remotes/origin/master

   This tells you that you are currently on the ``master`` branch, and
   that you also have a ``remote`` connection to ``origin/master``.
   What remote repository is ``remote/origin``? Try ``git remote -v`` to
   see the URLs for the remote.  They will point to your github fork.

   Now you want to connect to the upstream `SfePy github`_ repository, so
   you can merge in changes from trunk.

.. _linking-to-upstream:

Linking your repository to the upstream repo
--------------------------------------------

::

   cd sfepy
   git remote add upstream git://github.com/sfepy/sfepy.git

``upstream`` here is just the arbitrary name we're using to refer to the
main `SfePy`_ repository at `SfePy github`_.

Note that we've used ``git://`` for the URL rather than ``git@``.  The
``git://`` URL is read only.  This means we that we can't accidentally
(or deliberately) write to the upstream repo, and we are only going to
use it to merge into our own code.

Just for your own satisfaction, show yourself that you now have a new
'remote', with ``git remote -v show``, giving you something like::

   upstream	git://github.com/sfepy/sfepy.git (fetch)
   upstream	git://github.com/sfepy/sfepy.git (push)
   origin	git@github.com:your-user-name/sfepy.git (fetch)
   origin	git@github.com:your-user-name/sfepy.git (push)

.. include:: links.inc

