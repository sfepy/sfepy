.. highlight:: python
   :linenothreshold: 3

Developer Guide
===============

This section purports to document the *SfePy* internals. It is mainly useful for those who wish to develop *SfePy* and understand the inner workings of the code.

How to Contribute
-----------------

Read this section if you wish to contribute some work to the *SfePy* project. Contributions can be made in a variety of forms, not just code. Reporting bugs and contributing to the documentation, tutorials, and examples is in great need!

Below we describe

#. where to report or find current problems, issues, and suggestions of
   particular topics for additional development
#. what to do to apply changes/fixes
#. what to do after you made your changes/fixes

Reporting problems
^^^^^^^^^^^^^^^^^^

*Reporting a bug is the first way in which to contribute to an open source project*

We use the site http://sfepy.org to provide a hub for the developers to post
problems, ask questions, create wiki pages, etc. The address is currently just
an alias to the Google code site http://code.google.com/p/sfepy.

When you encounter a problem, try searching that site first - an answer may
already be posted in the `sfepy-devel
<http://groups.google.com/group/sfepy-devel>`_ mailing list (to which we
suggest you subscribe...), or the problem might have been added to the `Issues
<http://code.google.com/p/sfepy/issues/list>`_ web page. As is true in any open
source project, doing your homework by searching for existing known problems
greatly reduces the burden on the developers by eliminating duplicate issues.
If you find your problem already exists in the issue tracker, feel free to
gather more information and append it to the issue. In case the problem is not
there, create a new issue with proper labels for the issue type and priority,
and/or ask us using the mailing list.

**Note** A google account (e.g., gmail account) is needed to join the mailing
list and post comments to issues. It is, however, not needed to create a new
issue.

**Note** When reporting a problem, try to provide as much information as
possible concerning the version of *SfePy*, the OS / Linux distribution, and
the versions of *Python*, *NumPy* and *SciPy*, and other prerequisites.

Our persisting all star top priority issues include:

* missing docstrings in many functions/classes/modules
* incomplete documentation
* lowering the barrier for new users

  * e.g., through generation of additional tutorial material

So if you are a new user, please let us know what difficulties you have with
this documentation. We greatly welcome a variety of contributions not limited
to code only.

Making changes
^^^^^^^^^^^^^^

This step is simple, just keep in mind to use the latest development version of
the code from the `downloads tab
<http://code.google.com/p/sfepy/wiki/Downloads?tm=2>`_ at the developers'
site - the git repository (*not* the latest released sources).

We use `git <http://git-scm.com/>`_ to track source code, documentation,
examples, and other files related to the project.

It is not necessary to learn git in order to contribute to *SfePy* but we
strongly suggest you do so as soon as possible - it is an extremely useful tool
not just for writing code, but also for tracking revisions of articles,
Ph.D. theses, books, ... it will also look well in your CV :-) It is also much
easier for us to integrate changes that are in form of a nice git patch than in
another form.

Having said that, to download the latest snapshot, do either (with git):

- git clone git://git.sympy.org/sfepy.git

or (without git):

- click this link: http://git.sympy.org/?p=sfepy.git;a=snapshot;h=HEAD;sf=tgz

Then make the changes as you wish, following our `style guide
<http://code.google.com/p/sfepy/wiki/CodingStyle>`_.

**Note** Do not be afraid to experiment - git works with your *local* copy of
the repository, so it is not possible to damage the master repository. It is
always possible to re-clone a fresh copy, in case you do something that is
really bad.

Contributing changes
^^^^^^^^^^^^^^^^^^^^

Even if you do not use git, try to follow the spirit of :ref:`notes_patches`

Without git
"""""""""""

Without using git, send the modified files to the `sfepy-devel
<http://groups.google.com/group/sfepy-devel>`_ mailing list or attach them to
the corresponding issue at the `Issues
<http://code.google.com/p/sfepy/issues/list>`_ web page. Do not forget to
describe the changes properly.

With git
""""""""

**NOTE**: This section is somewhat superceded by the `Contributing to SfePy` documentation, which was shamelessly borrowed from the *Numpy* project. Read this from the links below.

.. toctree::
   :maxdepth: 2

   dev/index

With git there are some additional options. Before listing them, let us
describe a typical development session and the related git commands:

#. Either clone a fresh copy by::

     git clone git://git.sympy.org/sfepy.git

   or update your local repository::

     # look for changes at origin
     git fetch origin

     # difference between local and origin master branch
     git diff master origin/master

     # apply the changes to local master branch
     git pull origin master

#. Introduce yourself to git and make (optionally) some handy aliases
   either in ``.gitconfig`` in your home directory (global setting for
   all your git projects), or directly in ``.git/config`` in the repository::

     [user]
         email = mail@mail.org
         name = Name Surname

     [color]
         ui = auto
         interactive = true

     [alias]
         ci = commit
         di = diff --color-words
         st = status
         co = checkout

#. Change some file(s), and review the changes::

     # text diff
     git diff

     # use GUI to visualize of project history (all branches)
     gitk --all

#. Create one or more commits::

    # schedule some of the changed files for the next commit
    git add file1 file2 ...
    # an editor will pop up where you should describe the commit
    git commit

#. The commit(s) now reflect changes, but only in your *local* git
   repository. Then you must somehow allow others to see them. This can be done,
   for example, by sending a patch (or through the other option below). So
   create the patch(es)::

    # create patches for, e.g., the last two commits
    git format-patch HEAD~2

#. Send the patch(es) to the `sfepy-devel
   <http://groups.google.com/group/sfepy-devel>`_ mailing list or attach them
   to the corresponding issue at the `Issues
   <http://code.google.com/p/sfepy/issues/list>`_ web page.

#. If the patches are fine, they will appear in the master
   repository. Then synchronize your repository with the master:

   - either clone a fresh copy
   - or use the fetch, pull, merge or rebase commands. This may require
     a deeper git-fu in case of conflicts. For beginners, it is
     advisable to clone always a fresh copy if they see a conflict.

There is another option than submitting patches, however, useful when you wish
to get feedback on a larger set of changes. This option is to publish your
repository at a free git hosting web site like `Github <http://github.com/>`_
and let the other developers know about it. For example, Robert usually
publishes fixes to issues at http://github.com/rc/sfepy for review, before
pushing them to the main repository.

.. _notes_patches:

Notes on commits and patches
""""""""""""""""""""""""""""
- Follow our `style guide <http://code.google.com/p/sfepy/wiki/CodingStyle>`_.
- Do not use lines longer than 79 characters (exception: tables of
  values, e.g., quadratures).
- Write descriptive docstrings in correct style, see :ref:`docstrings`.
- There should be one patch for one topic - do not mix unrelated things in one
  patch. For example, when you add a new function, then notice a typo in
  docstring in a nearby function and correct it, create two patches: one fixing
  the docstring, the other adding the new function.
- The commit message and description should clearly state what the patch
  does. Try to follow the style of other commit messages. Some interesting
  notes can be found `here
  <http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html>`_,
  namely that the commit message is better to be written in the present tense:
  "fix bug" and not "fixed bug".

.. _docstrings:

Docstring standard
""""""""""""""""""

We use `sphinx <http://sphinx.pocoo.org>`_ with the `numpydoc
<http://pypi.python.org/pypi/numpydoc/0.3.1>`_ extension to generate
this documentation. Refer to the sphinx site for the possible markup
constructs.

Basically (with a little tweak), we try to follow the NumPy/SciPy
docstring standard as described in this `guide
<http://projects.scipy.org/numpy/wiki/CodingStyleGuidelines>`_. See also
the complete `example.py
<http://svn.scipy.org/svn/numpy/trunk/doc/example.py>`_. It is exaggerated
a bit to show all the possibilities. Use your common sense here - the
docstring should be sufficient for a new user to use the documented
object. A good way to remember the format
is to type::

    In [1]: import numpy as nm
    In [2]: nm.sin?

in `ipython`. The little tweak mentioned above is the starting newline::

    def function(arg1, arg2):
        """
	This is a function.

        Parameters
        ----------
        arg1 : array
            The coordinates of ...
        arg2 : int
            The dimension ...

        Returns
        -------
        out : array
           The resulting array of shape ....
        """

It seems visually better than::

    def function(arg1, arg2):
        """This is a function.

        Parameters
        ----------
        arg1 : array
            The coordinates of ...
        arg2 : int
            The dimension ...

        Returns
        -------
        out : array
           The resulting array of shape ....
        """

When using :math:`\mbox{\LaTeX}` in a docstring, use a raw string::

    def function():
        r"""
	This is a function with :math:`\mbox{\LaTeX}` math:
        :math:`\frac{1}{\pi}`.
	"""

to prevent Python from interpreting and consuming the backslashes in common
escape sequences like '\\n', '\\f' etc.

.. _how_to_regenerate_documentation:

How to Regenerate Documentation
-------------------------------

The following steps summarize how to regenerate this documentation.

#. Install `sphinx <http://sphinx.pocoo.org>`_ and `numpydoc
   <http://pypi.python.org/pypi/numpydoc/0.3.1>`_. Do not forget to set the
   path to numpydoc in site_cfg.py if it is not installed in a standard
   location for Python packages on your platform. A recent
   :math:`\mbox{\LaTeX}` distribution is required, too, for example `TeX Live
   <http://www.tug.org/texlive/>`_. Depending on your OS/platform, it can be in
   the form of one or several packages.

#. Edit the rst files in `doc/` directory using your favorite text editor - the
   ReST format is really simple, so nothing fancy is needed. Follow the
   existing files in `doc/`; for reference also check [1]_, [2]_ and [3]_.

   - When adding a new Python module, add a corresponding documentation file
     into `doc/src/sfepy/<path>`, where `<path>` should reflect the location of
     the module in `sfepy/`.

   - Figures belong to `doc/images`; subdirectories can be used.

#. (Re)generate the documentation (assuming GNU make is installed)::

    cd doc
    make html

#. View it (substitute your favorite browser)::

    firefox _build/html/index.html

.. [1] http://sphinx.pocoo.org/rest.html
.. [2] http://sphinx.pocoo.org/markup/index.html
.. [3] http://docutils.sourceforge.net/rst.html

How to Implement a New Term
---------------------------
*tentative documentation*

**Warning 1** Implementing a new term usually involves C.

**Warning 2** It is also more complicated than it should and could be.

We are slowly working to "fix" both problems.

Notes on terminology
^^^^^^^^^^^^^^^^^^^^

'Volume' refers to the whole domain (in space of dimension :math:`d`), while
'surface' to a subdomain of dimension :math:`d-1`, for example a part of the
domain boundary. So in 3D problems volume = volume, surface = surface, while in
2D volume = area, surface = curve.

Introduction
^^^^^^^^^^^^

A term in *SfePy* usually corresponds to a single integral term in (weak)
integral formulation of an equation. Both volume and surface integrals are
supported. There are three types of arguments a term can have:

  - *variables*, i.e. the unknown, test or parameter variables declared by the
    `variables` keyword, see :ref:`sec-problem-description-file`.
  - *materials*, corresponding to material and other parameters (functions)
     that are known, declared by the `materials` keyword
  - *user data* - anything, but user is responsible for passing them to the
     evaluation functions.

**The purpose of a term class** is:

  #. to extract the real data from its arguments
  #. to pass those data to a element matrix/residual evaluation function
     (usually in C)

So a term class basically transforms the arguments to a form suitable for the
actual computation.


Technical details
^^^^^^^^^^^^^^^^^
A term class should inherit from :class:`sfepy.terms.terms.Term` class, and
(optionally) also from classes in :mod:`sfepy.terms.terms_base`. Those classes
provide some common functionality, e.g. the `_call()` methods and data shape
setting. The simplest possible term class, however, just needs to have the
following attributes and methods:

- docstring (not really required per se, but we require it)
- `name` attribute - the name to be used in `equations`
- `arg_types` attribute - the types of arguments the term accepts
- `geometry` attribute - the kind of geometrical data the term needs (usually
  `Volume` or `Surface`)
- `__call__()` method - subclasses of `Term` either implement `__call__()` or
  plug in a proper `_call()` method that is called by the default
  implementation. It takes the following arguments::

      __call__(self, diff_var=None, chunk_size=None, **kwargs)

 - `diff_var` is either None (residual mode), or the name of the
   variable to differentiate with respect to (matrix mode)
 - `chunk_size` is the number of elements that should be processed in one
   call - `__call__()` is a generator that is called in the assembling loop
   in :mod:`sfepy.fem.evaluate`
 - `**kwargs` contains the arguments as described by `arg_types` and any
   other arguments passed by the caller

Let us show how a simple volume integral with the usual three arguments can
look like. Let us assume, that the variables are vector fields and the
evaluation function needs the FE base function. Then::

  class MyTerm(Term):
      r"""
      :Description:
      Some description.
      
      :Definition:
      .. math::
          \int_\Omega \dots
      """
      name = 'dw_my_term'
      arg_types = ('material', 'virtual', 'state')
      geometry = [(Volume, 'virtual')]

      def __call__(self, diff_var=None, chunk_size=None, **kwargs):
          mat, virtual, state = self.get_args(**kwargs)
          ap, vg = virtual.get_approximation(self.get_current_group(), 'Volume')
          n_el, n_qp, dim, n_ep = ap.get_v_data_shape(self.integral)

          if diff_var is None:
              shape = (chunk_size, 1, dim * n_ep, 1)

          elif diff_var == self.get_arg_name( 'state' ):
              shape = (chunk_size, 1, dim * n_ep, dim * n_ep)

          else:
              raise StopIteration

          bf = ap.get_base('v', 0, self.integral)
          for out, chunk in self.char_fun(chunk_size, shape):
              status = self.some_function(out, bf, other_args, chunk)

              yield out, chunk, status

Discussion:

- line 15: extract the three arguments from the argument list
- line 16: get the 'Volume' Approximation and VolumeGeometry instances
  corresponding to the term region and current element group given by the
  `self.get_current_group()` call
- line 17: get data shape:

  - `n_el` .. number of elements in the group
  - `n_qp` .. number of quadrature points in each element
  - `dim`  .. space dimension
  - `n_ep` .. number of element points (=nodes) of each element

- lines 19-26: determine data shape of the chunk, depending whether in residual
  or matrix mode.
- line 28: get volume base function corresponding to the integral used
- lines 29-32: evaluate for each chunk and yield the results

In practice, such a term would inherit also from
:class:`sfepy.terms.terms_base.VectorVector` - then it could look, for example,
like the :class:`sfepy.terms.termsMass.MassTerm`, i.e., just take care of
providing the correct arguments to the evaluation function (`self.function`
attribute)::

  class MassTerm( VectorVector, Term ):
      r"""
      :Description:
      Inertial forces term.
  
      :Definition:
      .. math::
          \int_{\Omega} \rho \ul{v} \cdot \frac{\ul{u} - \ul{u}_0}{\dt}
  
      :Arguments:
      material : :math:`\rho`,
      ts.dt : :math:`\dt`,
      parameter : :math:`\ul{u}_0`"""
      name = 'dw_mass'
      arg_types = ('ts', 'material', 'virtual', 'state', 'parameter')
      geometry = [(Volume, 'virtual')]
  
      def __init__(self, name, sign, **kwargs):
          Term.__init__(self, name, sign, function=terms.dw_mass, **kwargs)
  
      def get_fargs(self, diff_var=None, chunk_size=None, **kwargs):
          ts, mat, virtual, state, state0 = self.get_args(**kwargs)        
          ap, vg = virtual.get_approximation(self.get_current_group(), 'Volume')
  
          self.set_data_shape(ap)
          shape, mode = self.get_shape(diff_var, chunk_size)
  
          dvec = state() - state0()
          rhodt = mat / ts.dt
          bf = ap.get_base('v', 0, self.integral)
  
          fargs = (rhodt, dvec, 0, bf, vg, ap.econn)
          return fargs, shape, mode

This is just a very basic introduction to the topic of new term
implementation. Do not hesitate to ask the `sfepy-devel mailing list
<http://groups.google.com/group/sfepy-devel>`_, and look at the source code of
the already implemented terms.

Module Index
------------

sfepy.applications package
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 2

   src/sfepy/applications/application
   src/sfepy/applications/simple_app
   src/sfepy/applications/top_level

sfepy.base package
^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/base/base
   src/sfepy/base/conf
   src/sfepy/base/getch
   src/sfepy/base/ioutils
   src/sfepy/base/log
   src/sfepy/base/plotutils
   src/sfepy/base/progressbar
   src/sfepy/base/reader
   src/sfepy/base/tasks
   src/sfepy/base/testing

sfepy.fem package
^^^^^^^^^^^^^^^^^

WARNING: The code in the fem package is undergoing rapid change. It is best to
refer directly to the code base until the code stabilizes. The modules listed
below are already more or less updated.

.. toctree::
   :maxdepth: 2

   src/sfepy/fem/conditions
   src/sfepy/fem/dof_info
   src/sfepy/fem/equations
   src/sfepy/fem/evaluate
   src/sfepy/fem/functions
   src/sfepy/fem/geometry_element
   src/sfepy/fem/integrals
   src/sfepy/fem/mesh
   src/sfepy/fem/mesh_generators
   src/sfepy/fem/meshio
   src/sfepy/fem/periodic
   src/sfepy/fem/poly_spaces
   src/sfepy/fem/probes
   src/sfepy/fem/problemDef
   src/sfepy/fem/quadratures

sfepy.geom package
^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/geom/femlab
   src/sfepy/geom/geometry
   src/sfepy/geom/gmsh
   src/sfepy/geom/meshutils
   src/sfepy/geom/tetgen

sfepy.homogenization package
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/homogenization/coefficients
   src/sfepy/homogenization/coefs_base
   src/sfepy/homogenization/coefs_elastic
   src/sfepy/homogenization/convolutions
   src/sfepy/homogenization/engine
   src/sfepy/homogenization/phono
   src/sfepy/homogenization/recovery
   src/sfepy/homogenization/utils

sfepy.interactive package
^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/interactive

sfepy.linalg package
^^^^^^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 2

   src/sfepy/linalg/eigen
   src/sfepy/linalg/geometry
   src/sfepy/linalg/sparse
   src/sfepy/linalg/utils

sfepy.mechanics package
^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/mechanics/matcoefs
   src/sfepy/mechanics/tensors
   src/sfepy/mechanics/units

sfepy.physics package
^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/physics/dft

sfepy.postprocess package
^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/postprocess/dataset_manager
   src/sfepy/postprocess/sources
   src/sfepy/postprocess/time_history
   src/sfepy/postprocess/utils
   src/sfepy/postprocess/viewer

sfepy.solvers package
^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/solvers/eigen
   src/sfepy/solvers/generic
   src/sfepy/solvers/ls
   src/sfepy/solvers/nls
   src/sfepy/solvers/optimize
   src/sfepy/solvers/oseen
   src/sfepy/solvers/semismooth_newton
   src/sfepy/solvers/solvers
   src/sfepy/solvers/ts

sfepy.terms package
^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   terms_overview
   src/sfepy/terms/terms
   src/sfepy/terms/termsAcoustic
   src/sfepy/terms/termsAdjointNavierStokes
   src/sfepy/terms/termsBasic
   src/sfepy/terms/termsBiot
   src/sfepy/terms/termsElectric
   src/sfepy/terms/termsLaplace
   src/sfepy/terms/termsLinElasticity
   src/sfepy/terms/termsMass
   src/sfepy/terms/termsNavierStokes
   src/sfepy/terms/termsPiezo
   src/sfepy/terms/termsPoint
   src/sfepy/terms/termsSurface
   src/sfepy/terms/termsVolume
   src/sfepy/terms/terms_base
   src/sfepy/terms/terms_fibres
   src/sfepy/terms/terms_hyperelastic_base
   src/sfepy/terms/terms_hyperelastic_tl
   src/sfepy/terms/terms_hyperelastic_ul

sfepy.terms package - full inheritance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section is the same as the previous, but shows the full inheritance for
some of the terms classes.

.. toctree::
   :maxdepth: 2
   
   src/sfepy/terms/termsLinElasticity_full
