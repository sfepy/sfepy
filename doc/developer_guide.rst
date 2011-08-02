.. highlight:: python
   :linenothreshold: 3

Developer Guide
===============

.. contents:: Table of Contents
   :local:
   :backlinks: top

This section purports to document the *SfePy* internals. It is mainly useful
for those who wish to contribute to the development of  *SfePy* and understand
the inner workings of the code.

SfePy Directory Structure
-------------------------

Here we list and describe the directories that are in the main sfepy
directory.

.. list-table:: Top directory structure.
   :widths: 10 90
   :header-rows: 1

   * - name
     - description
   * - `build/`
     - directory created by the build process (generated)
   * - `doc/`
     - source files of this documentation
   * - `examples/`
     - example problem description files
   * - `meshes/`
     - finite element mesh files in various formats shared by the examples
   * - `output/`
     - default output directory for storing results of the examples
   * - `output-tests/`
     - output directory for tests
   * - `script/`
     - various small scripts (simple mesh generators, mesh format
       convertors etc.)
   * - `sfepy/`
     - the source code
   * - `tests/`
     - the tests run by `runTests.py`
   * - `tmp/`
     - directory for temporary files (generated)

New users/developers (after going through the :ref:`sec-tutorial`)
should explore the `examples/` directory. For developers, the principal
directory is `sfepy/`, which has the following contents:

.. list-table:: `sfepy/` directory structure.
   :widths: 10 80 10
   :header-rows: 1

   * - name
     - description
     - field-specific
   * - `applications/`
     - top level application classes (e.g. :class:`SimpleApp` that
       implements all that `simple.py` script does)
     -
   * - `base/`
     - common utilities and classes used by most of the other modules
     -
   * - `fem/`
     - the finite element core: modules taking care of boundary
       conditions, degrees of freedom, approximations, variables,
       equations, meshes, regions, quadratures, etc.
     -
   * - `geom/`
     - some utilities to interface with tetgen and triangle mesh generators
     -
   * - `homogenization/`
     - the homogenization engine and supporting modules - highly
       specialized code, one of the reasons of *SfePy* existence
     - *
   * - `interactive/`
     - setup of IPython-based shell `isfepy`
     -
   * - `linalg/`
     - linear algebra functions not covered by NumPy and SciPy
     -
   * - `mechanics/`
     - modules for (continuum) mechanics: elastic constant
       conversions, tensor, units utilities, etc.
     - *
   * - `optimize/`
     - modules for shape optimization based on free-form deformation
     - *
   * - `physics/`
     - small utilities for quantum physics (`schroedinger.py`)
     - *
   * - `postprocess/`
     - Mayavi-based post-processing modules (`postproc.py`)
     -
   * - `solvers/`
     - interface classes to various internal/external solvers (linear,
       nonlinear, eigenvalue, optimization, time stepping)
     -
   * - `terms/`
     - implementation of the terms (weak formulation integrals), see
       :ref:`term_overview`
     -

The directories in the "field-specific" column are mostly interesting
for specialists working in the respective fields.

The `fem/` is the heart of the code, while the `terms/` contains the
particular integral forms usable to build equations - new term writers
should look there.

How to Contribute
-----------------

Read this section if you wish to contribute some work to the *SfePy* project.
Contributions can be made in a variety of forms, not just code. Reporting bugs
and contributing to the documentation, tutorials, and examples is in great
need!

Below we describe

#. where to report or find current problems, issues, and suggestions of
   particular topics for additional development
#. what to do to apply changes/fixes
#. what to do after you made your changes/fixes

Reporting problems
^^^^^^^^^^^^^^^^^^

*Reporting a bug is the first way in which to contribute to an open source
project*

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

- git clone git://github.com/sfepy/sfepy.git

or (without git):

- click this link: http://github.com/sfepy/sfepy/tarball/master

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

.. toctree::
   :hidden:

   dev/gitwash/index

**Note**: This section will get quickly get you started using git and github.
For more in-depth reading about how these tools work with the *SfePy* source
code and the general git development, read :ref:`using-git`, which was adapted
from Matthew Brett's excellent `git tutorial
<http://github.com/matthew-brett/gitwash>`_.

With git there are some additional options for how to send changes to *SfePy*.
Before listing them, let us describe a typical development session and the
related git commands:

#. Either clone a fresh copy by::

     git clone git://github.com/sfepy/sfepy.git

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

**Warning** Implementing a new term usually involves C. As Cython is now
supported by our build system, it should not be that
difficult. Python-only terms are possible as well.

Notes on terminology
^^^^^^^^^^^^^^^^^^^^

*Volume* refers to the whole domain (in space of dimension :math:`d`), while
*surface* to a subdomain of dimension :math:`d-1`, for example a part of the
domain boundary. So in 3D problems volume = volume, surface = surface, while in
2D volume = area, surface = curve.

Introduction
^^^^^^^^^^^^

A term in *SfePy* usually corresponds to a single integral term in (weak)
integral formulation of an equation. Both volume and surface integrals are
supported. There are three types of arguments a term can have:

- *variables*, i.e. the unknown, test or parameter variables declared by the
  `variables` keyword, see :ref:`sec-problem-description-file`,
- *materials*, corresponding to material and other parameters (functions)
  that are known, declared by the `materials` keyword,
- *user data* - anything, but user is responsible for passing them to the
  evaluation functions.

Terms come in two flavors:

- standard terms are subclasses of :class:`sfepy.terms.terms.Term`
- *new* terms are subclasses of :class:`sfepy.terms.terms_new.NewTerm`

As new terms are now not much more than a highly experimental proof of
concept, we will focus on the standard terms here.

The purpose of a standard term class is to implement a (vectorized)
function that assembles the term contribution to residual/matrix and/or
evaluates the term integral in a group of elements simultaneously. Most
such functions are currently implemented in C, but some terms are pure
Python, vectorized using NumPy. A term with a C function needs to be
able to extract the real data from its arguments and then pass those
data to the C function.

Evaluation modes
^^^^^^^^^^^^^^^^

A term can support several evaluation modes:

- `'eval'` : evaluate the integral over a region, result has dimension
  like the quantity integrated;

- `'el_avg'` : element average - result is array of the quantity
  averaged in each element of a region - this is the mode for
  postprocessing;

- `'qp'` : quantity interpolated into quadrature points of each element
  in a region;

- `'weak'` : assemble either the vector or matrix depending on
  `diff_var` argument (the name of variable to differentiate with
  respect to).

Currently, not all terms support all the modes, one needs to look at the
sources. There are, however, certain naming conventions:

- `'dw_*'` terms support `'weak'` mode
- `'dq_*'` term support `'eval'` mode
- `'de_*'` term support `'el_avg'` mode

Actually most `'dq_*'`, `'de_*'`, `'di_*'`, `'d_'` terms support `'eval'`,
`'el_avg'` and `'qp'` modes.

Note that the naming prefixes are due to history when the `mode`
argument to `Term.evaluate()` was not available. Now they are mostly
redundant, but at least one has a notion what is the evaluation purpose
of each term. They may disappear after some more term
unification. easier_terms branch already resulted in a number of terms
disappearing.

Basic attributes
^^^^^^^^^^^^^^^^

A term class should inherit from :class:`sfepy.terms.terms.Term` base
class. The simplest possible term with volume integration and 'weak'
evaluation mode needs to have the following attributes and methods:

- docstring (not really required per se, but we require it);
- `name` attribute - the name to be used in `equations`;
- `arg_types` attribute - the types of arguments the term accepts;
- `integration` attribute, optional - the kind of integral the term
  implements, one of `'volume'` (the default, if not given), `'surface'` or
  `'surface_extra'`;
- `function()` static method - the assembling function;
- `get_fargs()` method - the method that takes term arguments and
  converts them to arguments for `function()`.

Argument types
""""""""""""""

The argument types can be ("[_*]" denotes an optional suffix):

- `'material[_*]'` for a material parameter, i.e. any function that can
  be can evaluated in quadrature points and that is not a variable.
- `'virtual'` for a virtual (test) variable (no value defined), `'weak'`
  evaluation mode;
- `'state[_*]'` for state (unknown) variables (have value), `'weak'`
  evaluation mode;
- `'parameter[_*]'` for parameter variables (have known value), any
  evaluation mode.

Only one `'virtual'` variable is allowed in a term.

Integration kinds
"""""""""""""""""

The integration kinds have the following meaning:

- `'volume'` for volume integral over a region that contains elements;
  uses volume element connectivity for assembling;
- `'surface'` for surface integral over a region that contains faces;
  uses surface face connectivity for assembling;
- `'surface_extra'` for surface integral over a region that contains
  faces; uses volume element connectivity for assembling - this is
  needed if full gradients of a variable are required on the boundary.

`function()`
""""""""""""

The `function()` static method has always the following arguments::

    out, *args

where `out` is the already preallocated output array (change it in
place!) and `*args` are any other arguments the function requires. These
function arguments have to be provided by the `get_fargs()` method. The
function returns zero `status` on success, nonzero on failure.

The `out` array has shape `(n_el, 1, n_row, n_col)`, where `n_el` is the
number of elements in a group and `n_row`, `n_col` are matrix dimensions
of the value on a single element.

`get_fargs()`
"""""""""""""

The `get_fargs()` method has always the same structure of arguments:

- positional arguments corresponding to `arg_types` attribute:

  - example for a typical weak term:

    - for::

        arg_types = ('material', 'virtual', 'state')

      the positional arguments are::

        material, virtual, state

- keyword arguments common to all terms::

    mode=None, term_mode=None, diff_var=None, **kwargs

  here:

  - `mode` is the actual evaluation mode, default is `'eval'`;
  - `term_mode` is an optional term sub-mode influencing what the term
    should return (example: `dw_tl_he_neohook` term has 'strain' and
    'stress' evaluation sub-modes);
  - `diff_var` is taken into account in the `'weak'` evaluation mode. It
    is either `None` (residual mode) or a name of variable with respect
    to differentiate to (matrix mode);
  - `**kwargs` are any other arguments that the term supports.

The `get_fargs()` method returns arguments for `function()`.

Example
^^^^^^^

Let us now discuss the implementation of a simple weak term
`dw_volume_integrate_w` defined as :math:`\int_\Omega c q`, where
:math:`c` is a weight (material parameter) and :math:`q` is a virtual
variable. This term is implemented as follows::

    class IntegrateVolumeOperatorWTerm(Term):
        r"""
        :Description:
        Volume integral of a test function weighted by a scalar function
        :math:`c`.


        :Definition:
        .. math::
            \int_\Omega c q

        :Arguments:
            material : :math:`c`,
            virtual  : :math:`q`
        """
        name = 'dw_volume_integrate_w'
        arg_types = ('material', 'virtual')

        @staticmethod
        def function(out, mat, bf, geo):
            bf_t = nm.tile(bf.transpose((0, 2, 1)), (out.shape[0], 1, 1, 1))
            bf_t = nm.ascontiguousarray(bf_t)
            status = geo.integrate(out, mat * bf_t)

            return status

        def get_fargs(self, mat, virtual,
                      mode=None, term_mode=None, diff_var=None, **kwargs):
            assert_(virtual.n_components == 1)
            geo, _ = self.get_mapping(virtual)

            return mat, geo.bf, geo

- lines 2-15: the docstring - always write one!
- line 16: the name of the term, that can be referred to in equations;
- line 17: the argument types - here the term takes a single material
  parameter, and a virtual variable;
- lines 19-25: the term function

  - its arguments are:

    - the output array `out`, already having the required shape,
    - the material coefficient (array) `mat` evaluated in physical
      quadrature points of all elements of an element group,
    - a base function (array) `bf` evaluated in the quadrature points of
      a reference element and
    - a reference element (geometry) mapping `geo`.

  - line 21: transpose the base function and tile it so that is has
    the correct shape - it is repeated for each element;
  - line 22: ensure C contiguous order;
  - line 23: perform numerical integration in C - `geo.integrate()`
    requires the C contiguous order;
  - line 25: return the status.

- lines 27-32: prepare arguments for the function above:

  - line 29: verify that the variable is scalar, as our implementation
    does not support vectors;
  - line 30: get reference element mapping corresponding to the virtual
    variable;
  - line 32: return the arguments for the function.

Concluding remarks
^^^^^^^^^^^^^^^^^^

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
   src/sfepy/fem/fields
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
   src/sfepy/fem/projections
   src/sfepy/fem/quadratures
   src/sfepy/fem/state
   src/sfepy/fem/variables

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
   src/sfepy/terms/terms_new

sfepy.terms package - full inheritance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section is the same as the previous, but shows the full inheritance for
some of the terms classes.

.. toctree::
   :maxdepth: 2
   
   src/sfepy/terms/termsLinElasticity_full
