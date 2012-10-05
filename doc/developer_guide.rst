.. highlight:: python
   :linenothreshold: 3

.. include:: links.inc

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
     - top level application classes (e.g. :class:`PDESolverApp` that
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
   * - `mesh/`
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

Short version: go to the main `SfePy`_ and follow the links given there.

When you encounter a problem, try searching that site first - an answer may
already be posted in the `SfePy mailing list`_ (to which we suggest you
subscribe...), or the problem might have been added to the `SfePy issues`_ web
page. As is true in any open source project, doing your homework by searching
for existing known problems greatly reduces the burden on the developers by
eliminating duplicate issues. If you find your problem already exists in the
issue tracker, feel free to gather more information and append it to the
issue. In case the problem is not there, create a new issue with proper labels
for the issue type and priority, and/or ask us using the mailing list.

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
the code from the `SfePy git repository`_ page.

We use `git`_ to track source code, documentation, examples, and other files
related to the project.

It is not necessary to learn git in order to contribute to *SfePy* but we
strongly suggest you do so as soon as possible - it is an extremely useful tool
not just for writing code, but also for tracking revisions of articles,
Ph.D. theses, books, ... it will also look well in your CV :-) It is also much
easier for us to integrate changes that are in form of a nice git patch than in
another form.

Having said that, to download the latest snapshot, do either (with git):

- git clone git://github.com/sfepy/sfepy.git

or (without git):

- use the `SfePy tarball`_ link

Then make the changes as you wish, following our :ref:`coding_style`.

**Note** Do not be afraid to experiment - git works with your *local* copy of
the repository, so it is not possible to damage the master repository. It is
always possible to re-clone a fresh copy, in case you do something that is
really bad.

.. _coding_style:

Coding style
^^^^^^^^^^^^

All the code in SfePy should try to adhere to python style guidelines, see
`PEP-0008`_.

There are some additional recommendations:

- Prefer whole words to abbreviations in public APIs - there is completion
  after all. If some abbreviation is needed (*really* too long name), try to
  make it as comprehensible as possible. Also check the code for similar
  names - try to name things consistently with the existing code. Examples:

  - yes: ``equation``, ``transform_variables()``, ``filename``
  - rather not: ``eq``, ``transvar()``, ``fname``

- Functions have usually form ``<action>_<subject>()`` e.g.: ``save_data()``,
  ``transform_variables()``, do not use ``data_save()``,
  ``variable_transform()`` etc.
- Variables like ``V``, ``c``, ``A``, ``b``, ``x`` should be tolerated only
  locally when expressing mathematical ideas.

Really minor recommendations:

- Avoid single letter names, if you can:

  - not even for loop variables - use e.g. ir, ic, ... instead of i, j for rows
    and columns
  - not even in generators, as they "leak" (this is fixed in Python 3.x)

These are recommendations only, we will not refuse code just on the ground that
it uses slightly different formatting, as long as it follows the PEP.

Note: some old parts of the code might not follow the PEP, yet. We fix them
progressively as we update the code.

Contributing changes
^^^^^^^^^^^^^^^^^^^^

Even if you do not use git, try to follow the spirit of :ref:`notes_patches`

Without git
"""""""""""

Without using git, send the modified files to the `SfePy mailing list`_ or
attach them using gist to the corresponding issue at the `Issues`_ web page. Do
not forget to describe the changes properly.

With git
""""""""

.. toctree::
   :hidden:

   dev/gitwash/index

**Note**: This section will get quickly get you started using git and github.
For more in-depth reading about how these tools work with the *SfePy* source
code and the general git development, read :ref:`using-git`, which was adapted
from Matthew Brett's excellent `gitwash`_ git tutorial.

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

#. Send the patch(es) to the `SfePy mailing list`_ or attach them
   to the corresponding issue at the `Issues`_ web page.

#. If the patches are fine, they will appear in the master
   repository. Then synchronize your repository with the master:

   - either clone a fresh copy
   - or use the fetch, pull, merge or rebase commands. This may require
     a deeper git-fu in case of conflicts. For beginners, it is
     advisable to clone always a fresh copy if they see a conflict.

There is another option than submitting patches, however, useful when you wish
to get feedback on a larger set of changes. This option is to publish your
repository using `Github`_ and let the other developers know about it - follow
the instructions in :ref:`git-development` of :ref:`using-git`.

.. _notes_patches:

Notes on commits and patches
""""""""""""""""""""""""""""
- Follow our :ref:`coding_style`.
- Do not use lines longer than 79 characters (exception: tables of
  values, e.g., quadratures).
- Write descriptive docstrings in correct style, see :ref:`docstrings`.
- There should be one patch for one topic - do not mix unrelated things in one
  patch. For example, when you add a new function, then notice a typo in
  docstring in a nearby function and correct it, create two patches: one fixing
  the docstring, the other adding the new function.
- The commit message and description should clearly state what the patch
  does. Try to follow the style of other commit messages. Some interesting
  notes can be found at `tbaggery.com`_, namely that the commit message is
  better to be written in the present tense: "fix bug" and not "fixed bug".

.. _docstrings:

Docstring standard
""""""""""""""""""

We use `sphinx`_ with the `numpydoc`_ extension to generate this
documentation. Refer to the sphinx site for the possible markup constructs.

Basically (with a little tweak), we try to follow the NumPy/SciPy docstring
standard as described in `NumPy documentation guide`_. See also the complete
`docstring example`_. It is exaggerated a bit to show all the
possibilities. Use your common sense here - the docstring should be sufficient
for a new user to use the documented object. A good way to remember the format
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

#. Install `sphinx`_ and `numpydoc`_. Do not forget to set the path to numpydoc
   in site_cfg.py if it is not installed in a standard location for Python
   packages on your platform. A recent :math:`\mbox{\LaTeX}` distribution is
   required, too, for example `TeX Live`_. Depending on your OS/platform, it
   can be in the form of one or several packages.

#. Edit the rst files in `doc/` directory using your favorite text editor - the
   ReST format is really simple, so nothing fancy is needed. Follow the
   existing files in `doc/`; for reference also check `reStructuredText
   Primer`_, `Sphinx Markup Constructs`_ and `docutils reStructuredText`_.

   - When adding a new Python module, add a corresponding documentation file
     into `doc/src/sfepy/<path>`, where `<path>` should reflect the location of
     the module in `sfepy/`.

   - Figures belong to `doc/images`; subdirectories can be used.

#. (Re)generate the documentation (assuming GNU make is installed)::

    cd doc
    make html

#. View it (substitute your favorite browser)::

    firefox _build/html/index.html

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

Not all terms support all the modes, one usually needs to look at the
sources. There are, however, certain naming conventions:

- `'dw_*'` terms support `'weak'` mode
- `'dq_*'` terms support `'qp'` mode
- `'d_*'`, `'di_*'` terms support `'eval'` mode
- `'ev_*'` terms support `'eval'`, `'el_avg'` and `'qp'` modes

Note that the naming prefixes are due to history when the `mode` argument to
`Term.evaluate()` was not available. Now they are often redundant, but at least
one has a notion what is the evaluation purpose of each term. They may
disappear after some more term unification - "easier_terms" branch already
resulted in a number of terms disappearing.

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
            bf_t = nm.tile(bf.transpose((0, 1, 3, 2)), (out.shape[0], 1, 1, 1))
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
implementation. Do not hesitate to ask the `SfePy mailing list`_, and look at
the source code of the already implemented terms.

How To Make a Release
---------------------

.. toctree::
   :maxdepth: 2

   release_tasks.rst

.. _using-git:

Working with *SfePy* source code
--------------------------------

This section was adapted from Matthew Brett's excellent `gitwash`_ git
tutorial. It complements the above sections and details several aspects of
working with Git and Github.

It can be updated by running::

   $ curl -O https://raw.github.com/matthew-brett/gitwash/master/gitwash_dumper.py
   $ python gitwash_dumper.py doc/dev SfePy --repo-name=sfepy --github-user=sfepy --project-url=http://sfepy.org --project-ml-url=http://groups.google.com/group/sfepy-devel

in the SfePy source directory. Do not forget to delete the section title
in `doc/dev/gitwash/index.rst`, as it is already here.

.. toctree::
   :maxdepth: 2

   dev/gitwash/index

Module Index
------------

Main scripts
^^^^^^^^^^^^
.. toctree::
   :maxdepth: 2

   src/extractor
   src/homogen
   src/phonon
   src/postproc
   src/probe
   src/runTests
   src/schroedinger
   src/shaper
   src/simple

Utility scripts
^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/build_helpers
   src/findSurf
   src/genPerMesh
   src/plotPerfusionCoefs
   src/test_install

   src/script/blockgen
   src/script/config
   src/script/convert_mesh
   src/script/cylindergen
   src/script/edit_identifiers
   src/script/evalForms
   src/script/eval_tl_forms
   src/script/gen_gallery
   src/script/gen_lobatto_pyx
   src/script/gen_term_table
   src/script/plot_condition_numbers
   src/script/save_basis
   src/script/show_authors
   src/script/sync_module_docs

sfepy package
^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 2

   src/sfepy/config
   src/sfepy/version

sfepy.applications package
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 2

   src/sfepy/applications/application
   src/sfepy/applications/pde_solver_app

sfepy.base package
^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/base/base
   src/sfepy/base/compat
   src/sfepy/base/conf
   src/sfepy/base/getch
   src/sfepy/base/goptions
   src/sfepy/base/ioutils
   src/sfepy/base/log
   src/sfepy/base/parse_conf
   src/sfepy/base/plotutils
   src/sfepy/base/progressbar
   src/sfepy/base/reader
   src/sfepy/base/tasks
   src/sfepy/base/testing

sfepy.fem package
^^^^^^^^^^^^^^^^^

WARNING: The code in the fem package is undergoing rapid change. It is best to
refer directly to the code base until the code stabilizes.

.. toctree::
   :maxdepth: 2

   src/sfepy/fem/conditions
   src/sfepy/fem/dof_info
   src/sfepy/fem/domain
   src/sfepy/fem/equations
   src/sfepy/fem/evaluate
   src/sfepy/fem/evaluate_variable
   src/sfepy/fem/facets
   src/sfepy/fem/fe_surface
   src/sfepy/fem/fea
   src/sfepy/fem/fields_base
   src/sfepy/fem/fields_hierarchic
   src/sfepy/fem/fields_nodal
   src/sfepy/fem/functions
   src/sfepy/fem/geometry_element
   src/sfepy/fem/global_interp
   src/sfepy/fem/history
   src/sfepy/fem/linearizer
   src/sfepy/fem/integrals
   src/sfepy/fem/mappings
   src/sfepy/fem/mass_operator
   src/sfepy/fem/materials
   src/sfepy/fem/mesh
   src/sfepy/fem/meshio
   src/sfepy/fem/parseEq
   src/sfepy/fem/parseReg
   src/sfepy/fem/periodic
   src/sfepy/fem/poly_spaces
   src/sfepy/fem/probes
   src/sfepy/fem/problemDef
   src/sfepy/fem/projections
   src/sfepy/fem/quadratures
   src/sfepy/fem/refine
   src/sfepy/fem/region
   src/sfepy/fem/simplex_cubature
   src/sfepy/fem/state
   src/sfepy/fem/utils
   src/sfepy/fem/variables
   src/sfepy/fem/extmods/_fmfield
   src/sfepy/fem/extmods/assemble
   src/sfepy/fem/extmods/bases
   src/sfepy/fem/extmods/lobatto
   src/sfepy/fem/extmods/lobatto_template
   src/sfepy/fem/extmods/mappings
   src/sfepy/fem/extmods/mesh

sfepy.mesh package
^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/mesh/femlab
   src/sfepy/mesh/geom_tools
   src/sfepy/mesh/mesh_generators
   src/sfepy/mesh/mesh_tools
   src/sfepy/mesh/meshutils
   src/sfepy/mesh/splinebox

sfepy.homogenization package
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/homogenization/band_gaps_app
   src/sfepy/homogenization/coefficients
   src/sfepy/homogenization/coefs_base
   src/sfepy/homogenization/coefs_elastic
   src/sfepy/homogenization/coefs_perfusion
   src/sfepy/homogenization/coefs_phononic
   src/sfepy/homogenization/convolutions
   src/sfepy/homogenization/engine
   src/sfepy/homogenization/homogen_app
   src/sfepy/homogenization/micmac
   src/sfepy/homogenization/recovery
   src/sfepy/homogenization/utils

sfepy.interactive package
^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/interactive/session

sfepy.linalg package
^^^^^^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 2

   src/sfepy/linalg/eigen
   src/sfepy/linalg/geometry
   src/sfepy/linalg/sparse
   src/sfepy/linalg/utils
   src/sfepy/linalg/extmods/crcm

sfepy.mechanics package
^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/mechanics/elastic_constants
   src/sfepy/mechanics/friction
   src/sfepy/mechanics/matcoefs
   src/sfepy/mechanics/membranes
   src/sfepy/mechanics/tensors
   src/sfepy/mechanics/units

sfepy.optimize package
^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/optimize/freeFormDef
   src/sfepy/optimize/shapeOptim

sfepy.physics package
^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/physics/energy
   src/sfepy/physics/potentials
   src/sfepy/physics/radial_mesh
   src/sfepy/physics/schroedinger_app

sfepy.postprocess package
^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/postprocess/dataset_manager
   src/sfepy/postprocess/domain_specific
   src/sfepy/postprocess/plot_dofs
   src/sfepy/postprocess/plot_facets.rst
   src/sfepy/postprocess/sources
   src/sfepy/postprocess/time_history
   src/sfepy/postprocess/utils
   src/sfepy/postprocess/viewer

sfepy.solvers package
^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/solvers/eigen
   src/sfepy/solvers/ls
   src/sfepy/solvers/nls
   src/sfepy/solvers/optimize
   src/sfepy/solvers/oseen
   src/sfepy/solvers/petsc_worker
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
   src/sfepy/terms/termsNavierStokes
   src/sfepy/terms/termsPiezo
   src/sfepy/terms/termsPoint
   src/sfepy/terms/termsSurface
   src/sfepy/terms/termsVolume
   src/sfepy/terms/terms_constraints
   src/sfepy/terms/terms_dot
   src/sfepy/terms/terms_fibres
   src/sfepy/terms/terms_hyperelastic_base
   src/sfepy/terms/terms_hyperelastic_tl
   src/sfepy/terms/terms_hyperelastic_ul
   src/sfepy/terms/terms_membrane
   src/sfepy/terms/terms_new
   src/sfepy/terms/terms_th
   src/sfepy/terms/utils

   src/sfepy/terms/extmods/terms
