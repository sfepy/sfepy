.. highlight:: python
   :linenothreshold: 3

.. include:: links.inc

.. _developer_guide:

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
     - the tests run by `run_tests.py`
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
   * - `discrete/`
     - general classes and modules for describing a discrete problem, taking
       care of boundary conditions, degrees of freedom, approximations,
       variables, equations, meshes, regions, quadratures, etc.

       Discretization-specific classes are in subdirectories:

       - `common/` - common parent classes for discretization-specific classes
       - `fem/` - finite element specific classes
       - `iga/` - isogeometric analysis specific classes
     -
   * - `mesh/`
     - some utilities to interface with tetgen and triangle mesh generators
     -
   * - `homogenization/`
     - the homogenization engine and supporting modules - highly
       specialized code, one of the reasons of *SfePy* existence
     - *
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

.. _how_to_contribute:

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

A term can support several evaluation modes, as described in
:ref:`term_evaluation`.

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
  be can evaluated in quadrature points and that is not a variable;
- `'opt_material[_*]'` for an optional material parameter, that can be left
  out - there can be only one in a term and it must be the first argument;
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

Additional attributes
^^^^^^^^^^^^^^^^^^^^^

These attributes are used mostly in connection with the
`tests/test_term_call_modes.py` test for automatic testing of term calls.

- `arg_shapes` attribute - the possible shapes of term arguments;
- `geometries` attribute - the list of reference element geometries that the
  term supports;
- `mode` attribute - the default evaluation mode.

Argument shapes
"""""""""""""""

The argument shapes are specified using a dict of the following form::

    arg_shapes = {'material' : 'D, D', 'virtual' : (1, 'state'),
                  'state' : 1, 'parameter_1' : 1, 'parameter_2' : 1}

The keys are the argument types listed in the `arg_types` attribute, for
example::

    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'parameter_1', 'parameter_2'))

The values are the shapes containing either integers, or 'D' (for space
dimension) or 'S' (symmetric storage size corresponding to the space
dimension).  For materials, the shape is a string `'nr, nc'` or a single value,
denoting a special-valued term, or `None` denoting an optional material that is
left out. For state and parameter variables, the shape is a single value. For
virtual variables, the shape is a tuple of a single shape value and a
name of the corresponding state variable; the name can be `None`.

When several alternatives are possible, a list of dicts can be used. For
convenience, only the shapes of arguments that change w.r.t. a previous dict
need to be included, as the values of the other shapes are taken from the
previous dict. For example, the following corresponds to a case, where an
optional material has either the shape (1, 1) in each point, or is left out::

    arg_types = ('opt_material', 'parameter')
    arg_shapes = [{'opt_material' : '1, 1', 'parameter' : 1},
                  {'opt_material' : None}]

Geometries
""""""""""

The default that most terms use is a list of all the geometries::

    geometries = ['2_3', '2_4', '3_4', '3_8']

In that case, the attribute needs not to be define explicitly.

Example
^^^^^^^

Let us now discuss the implementation of a simple weak term
`dw_volume_integrate` defined as :math:`\int_\Omega c q`, where :math:`c` is a
weight (material parameter) and :math:`q` is a virtual variable. This term is
implemented as follows::

    class IntegrateVolumeOperatorTerm(Term):
        r"""
        Volume integral of a test function weighted by a scalar function
        :math:`c`.

        :Definition:

        .. math::
            \int_\Omega q \mbox{ or } \int_\Omega c q

        :Arguments:
            - material : :math:`c` (optional)
            - virtual  : :math:`q`
        """
        name = 'dw_volume_integrate'
        arg_types = ('opt_material', 'virtual')
        arg_shapes = [{'opt_material' : '1, 1', 'virtual' : (1, None)},
                      {'opt_material' : None}]

        @staticmethod
        def function(out, material, bf, geo):
            bf_t = nm.tile(bf.transpose((0, 1, 3, 2)), (out.shape[0], 1, 1, 1))
            bf_t = nm.ascontiguousarray(bf_t)
            if material is not None:
                status = geo.integrate(out, material * bf_t)
            else:
                status = geo.integrate(out, bf_t)
            return status

        def get_fargs(self, material, virtual,
                      mode=None, term_mode=None, diff_var=None, **kwargs):
            assert_(virtual.n_components == 1)
            geo, _ = self.get_mapping(virtual)

            return material, geo.bf, geo

- lines 2-14: the docstring - always write one!
- line 15: the name of the term, that can be referred to in equations;
- line 16: the argument types - here the term takes a single material
  parameter, and a virtual variable;
- lines 17-18: the possible argument shapes
- lines 20-28: the term function

  - its arguments are:

    - the output array `out`, already having the required shape,
    - the material coefficient (array) `mat` evaluated in physical
      quadrature points of all elements of an element group,
    - a base function (array) `bf` evaluated in the quadrature points of
      a reference element and
    - a reference element (geometry) mapping `geo`.

  - line 22: transpose the base function and tile it so that is has
    the correct shape - it is repeated for each element;
  - line 23: ensure C contiguous order;
  - lines 24-27: perform numerical integration in C - `geo.integrate()`
    requires the C contiguous order;
  - line 28: return the status.

- lines 30-35: prepare arguments for the function above:

  - line 32: verify that the variable is scalar, as our implementation
    does not support vectors;
  - line 33: get reference element mapping corresponding to the virtual
    variable;
  - line 35: return the arguments for the function.

Concluding remarks
^^^^^^^^^^^^^^^^^^

This is just a very basic introduction to the topic of new term
implementation. Do not hesitate to ask the `SfePy mailing list`_, and look at
the source code of the already implemented terms.

How To Make a Release
---------------------

.. toctree::
   :maxdepth: 2

   release_tasks

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
   src/run_tests
   src/schroedinger
   src/shaper
   src/simple

Utility scripts
^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/build_helpers
   src/test_install

   src/script/blockgen
   src/script/convert_mesh
   src/script/cylindergen
   src/script/edit_identifiers
   src/script/eval_ns_forms
   src/script/eval_tl_forms
   src/script/extract_surface
   src/script/gen_gallery
   src/script/gen_iga_patch
   src/script/gen_lobatto1d_c
   src/script/gen_mesh_prev
   src/script/gen_term_table
   src/script/plot_condition_numbers
   src/script/plot_logs
   src/script/plot_mesh
   src/script/plot_quadratures
   src/script/plot_times
   src/script/save_basis
   src/script/show_authors
   src/script/show_terms_use
   src/script/sync_module_docs
   src/script/tile_periodic_mesh

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
   src/sfepy/base/log_plotter
   src/sfepy/base/mem_usage
   src/sfepy/base/parse_conf
   src/sfepy/base/plotutils
   src/sfepy/base/reader
   src/sfepy/base/resolve_deps
   src/sfepy/base/testing

sfepy.discrete package
^^^^^^^^^^^^^^^^^^^^^^

This package implements various PDE discretization schemes (FEM or IGA).

.. toctree::
   :maxdepth: 2

   src/sfepy/discrete/conditions
   src/sfepy/discrete/equations
   src/sfepy/discrete/evaluate
   src/sfepy/discrete/evaluate_variable
   src/sfepy/discrete/functions
   src/sfepy/discrete/integrals
   src/sfepy/discrete/mass_operator
   src/sfepy/discrete/materials
   src/sfepy/discrete/parse_equations
   src/sfepy/discrete/parse_regions
   src/sfepy/discrete/probes
   src/sfepy/discrete/problem
   src/sfepy/discrete/projections
   src/sfepy/discrete/quadratures
   src/sfepy/discrete/simplex_cubature
   src/sfepy/discrete/state
   src/sfepy/discrete/variables

sfepy.discrete.common sub-package
"""""""""""""""""""""""""""""""""

Common lower-level code and parent classes for FEM and IGA.

.. toctree::
   :maxdepth: 2

   src/sfepy/discrete/common/dof_info
   src/sfepy/discrete/common/domain
   src/sfepy/discrete/common/fields
   src/sfepy/discrete/common/mappings
   src/sfepy/discrete/common/region

sfepy.discrete.fem sub-package
""""""""""""""""""""""""""""""

.. toctree::
   :maxdepth: 2

   src/sfepy/discrete/fem/domain
   src/sfepy/discrete/fem/extmods/_fmfield
   src/sfepy/discrete/fem/extmods/_geommech
   src/sfepy/discrete/fem/extmods/assemble
   src/sfepy/discrete/fem/extmods/bases
   src/sfepy/discrete/fem/extmods/cmesh
   src/sfepy/discrete/fem/extmods/lobatto_bases
   src/sfepy/discrete/fem/extmods/mappings
   src/sfepy/discrete/fem/facets
   src/sfepy/discrete/fem/fe_surface
   src/sfepy/discrete/fem/fea
   src/sfepy/discrete/fem/fields_base
   src/sfepy/discrete/fem/fields_hierarchic
   src/sfepy/discrete/fem/fields_nodal
   src/sfepy/discrete/fem/geometry_element
   src/sfepy/discrete/fem/global_interp
   src/sfepy/discrete/fem/history
   src/sfepy/discrete/fem/lcbc_operators
   src/sfepy/discrete/fem/linearizer
   src/sfepy/discrete/fem/mappings
   src/sfepy/discrete/fem/mesh
   src/sfepy/discrete/fem/meshio
   src/sfepy/discrete/fem/periodic
   src/sfepy/discrete/fem/poly_spaces
   src/sfepy/discrete/fem/refine
   src/sfepy/discrete/fem/utils

sfepy.discrete.iga sub-package
""""""""""""""""""""""""""""""

.. toctree::
   :maxdepth: 2

   src/sfepy/discrete/iga/domain
   src/sfepy/discrete/iga/domain_generators
   src/sfepy/discrete/iga/extmods/igac
   src/sfepy/discrete/iga/fields
   src/sfepy/discrete/iga/iga
   src/sfepy/discrete/iga/io
   src/sfepy/discrete/iga/mappings
   src/sfepy/discrete/iga/plot_nurbs
   src/sfepy/discrete/iga/utils

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

sfepy.linalg package
^^^^^^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 2

   src/sfepy/linalg/check_derivatives
   src/sfepy/linalg/eigen
   src/sfepy/linalg/geometry
   src/sfepy/linalg/sparse
   src/sfepy/linalg/utils
   src/sfepy/linalg/extmods/crcm

sfepy.mechanics package
^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/mechanics/contact_bodies
   src/sfepy/mechanics/elastic_constants
   src/sfepy/mechanics/matcoefs
   src/sfepy/mechanics/membranes
   src/sfepy/mechanics/tensors
   src/sfepy/mechanics/units

sfepy.mesh package
^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/mesh/geom_tools
   src/sfepy/mesh/mesh_generators
   src/sfepy/mesh/mesh_tools
   src/sfepy/mesh/splinebox

sfepy.optimize package
^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/optimize/free_form_def
   src/sfepy/optimize/shape_optim

sfepy.physics package
^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/physics/energy
   src/sfepy/physics/potentials
   src/sfepy/physics/schroedinger_app

sfepy.postprocess package
^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/postprocess/dataset_manager
   src/sfepy/postprocess/domain_specific
   src/sfepy/postprocess/plot_cmesh
   src/sfepy/postprocess/plot_dofs
   src/sfepy/postprocess/plot_facets
   src/sfepy/postprocess/plot_quadrature
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
   src/sfepy/solvers/ts_solvers

sfepy.terms package
^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   terms_overview
   src/sfepy/terms/terms
   src/sfepy/terms/terms_acoustic
   src/sfepy/terms/terms_adj_navier_stokes
   src/sfepy/terms/terms_basic
   src/sfepy/terms/terms_biot
   src/sfepy/terms/terms_constraints
   src/sfepy/terms/terms_diffusion
   src/sfepy/terms/terms_dot
   src/sfepy/terms/terms_elastic
   src/sfepy/terms/terms_electric
   src/sfepy/terms/terms_fibres
   src/sfepy/terms/terms_hyperelastic_base
   src/sfepy/terms/terms_hyperelastic_tl
   src/sfepy/terms/terms_hyperelastic_ul
   src/sfepy/terms/terms_membrane
   src/sfepy/terms/terms_navier_stokes
   src/sfepy/terms/terms_new
   src/sfepy/terms/terms_piezo
   src/sfepy/terms/terms_point
   src/sfepy/terms/terms_surface
   src/sfepy/terms/terms_th
   src/sfepy/terms/terms_volume
   src/sfepy/terms/utils

   src/sfepy/terms/extmods/terms
