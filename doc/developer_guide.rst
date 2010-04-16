.. highlight:: python
   :linenothreshold: 3

Developer Guide
===============

This section purports to document the *SfePy* internals. It is mainly useful for those who wish to develop *SfePy* and understand the inner workings of the code.

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
  #. to pass those data to a element matrix/rezidual evaluation function
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
  `Volume` of `Surface`)
- `__call__()` method - subclasses of `Term` either implement `__call__()` or
  plug in a proper `_call()` method that is called by the default
  implementation. It takes the following arguments::

      __call__(self, diff_var=None, chunk_size=None, **kwargs)

 - `diff_var` is either None (rezidual mode), or the name of the
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
          n_el, n_qp, dim, n_ep = ap.get_v_data_shape(self.integral_name)

          if diff_var is None:
              shape = (chunk_size, 1, dim * n_ep, 1)

          elif diff_var == self.get_arg_name( 'state' ):
              shape = (chunk_size, 1, dim * n_ep, dim * n_ep)

          else:
              raise StopIteration

          bf = ap.get_base('v', 0, self.integral_name)
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

- lines 19-26: determine data shape of the chunk, depending whether in rezidual
  or matrix mode.
- line 28: get volume base function corresponding to the integral used
- lines 29-32: evaluate for each chunk and yield the results

In practice, such a term would inherit also from
:class:`sfepy.terms.terms_base.VectorVector` - then it could look, for example,
like the :class:`sfepy.terms.termsMass.MassTerm`, i.e. take care just about
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
          bf = ap.get_base('v', 0, self.integral_name)
  
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
   src/sfepy/base/la
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

   src/sfepy/fem/functions
   src/sfepy/fem/geometry_element
   src/sfepy/fem/mesh
   src/sfepy/fem/mesh_generators
   src/sfepy/fem/meshio
   src/sfepy/fem/periodic
   src/sfepy/fem/poly_spaces
   src/sfepy/fem/probes

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
   src/sfepy/homogenization/coefs_piezo
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
