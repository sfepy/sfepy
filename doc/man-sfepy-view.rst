sfepy-view
==========

synopsis
--------

sfepy-view [-h] [-f field_spec [field_spec ...]]
           [--fields-map map [map ...]] [-s step] [-l] [-i ISOSURFACES]
           [-e] [-w field] [--factor factor] [--opacity opacity]
           [--color-map cmap] [--axes-options options [options ...]]
           [--no-axes] [--grid-vector1 grid_vector1]
           [--grid-vector2 grid_vector2] [--max-plots MAX_PLOTS]
           [--no-labels] [--label-position position] [--no-scalar-bars]
           [--scalar-bar-size size] [--scalar-bar-position position]
           [-v position] [--camera-position camera_position]
           [--window-size window_size] [-a output_file] [-r rate]
           [-o output_file] [--off-screen] [-2]
           filenames [filenames ...]

description
-----------

This is a script for quick VTK-based visualizations of finite element
computations results.

In the examples below it is supposed that sfepy is installed. When using the
in-place build, replace ``sfepy-view`` by ``python3 sfepy/scripts/resview.py``.

Examples
''''''''

The examples assume that
``python -c "import sfepy; sfepy.test('--output-dir=output-tests')"``
has been run successfully and the resulting data files are present.

- View data in output-tests/test_navier_stokes.vtk::

    sfepy-view output-tests/navier_stokes-navier_stokes.vtk

- Customize the above output:
  plot0: field "p", switch on edges,
  plot1: field "u", surface with opacity 0.4, glyphs scaled by factor 2e-2::

    sfepy-view output-tests/navier_stokes-navier_stokes.vtk -f p:e:p0 u:o.4:p1 u:g:f2e-2:p1

- As above, but glyphs are scaled by the factor determined automatically as
  20% of the minimum bounding box size::

    sfepy-view output-tests/navier_stokes-navier_stokes.vtk -f p:e:p0 u:o.4:p1 u:g:f10%:p1

- View data and take a screenshot::

    sfepy-view output-tests/diffusion-poisson.vtk -o image.png

- Take a screenshot without a window popping up::

    sfepy-view output-tests/diffusion-poisson.vtk -o image.png --off-screen

- Create animation from output-tests/diffusion-time_poisson.*.vtk::

    sfepy-view output-tests/diffusion-time_poisson.*.vtk -a mov.mp4

- Create animation from output-tests/test_hyperelastic.*.vtk,
  set frame rate to 3, plot displacements and mooney_rivlin_stress::

    sfepy-view output-tests/test_hyperelastic_TL.*.vtk -f u:wu:e:p0 mooney_rivlin_stress:p1 -a mov.mp4 -r 3

positional arguments
--------------------

| filenames             results/mesh files

optional arguments
------------------

  -h, --help            show this help message and exit
  -f field_spec [field_spec ...], --fields field_spec [field_spec ...]
                        fields to plot, options separated by ":" are possible:
                        "cX" - plot only Xth field component; "e" - print
                        edges; "fX" - scale factor for warp/glyphs, see
                        --factor option; "g - glyphs (for vector fields only),
                        scale by factor; "iX" - plot X isosurfaces; "tX" -
                        plot X streamlines, gradient employed for scalar
                        fields; "mX" - plot cells with mat_id=X; "oX" - set
                        opacity to X; "pX" - plot in slot X; "r" - recalculate
                        cell data to point data; "sX" - plot data in step X;
                        "vX" - plotting style: s=surface, w=wireframe,
                        p=points; "wX" - warp mesh by vector field X, scale by
                        factor
  --fields-map map [map ...]
                      map fields and cell groups, e.g. 1:u1,p1 2:u2,p2
-s step, --step step  select data in a given time step
-l, --outline         plot mesh outline
-i ISOSURFACES, --isosurfaces ISOSURFACES
                      plot isosurfaces [default: 0]
-e, --edges           plot cell edges
-w field, --warp field
                      warp mesh by vector field
--factor factor       scaling factor for mesh warp and glyphs. Append "%" to
                      scale relatively to the minimum bounding box size.
--opacity opacity     set opacity [default: 1.0]
--color-map cmap      set color_map, e.g. hot, cool, bone, etc. [default:
                      viridis]
--axes-options options [options ...]
                      options for directional axes, e.g. xlabel="z1"
                      ylabel="z2", zlabel="z3"
--no-axes             hide orientation axes
--grid-vector1 grid_vector1
                      define positions of plots along grid axis 1 [default:
                      "0, 0, 1.6"]
--grid-vector2 grid_vector2
                      define positions of plots along grid axis 2 [default:
                      "0, 1.6, 0"]
--max-plots MAX_PLOTS
                      maximum number of plots along grid axis 1 [default: 4]
--no-labels           hide plot labels
--label-position position
                      define position of plot labels [default: "-1, -1, 0,
                      0.2"]
--no-scalar-bars      hide scalar bars
--scalar-bar-size size
                      define size of scalar bars [default: "0.15, 0.05"]
--scalar-bar-position position
                      define position of scalar bars [default: "0.8, 0.02,
                      0, 1.5"]
-v position, --view position
                      camera azimuth, elevation angles, and optionally zoom
                      factor [default: "225,75,0.9"]
--camera-position camera_position
                      define camera position
--window-size window_size
                      define size of plotting window
-a output_file, --animation output_file
                      create animation, mp4 file type supported
-r rate, --frame-rate rate
                      set framerate for animation
-o output_file, --screenshot output_file
                      save screenshot to file
--off-screen          off screen plots, e.g. when screenshotting
-2, --2d-view         2d view of XY plane
