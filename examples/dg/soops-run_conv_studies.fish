#!/usr/bin/fish
set -x PYTHONPATH "/home/tomas/Python_projects/sfepy/"

set elout "/media/tomas/Elements/Data/soops_dg"

soops-run -o $elout/example_dg_advection2D "problem_file='advection/example_dg_advection2D',
                                            --cfl=[0.1, .5], --adflux=[0.0, 0.5, 1.0] , --limit=[@defined, @undefined],
                                            mesh=['mesh/mesh_tensr_2D_01_4.vtk'], output_dir='$elout/example_dg_advection2D/%s'" run_dg_conv_study.py

soops-run -o $elout/example_dg_advection1D "problem_file='advection/example_dg_advection1D',
                                            --cfl=[0.1, .5], --adflux=[0.0, 0.5, 1.0] , --limit=[@defined, @undefined],
                                            mesh=['mesh/mesh_tensr_1D_01_2.vtk'], output_dir='$elout/example_dg_advection1D/%s'" run_dg_conv_study.py

soops-run -o $elout/example_dg_quarteroni1 "problem_file='advection/example_dg_quarteroni1',
                                            --flux=[0.0, 0.5, 1.0], --cw=[1, 10, 100, 1e3, 1e4, 1e5], --diffcoef=[1e-3, 1e-2, .1, 1],
                                            mesh=['mesh/mesh_tensr_2D_01_4.vtk', 'mesh/mesh_simp_2D_01_4.vtk'], output_dir='$elout/example_dg_quarteroni1/%s'" run_dg_conv_study.py

soops-run -o $elout/example_dg_quarteroni2 "problem_file='advection/example_dg_quarteroni2',
                                            --flux=[0.0, 0.5, 1.0], --cw=[1, 10, 100, 1e3, 1e4, 1e5], --diffcoef=[1e-3, 1e-2, .1, 1],
                                            mesh=['mesh/mesh_tensr_2D_01_4.vtk', 'mesh/mesh_simp_2D_01_4.vtk'], output_dir='$elout/example_dg_quarteroni2/%s'" run_dg_conv_study.py
soops-run -o $elout/example_dg_quarteroni3 "problem_file='advection/example_dg_quarteroni3',
                                            --flux=[0.0, 0.5, 1.0], --cw=[1, 10, 100, 1e3, 1e4, 1e5], --diffcoef=[1e-3, 1e-2, .1, 1],
                                            mesh=['mesh/mesh_tensr_2D_01_4.vtk', 'mesh/mesh_simp_2D_01_4.vtk'], output_dir='$elout/example_dg_quarteroni3/%s'" run_dg_conv_study.py

soops-run -o $elout/example_dg_quarteroni4 "problem_file='advection/example_dg_quarteroni1',
                                            --flux=[0.0, 0.5, 1.0], --cw=[1, 10, 100, 1e3, 1e4, 1e5], --diffcoef=[1e-3, 1e-2, .1, 1],
                                            mesh=['mesh/mesh_tensr_2D_01_4.vtk', 'mesh/mesh_simp_2D_01_4.vtk'], output_dir='$elout/example_dg_quarteroni4/%s'" run_dg_conv_study.py

soops-run -o $elout/example_dg_burgess1D_hesthaven "problem_file='burgess/example_dg_burgess1D_hesthaven',
                                                    --adflux=[0.0, 0.5, 1.0], --limit=[@defined, @undefined], --cfl=[1e-3, 1e-2],
                                                    --cw=[1, 10, 100, 1e3, 1e4, 1e5], --diffcoef=[1e-3, 1e-2],
                                                    mesh=['mesh/mesh_tensr_1D_11_2.vtk'], output_dir='$elout/example_dg_burgess1D_hesthaven/%s'" run_dg_conv_study.py

soops-run -o $elout/example_dg_burgess2D_kucera "problem_file='burgess/example_dg_burgess2D_kucera',
                                                 --adflux=[0.0, 0.5, 1.0], --limit=[@defined, @undefined], --dt=[1e-5],
                                                 --cw=[1, 10, 100, 1e3, 1e4, 1e5], --diffcoef=[1e-2],
                                                 mesh=['mesh/mesh_tensr_2D_11_4.vtk', 'mesh/mesh_simp_2D_11_4.vtk'], output_dir='$elout/example_dg_burgess2D_kucera/%s'" run_dg_conv_study.py

soops-run -o $elout/example_dg_diffusion1D_hesthaven "problem_file='diffusion/example_dg_diffusion1D_hesthaven',
                                                      --cfl=[1e-3, 1e-2], --cw=[1, 10, 100, 1e3, 1e4, 1e5], --diffcoef=[1e-2],
                                                      mesh=['mesh/mesh_tensr_1D_0-2pi_100.vtk'], --refines='[0,1,2]', output_dir='$elout/example_dg_diffusion1D_hesthaven/%s'" run_dg_conv_study.py

soops-run -o $elout/example_dg_diffusion2D_hartmann "problem_file='diffusion/example_dg_diffusion2D_hartmann',
                                                     --cw=[1, 10, 100, 1e3, 1e4, 1e5], --diffcoef=[1e-2],
                                                     mesh=['mesh/mesh_tensr_2D_01_4.vtk', 'mesh/mesh_simp_2D_01_4.vtk'], output_dir='$elout/example_dg_diffusion2D_hartmann/%s'" run_dg_conv_study.py

soops-run -o $elout/example_dg_diffusion2D_qart "problem_file='diffusion/example_dg_diffusion2D_qart',
                                                 --cw=[1, 10, 100, 1e3, 1e4, 1e5], --diffcoef=[1e-2],
                                                 mesh=['mesh/mesh_tensr_2D_01_4.vtk', 'mesh/mesh_simp_2D_01_4.vtk'], output_dir='$elout/example_dg_diffusion2D_qart/%s'" run_dg_conv_study.py

soops-run -o $elout/example_dg_laplace2D "problem_file='diffusion/example_dg_laplace2D',
                                          --cw=[1, 10, 100, 1e3, 1e4, 1e5],
                                          mesh=['mesh/mesh_tensr_2D_01_4.vtk', 'mesh/mesh_simp_2D_01_4.vtk'], output_dir='$elout/example_dg_laplace2D/%s'" run_dg_conv_study.py

