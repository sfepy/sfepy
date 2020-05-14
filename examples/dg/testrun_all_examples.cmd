
REM ADVECTION PARAMETRIZED
python .\run_dg_example.py .\advection\example_dg_quarteroni1.py
python .\run_dg_example.py .\advection\example_dg_quarteroni2.py
python .\run_dg_example.py .\advection\example_dg_quarteroni3.py
python .\run_dg_example.py .\advection\example_dg_quarteroni4.py

python .\run_dg_example.py .\advection\example_dg_advection2D.py
python .\run_dg_example.py .\advection\example_dg_advection1D.py


REM BURGESS PARAMETRIZED
python .\run_dg_example.py .\burgess\example_dg_burgess1D_hesthaven.py -dp
python .\run_dg_example.py .\burgess\example_dg_burgess2D_kucera.py


REM DIFFUSION PARAMETRIZED
python .\run_dg_example.py .\diffusion\example_dg_diffusion1D.py
python .\run_dg_example.py .\diffusion\example_dg_diffusion1D_hesthaven.py -dp
python .\run_dg_example.py .\diffusion\example_dg_diffusion2D_hartmann.py
python .\run_dg_example.py .\diffusion\example_dg_diffusion2D_qart.py
python .\run_dg_example.py .\diffusion\example_dg_laplace2D.py