
REM ADVECTION PARAMETRIZED
python .\simple.py .\examples\dg\advection\example_dg_quarteroni1.py  -o .\examples\dg\outputs\output\quarteroni1_sol
python .\simple.py .\examples\dg\advection\example_dg_quarteroni2.py  -o .\examples\dg\outputs\output\quarteroni2_sol
python .\simple.py .\examples\dg\advection\example_dg_quarteroni3.py  -o .\examples\dg\outputs\output\quarteroni3_sol
python .\simple.py .\examples\dg\advection\example_dg_quarteroni4.py  -o .\examples\dg\outputs\output\quarteroni4_sol

python .\simple.py .\examples\dg\advection\example_dg_advection2D.py  -o .\examples\dg\outputs\output\advection2D_sol
python .\simple.py .\examples\dg\advection\example_dg_advection1D.py  -o .\examples\dg\outputs\output\advection1D_sol

REM BURGESS PARAMETRIZED
python .\simple.py .\examples\dg\burgess\example_dg_burgess1D_hesthaven.py  -o .\examples\dg\outputs\output\burgess1D_hesthaven_sol
python .\simple.py .\examples\dg\burgess\example_dg_burgess2D_kucera.py  -o .\examples\dg\outputs\output\burgess2D_kucera_sol

REM DIFFUSION PARAMETRIZED
python .\simple.py .\examples\dg\diffusion\example_dg_diffusion1D.py  -o .\examples\dg\outputs\output\diffusion1D_sol
python .\simple.py .\examples\dg\diffusion\example_dg_diffusion1D_hesthaven.py  -o .\examples\dg\outputs\output\diffusion1D_hesthaven_sol
python .\simple.py .\examples\dg\diffusion\example_dg_diffusion2D_hartmann.py  -o .\examples\dg\outputs\output\diffusion2D_hartmann_sol
python .\simple.py .\examples\dg\diffusion\example_dg_diffusion2D_qart.py  -o .\examples\dg\outputs\output\diffusion2D_qart_sol
python .\simple.py .\examples\dg\diffusion\example_dg_laplace2D.py  -o .\examples\dg\outputs\output\laplace2D_sol