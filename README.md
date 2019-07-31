# Fitz Hugh Nagumo

This small project to test various ODE integrators in C++. In particular, we use the Fitz Hugh Nagumo model as
toy model. 

The numerical solution of the Fitz Hugh Nagumo model was previously shown to be quite sensitive with the default Matlab
ODE integrations here `https://doi.org/10.1016/j.cnsns.2014.07.004`.

This project uses the parsing functionality from the repository `https://github.com/gregorweiss/argparse` and the numpy IO
 utilities from `https://github.com/llohse/libnpy`.