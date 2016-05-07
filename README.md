# epilevy
Generalized random walks coupled to disease dynamics for epidemic and outbreak scenarios, plus surveillance tools based on compressed sensing.

Open source code but please write directly to Kyle Gustafson <kgustafson@idmod.org> before using this code for research or industry purposes.

Long-term, this should be a package capable of producing a variety of random walk processes suitable for modeling human mobility relevant to epidemiology modeling. Some inspiration is taken from "Levy walks" review article by Zaburdaev et al, Rev Mod Phys 2015. Code is based on my previous work http://dx.doi.org/10.1063/1.3690097. The first incarnation here will be only the Levy walks defined by a delta function coupling of step size to step duration, with no waiting time between steps.

Initial files:

Random walk generation
levy_batch : computing a parameter scan through a variety of Levy walks
levywalk_gen.m : generating the sequence of random numbers as steps in a Levy walk
dispLW_dist.m : converting steps into 2D displacements measured in Cartesian and polar coordinates

Feature extraction:
first_pass.m : compute first passage time for random walkers to reach a certain distance from the initial average position
disp_slope.m : compute the slope of the variance curve for the ensemble of walkers
cuml_pdf.m : compute the cumulative pdf of step sizes and find the slope of the assumed power law

Post-processing
levy_process.m : batch processing of Levy walks using the three feature extraction scripts
features_Levy.m : place features into a "big data" matrix for PCA/LDA/SVD analysis
