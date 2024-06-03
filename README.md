# 1. Alanine dipeptide

  In this example, wo gonna to reweight the dihedral psi of alanine dipeptide to
generated two-dimensional fes based on phi and psi. firstly, we use Gromacs to run 
biased MD simulation of total timelength of 8 ns , note that Gromacs should be properly patched with PLUMED code.

`gmx mdrun -v -deffnm md -pin on -nt 6 -nsteps 4000000 --plumed plumed.dat`

  Once the simulation finished, we can obtain the FES based on dihedral phi, using the sum_hills kit in the PLUMED module.

`plumed sum_hills --hills HILLS`

  the default output file named "fes.dat" will be produce under current folder, then we can plot FES using Gnuplot software,
just type gnuplot in bash and we can inter into gnuplot terminal.

`gnuplot> plot "fes.dat" u 1:2 w l lw 2 not`

`gnuplot> set xlabel "Phr / rad"`

`gnuplot> set ylabel "Free energy / kj/mol"`

`gnuplot> replot`

![fes_sumhills.png](examples%2Falanine_dipeptide%2Ffes_sumhills.png)

  Because the presence of Metadynamics bias potential, the statistical weight of other degree of freedom (or other CVs) are 
already altered, so we can not calculate the histogram of other variables to derive the multi-dimensional FES directly, this
is the knowledge of why we need reweighting technique to calculate the unbiased histogram of other variables. Let's reweight 
  dihedral psi using MTD_reweight.py script, the prerequisite is that the python environment needs to be properly set, 
  and additional packages like numpy and scipy are also installed (using conda environmental is more recommended!), 
  running below command:

`python MTD_reweight.py --rcv 2 3 --bias 5 --smooth --sigma 0.1`

"--rcv" argument defines which column(s) of CVs to be used in reweighting, "--bias" means which column(s) of bias will be 
taken into account, e.g., MTD bias, wall bias, etc. --smooth" argument means that we will use gaussian function to smooth 
the image, and the sigma parameter determines the level of smoothing. all the grid dimension will be smoothed at the same 
level if single sigma parameter exists, while different smoothing level will be applied to each corresponding grid dimension 
if multiple sigma exist. The smoothed 2d-fes based on dihedral phi and psi is shown in below:

![fes_phi_psi.png](examples%2Falanine_dipeptide%2Ffes_phi_psi.png)


# 2. The dissociation of NaCl in water

  In this example, we will use reweighting method to discount the effect of wall potential. To limit the exploration of 
uninterested region of the distance Na-Cl beyond than 0.6 nm, an upper wall bias with harmonic format is thus introduced, as 
  show in below plumed input:

``UPPER_WALLS ...
 ARG=d1
 AT=0.6
 KAPPA=2000.0
 EXP=2
 EPS=1
 OFFSET=0.
 LABEL=uwall
... UPPER_WALLS``

it can be seen clearly that the obtained final FES profile goes up steeply with Na-Cl distance larger than 0.6 nm due to 
the restriction of wall bias, but once the wall bias is taken into account when performing reweighting using both the 
metadynamics bias and the upper wall bias, the effect of wall restriction is removed.

`python MTD_reweight.py --rcv 2 --bias 5 7 --smooth --sigma 0.1`

![reweighted_fes.png](examples%2FNaCl_dissociation_in_water%2Freweighted_fes.png)







