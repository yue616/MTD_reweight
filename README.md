# 1. Alanine dipeptide

  In this example, we will reweight the dihedral angle psi of alanine dipeptide to generate a two-dimensional free energy 
  surface (FES) based on the dihedral angles phi and psi. First, we use GROMACS to run a biased molecular dynamics (MD) 
  simulation for a total duration of 8 ns. Note that GROMACS must be properly patched with the PLUMED code to perform this task.

`gmx mdrun -v -deffnm md -pin on -nt 6 -nsteps 4000000 --plumed plumed.dat`

  Once the simulation is finished, we can obtain the FES based on dihedral phi, using the sum_hills kit in the PLUMED module.

`plumed sum_hills --hills HILLS`

  The default output file named "fes.dat" will be produced in the current folder. To plot the FES, you can use Gnuplot software. 
  Simply type gnuplot in the bash terminal to enter the Gnuplot command interface.

`gnuplot> plot "fes.dat" u 1:2 w l lw 2 not`

`gnuplot> set xlabel "Phr / rad"`

`gnuplot> set ylabel "Free energy / kj/mol"`

`gnuplot> replot`

![fes_sumhills.png](examples%2Falanine_dipeptide%2Ffes_sumhills.png)

  Due to the presence of the Metadynamics bias potential, the statistical weights of other degrees of freedom (or other 
  collective variables, CVs) are altered. Therefore, we cannot directly calculate the histogram of other variables to 
  derive the multi-dimensional FES. This is why the reweighting technique is necessary to calculate the unbiased histogram 
  of other variables.

  To reweight the dihedral angle psi, we will use the MTD_reweight.py script. Ensure that your Python environment is properly 
  set up, with necessary packages like NumPy and SciPy installed. Using a Conda environment is highly recommended for 
  this purpose. Run the following command:

`python MTD_reweight.py --rcv 2 3 --bias 5 --smooth --sigma 0.1`

The --rcv argument defines which column(s) of collective variables (CVs) will be used in reweighting. The --bias argument 
specifies which column(s) of bias will be considered, such as Metadynamics bias or wall bias. The --smooth argument indicates 
that a Gaussian function will be used to smooth the image, with the sigma parameter determining the level of smoothing. 
If a single sigma parameter is provided, all grid dimensions will be smoothed at the same level. If multiple sigma values 
are provided, different smoothing levels will be applied to each corresponding grid dimension, note that the number of sigma
and the reweighted CVs should be equivalent.

The smoothed 2D free energy surface (FES) based on the dihedral angles psi and phi is shown below:

![fes_phi_psi.png](examples%2Falanine_dipeptide%2Ffes_phi_psi.png)


# 2. The dissociation of NaCl in water

  In this example, we will use the reweighting method to discount the effect of the wall potential. To limit the exploration 
of the uninteresting region where the distance between Na and Cl exceeds 0.6 nm, an upper wall bias with a harmonic format 
is introduced. The PLUMED input for this setup is shown below:

``UPPER_WALLS ...
 ARG=d1
 AT=0.6
 KAPPA=2000.0
 EXP=2
 EPS=1
 OFFSET=0.
 LABEL=uwall
... UPPER_WALLS``

It is evident that the final obtained FES profile rises steeply when the Na-Cl distance exceeds 0.6 nm due to the restriction 
imposed by the wall bias. However, when performing reweighting using both the Metadynamics bias and the upper wall bias, 
the effect of the wall restriction is removed.  

`python MTD_reweight.py --rcv 2 --bias 5 7 --smooth --sigma 0.1`

![reweighted_fes.png](examples%2FNaCl_dissociation_in_water%2Freweighted_fes.png)







