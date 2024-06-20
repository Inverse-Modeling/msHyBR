# Model Selection Hybrid Regularization

Here we provide MATLAB implementations of the model selection hybrid regularization method (msHybr), as well as two example scripts featuring a simple 1D deblurring problem and a synthetic atmospheric modeling problem.

## Project Description

Hybrid regularization methods have been proposed as effective approaches for solving large-scale ill-posed inverse problems. In particular, the model selection hybrid regularization method (msHybr), is designed to incorporate predictor variables to help estimate natural processes or parameters of interest from observed data.

This method follows a one-step approach, where model selection for a subset of relevant predictior variables and the reconstruction of the object of interest are done simultaneusly using flexible Krylov subspace methods for efficient optimization. 
 
## Installation 
### Software language

       MATLAB 9.14 (R2023a)
       For those without access to MATLAB, Octave provides an alternative platform.  
       Note that these codes have not been tested in Octave. 

### Requirements
The MainDrivers require the following package:

    "IR tools: A MATLAB Package of Iterative Regularization"
    by Silvia Gazzola, Per Christian Hansen and James G. Nagy
    https://github.com/jnagy1/IRtools.git

and require the data sets:
    
    "Geostatistical inverse modeling with large atmospheric data: 
    data files for a case study from OCO-2"
    by Scot M. Miller, Arvind K. Saibaba, Michael E. Trudeau, 
    Marikate E. Mountain, and Arlyn E. Andrews
    https://doi.org/10.5281/zenodo.3241466

     - Put the matrix 'H_all_OCO2.mat' in the folder "Data/" 

## How to Use
For Experiment 1 in [1], the main drivers include:
    
    main_deblurring_problem.m       Generate synthetic problem 
                                    Run msHyBR
                                    Run 2 steps process: forward_selection + geostatistical_inversion
                                    (forward selection is done using the 'BIC', 'AIC' and 'F-test' methods)
                                    Run 2 steps process: exhaustive_selection + geostatistical_inversion
                                    (exhaustive selection is done using the 'BIC' and 'AIC' methods)
                                    
    plots_deblurring_problem.m      Plot the results from 'main_deblurring_problem.m' including:
                                    - Basis vectors in the problem formulation
                                    - Solutions Reconstructions
                                    - Coefficient Reconstructions
                                    - Relative error norm histories
                                    - Print out of 'Binary classifier metrics'
                                    Save figures in the folder 'Figures'
                                    * Requires running 'main_deblurring_problem.m'

To run these codes, open a MATLAB command window, and type 
     
     >> main_deblurring_problem [press enter/return]
     >> plots_deblurring_problem [press enter/return]



For Experiment 2 in [1], the main drivers include:

    generate_synthetic_data.m		Generate synthetic data with the same distribution than the one corresponding to
						Experiment 2 in [1]. (Note that these do not generate fully reproducible codes)
   
    main_atmospheric_synthetic_problem.m	Load synthetic problem 
						(Optionally, generate new synthetic data instead running 'generate_synthetic_data.m')
                                    	Run msHyBR
                                    	Run 2 steps process: forward_selection + geostatistical_inversion
                                    	(forward selection is done using the 'BIC' method)
                                    
    plots_atmospheric_synthetic_problem.m       Plot the results from 'main_atmospheric_synthetic_problem.m' including:
                                    		- Representation of the donation model in the problem formulation
                                    		- Solutions Reconstructions
                                   			- Coefficient Reconstructions
                                    		- Relative error norm histories
                                    		- Representation of the selected basis vectors
                                    		Save figures in the folder 'Figures'
                                    		* Requires running 'main_atmospheric_synthetic_problem.m'

To run these codes, open a MATLAB command window, and type 
     
     >> main_atmospheric_synthetic_problem [press enter/return]
     >> plots_atmospheric_synthetic_problem [press enter/return]
     
### Contributors
        Malena Sabaté Landman, 
        Department of Mathematics, Emory University

        Julianne Chung, 
        Department of Mathematics, Emory University
        
        Jiahua Jiang,
        School of Mathematics, University of Birmingham
        
        Scot M. Miller, 
        Department of Environmental Health and Engineering, Johns Hopkins University
        
        Arvind K. Saibaba, 
        Department of Mathematics, North Carolina State University
	
## Licensing

If you use this codes, you *must* cite the original authors:

       [1] Sabaté Landman et al. "A Joint Reconstruction and Model Selection Approach for Large Scale Inverse Modeling". 2024.


[MIT](LICENSE)

## Acknowledgement

This work was partially supported by the National Science Foundation under grants DMS-2208294, DMS-2341843, DMS-2026830, and DMS-2026835. Any opinions, findings, conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
