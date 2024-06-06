  ***************************************************************************
  * All the software  contained in this library  is protected by copyright. *
  * Permission  to use, copy, modify, and  distribute this software for any *
  * purpose without fee is hereby granted, provided that this entire notice *
  * is included  in all copies  of any software which is or includes a copy *
  * or modification  of this software  and in all copies  of the supporting *
  * documentation for such software.                                        *
  ***************************************************************************
  * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED *
  * WARRANTY. IN NO EVENT, NEITHER  THE AUTHORS, NOR THE PUBLISHER, NOR ANY *
  * MEMBER  OF THE EDITORIAL BOARD OF  THE JOURNAL  "NUMERICAL ALGORITHMS", *
  * NOR ITS EDITOR-IN-CHIEF, BE  LIABLE FOR ANY ERROR  IN THE SOFTWARE, ANY *
  * MISUSE  OF IT  OR ANY DAMAGE ARISING OUT OF ITS USE. THE ENTIRE RISK OF *
  * USING THE SOFTWARE LIES WITH THE PARTY DOING SO.                        *
  ***************************************************************************
  * ANY USE  OF THE SOFTWARE  CONSTITUTES  ACCEPTANCE  OF THE TERMS  OF THE *
  * ABOVE STATEMENT AND OF THE ACCOMPANYING FILE LICENSE.txt.               *
  ***************************************************************************
  
AUTHORS:
        
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
   

REFERENCE:

       "A Joint Reconstruction and Model Selection Approach for Large Scale Inverse Modeling". 2024.


SOFTWARE LANGUAGE:

       MATLAB 9.14 (R2023a)
       For those without access to MATLAB, Octave provides an alternative platform.  Note that these codes have not been tested in Octave. 

SOFTWARE:

For Experiment 1, the main drivers include:
    main_deblurring_problem.m       Generate synthetic problem 
                                    Run msHyBR
                                    Run 2 steps process: forward_selection + geostatistical_inversion
                                    (forward selection is done using the 'BIC', 'AIC' and 'F-test' methods)
                                    Run 2 steps process: exhaustive_selection + geostatistical_inversion
                                    (exhaustive selection is done using the 'BIC' and 'AIC' methods)
                                    
    plots_deblurring_problem.m      Plot the results from 'main_deblurring_problem.m' including:
                                    - Basis vectors in the problem formulation
                                    - Solutions Reconstrucions
                                    - Coefficient Reconstrucions
                                    - Relative error norm histories
                                    - Print out of 'Binary classifier metrics'
                                    Save figures in the folder 'Figures'
                                    * Requires running 'main_deblurring_problem.m'
To run these codes, open a MATLAB command window, and type 
        >> main_deblurring_problem [press enter/return]
        >> plots_deblurring_problem [press enter/return]

