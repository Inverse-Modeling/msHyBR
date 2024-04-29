
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
       

SOFTWARE:

    main_deblurring_problem.m       Generate synthetic problem 
                                    Run msHyBR
                                    Run 2 steps process: forward_selection + geostatistical_inversion
                                    (forward selection is done using the 'BIC', 'AIC' and 'F-test' methods)
                                    Run 2 steps process: exhaustive_selection + geostatistical_inversion
                                    (exhaustive selection is done using the 'BIC' and 'AIC' methods)
                                    
    plots_deblurring_problem.m      Plot the results from 'main_deblurring_problem.m'including:
                                    - Basis vectors in the problem formulation
                                    - Solutions Reconstrucions
                                    - Coefficient Reconstrucions
                                    - Relative error norm histories
                                    - Print out of 'Binary classifier metrics'
                                    Save figures in the folder 'Figures'
                                    * Requieres running 'main_deblurring_problem.m'

