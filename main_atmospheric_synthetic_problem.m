
% This script concerns a synthetic inverse modeling problem where the true 
% solution features distinct blocks of emissions in different regions of 
% North America (i.e., using a zonation model)

%   - Load synthetic problem
%    (Optionally, generate new synthetic data instead running 'generate_synthetic_data.m')
%   - Run msHyBR
%   - Run 2 steps process: forward_selection + geostatistical_inversion
%     (forward selection is done using the 'BIC' method)

% Latest Update by M. Sabaté Landman 06/2024

clear all
close all

% Set paths
addpath './Data'
addpath './Auxiliary'
addpath './msHyBR'
addpath './genHyBR-master'
addpath './genHyBR-master/genHyBR'
addpath './genHyBR-master/toeplitz'

% Load system matrix 
% Load OCO2 dataset
load('Data/H_all_OCO2.mat', 'H');
load('Data/land_mask_OCO2.mat','land');
lon = -179.5:1:-10.5; ny = length(lon);
lat = 10.5:1:79.5; nx = length(lat); 

% Load synthetic solution
load('Data/Ex_2_true_solution');

% Alternatively, uncomment the following line to build a similar
% experiment. Warning: the generated data will not be fully reproducible.
% generate_data_Ex2

%% Construct forward model 
% Modify H to include sampling
n = length(land);
Samp = sparse(1:n,double(land),ones(size(land)),length(land),nx*ny);
H = H*Samp; 
b = H*s_true;

% Add noise to the data
noise  = randn(size(b(:)));
sigma = 0.04*(norm(b)/norm(noise));
N = sigma*noise;
bn = b(:) + N(:);

%% Solve the inverse problems

% Set parameters to solve the inverse problem
R = sigma^2*speye(size(H,1));
x_true = reshape(s_true,nx,ny);

% Select different Q that the one used in the forward model
% to avoid commiting the inverse crime
xmin = [0 0];           % Coordinates of left corner
xmax = [1 1];           % Coordinates of right corner
nvec = [nx, ny];        % Number of points. Here this is a 256 x 256 grid

nu = .5; ell = .5;          % Matern parameters
k = @(r) matern(r,nu,ell);  % Matern kernel
theta = [1.0 1.0];          % Additional parameters governing length scales. 

% Build row/column of the Toeplitz matrix
Qr = createrow(xmin,xmax,nvec,k,theta);
Qfun = @(x)toeplitzproduct(x, Qr, nvec);
Q = funMat(Qfun,Qfun, nvec);

%% Reconstruction using the new solver
thr = 10^(-3);
maxit = 50;
solver = 'tikhonov';
regpar = 'dp';
input = mshybr_set('InSolv', solver, 'RegPar',regpar, 'Iter', maxit,'x_true',x_true, 'nLevel',sigma,'thr',thr,'mask',land,'FlatTol', 10^-6);
disp('mshybr solver is running')
tic
[x_mshybr, beta_mshybr, output_mshybr] = msHyBR(H, bn, Q, R, X, input);
toc
e_sdhybr = output_mshybr.Enrm;
iteration_stop_sdhybr = output_mshybr.iterations;

%% Reconstruction using the two-step approach
disp('Two-step solver is running')
tic
disp('Starting subset selection')
% Add the true scaling parameter to the covariance Q 
% (this method does not allow for automatic tuning of this parameter)
nu = .5; ell = .5;
k = @(r) matern(r,nu,ell);
% Build row/column of the Toeplitz matrix
Qr = createrow(xmin,xmax,nvec,k,theta);
Qr = Qr*(0.3^2);
Qfun = @(x)toeplitzproduct(x, Qr, nvec); 
Q2 = funMat(Qfun,Qfun, nvec);
[X_selected,info_selected] = forward_selection_highdim(X, H, Q2, R, bn, 'BIC');

% Run solver given the selected covariates
% The used solver is genGKmean from [2]
disp('Starting reconstruction')
nbeta = size(info_selected.indices_selected,2);
gamma = 5; % Set parameter 
Qbeta = (1/gamma^2) * eye(nbeta); 
QgenGKmean = QtilMat(Qbeta,Q,X_selected);
HgenGKmean = matvecH_aug(1,size(info_selected.indices_selected,2),H);
input = HyBR_lsmrset('InSolv', 'tikhonov','RegPar','dp','Reorth','on','Iter', maxit,'x_true',x_true,'nLevel',sigma);
[xbeta_genGKmean, output_genGKmean] = genHyBR_Separ(HgenGKmean, bn(:), QgenGKmean, R, land, input);
toc
x_twostep=xbeta_genGKmean(1:11900);
beta_twostep=xbeta_genGKmean(11901:end);
save('atmospheric_synthetic_results')