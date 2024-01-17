clear all

% Add paths
addpath('Auxiliary')
addpath('msHyBR')

% Set seeed for reproducibility
rng(21) 

%% Build problem 
% The predictor variables are 5 polynomial basis vectors on an equispaced
% grid and two mirrored Heavyside functions on the same grid
n_cols = 7;
n = 100;
domain = linspace(0,1,n);
X = zeros(n,n_cols);
for i=1:n_cols-2
    X(:,i)=chebyshevT(i-1,domain);
end
% 
X(1:50,n_cols-1) = ones(50,1);
X(51:end,n_cols) = ones(50,1);

% True solution has a sparse representation in that dictionary
% Coefficients follow a Bernoulli distribution
p = 0.2;
beta_true = [1, 0, 0,1,0,1,0]'; 

x_true = X*beta_true;

% Build solution covariance - Typical Matern covariance kernel
[x1,x2] = meshgrid(domain);
matern_R = abs(x1-x2);
Q = 0.01*matern(matern_R,0.5,1); 

% Add model noise with covariance Q
cholQ = chol(Q);
x_mod_noisy = x_true + cholQ*randn(n,1); 

%% Build measurements
sigma_blur = 1;
H = blur_1d(n,sigma_blur);
b = H*x_mod_noisy;

% Build noise covariance and add measument noise (Gaussian iid)
sigma_r = 0.01;
R = sigma_r*eye(n);
bn = b + sqrt(sigma_r)*randn(n,1);  

%% Compute solutions % Now with inverse crime

% Naive solution
x_naive = H\bn;

% New solver - msHyBR
maxit = 50;
solver = 'tikhonov';
regpar = 'optimal';
input = mshybr_set('InSolv', solver, 'RegPar',regpar, 'Iter', maxit,'x_true',x_true, 'nLevel',sqrt(sigma_r));
[x_ss, beta_ss, output_ss] = msHyBR(H, bn, Q, R, X, input);


%% COMPARISON: Forward selection + inverse solve
[Xhat_forward_VT,info_forward_VT] = forward_selection(X, H, Q, R, bn, 'VT');
[Xhat_forward_BIC,info_forward_BIC] = forward_selection(X, H, Q, R, bn, 'BIC');
[Xhat_forward_AIC,info_forward_AIC] = forward_selection(X, H, Q, R, bn, 'AIC');

[s_forward_VT,beta_forward_VT] = geostatistical_inversion(Xhat_forward_VT, H, Q, R, bn);
[s_forward_BIC,beta_forward_BIC] = geostatistical_inversion(Xhat_forward_BIC, H, Q, R, bn);
[s_forward_AIC,beta_forward_AIC] = geostatistical_inversion(Xhat_forward_AIC, H, Q, R, bn);

beta_forward_VT_for_comparison = zeros(size(beta_true));
beta_forward_VT_for_comparison(info_forward_VT.indices_selected)=beta_forward_VT;
beta_forward_BIC_for_comparison = zeros(size(beta_true));
beta_forward_BIC_for_comparison(info_forward_BIC.indices_selected)=beta_forward_BIC;
beta_forward_AIC_for_comparison = zeros(size(beta_true));
beta_forward_AIC_for_comparison(info_forward_AIC.indices_selected)=beta_forward_AIC;


%% COMPARISON: Exhaustive selection + inverse solve
% %Compute all possible combinations of columns and measure WSS. 
[Xhat_exhaustive_BIC,info_exhaustive_BIC] = exhaustive_selection(X, H, Q, R, bn, 'BIC');
[Xhat_exhaustive_AIC,info_exhaustive_AIC] = exhaustive_selection(X, H, Q, R, bn, 'AIC');

[s_exhaustive_BIC,beta_exhaustive_BIC] = geostatistical_inversion(Xhat_exhaustive_BIC, H, Q, R, bn);
[s_exhaustive_AIC,beta_exhaustive_AIC] = geostatistical_inversion(Xhat_exhaustive_AIC, H, Q, R, bn);

beta_exhaustive_BIC_for_comparison = zeros(size(beta_true));
beta_exhaustive_BIC_for_comparison(info_exhaustive_BIC.indices_selected)=beta_exhaustive_BIC;
beta_exhaustive_AIC_for_comparison = zeros(size(beta_true));
beta_exhaustive_AIC_for_comparison(info_exhaustive_AIC.indices_selected)=beta_exhaustive_AIC;

plots_deblurring_problem

function A = blur_1d(N,sigma)
z = [exp(-((0:N-1).^2)/(2*sigma^2))];
A = toeplitz(z);
A = sparse(A);
A = (1/((2*pi*sigma^2).^0.5))*A;
end