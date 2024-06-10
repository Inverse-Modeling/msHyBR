
% The function call eigs is not fully reproducible, so the data used to
% generate the experiments is give as a data file. However, similar
% experiments can be generated using the code above.
warning('the following code does not produce universally reproducible data')

rng(0,'twister')

% Build the stochastic component ($\zeta$ in Eq(7) in [1])
nx = length(lat); ny = length(lon);
xmin = [0 0];           % Coordinates of left corner
xmax = [1 1];           % Coordinates of right corner
nvec = [nx, ny];        % Number of points. Here this is a 256 x 256 grid

% Build covariance 
nu = 2.5; ell = .05; 
k = @(r) matern(r,nu,ell);

% Additional parameters governing length scales. Can be used to model
theta = [1.0 1.0];      %For now set them as isotropic

% The covariance matrix is constructed exploiting its toeplix structure: 
% only a row/column of the Toeplitz matrix needs to be explicitely given 
Qr = createrow(xmin,xmax,nvec,k,theta);
Qfun = @(x)toeplitzproduct(x, Qr, nvec);
[V,D] = eigs(Qfun,prod(nvec),20);
d = diag(D);
xt = 0.3*V*sqrt(D)*randn(size(d,1),1); 

% Build mean, ($X \beta$ in Eq(7) in [1])
% Using a zonation model for the mean and the modeling assumption that
% just a small amount of (hypothetical) states have gas concentration

Build X
X_reduced=[];   % matrix storing only the non-zero columns 
X=[];           % full zonation model including the sea
for j=1:10:170
    for i= 1:7:70
        x_vec = zeros(ny,nx);
        x_vec_aux = zeros(ny,nx);
        x_vec_aux(j:j+9,i:i+6)=1;
        x_vec(land) = x_vec_aux(land);
        if sum(x_vec)==0
            X = [X,x_vec(:)];
        else
            X_reduced = [X_reduced,x_vec(:)];
            X = [X,x_vec(:)];
        end
    end
end

% Note that the non zero coefficients of beta are forced by construction to
correspond to areas of land. 
beta_reduced = zeros(size(X_reduced,2),1);
ind1 = [14:16,22:24];
beta_reduced(ind1)=1;
ind2 = [31,41,51,60];
beta_reduced(ind2)=0.5;

beta_true=zeros(size(X,2),1);
beta_true(sum(X)~=0) = beta_reduced;

mean_true = X_reduced*beta_reduced;
s_true = stochastic_true + mean_true;