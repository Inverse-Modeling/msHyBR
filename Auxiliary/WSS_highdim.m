function out = WSS_highdim(X, H, Q, R, z)

HX = H*X;

% Convert psi into a function
psi = @(x) H*(Q*(H'*x))+R*x;

% Build the augmented system
hSPS = @(x,trans) [psi(x(1:size(z,1))) + HX*x(size(z,1)+1:end); HX'*x(1:size(z,1))];
rhs = [z; zeros(size(X,2),1)];

% This is an ill-posed problem, so we use explicit regularization 
x = IRhybrid_gmres(hSPS,rhs,300);
out   = z'* x(1:size(z,1));
