function [s,beta] = geostatistical_inversion(X, H, Q, R, z)
%
%   Computes the closed for solution of the geostatistical 
%   inverse problems. 
%
%   Note: this code only works for small problems.
%
%   Input:
%          H - syste matrix H
%       R, Q - (noise and solution) covariance matrices
%          X - a fixed matrix that includes covariates
%          z - right-hand side
% Output:
%          s - estimated solution
%       beta - estimated parameters

phi = H*Q*H'+R;

HX = H*X;
HQ = H*Q;

aux = [phi, HX; HX', zeros(size(HX,2))]\[HQ;X'];

Lambda = aux(1:size(HX,1),:);

s = Lambda*z;
beta = X'*(Q\s);
aux2 = X'*(Q\X);
beta = aux2\beta;