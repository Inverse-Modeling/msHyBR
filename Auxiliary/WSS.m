function [out,beta] = WSS(X, H, Q, R, z)
%
%   Computes the weighted sum of squares (WSS). 
%   Note: this code only works for small problems.
%
%   Input:
%          H - syste matrix H
%       R, Q - covariance matrices
%          X - a fixed matrix that includes covariates
%          z - right-hand side
% Output:
%        out - WWS
%       beta - test parameter

phi = H*Q*H'+R;

HX = H*X;

psiinv = inv(phi);
A = z' * psiinv * z;
B = z' * (psiinv * HX);
C = HX' * psiinv * HX;
beta = (C \ B');
out   = A - B * beta;

