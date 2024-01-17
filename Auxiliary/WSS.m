function [out,beta] = WSS(X, H, Q, R, z)

phi = H*Q*H'+R;

HX = H*X;

% inner_braket = HX'*(phi\HX);
% phi_inv_z = phi\z;
% out = z'* phi_inv_z - z'*(phi\(HX*(inner_braket\(HX'*phi_inv_z))));

psiinv = inv(phi);
A = z' * psiinv * z;
B = z' * (psiinv * HX);
C = HX' * psiinv * HX;
beta = (C \ B');
out   = A - B * beta;

