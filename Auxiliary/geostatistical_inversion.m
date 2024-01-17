function [s,beta] = geostatistical_inversion(X, H, Q, R, z)

phi = H*Q*H'+R;

HX = H*X;
HQ = H*Q;

aux = [phi, HX; HX', zeros(size(HX,2))]\[HQ;X'];

Lambda = aux(1:size(HX,1),:);

s = Lambda*z;
beta = X'*(Q\s);
aux2 = X'*(Q\X);
beta = aux2\beta;