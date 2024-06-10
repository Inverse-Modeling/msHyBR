function ERR = TikOPT_land(lambda, V, bhat, s, QV, xtrue, land)
%
%    ERR = TikOPT(lambda, V, bhat, s, QV, xtrue)
%
%  This function evaluates the error function for Tikhonov
%  regularization.  
%
%  Input:  lambda -  regularization parameter
%               V - singular vectors of subproblem
%            bhat -  vector U'*b, where U = left singular vectors
%               s -  vector containing the singular values
%           QV - Q*basis vectors from GK bidiagonalization
%           xtrue - true image
%
%  Output:  ERR = ||xtrue - V*y_k|| where y_k is the Tikhonov solution to
%  the subproblem
%
%  J.Chung and K. Palmer 3/2014

D = abs(s).^2 + lambda^2;
bhat = conj(s) .* bhat(1:length(s));
xhat = bhat ./ D;
y = V * xhat;
p = 2;

QVy = QV*y;
ERR = norm(xtrue(land) - QVy(land), p);

