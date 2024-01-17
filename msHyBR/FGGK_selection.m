function [U, M, T, V, z, lv] = FGGK_selection(A, R, Q, U, M, T, V, L, X)
%
% [U, M, T, V, z] = FGGK(A, R, Q, U, M, T, V, L)
%
%  Performs one step of the flexible, generalized Golub Kahan process, with
%  reorthogonalization.
%
% Input:
%          A - matrix A and A-transpose functions
%       R, Q - covariance matrices
%       U, V - accumulation of orthonormal vectors
%          M - k x k upper triangular matrix
%          T - (k+1) x k upper Hessenberg matrix
%          L - Changing preconditioner matrix
%          X - a fixed matrix that includes covariates
% Output:
%       U, V - updated matrix or orthonormal columns
%       M, T - updated matrices
%          z - new z vector
%        lv - new flexible hybr basis
%

% Get the previous dimension
% written by J. Jiang 10/2022

k = size(T,2)+1;

if k == 1
    v = A'*(R\U(:,k)); % initial v vector
    v = v(:);
else
    v = A'*(R\U(:,k)); % get net v vector
    v = v(:);
    for j = 1:k-1 % reorthogonalize against previous vectors
        M(j,k)=V(:,j)'*(Q*v);
        v = v - M(j,k)*V(:,j);
    end
end
M(k,k) = normM(v,Q);
v = v / M(k,k);

lv = L*(X'*v); % generate next flexible vector
z = [v; lv]; % augment

u = A*(Q*v) + A*(X*lv);
% u = A_times_vec(A, Q*v) + A_times_vec(A, lv);
u = u(:);

for j = 1:k % reorthogonalize against previous vectors
    T(j,k) = U(:,j)'*(R\u);
    u = u - T(j,k)*U(:,j);
end
T(k+1,k) = normM(u,@(x)R\x);
u = u / T(k+1,k);

U = [U, u]; % append 
V = [V, v]; % append

%% Subfunctions 
    function nrm = normM(v, M)
        if isa(M, 'function_handle')
            Mv = M(v);
        else
            Mv = M*v;
        end
        nrm = sqrt(v'*Mv);
    end
end
