function [x_out, beta_out, output] = msHyBR(A, b, Q, R, X, options)
%
% msHyBR is a subset-selection method used for solving large-scale, 
% ill-posed inverse problems of the form:
%               b = A*s + noise
% where s = X \beta + \psi.
%
% The mshybr method combines a flexible, generalized Golub-Kahan (FGGK)
% projection method with regularization on the projected problem to tackle
% problems where the regularization term is of the form,
%       ....
% where \lambda and \alpha are regularization parameters.
%
% Inputs:
%     A : either (a) a full or sparse matrix
%                (b) a matrix object that performs mat*vec and mat'*vec
%         Note: If A is a function handle, create an object called funMat
%         e.g., A = funMat(Afun, Atfun) where Afun and Atfun are function
%                                       handles for A*vec and A'*vec
%     b : rhs vector
%  Q, R : covariance matrices, Q is either (a) a full or sparse matrix or
%                (b) a matrix object that performs matrix*vector operations
%  X: a fixed matrix that includes covariates (e.g., environmental data or activity data) or a bottom-up inventory/flux model
%  thr : threshhold for the regularization matrix
%         Note: If Q is a function handle, create an object called funMat
%         e.g., Q = funMat(Qfun, Qtfun) where Qfun and Qtfun are function
%                                       handles for Q*vec and Q'*vec
%  mask :  mask used to compute relative reconstruction errors and optimal
%                   regularization parameter
% options : structure with the following fields (optional)
%         InSolv - solver for the inner problem: [none | TSVD | {Tikhonov}]
%         RegPar - a value or method to find the regularization parameter:
%                       [non-negative scalar | DP | GCV | {WGCV} | optimal]
%                   Note: 'optimal' requires x_true
%         nLevel - if RegPar is 'DP', then nLevel represents the noise level and must be
%                       [non-negative scalar | {est}]
%          Omega - if RegPar is 'WGCV', then omega must be
%                       [non-negative scalar | {adapt}]
%           Iter - maximum number of GK iterations:
%                       [ positive integer | {min(m,n,100)} ]
%         Reorth - reorthogonalize Lanczos subspaces: [on | {off}]
%         x_true - True solution : [ array | {off} ]
%                Returns error norms with respect to x_true at each iteration
%                and is used to compute 'optimal' regularization parameters
%         BegReg - Begin regularization after this iteration:
%                   [ positive integer | {2} ]
%             Vx - extra space needed for finding optimal reg. parameters
%        FlatTol - Tolerance for detecting flatness in the GCV curve as a
%                    stopping criteria
%                   [ non-negative scalar | {10^-6}]
%         MinTol - Window of iterations for detecting a minimum of the GCV curve
%                    as a stopping criteria
%                   [ positive integer | {3}]
%         ResTol - Residual tolerance for stopping the LBD iterations,
%                    similar to the stopping criteria from [1]: [atol, btol]
%                   [non-negative scalar  | {[10^-6, 10^-6]}]
%            thr - Threshhold used for flexible methods
%                   [ non-negative scalar | {10^-3}]
%           mask - mask : [ array | {off} ]
%                mask used to compute 'optimal' regularization parameters
%                and relative reconstruction errors (e.g., a land mask)
%
%       Note: options is a structure created using the function 'HyBRset'
%               (see 'HyBRset' for more details)
%
% Outputs:
%     x_out : computed solution
%     beta_out - coefficient beta
%     output : structure with the following fields:
%         iterations - stopping iteration (options.Iter | GCV-determined)
%         GCVstop - GCV curve used to find stopping iteration
%         Enrm - relative error norms (requires x_true)
%         Rnrm - relative residual norms
%         Xnrm - relative solution norms
%         U,QV - U and V are genGK basis vectors and Q is from the prior
%         B - bidiagonal matrix from genGK
%         flag - a flag that describes the output/stopping condition:
%                       1 - flat GCV curve
%                       2 - min of GCV curve (within window of MinTol its)
%                       3 - performed max number of iterations
%                       4 - achieved residual tolerance
%         alpha - regularization parameter at (output.iterations) its
%         Alpha - vector of all regularization parameters computed
%       
% References:
%   [1] to appear.
%
% J. Jiang, modified 10/2022
% M. Sabaté Landman modified 12/2023
%
% Based on codes from: 
% J. Chung and J. Nagy 3/2007
% J. Chung and A. Saibaba, modified 1/2016
% J. Jiang, modified 5/2020

%% Initialization for sdhybr
defaultopt = struct('InSolv','tikhonov','RegPar','wgcv','nLevel', 'est',...
    'Omega', 'adapt', 'Iter', [] , 'Reorth', 'off', 'x_true', 'off', 'BegReg', 2,...
    'Vx' , [], 'FlatTol', 10^-6, 'MinTol', 10, 'ResTol', [10^-6, 10^-6],'thr',1e-3,'mask','off');

% If input is 'defaults,' return the default options in x_out
if nargin==1 && nargout <= 1 && isequal(A,'defaults')
    x_out = defaultopt;
    return;
end

% Check for acceptable number of input arguments
if nargin < 4
    error('sdhybr: Not Enough Inputs')
elseif nargin < 5
    options = [];
end
if isempty(options)
    options = defaultopt;
end

% Get options:
[m,n] = size(A);
defaultopt.Iter = min([m, n, 100]);
% options = sdhybr_get(defaultopt, options);
options = mshybr_set(defaultopt, options);

solver = mshybr_get(options,'InSolv',[],'fast');
regPar = mshybr_get(options,'RegPar',[],'fast');
nLevel = mshybr_get(options,'nLevel',[],'fast');
omega = mshybr_get(options,'Omega',[],'fast');
maxiter = mshybr_get(options,'Iter',[],'fast');
x_true = mshybr_get(options,'x_true',[],'fast');
regstart = mshybr_get(options,'BegReg',[],'fast');
degflat = mshybr_get(options,'FlatTol',[],'fast');
mintol = mshybr_get(options,'MinTol',[],'fast');
restol = mshybr_get(options,'ResTol',[],'fast');
thr = mshybr_get(options,'thr',[],'fast');
mask = mshybr_get(options,'mask',[],'fast');

adaptWGCV = strcmp(regPar, {'wgcv'}) && strcmp(omega, {'adapt'});
notrue = strcmp(x_true,{'off'});
nomask = strcmp(mask,{'off'});

if (strcmpi(regPar,'dp') ||strcmpi(regPar,'upre')) && strcmpi(nLevel,'est')
    % Estimate the noise level using finest wavelet coefficients
    if size(b,2) ==1 %1D
        [cA, cD] = dwt(b,'db1');
        nLevel = median(abs(cD(:)))/.67
        options = sdhybr_get(options, 'nLevel', nLevel);
    else %2D
        [cA2,cH2,cV2,cD2] = dwt2(b,'db1');
        nLevel = median(abs(cD2(:)))/.67
        options = sdhybr_get(options, 'nLevel', nLevel);
    end
end
%--------------------------------------------
%  The following is needed for RestoreTools:
if isa(A, 'psfMatrix')
    bSize = size(b);
    b = b(:);
    A.imsize = bSize;
    if ~notrue
        x_true = x_true(:);
    end
end
%  End of new stuff needed for RestoreTools
%--------------------------------------------

if ~notrue
    if nomask
        x_truemask = x_true;
    else
        x_truemask = x_true(mask);
    end
    nrmtruemask = norm(x_truemask(:));
end

% Set-up output parameters:
outputparams = nargout>1;
if outputparams
    output.iterations = maxiter;
    output.GCVstop = [];
    output.Enrm = ones(maxiter,1);
    output.Rnrm = ones(maxiter,1);
    output.Xnrm = ones(maxiter,1);
    output.U = [];
    output.QV = [];
    output.B = [];
    output.flag = 3;
    output.alpha = 0;
    output.Alpha = [];
    output.lambda = 0;
    output.Lambda = [];
end


%% Initialization for FGGK
beta = normM(b,@(x)R\x);
U = (1 / beta)*b;
N = size(A,2);
% L = speye(N);
L = 1;
T = [];
M = [];
Z = [];
V_flex = [];
LV = [];

%% Main Code Begins Here
B = []; V = []; QV = []; GCV = []; Omega= []; x_out = []; Alpha = []; Lambda = [];
RegPar = 'none'; terminate = 1; warning = 1; norma = 0; normr = beta; iterations_save = 1; diff_x = [];diff_res = [];
h = waitbar(0, 'Beginning iterations: please wait ...'); G_ob = [];
for i = 1:maxiter+1 %Iteration (i=1) is just an initialization
    vector = zeros(i+1,1); vector(1) = 1; vector = beta*vector; % beta*e1
    
    % the Krylov subspace is expanded at each step with FGGK
    [U, T, M, V_flex, z, lv] = FGGK_selection(A, R, Q, U, T, M, V_flex, L, X); 
    Z = [Z, z]; % Extend the solution subspace
    
    %% generate matrix requested by projected problem
    
 
    if  i==1
        [Qs,Rs] = qr(lv,0);
    else
        [Qs, Rs] = qr_lastcol(Qs,Rs,lv);
    end
    
    if i >= 2 %Begin GK iterations
        
        
        if i >= regstart %Begin to regularize projected problem
            RegPar = regPar;
        end
        
        switch RegPar
            case 'optimal'
                p0 = [-0.5 -0.5];
                if isfield(options,'mu')
                    mu = options.mu;
                    errhan = @(p)sdecomOPT(p,G,Rs,QV,Z,beta,x_true-mu);
                else
                    errhan = @(p)sdOPT(p,M,Rs,L,Q,Z,beta,x_true,X,mask);
                    %errhan = @(p)sdOPT_s2(p,M,Rs,L,Q,Z,beta,[1,0,0,1,0,1,0],X,mask);
                end
            case 'upre'
                p0 = [-0.5 -0.5];
                sigma = nLevel;
                errhan = @(p)sdUPRE(p,M,Rs,beta,sigma);
            case 'wgcv'
                p0 = [-0.5 -0.5];
                omega = size(M,2)/m;
                errhan = @(p)sdWGCV(p,M,Rs,beta,omega);
            case 'dp'
                p0 = [-0.5 -0.5];
                sigma = nLevel;
                errhan = @(p)sdDP(p,M,Rs,beta,sigma,m);
                
            otherwise
                error('sdhybr error: No inner solver!')
        end
        
        p = fminunc(errhan,p0);
        % p0 = [5 100];
        % lb = [0,0];
        % ub = [10,100];
        % p = fmincon(errhan,p0,[],[],[],[],lb,ub)
        lambda = p(1);
        alpha = p(2);
       
        %     p
        IR = lambda^2*eye(i) + alpha^2*(Rs'*Rs);
        GIR = M'*M + IR;
        C = GIR\M';
        f = C*vector;
        yx = Z*f;
        
        %s2 = L*yx(N+1:end);
        s2 = yx(N+1:end);
        
        % update L (flexible preconditioner)
        dxj = 2*sqrt(s2(:).^2+thr).^(0.5);
        %dxj(dxj<thr) = eps; % why don't we set this to thr?
        N_d = length(dxj);
        L = spdiags(sqrt(dxj),0:0,N_d,N_d);
        
        % compute x 
        s1 = yx(1:N);
        beta_out = s2;
        x = Q*s1 + X*beta_out;
        
        
        r = b(:) - A*x(:);
        normr = norm(r(:));

        
        if outputparams
            if ~notrue
                if nomask
                   xmask = x;
                else
                   xmask = x(mask);
                end
                output.Enrm(i,1) = norm(xmask(:)-x_truemask(:))/nrmtruemask;
                
                output.s1(:,i) = s1(:);
                output.s2(:,i) = s2(:);
            end
            output.Rnrm(i,1) = normr;
            output.Xnrm(i,1) = norm(x(:));
        end
        
        %    norma = norm([norma B(i,i) B(i+1,i)]);
        
        if isa(A,'function_handle')
            normar = norm(A(r,'transp'));
        else
            normar = norm(A'*r);
        end
        normx = norm(x(:));
        
        % Compute the GCV value used to find the stopping criteria

        GCV(i-1) = sdWGCV(p,M,Rs,beta,1);
        
        if i > 2 && terminate
            %%-------- If GCV curve is flat, we stop -----------------------
            if abs((GCV(i-1)-GCV(i-2)))/GCV(regstart-1) < degflat
                if notrue %Set all the output parameters and return
                    if outputparams
                        output.U = U;
                        output.V = V;
                        output.QV = QV;
                        output.B = B;
                        output.GCVstop = GCV(:);
                        output.iterations = i-1;
                        output.flag = 1;
                        output.alpha = alpha; % Reg Parameter at the (i-1)st iteration
                        %output.Alpha = Alpha(1:i-1); % Reg Parameters
                    end
                    close(h)
                    
                    return;
                else % Flat GCV curve means stop, but continue since have x_true
                    if outputparams
                        iterations_save = i-1;
                        output.iterations = i-1; % GCV says stop at (i-1)st iteration
                        output.flag = 1;
                        output.alpha = alpha; % Reg Parameter at the (i-1)st iteration
                    end
                end
                terminate = 0; % Solution is already found!
                
                %%--- Have warning : Avoid bumps in the GCV curve by using a
                %    window of (mintol+1) iterations --------------------
            elseif warning && length(GCV) > iterations_save + mintol %Passed window
                if GCV(iterations_save) < GCV(iterations_save+1:end)
                    % We should have stopped at iterations_save.
                    if notrue %Set all the output parameters and return
                        if outputparams
                            output.U = U;
                            output.V = V;
                            output.QV = QV;
                            output.B = B;
                            output.GCVstop = GCV(:);
                            output.iterations = iterations_save;
                            output.flag = 2;
                            output.alpha = alpha_save;
                            output.Alpha = Alpha(1:iterations_save); % Reg Parameters
                        end
                        close(h)
                        return;
                    else % GCV says stop at iterations_save, but continue since have x_true
                        if outputparams
                            output.iterations = iterations_save;
                            output.flag = 2;
                            output.alpha = alpha_save;
                        end
                    end
                    terminate = 0; % Solution is already found!
                    
                else % It was just a bump... keep going
                    warning = 0;
                    iterations_save = maxiter;
                    alpha_save = 0;
                end
                
                %% ----- No warning yet: Check GCV function---------------------
            elseif ~warning
                if GCV(i-2) < GCV(i-1) %Potential minimum reached.
                    warning = 1;
                    iterations_save = i-1;
                    alpha_save = alpha;
                end
            end
        end
 
    else
        
        f = M \ vector;
        yx = Z*f;
        %     x = Q*(V*f);
        s1 = yx(1:N);
        s2 = yx(N+1:end);
        beta_out = s2;
        x = Q*s1 + X*beta_out;

        
        if outputparams
            if ~notrue
                %output.Enrm(i,1) = norm(x(:)-x_true(:))/nrmtrue;
                if nomask
                    xmask = x;
                else
                    xmask = x(mask);
                end
                output.Enrm(i,1) = norm(xmask(:)-x_truemask(:))/nrmtruemask;
                
                output.s1(:,i) = s1(:);
                output.s2(:,i) = s2(:);
            end
            output.Rnrm(i,1) = normr;
            output.Xnrm(i,1) = norm(x(:));
        end
    end
    waitbar(i/(maxiter+1), h)
end
close(h)

if isempty(x_out) % GCV did not stop the process, so we reached max. iterations
    x_out = x;
end

%--------------------------------------------
%  The following is needed for RestoreTools:
%
if isa(A, 'psfMatrix')
    x_out = reshape(x, bSize);
end
%
%  End of new stuff needed for RestoreTools
%--------------------------------------------

if outputparams
    output.U = U;
    output.V = V;
    output.QV = QV;
    output.B = B;
    output.GCVstop = GCV(:);
    output.G_ob = G_ob;
    %output.diffx = diff_x(:);
    %output.diffres = diff_res(:);
    %   output.s1 = s1;
    %   output.s2 = s2;
    %  output.coef = yz1;
    if output.alpha == 0
        output.alpha = alpha;
    end
    output.Alpha = Alpha;
    if output.lambda == 0
        output.lambda = lambda;
    end
    output.Lambda = Lambda;
    %output.iterations = iterations_save;
end
end

%% -----------------------SUBFUNCTION---------------------------------------
function omega = findomega(bhat, s, insolv)
%
%   omega = findomega(bhat, s, insolv)
%
%  This function computes a value for the omega parameter.
%
%  The method: Assume the 'optimal' regularization parameter to be the
%  smallest singular value.  Then we take the derivative of the GCV
%  function with respect to alpha, evaluate it at alpha_opt, set the
%  derivative equal to zero and then solve for omega.
%
%  Input:   bhat -  vector U'*b, where U = left singular vectors
%              s -  vector containing the singular values
%         insolv -  inner solver method for HyBR
%
%  Output:     omega - computed value for the omega parameter.

%
%   First assume the 'optimal' regularization parameter to be the smallest
%   singular value.
%

%
% Compute the needed elements for the function.
%
m = length(bhat);
n = length(s);
switch insolv
    case 'tsvd'
        k_opt = n;
        omega = (m*bhat(k_opt)^2) / (k_opt*bhat(k_opt)^2 + 2*bhat(k_opt+1)^2);
        %  omega = ((m/2)*(bhat(k_opt)^2 + bhat(k_opt+1)^2)) / ((k_opt/2)*(bhat(k_opt)^2 +bhat(k_opt+1)^2) + 2*bhat(k_opt+1)^2)
    case 'tikhonov'
        t0 = sum(abs(bhat(n+1:m)).^2);
        alpha = s(end);
        s2 = abs(s) .^ 2;
        alpha2 = alpha^2;
        
        tt = 1 ./ (s2 + alpha2);
        
        t1 = sum(s2 .* tt);
        t2 = abs(bhat(1:n).*alpha.*s) .^2;
        t3 = sum(t2 .* abs((tt.^3)));
        
        t4 = sum((s.*tt) .^2);
        t5 = sum((abs(alpha2*bhat(1:n).*tt)).^2);
        
        v1 = abs(bhat(1:n).*s).^2;
        v2 = sum(v1.* abs((tt.^3)));
        
        %
        % Now compute omega.
        %
        omega = (m*alpha2*v2)/(t1*t3 + t4*(t5 + t0));
        
    otherwise
        error('Unknown solver');
end
end

%% ---------------SUBFUNCTION ---------------------------------------


function nrm = normM(v, M)
if isa(M, 'function_handle')
    Mv = M(v);
else
    Mv = M*v;
end
nrm = sqrt(v'*Mv);
end

function [Q_out,R_out] = qr_lastcol(Q_in, R_in, w)
QTw = Q_in'*w;
QQTw = Q_in*QTw;
new_col = w - QQTw;
beta = norm(new_col);
Q_out = [Q_in, new_col/beta];
k = size(R_in,2);
R_out = [R_in, QTw; zeros(1,k), beta];
end

function err = sdOPT(p,M,R,L,Q,Z,beta,xtrue,X,mask)
%
% err = sdOPT(p,M,R,L,Q,Z,beta,xtrue,mask)
%
% Input:
%
%
% Output:
%
%
%
%

m = size(M,2);
lambda = p(1);
alpha = p(2);
IR = lambda^2*eye(m) + alpha^2*(R'*R);
GIR = M'*M + IR;
vector = zeros(m+1,1);
vector(1) = 1;
vector = beta*vector;
Gtb = M'*vector;
f = GIR\Gtb;
yx = Z*f;
n = length(xtrue(:));
s1 = yx(1:n);
s2 = yx(n+1:end);
beta_out = s2;
s = Q*s1 + X*beta_out;
%s = s1+s2;

if strcmp(mask,{'off'})
    err = norm(s - xtrue);
else
    err = norm(s(mask) - xtrue(mask));
end

end

function err = sdWGCV(p,M,R,beta,omega)
%
% err = sdWGCV(p,M,R,beta,omega)
%
% Input:
%
%
% Output:
%
%
%
%
m = size(M,2); n = size(M,1);
lambda = p(1);
alpha = p(2);
IR = lambda^2*eye(m) + alpha^2*(R'*R);
GIR = M'*M + IR;
vector = zeros(m+1,1);
vector(1) = 1;
vector = beta*vector;
Gtb = M'*vector;
f = GIR\Gtb;
part_r = M*f - vector;
whole_r = part_r'*part_r;

MC = GIR\(M'*M);
tr  = trace(MC);
err = whole_r/(n - omega*tr)^2;
% yx = Z*f;
% n = size(Z,1)/2;
% s1 = Q*yx(1:n);
% s2 = L*yx(n+1:end);
% s = s1+s2;
% err = norm(s(land) - xtrue(land));
end

function err = sdUPRE(p,M,R,beta,sigma)
%
% err = sdUPRE(p,M,R,beta,sigma)
%
% Input:
%
%
% Output:
%
%
%
%
m = size(M,2);
lambda = p(1);
alpha = p(2);
IR = lambda^2*eye(m) + alpha^2*(R'*R);
GIR = M'*M + IR;
vector = zeros(m+1,1);
vector(1) = 1;
vector = beta*vector;
Gtb = M'*vector;
f = GIR\Gtb;
part_r = M*f - vector;
whole_r = part_r'*part_r;

MC = GIR\(M'*M);
tr  = trace(MC);
err = whole_r/m + 2*sigma*tr/m - sigma^2;
% yx = Z*f;
% n = size(Z,1)/2;
% s1 = Q*yx(1:n);
% s2 = L*yx(n+1:end);
% s = s1+s2;
% err = norm(s(land) - xtrue(land));
end

function err = sdDP(p,M,R,beta,sigma,N)
%
% err = sdDP(p,M,R,beta,omega)
%
% Input:
%
%
% Output:
%
%
%
%

tau = 1.1;
m = size(M,2);
lambda = p(1);
alpha = p(2);
IR = lambda^2*eye(m) + alpha^2*(R'*R);
GIR = M'*M + IR;
vector = zeros(m+1,1);
vector(1) = 1;
vector = beta*vector;
Gtb = M'*vector;
f = GIR\Gtb;
part_r = M*f - vector;
whole_r = part_r'*part_r;

% MALENA
err = abs(whole_r - tau*N);

% MC = GIR\(M'*M);
% tr  = trace(MC);

%%% err = abs(whole_r - tau^2*sigma^2*N);

% yx = Z*f;
% n = size(Z,1)/2;
% s1 = Q*yx(1:n);
% s2 = L*yx(n+1:end);
% s = s1+s2;
% err = norm(s(land) - xtrue(land));


%err = abs(norm(part_r)-tau*sigma);
end