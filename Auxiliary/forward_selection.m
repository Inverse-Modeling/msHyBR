function [Xhat,info] = forward_selection(X, H, Q, R, z, method)
%
%   Model selection method to determine which covariates 
%   (stored in the columns of X) are relevant for a given
%   inverse problem.
%
%   Note: Forwards selection implies that we start with 
%   0 covariates and we add one at each iteration until 
%   some stopping criterion is satisfied.
%
%   Input:
%          H - syste matrix H
%       R, Q - (noise and solution) covariance matrices
%          X - a fixed matrix that includes possible covariates
%          z - right-hand side
%     method - model selection criterion: 'BIC', 'AIC' or 'VT' 
%
%   Output:
%       Xhat - Matrix whose columns are the selected relevant 
%       covariates (in order)
%      info.indices_selected - indeces corresponding to the solected
%      covariates in the original matrix X.
%


[n,p] = size(X);

% We keep track of the indices with the following variables:
indices_left = 1:p;
indices_selected = [];

if strcmp(method,'VT')
    %%%% Compute the best x0
    for i = 1:p
        ws(i) = WSS(X(:,i), H, Q, R, z);
    end
    [val, ind] = min(ws);
    indices_left(indices_left==ind)=[];
    indices_selected = [indices_selected, ind];
    ws_0 = val;
    ws = [];
    beta = [];
    %%%% Propose all xk options and choose one
    for k = 1:p-1 % add one column of Hhat at a time
        for j = 1:size(indices_left,2) % test all columns of X that are left
            [ws_help, beta_help ] = WSS(X(:,[indices_selected,indices_left(j)]), H, Q, R, z); 
            beta(:,j)= beta_help;
            ws(j) = ws_help;
            vt(j) = nu(ws_0,ws(j),n,1,size(indices_selected,2));
        end
        % F distribution with q and n - p - q degrees of freedom
        critical_value = finv(0.95,1,n - size(indices_selected,2) - 1);
        [val, ind] = max(vt);
        if val < critical_value % stop the iterations if no significative vectors are left
            break
        end
        ws_0 = ws(ind);
        indices_selected = [indices_selected, indices_left(ind)];
        indices_left(indices_left==indices_left(ind))=[];
        vt = [];
        ws = [];
        beta =[];
        info.beta = beta_help;
    end
    Xhat = X(:,indices_selected);
    info.indices_selected=indices_selected;
elseif strcmp(method,'BIC')
    %%%%
    psi = H*Q*H'+R;
    detpsi = det(psi);
    logdetpsi = log(detpsi);
    val_0 = Inf;
    %%%% Propose all xk options and choose one
    for k = 1:p 
        BIC_0 = logdetpsi+k*log(n);
        for j = 1:size(indices_left,2)
            ws(j) = WSS(X(:,[indices_selected,indices_left(j)]), H, Q, R, z);
            bic(j) = ws(j) + BIC_0;
        end

        [val, ind] = min(bic);
        if val > val_0 % stop the iterations if no significative vectors are left
            Xhat = X(:,indices_selected);
            info.indices_selected=indices_selected;
            break
        end

        val_0 = val;
        indices_selected = [indices_selected, indices_left(ind)]
        indices_left(indices_left==indices_left(ind))=[]
        bic = [];
        ws = [];
    end
    Xhat = X(:,indices_selected);
    info.indices_selected=indices_selected;
elseif strcmp(method,'AIC')
    %%%%
    psi = H*Q*H'+R;
    detpsi = det(psi);
    logdetpsi = log(detpsi);
    val_0 = Inf;
    %%%% Propose all xk options and choose one
    for k = 1:p 
        AIC_0 = logdetpsi+k*2;
        for j = 1:size(indices_left,2)
            ws(j) = WSS(X(:,[indices_selected,indices_left(j)]), H, Q, R, z);
            aic(j) = ws(j) + AIC_0;
        end

        [val, ind] = min(aic);
        if val > val_0 % stop the iterations if no significative vectors are left
            Hhat = H(:,indices_selected);
            break
        end

        val_0 = val;
        indices_selected = [indices_selected, indices_left(ind)]
        indices_left(indices_left==indices_left(ind))=[]
        aic = [];
        ws = [];
    end
    Xhat = X(:,indices_selected);
    info.indices_selected=indices_selected;
end
end

function p_out = nu(w0,w1,n,q,p)
% for forward selection q = 1;
p_out = ((w0-w1)*(n-p-q))/(q*w1);
end