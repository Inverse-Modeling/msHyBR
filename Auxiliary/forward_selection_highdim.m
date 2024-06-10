function [Xhat,info] = forward_selection_highdim(X, H, Q, R, z, method)

% Input:
% X of size (n,p)
% method = 'BIC', 'AIC' or 'VT' 
% Output:
% hat(X) = [x0 , x1, ...]

p = size(X,2);
n = size(R,2);

% Forwards selection implies that we start with 0 columns in Hhat and
% we add one colum at each iteration. We keep track of the indices with
% the following variables:
indices_left = 1:p;
indices_selected = [];

if strcmp(method,'VT')
    %%%% Compute the best x0
    for i = 1:p
        ws(i) = WSS_highdim(X(:,i), H, Q, R, z);
    end
    [val, ind] = min(ws);
    indices_left(indices_left==ind)=[];
    indices_selected = [indices_selected, ind];
    ws_0 = val;
    ws = [];
    %%%% Propose all xk options and choose one
    for k = 1:p-1 % add one column of Hhat at a time
        for j = 1:size(indices_left,2) % test all columns of X that are left
            ws_help = WSS_highdim(X(:,[indices_selected,indices_left(j)]), H, Q, R, z); 
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
    end
    Xhat = X(:,indices_selected);
    info.indices_selected=indices_selected;
elseif strcmp(method,'BIC')
    % psi = H*Q*H'+R;
    % detpsi = det(psi);
    % logdetpsi = log(detpsi);
    logdetpsi = 0;  % The same value applies for all subsets, so we don't consider this.
    val_0 = Inf;
    %%%% Propose all xk options and choose one
    for k = 1:p 
        BIC_0 = logdetpsi+k*log(n);
        for j = 1:size(indices_left,2)
            ws = WSS_highdim(X(:,[indices_selected,indices_left(j)]), H, Q, R, z);
            ws_vec(j) = ws;
            bic(j) = ws_vec(j) + BIC_0;
        end
        [val, ind] = min(bic);
        if val > val_0 % stop the iterations if no significative vectors are left
            break
        end
        val_0 = val;
        indices_selected = [indices_selected, indices_left(ind)];
        indices_left(indices_left==indices_left(ind))=[];
        bic = [];
        ws_vec = [];
    end
    Xhat = X(:,indices_selected);
    info.indices_selected=indices_selected;
    
elseif strcmp(method,'AIC')
    % psi = H*Q*H'+R;
    % detpsi = det(psi);
    % logdetpsi = log(detpsi);
    logdetpsi = 0;  % The same value applies for all subsets, so we don't consider this.
    val_0 = Inf;
    %%%% Propose all xk options and choose one
    for k = 1:p 
        AIC_0 = logdetpsi+k*2;
        for j = 1:size(indices_left,2)
            ws= WSS_highdim(X(:,[indices_selected,indices_left(j)]), H, Q, R, z);
            ws_vec(j) = ws;
            aic(j) = ws + AIC_0;
        end

        [val, ind] = min(aic);
        if val > val_0 % stop the iterations if no significative vectors are left
            break
        end
        val_0 = val;
        indices_selected = [indices_selected, indices_left(ind)];
        indices_left(indices_left==indices_left(ind))=[];
        aic = [];
        ws_vec = [];
        info.beta= beta_help;
    end
    Xhat = X(:,indices_selected);
    info.indices_selected=indices_selected;
end