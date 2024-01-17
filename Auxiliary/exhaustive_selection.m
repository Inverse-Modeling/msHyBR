function [Xhat,info] = exhaustive_selection(X, H, Q, R, z, method)

% Input:
% X of size (n,p)
% method = 'BIC', 'AIC' 
% Output:
% hat(X) = [x0 , x1, ...]

[n,p] = size(X);

results_WSS = cell(2^p-1,3);
counter = 1;
results_WSS{counter,1} = 'WSS';
results_WSS{counter,2} = 'lnQ';
if strcmp(method,'BIC')
    results_WSS{counter,3} = 'BIC';
elseif strcmp(method,'AIC')
    results_WSS{counter,3} = 'AIC';
end
results_WSS{counter,4} = 'subset indices';

counter = counter + 1;

psi=H*Q*H'+R;
detpsi = det(psi); % This is basically zero!

for k = 1:p
    subsets = combnk(1:p, k);
    for i = 1:size(subsets,1)
        results_WSS{counter,1} = WSS(X(:,subsets(i,:)), H, Q, R, z);
        results_WSS{counter,2} = log(detpsi);
        if strcmp(method,'BIC')
            results_WSS{counter,3} = results_WSS{counter,2}+results_WSS{counter,1} + k*log(n);
        elseif strcmp(method,'AIC')
            results_WSS{counter,3} = results_WSS{counter,2}+results_WSS{counter,1} + 2*k;
        end
        results_WSS{counter,4} = subsets(i,:);
        counter = counter + 1;
    end
end

[val_ind,min_ind] = min(cell2mat(results_WSS(2:end,3)));
indices_selected=cell2mat(results_WSS(min_ind+1,4));

Xhat = X(:,indices_selected);
info.indices_selected=indices_selected;
info.WSStable = results_WSS;
