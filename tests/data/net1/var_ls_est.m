function [ mu, para] = var_ls_est( panel, lb_pos)
%VAR_LS_EST estimate the parameter of a var process using least square 
%   Detailed explanation goes here

mu = mean(panel, 2);
tl = size(panel, 2);

panel = panel - repmat(mu, 1, tl);

p = length(lb_pos);
dim = size(panel, 1);
lb = max(lb_pos);
T = tl - lb;

X = zeros(p * dim,  T);
for ii = 1:p
    X(dim*(ii-1)+1:dim*ii, :) = ...
        panel(:, lb-lb_pos(ii)+1:lb-lb_pos(ii)+T);
end

y0 = panel(:, lb+1:end);
y0 = y0(:);

para_all = regress(y0, kron(X', eye(dim)));

para = cell(1, p);
for ii = 1:p
    para{ii} = reshape(para_all(dim*dim*(ii-1)+1:dim*dim*ii), dim, dim);
end

end

