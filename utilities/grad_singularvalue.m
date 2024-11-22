function grad_sigma_k = grad_singularvalue(M,k,varargin)
%GRAD_SINGULARVALUE compute the gradient of the kth singular value of M w.r.t. M.

flagPSD = false;
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'psd'
            flagPSD = true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

[m,n,p] = size(M);
grad_sigma_k = nan(m*n,p);
if p>1 && isscalar(k)
    k = k.*ones(p,1);
end

for i_p = 1:p
    [U,~,V] = svd(M(:,:,i_p));
    k_ip = k(i_p);
    if k_ip <= 0
        J_sigma_k = zeros(m,n);
    else
        J_sigma_k = nan(m,n);
        for i = 1:m
            for j = 1:n
                if ~flagPSD
                    J_sigma_k(i,j) = U(i,k_ip)*V(j,k_ip);
                else
                    J_sigma_k(i,j) = U(i,k_ip)*U(j,k_ip);
                end
            end
        end
    end
    grad_sigma_k(:,i_p) = vec(J_sigma_k);
end

end

