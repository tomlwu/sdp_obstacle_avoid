function lambda_1 = max_eigen(M)
%MAX_EIGEN Find largest eigenvalues of a set of matrices M(m,n,d)

lambda_1 = zeros(size(M,3),1);
for i = 1:size(M,3)
    lambda = sort(eig(M(:,:,i)),'descend');
    lambda_1(i) = max(lambda);
end

end

