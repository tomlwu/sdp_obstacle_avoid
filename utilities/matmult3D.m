function M = matmult3D(A,B)
%MATMULT3D multiply 3D matrices A and B such that M(:,:,i) = A(:,:,i)*B*(:,:,i)
%credit to https://www.mathworks.com/matlabcentral/answers/10161-3d-matrix-multiplication
Ap = permute(A,[2,1,4,3]); % (c x a x 1 x Z)
Bp = permute(B,[1,4,2,3]); % (c x 1 x b x Z)
M = Ap .* Bp;              % (c x a x b x Z)
M = sum(M,1);              % (1 x a x b x Z)
M = permute(M,[2,3,4,1]);  % (a x b x Z)
end

