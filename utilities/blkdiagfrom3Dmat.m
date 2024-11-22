function A = blkdiagfrom3Dmat(A3d)
%BLKDIAGFROM3DMAT given a [d1xd2xN] matrix A3d construct a [d1Nxd2N] matrix A 
% whose d1xd2 diagnal blocks are each layer of A3D and zeros elsewhere.

N = size(A3d,3);
d1 = size(A3d,1);
d2 = size(A3d,2);
A = zeros(d1*N,d2*N);

for i = 1:N
    A((i-1)*d1+1:i*d1,(i-1)*d2+1:i*d2)=A3d(:,:,i);
end
end

