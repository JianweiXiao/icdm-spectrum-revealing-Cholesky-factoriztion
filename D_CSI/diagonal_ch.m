function [L,P] = diagonal_ch(A,r)

% return L factor and column pivoting vector P
% input: positive definite matrix A, the incomplete cholesky step r
% output: L factor and column pivoting vector P. LL' = A(P,P)
n = size(A,1); % n is the size of A
P = 1:n; % initialize permutation vector P
L = zeros(n,n);
for i = 1:r
    [pivot,index] = max(diag(A(i:end,i:end)));
    index = index+i-1;
    P([i,index]) = P([index,i]);
    A([i,index],:) = A([index,i],:);
    A(:,[i,index]) = A(:,[index,i]);
    L([i,index],:) = L([index,i],:);
    L(i,i) = sqrt(pivot);
    L(i+1:end,i) = A(i+1:end,i) / L(i,i);
    A(i+1:end,i+1:end) = A(i+1:end,i+1:end) - L(i+1:end,i) * L(i+1:end,i)'; % update the trailing matrix. 
end
end
    
    