function [G, P] = my_srch(A,r,nb,np)

% A is kernel matrix
% r is target rank
% nb is blocksize
% np is oversampling

n = size(A,1);
P = 1:n; % permutation
W = randn(np,n);
RP = W * A; % apply random projection to A

icol = 1; % denote the process of the loop
while icol <= r % everytime we do a loop. We may add nb columns.
    %
    % column pivoting
    %
    nb = min(nb, r-icol+1); % dealing with the boundary case
    %         if mod(icol,100) == 0 | (mod(icol,100) > mod(icol+nb-1,100))
    %             fprintf('row %6d\n',icol)
    %         end
    [~,~,pivc] = lu(RP(:,icol:end)','vector'); % use RP to decide the column pivoting.
    pivc = pivc(:); % change pivc to a column vector. Denote which columns we move to the front of trailing matrix.
    %%%%%%%%%%%%%%%%%%%
    %pivc = fliplr(pivc);
    %%%%%%%%%%%%%%%%%%%
    msk = find(abs(pivc-(1:length(pivc))')>0); % msk denote all the columns that are not previous columns. i.e. those columns we interchange.
    A(:,msk+icol-1) = A(:,icol+pivc(msk)-1); % interchange the columns in A, PC, PR and RP. RP = W*A, then RP*P=W*P*(P'*A*P). Then
    P(msk+icol-1)  = P(pivc(msk)+icol-1); % A <- P'*A*P, W<-W*P, RP<-RP*P, PC <- PC*P, PR<-PR*P(P'*PR'=PR')
    RP(:,msk+icol-1)= RP(:,pivc(msk)+icol-1);
    A(msk+icol-1,:) = A(icol+pivc(msk)-1,:);
    W(:,msk+icol-1) = W(:,pivc(msk)+icol-1);

    A(icol:n, icol:icol+nb-1) = A(icol:n, icol:icol+nb-1) - A(icol:n, 1:icol-1) * A(icol:icol+nb-1, 1:icol-1)';
    A(icol:icol-1+nb,icol:icol-1+nb) = chol(A(icol:icol-1+nb,icol:icol-1+nb), 'lower');
    
    %A(icol:icol-1+nb,icol:icol-1+nb) = L; % rewrite L in the upper left part of trailing matrix
    A(icol+nb:n,icol:icol-1+nb) = A(icol+nb:n,icol:icol-1+nb) / A(icol:icol-1+nb,icol:icol-1+nb)'; % update the lower left part of trailing mtraix

    %
    % update RP
    if icol+nb-1 < r
        RP(:,icol+nb:end) = RP(:,icol+nb:end) - ...
            (W(:, icol:icol+nb-1) * A(icol:icol+nb-1, icol:icol+nb-1) +  W(:, icol+nb:n)* A(icol+nb:n,icol:icol+nb-1)) * ...
            A(icol+nb:n, icol:icol+nb-1)';
    end
    icol = icol + nb; % increase icol by nb. Record our process.
end
G = tril(A); % compute L factor
G = G(1:n,1:r);
end


