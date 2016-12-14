function [G,P,Q,R,error1] = my_csi_gaussian_srch(x,alpha,y,m,centering)

n = size(y,1);
%kappa = 0.1;
%y = kappa * y;
%if centering, y = y - 1/n * repmat(sum(y,1),n,1); end % center Y, side information

%training_X = [x,y];
%training_X = x;

%K = exp( - alpha * sqdist(training_X',training_X') );
nb = 50;
np = nb+5;

K = exp( - alpha * sqdist(x,x) );
[G, P] = my_srch(K,m,nb,np);

%P = 1:n;
%K = K(P,P);
%G = my_partial_cholesky(K,m);

temp = G - 1/n * repmat(sum(G,1),n,1);
[Q,R] = qr(temp);
R = R(1:m,1:m);

% error1 : tr(K-G*G')/tr(K) at each step of the decomposition
% error2 : ||Y-Q*Q'*Y||_F^2 / ||Y||_F^2 at each step of the decomposition
error1 = zeros(1,m+1);
%error2 = zeros(1,m+1);
trace_K = trace(K);
temp = trace_K;
%temp_y = y;
%norm_YF2 = sumsqr(y);
error1(1) = 1.0;
%error2(1) = 1.0;
for i = 1:m
    temp = temp - sumsqr(G(i:end,i));
    error1(i+1) = temp / trace_K;
    %temp_y = temp_y - Q(:,i) * (Q(:,i)' * y);
    %error2(i+1) = sumsqr(temp_y) / norm_YF2;
end
end