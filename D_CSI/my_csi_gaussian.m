function [G,P, Q, R] = my_csi_gaussian(x,alpha,Y,m,centering)
% used for qrcp
% INPUT
% x  : input data nvar x n
% alpha : kernel parameter of the RBF kernel
% Y  : target vector n x d
% m  : maximal rank
% centering : 1 if centering, 0 otherwise (suggested: 1)
%
% OUTPUT
% G : Cholesky decomposition -> K(P,P) is approximated by G*G' (n by p)
% P : permutation matrix (1 by n)
% Q,R : QR decomposition of G (or center(G) if centering)
%
% Copyright (c) Jianwei Xiao, 2016.2.

% K = exp( - alpha * sqdist(x,x) ) + eye(size(x,2)) * 1e-6;
K = exp( - alpha * sqdist(x,x) );
[G, P, Q, R] = my_qrcp_rankversion(K, Y, m, centering);