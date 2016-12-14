% compare methods
data = load('new_train_images.mat');
data = data.train_images;
label = load('new_train_labels.mat');
label = label.train_labels;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 3000; % number of training samples
m = 200; % target rank
n_test = 3000; % number of test samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 10; % number of labels
centering = 1; % center the data matrix

train_data = data(1:n,:);
train_label = label(1:n,:);
test_data = data(n+1:n+n_test,:);
test_label = label(n+1:n+n_test,:);
x = train_data';
xtest = test_data';

% one-hot encoding
y = zeros(d,n); % y stores side information of training data
for i = 1:n
    y(train_label(i)+1,i) = 1;
end

ytest = zeros(d,n_test); % ytest stores side information of testing data
for i = 1:n_test
    ytest(test_label(i)+1,i) = 1;
end

% LS SVM parameters
alpha = 0.5;
kappa = 1e-4;

% CSI decomposition with two look ahead parameters
tic
[G,P,Q,R,error1,error2] = csi_gaussian(x,alpha,y',m,centering,.99,40,1e-8);
csi_time_40=toc

tic
[GL,PL,QL,RL,error1L,error2L] = csi_gaussian(x,alpha,y',m,centering,.99,80,1e-8);
csi_time_80=toc

% regular Cholesky without look_ahead
tic
[GC,PC,QC,RC,error1C,error2C] = csi_gaussian(x,alpha,y',m,centering,0,0,1e-8);
csi_no_lookahead_time=toc

% Run cholesky with random projection on it to find pivots
tic
[Gsrch,Psrch,Qsrch,Rsrch,error1srch] = my_csi_gaussian_srch(x,alpha,y',m,centering);
csi_srch_time=toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P = randperm(n);
% PL = randperm(n);
% PC = randperm(n);
% Psrch = randperm(n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% compute test set accuracies for all steps of the decomposition
% using an LS SVM
n = size(G,1);
ntest = size(ytest,2);
k = size(G,2);
L = exp( - alpha * sqdist(xtest,x(:,P(1:k))));
testing_errors=zeros(k,1);
Minv=[];
for i=1:k
    r = R(1:i,i);
    Minv = [ Minv zeros(i-1,1); zeros(1,i-1) 1/n/kappa];
    Minv = Minv - ( Minv * r ) * ( r' * Minv ) * inv( 1 + r' * Minv * r ) ;
    Gtalpha = R(1:i,1:i)' * ( Minv * ( Q(:,1:i)' * y(:,P)' ) );
    b =  - mean( G(:,1:i) * ( Gtalpha ) ) + mean(y');
    beta = ( G(1:i,1:i)' ) \ Gtalpha;
    Zhat = L(:, 1:i )*beta + ones(ntest,1)*b;
    distances=zeros(ntest,d);
    for j=1:d
        delta=zeros(1,d); delta(j)=1;
        distances(:,j)= sum( ( Zhat - repmat(delta,ntest,1) ).^2 ,2 );
    end
    [a,Zpred] = min(distances,[],2);
    Ztest = test_label'+1;
    testing_errors(i)=length(find(Ztest-Zpred'~=0))/ntest;
end

% compte test set accuracies for all steps of the decomposition
% using an LS SVM
n = size(GL,1);
ntest = size(ytest,2);
k = size(GL,2);
L = exp( - alpha * sqdist(xtest,x(:,PL(1:k))));
testing_errorsL=zeros(k,1);
Minv=[];
for i=1:k
    r = RL(1:i,i);
    Minv = [ Minv zeros(i-1,1); zeros(1,i-1) 1/n/kappa];
    Minv = Minv - ( Minv * r ) * ( r' * Minv ) * inv( 1 + r' * Minv * r ) ;
    Gtalpha = RL(1:i,1:i)' * ( Minv * ( QL(:,1:i)' * y(:,PL)' ) );
    b =  - mean( GL(:,1:i) * ( Gtalpha ) ) + mean(y');
    beta = ( GL(1:i,1:i)' ) \ Gtalpha;
    Zhat = L(:, 1:i )*beta + ones(ntest,1)*b;
    distances=zeros(ntest,2);
    for j=1:d
        delta=zeros(1,d); delta(j)=1;
        distances(:,j)= sum( ( Zhat - repmat(delta,ntest,1) ).^2 ,2 );
    end
    [a,Zpred] = min(distances,[],2);
    Ztest = test_label'+1;
    testing_errorsL(i)=length(find(Ztest-Zpred'~=0))/ntest;
end

% compte test set accuracies for all steps of the decomposition
% using an LS SVM for the Cholesky decompostion with no look ahead
n = size(GC,1);
ntest = size(ytest,2);
k = size(GC,2);
L = exp( - alpha * sqdist(xtest,x(:,PC(1:k))));
testing_errorsC=zeros(k,1);
Minv=[];
for i=1:k
    r = RC(1:i,i);
    Minv = [ Minv zeros(i-1,1); zeros(1,i-1) 1/n/kappa];
    Minv = Minv - ( Minv * r ) * ( r' * Minv ) * inv( 1 + r' * Minv * r );
    Gtalpha = RC(1:i,1:i)' * ( Minv * ( QC(:,1:i)' * y(:,PC)' ) );
    b =  - mean( GC(:,1:i) * ( Gtalpha ) ) + mean(y');
    beta = ( GC(1:i,1:i)' ) \ Gtalpha;
    Zhat = L(:, 1:i )*beta + ones(ntest,1)*b;
    distances=zeros(ntest,2);
    for j=1:d
        delta=zeros(1,d); delta(j)=1;
        distances(:,j)= sum( ( Zhat - repmat(delta,ntest,1) ).^2 ,2 );
    end
    [a,Zpred] = min(distances,[],2);
    Ztest = test_label'+1;
    testing_errorsC(i)=length(find(Ztest-Zpred'~=0))/ntest;
end

% using srch
n = size(Gsrch,1);
ntest = size(ytest,2);
k = size(Gsrch,2);
L = exp( - alpha * sqdist(xtest,x(:,Psrch(1:k))));
testing_errorssrch=zeros(k,1);
Minv=[];
for i=1:k
    r = Rsrch(1:i,i);
    Minv = [ Minv zeros(i-1,1); zeros(1,i-1) 1/n/kappa];
    Minv = Minv - ( Minv * r ) * ( r' * Minv ) * inv( 1 + r' * Minv * r ) ;
    Gtalpha = Rsrch(1:i,1:i)' * ( Minv * ( Qsrch(:,1:i)' * y(:,Psrch)' ) );
    b =  - mean( Gsrch(:,1:i) * ( Gtalpha ) ) + mean(y');
    beta = ( Gsrch(1:i,1:i)' ) \ Gtalpha;
    Zhat = L(:, 1:i )*beta + ones(ntest,1)*b;
    distances=zeros(ntest,2);
    for j=1:d
        delta=zeros(1,d); delta(j)=1;
        distances(:,j)= sum( ( Zhat - repmat(delta,ntest,1) ).^2 ,2 );
    end
    [a,Zpred] = min(distances,[],2);
    Ztest = test_label'+1;
    testing_errorssrch(i)=length(find(Ztest-Zpred'~=0))/ntest;
end

figure(1)
plot(error1,'b','LineWidth',1.5); hold on;
plot(error1L,'g','LineWidth',1.5); hold on;
plot(error1C,'r','LineWidth',1.5); hold on;
plot(error1srch,'m','LineWidth',1.5); hold off
leg=legend('CSI-40','CSI-80','Cholesky','SRCH');
set(leg,'FontSize',20);
xlabel('number of pivoting training samples used','FontSize',20);
ylabel('approximation error','FontSize',20);
title('approximation of kernel matrix','FontSize',20);

figure(2)
plot(testing_errors,'b','LineWidth',1.5); hold on
plot(testing_errorsL,'g','LineWidth',1.5); hold on
plot(testing_errorsC,'r','LineWidth',1.5); hold on
plot(testing_errorssrch,'m','LineWidth',1.5); hold off
leg=legend('CSI-40','CSI-80','Cholesky','SRCH');
set(leg,'FontSize',20);
xlabel('number of pivoting training samples used','FontSize',20);
ylabel('prediction error','FontSize',20);
title('prediction error comparison','FontSize',20)

[csi_time_40, csi_time_80, csi_no_lookahead_time, csi_srch_time]

