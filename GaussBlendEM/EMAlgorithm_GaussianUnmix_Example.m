% Test EM: Guassian unmixing with EM algorithm illustration script

% Main funciton: UnmixGaussEM
% Auxilary function: CalculateElliplse (for plotting)

% JIN CHIY
% version 2024-07


clc, clear

%% Three centers of 2-D Gaussian pdf
Mu = [0, 0; 
      3, 3
      2,-4];
%% Three covariance matrix of 2-D Gaussian pdf
Sig1 = [0.9,0;0,0.9];
Sig2 = [1,0.5;0.5,1];
Sig3 = [1,-0.7;-0.7,1];
Sig = cat(3, Sig1, Sig2, Sig3);

%% Number of samples
N = 1000;

%% Mixing proportion [0.5, 0.3, 0.2]
s0 = rand(N,1);
s(s0<0.5) = 1;
s(s0>0.5&s0<0.8) = 2;
s(s0>=0.8)=3;

%% Generate data following the mixed Gaussian distribution

for  i = 1 : N
    x(i,:) =  mvnrnd(Mu(s(i),:), Sig(:,:,s(i)));
end

%% Initializaiton of parameters
a0 = ones(3,1)/3;
MuE0 = [0.5, 0.5;
        1, 1;
        1, -1];
SigE0(:,:,1) = eye(2);
SigE0(:,:,2) = eye(2);
SigE0(:,:,3) = eye(2);

[a, MuE, SigE,Lh]=UnmixGaussEM(x,a0, MuE0, SigE0, 100);

%% Illustration
clr = ['rgm'];

% Original
figure, plot(x(:,1),x(:,2),'.'), title('Original Data - circle: centers, ellipses: 3-sigma region','fontsize',15)
hold on
for l = 1 : 3
   plot(Mu(l,1),Mu(l,2),[clr(l),'o'],'linewidth',2)    
   [U,D] = eig(Sig(:,:,l));
   d = diag(D);
   [ds, idx] = sort(d,'descend');
   ax = 3*sqrt(ds(1));   % 3 \sigma reigon
   bx = 3*sqrt(ds(2));
   tan1 = U(2,idx(1))/U(1,idx(1));
   [EX, EY]=calculateEllipse(Mu(l,1),Mu(l,2), ax, bx, 90-atan(tan1)/pi*180 , 360);
   plot(EX,EY,clr(l));
end
set(gcf, 'color',[1,1,1])

% Estimated
figure, plot(x(:,1),x(:,2),'.'), title('Estimated - circle: centers, ellipses: 3-sigma region','fontsize',15)
hold on
for l = 1 : 3
   plot(Mu(l,1),Mu(l,2),[clr(l),'o'],'linewidth',2)    
   [U,D] = eig(SigE(:,:,l,end));
   d = diag(D);
   [ds, idx] = sort(d,'descend');
   ax = 3*sqrt(ds(1));   % 3 \sigma reigon
   bx = 3*sqrt(ds(2));
   tan1 = U(2,idx(1))/U(1,idx(1));
   [EX, EY]=calculateEllipse(Mu(l,1),Mu(l,2), ax, bx, 90-atan(tan1)/pi*180 , 360);
   plot(EX,EY,clr(l));
   plot(squeeze(MuE(l,1,:)),squeeze(MuE(l,2,:)),['.-',clr(l)])
end
set(gcf, 'color',[1,1,1])

figure,
subplot(121),plot(a','linewidth',2), xlabel('Iterations'), ylabel('Mixing proportion');
subplot(122),plot(Lh,'linewidth',2), xlabel('Iterations'), ylabel('likelihood');
set(gcf, 'color',[1,1,1])

display('Estimated portion:'), display(a(:,end))
display('Estimated centers:'), display(MuE(:,:,end))
display('Estimated covariance:'), display(SigE(:,:,:,end))
