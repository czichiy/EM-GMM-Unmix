
%# Key Points:
%#	•	Input Parameters:
%#	•	x: Data matrix of size (N*D).
%#	•	a0: Initial guess of mixture proportions, size (L*1).
%#	•	MuE0: Initial guess of means, size (L*D).
%#	•	SigE0: Initial guess of covariance matrices, size (D*D*L).
%#	•	Iter_stop: Stopping criterion for iterations (either as a maximum number or a convergence threshold).
%#	•	Output:
%#	•	a: Trajectory of estimated mixture proportions over iterations.
%#	•	MuE: Trajectory of estimated means over iterations.
%#	•	SigE: Trajectory of estimated covariances over iterations.
%#	•	Lh: Trajectory of log-likelihood values over iterations.
%#	•	Iterations:
%#	•	The algorithm iterates to maximize the log-likelihood function using an Expectation-Maximization approach.
%#	•	It computes responsibilities (py_x), updates parameters (a, MuE, SigE), and calculates the log-likelihood (Lh).





function [a, MuE, SigE,Lh]=UnmixGaussEM(x,a0, MuE0, SigE0, Iter_stop)
% EM Algorithm for unmix of Gaussian processes
% --- Input ---
% x:     (N*D) input random variables data from mixed Gaussian processes
% a0:    (L*1) initial guess of mixture portions
% MuE0:  (L*D) initial guess of means
% SigE0: (D*D*L) inital guess of covariance matrices
% Iter_max: (scaler):
%         If Iter_stop > 1:  number of maximum iteration (by default : 100)
%         If 0< Iter_stop < 1: The ratio over which likelihood donot
%         increase, i.e., (Lh(new)-Lh(old))/Lh(old) < Iter_stop;

% --- Output ---
% a:     (L*Iter_number)   Trajactory of poriton estimation
% MuE:   (L*D*Iter_number) Trajactory of mean estimation
% SigE:  (D*D*Iter_number) Trajactory of covariance estimation
% Lh:    (Iter_number)     Trajactory of likelihood

% Implemented by
% JIN CHIY
% version 2024-07

%% Iteration stop criterion
Iter_p = 0.001;
Iter_max = 100;
if nargin < 5
    Iter_stop = Iter_max;
elseif nargin == 5
    if Iter_stop > 1
        Iter_max = Iter_stop;
    else
        Iter_p = Iter_stop;
    end
end

L = length(a0);
[N,D] = size(x);

%% Initialization
a(:,1) = a0;
MuE(:,:,1) = MuE0;
SigE(:,:,:,1) = SigE0;
Lh(1) = -Inf; % Inf is the standard symbol for representing infinity, while inf is a less conventional notation.

%% Iterations
for t = 2 : Iter_max
    % E-step: Compute responsibilities
    for l = 1 : L
       pyx(:,l) = mvnpdf(x,MuE(l,:,t-1),SigE(:,:,l,t-1));
       py_x(:,l) = a(l,t-1)*pyx(:,l);
    end
       py_x = bsxfun(@rdivide, py_x, sum(py_x, 2)); % Improve computational efficiency by using bsxfun to avoid generating unnecessary ones matrices.
       a(:,t) = mean(py_x)';
 
    for l = 1 : L
    % M-step: Update parameters
    for l = 1:L
        MuE(l,:,t) = sum((py_x(:,l)*ones(1,size(x,2))).*x)/sum(py_x(:,l));
        SigE(:,:,l,t) = (x-ones(N,1)*MuE(l,:,t))'* diag(py_x(:,l))*(x-ones(N,1)*MuE(l,:,t))/sum(py_x(:,l));
    end
    % Compute log-likelihood
    Lh(t) =sum(sum((ones(N,1)* log(a(:,t))').*py_x + log(py_x).*py_x));
    if (Lh(t)-Lh(t-1))/abs(Lh(t-1)) < Iter_p && Iter_stop <1
        break;
    end
end

Lh(1) =Lh(2);
