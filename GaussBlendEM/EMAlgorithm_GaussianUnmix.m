
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
