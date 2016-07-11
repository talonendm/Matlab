% KF_opt - using discrete Kalman filter for obtaining 
% multiple parameters estimations dependence on time based upon an observation.
%
% INTRODUCTION
% The problem of applying Kalman filtering for observation over various
% systems is highly interesting. The possibility of using KF for arbitrary
% noise estimations deserve proper attention. 
% In this Matlab file is considered application of optimal Kalman filter
% for estimation of system with some parameters.
%
% PURPOSE
% The purpose of model is to obtain system parameter variance dependence on
% time. Expected that on each iteration of a Kalman filter system parameter
% variance decreases. Assume that observations occur at fixed discrete time
% intervals. Suppose system with state vector varying on linear rule like: x = Ax + Bu
% 
% NOTATION
% A = state transition matrix.
% P = covariance of the state vector estimate. 
% B = input matrix.
% Q = process noise covariance. 
% R = measurement noise covariance.
% H = observation matrix.
%
% OPERATION
% Covariance of the state vector estimate P calculated on each filter
% iteration. Diagonal elements of P related with estimate variance. So 
% dependence of system parameters estimates variations on time is obtained.

% System parameters:
alpha=1;
beta=1; 
sigma=0.04;
sigma0=0.01;
sigma1=0.02;
sigma2=0.03;

% Define system dimensions
n=6;
A = [1 0.5 0 0 0 0; 0 1 0 0 0 0; 0 0 0.5 0 0 0; 0 0 0 0.5 0.5 0; 0 0 0 -0.5 0.823 0; 0 0 0 0 0 1];     % matrix
lambda = eig(A); % Matrix A eigenvalues
% Define system input control functions:
B = [0 0; 0 0; 0.02 0; 0 0.02; 0 -0.01; 0 0];     % matrix
% Define a process noise:
Q = [1 0; 0 1]; % normalized variance matrix
% Define a measurement error
H = [1 0 1 -1 0 -1]; % normalized measurement vector
R = 0.2; % Error variance
% Covariation matrix of estimation errors (Initial state of matrix)
P = [0.25 0 0 0 0 0; 0 sigma0^2 0 0 0 0; 0 0 sigma1^2 0 0 0; 0 0 0 sigma^2 0 0; 0 0 0 0 ((sigma^2)*(alpha^2+beta^2)) 0; 0 0 0 0 0 sigma2^2];     
% Diagonal contains squares of standard deviation variations
% 5-th diagonal element embodies complex error function
%
m=[]; %Vector representing diagonal elements on each step
c=[]; %Matrix representing time dependence of system variations
k=100;

for i=1:k
   m=sqrt(diag(P));
   S = A * P * A' + B*Q*B';
   % Compute Kalman gain factor:
   K = S*H'*inv(H*S*H'+R);
   % Identity matrix 
   I=eye(n);
   P=(I-K*H)*S;
   c(i,:) = m;
end 
c

figure
 hold on
 grid on
 % plot error dependencies on time:
 hz=plot(c(1:end,1),'r.');
 %hz=plot(c(1:end,2),'g.');
 %hz=plot(c(1:end,3),'b.');
 %hz=plot(c(1:end,4),'c.');
 hz=plot(c(1:end,5),'m.');
 %hz=plot(c(1:end,6),'y.');
 title('Covariance Matrix of Diagonal Components Dependence on Time')
 hold off

