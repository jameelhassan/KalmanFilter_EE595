clc; clear all;
%%% Simulating a particle moving in a plane using Kalman filter estimate
%%% INPUTS
%dur - duration of simulation in secs
%dt - Time step in secs
%sigma_acc - SD of acceleration ms-2
%sigma_gps - SD of GPS measurement m
dur = 30;
dt = 0.1;
sigma_acc = 0.5;
sigma_dis = 2;

%System Model
A = [1 dt 0 0; 0 1 0 0; 0 0 1 dt; 0 0 0 1];  %Transition/System matrix
B = [dt^2/2 0; dt 0; 0 dt^2/2; 0 dt];       %Input matrix
H = [1 0 0 0; 0 0 1 0];                     %Output matrix
x = [0;0;0;0];          %Initial x
xhat = x;               %Initial extimate

R = diag([sigma_dis^2 sigma_dis^2]);        %Measurement noise cov mat
Q = sigma_acc^2*[dt^4/4 dt^3/2 0 0; dt^3/2 dt^2 0 0; 0 0 dt^4/4 dt^3/2;...
    0 0 dt^3/2 dt^2];                       %Process noise cov mat
P = Q;

for t=1:dt:dur
    %Input forcing function
    u = [1+sin(t); 2+5*sin(t)];
    
    %Process noise model
    processnoise = sigma_acc*[dt^2/2*randn 0; dt*randn 0; 0 dt^2/2*randn; 0 dt*randn];
    x = A*x +B*u + processnoise;
    
    %Measurement noise model
    measnoise = sigma_dis*[randn; randn];
    y = H*x + measnoise;
    
    %%% Kalman Filter equations %%%
    
    %% Prediction step
    
    % A priori estimate
    xhat = A*xhat + B*u;   
    % State/Error covariance matrix update
    P = A*P*A' + Q;         
    
    %% Update step
    
    %inverse matrix for Kalman gain
    invmat = H*P*H' + R;    
    %Kalman gain
    K = P*H'*inv(invmat); 
    %Error vec for estimate update
    err_vec = y - H*xhat;   
    %Estimate update
    xhat = xhat + K*(err_vec);
    %Covariance of state matrix
    n = length(K*H);
    P = (eye(n) - K*H)*P;
    
    
    
end
    