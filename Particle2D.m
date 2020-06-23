clc; clear all;
%%% Simulating a particle moving in a plane using Kalman filter estimate
%%% INPUTS
%dur - duration of simulation in secs
%dt - Time step in secs
%sigma_acc - SD of acceleration ms-2
%sigma_gps - SD of GPS measurement m
dur = 10;
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
P = Q;      %State error cov matrix

%Parameters for plotting
pos = [];
posmeas=[];
poshat=[];
vel=[];
velmeas=[];
velhat=[];

for t=0:dt:dur
    %Input forcing function
    u = [1+sin(t); 2+5*sin(t)];
    
    %Process noise model
    processnoise = sigma_acc*[dt^2/2*randn; dt*randn; dt^2/2*randn; dt*randn];
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
    
    %Parameters for plotting
    pos = [pos; x(1), x(3)];
    posmeas = [posmeas; y(1), y(2)];
    poshat = [poshat; xhat(1), xhat(3)];
    vel = [vel; x(2), x(4)];
    velhat = [velhat; xhat(2), xhat(4)];
    
end

%% Plots
t = 0:dt:dur;
figure()
plot(t,vel(:,1), t, velhat(:,1));
legend('Actual velocity', 'Estimated velocity')
xlabel('Time/ secs'); ylabel('Velocity/ cm s^-1')
title('Variation of actual and estimated x-axis velocity');
grid on;
figure()
plot(t,vel(:,2), t, velhat(:,2))
legend('Actual velocity', 'Estimated velocity')
xlabel('Time/ secs'); ylabel('Velocity/ cm s^-1')
title('Variation of actual and estimated y-axis velocity');
grid on;

figure()
plot(pos(:,1), pos(:,2), poshat(:,1), poshat(:,2),'g','LineWidth',1.5)
hold on;
scatter(posmeas(:,1), posmeas(:,2), 5,'r','filled')
legend('Actual position', 'Estimated position', 'Measured position')
xlabel('x-axis position/ cm'); ylabel('y-axis position/ cm')
title('Variation of position of particle');
grid on;
    