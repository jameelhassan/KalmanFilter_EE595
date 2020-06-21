clc; clear all;
%%% Ackermann steering model of robot simulation
%%% Inputs 
% dur - Duration of simulation (secs)
% dt - time step (secs)
% phi - angle of turn (radians)
%L = length of robot (cm)
dur = 30;
dt = 0.1;
phi = 0.1;
L = 5;

%Noise parameters
measnois_dis = 4;      %Measurement noise of distances
measnois_ang = 0.1;     %Measurement noise of angle in radians
velnois =2;          %Velocity input noise
u = 10;     %Velocity of 10 cm/s

%System model
x = [0;0;0.2];    %Initial condition
xhat = [0;0;0]; %Initial estimatd
A = eye(3);     %System matrix
B = [dt*cos(x(3)); dt*sin(x(3)); tan(phi)*dt/L];    %Input matrix
C = eye(3);     %Output matrix

R = diag([measnois_dis^2, measnois_dis^2, measnois_ang^2]);  %Measurement noise cov mat
Q = velnois^2 * [dt^2, dt^2, tan(phi)*dt^2/L; dt^2, dt^2, tan(phi)*dt^2/L;...
    tan(phi)*dt^2/L, tan(phi)*dt^2/L, (tan(phi)*dt/L)^2];     %Process noise cov mat
P = Q;      % State/Error covariance matrix

%Plotting parameters
pos = [];
poshat = [];
posmeas = [];

for t = 1:dt:dur
    B = [dt*cos(x(3)); dt*sin(x(3)); tan(phi)*dt/L];
    %Process noisy model
    processnois = velnois* [dt*cos(x(3))*randn; dt*sin(x(3))*randn; dt*tan(phi)*randn/L];
    x = A*x + B*u + processnois;
    
    %Measurement noisy model
    measnois = [measnois_dis*randn; measnois_dis*randn; measnois_ang*randn];
    y = C*x + measnois;
    
    %%% Kalman Filter equations %%%
    
    %% Prediction step
    
    % A priori estimate
    xhat = A*xhat + B*u;   
    % State/Error covariance matrix update
    P = A*P*A' + Q;         
    
    %% Update step
    
    %inverse matrix for Kalman gain
    invmat = C*P*C' + R;    
    %Kalman gain
    K = P*C'*inv(invmat); 
    %Error vec for estimate update
    err_vec = y - C*xhat;   
    %Estimate update
    xhat = xhat + K*(err_vec);
    %Covariance of state matrix
    n = length(K*C);
    P = (eye(n) - K*C)*P;
    
    
    % State updates for plotting
    pos = [pos; x(1), x(2), x(3)];
    poshat = [poshat; xhat(1), xhat(2), xhat(3)];
    posmeas = [posmeas; y(1), y(2), y(3)];
    
end

t = 1:dt:dur;
plot(poshat(:,1),poshat(:,2),'g-', 'Linewidth', 1.5)
hold on
plot(pos(:,1),pos(:,2),'color',[0.8500, .250, 0],'Linewidth', 1.5)
scatter(posmeas(:,1), posmeas(:,2),10,'m','filled')