clc;
clear all;

%Inputs
%dt - time step (secs)
%dur - duration (secs)
dt = 0.1;
dur = 8;
measnois = 10;      %Measurement noise in postion
accelnois = 1;      %Noise in acceleration

%Initial condition
u = 10; %(cm/s^2)

%System model
A = [1 dt; 0 1];        %Transition matrix
B = [dt^2/2; dt];       %Input matrix
C = [1 0];              %Outut matrix
x = [0; 0];             %Initial state
xhat = x + randn;       %Initial estimate

Q = accelnois^2*[dt^4/4 dt^3/2; dt^3/2 dt^2];   %Process noise cov mat
R = measnois^2;                                 %Measurement noise cov mat
P = Q;                                          %State cov matrix

%Matrices holding past data
pos = [];       %True position
poshat = [];    %Estimated position
posmeas = [];   %Measured position
vel = [];       %Actual velocity
velhat = [];    %Estimated velocity
% kmat = [];

%Kalman loop
for t = 0: dt : dur
    %Process model with noise
    processnois = accelnois * [(dt^2/2)*randn; dt*randn];
    x = A*x + B*u + processnois;

    %Measurement model with noise
    measnoise = measnois * randn;
    y = C*x + measnoise;
    
    %%% Kalman Filter equations %%%%
    
    %Predict step
    xhat = A*xhat + B*u;
    P = A*P*A' + Q; 
    
    %Update step
    invmat = C*P*C' + R;
    K = P*C'*inv(invmat);
    x = xhat + K*(y - C*xhat);
    n = length(K*C);
    P = (eye(n)-K*C)*P;
    
    %Plotting parameters
    pos = [pos x(1)];
    posmeas = [posmeas y];
    poshat = [poshat xhat(1)];
    vel = [vel x(2)];
    velhat = [velhat xhat(2)];
%     kmat = [kmat K];
end
    
t = 0: dt : dur;
figure();
plot(t,pos,'b.-', t,posmeas,'r.-', t,poshat,'g.-');
legend('True position', 'Measured position', 'Estimated position');
grid;
xlabel('Time/sec');
ylabel('Position/cm');
title('Robot position (True, Measured and Estimated)');

figure();
plot(t, pos-posmeas,'r.-', t, pos-poshat,'g.-');
legend('Measured position error', 'Estimated position error');
grid;
xlabel('Time/ sec');
ylabel('Position/ cm');
title('Robot position errors');

figure();
plot(t,vel,'b.-', t,velhat,'g.-');
legend('True velocity', 'Estimated velocity');
grid;
xlabel('Time/ sec');
ylabel('Velocity/ cm s^-1');
title('Robot Velocity (True and Estimated)');

figure();
plot(t, vel-velhat);
grid;
xlabel('Time/ sec');
ylabel('Velocity/ cm s^-1');
title('Robot Velocity error');


