clear all
close all
clc

%% Parameter Definition
% Daniel Salgado Varela - 5321263
D_a1 = 2;
D_a2 = 6;
D_a3 = 3;

% Rohan Chandrashekar - 5238382
D_b1 = 3;
D_b2 = 8;
D_b3 = 2;

%Parameter Definition
E1 = D_a1 + D_b1;
E2 = D_a2 + D_b2;
E3 = D_a3 + D_b3;
%% Loading Data
measurements = table2array(readtable('measurements_physical.csv', 'HeaderLines',1));

delta_t = 3600;
q_dot_occ = measurements(:,1);  % kW
q_dot_ac = measurements(:,2);  % kW
q_dot_vent = measurements(:,3);  % kW
q_dot_solar = measurements(:,4);  % kW
T_amb = measurements(:,5);  % ºC
T_b = measurements(:,6);  % ºC
Phi = measurements(:,7);  % Eur/KWh

%% Task 3
N = 2159;
Y = zeros(N,1);
phi = zeros(N,3);
for i = 1:2159
Y(i) = T_b(i+1) - T_b(i);
phi(i,:) = delta_t * [q_dot_solar(i) q_dot_occ(i)+q_dot_ac(i)-q_dot_vent(i) T_amb(i)-T_b(i)];
end

H3 = 2 * (phi)' * (phi);
c3 = (-2 * Y' * phi)';

[x3,fval3,exitflag3] = quadprog(H3,c3);
a1 = x3(1);
a2 = x3(2);
a3 = x3(3);

A = 1-a3 * delta_t;
B = delta_t * [a1 a2 a2 -a2 a3];

sum = 0;
for i = 1:2159
    s1 = A * T_b(i);
    s2 = B*[q_dot_solar(i); q_dot_occ(i) ; q_dot_ac(i) ; q_dot_vent(i); T_amb(i)];
    sum = sum + (T_b(i+1) - (s1 + s2))^2;
end


%% Task 4
N = 2160;
T_b1= 22.43;
q_ac_max = 100;
T_min = 15;
T_max = 28;
T_ref = 22;

% Auxiliary Variables
D = (E2 +1) / (10 * (1- a3 * delta_t)^2);
F = (a3 * delta_t - 1)*T_ref - delta_t *(a1 * q_dot_solar ...
    + a2 * (q_dot_occ - q_dot_vent) +a3 * T_amb);
G = (Phi - 2*D*F*a2)*delta_t;

c = [2*D*F ; G];

H = 2* D * [eye(N), -delta_t * a2 * eye(N);
        -delta_t * a2 * eye(N), (delta_t * a2)^2 * eye(N)];

% Equality Constraints
Aeq = zeros(N, 2*N);
beq = zeros(N,1);
Aeq(1,1) = 1;
Aeq(1,N+1) = - a2*delta_t;
beq(1) = (1-a3*delta_t)*T_b1 + delta_t * (a1*q_dot_solar(1) + a2*(q_dot_occ(1) - q_dot_vent(1)) + a3 * T_amb(1));
for i = 2:N
    Aeq(i,i-1) = a3*delta_t - 1;
    Aeq(i,i) = 1;
    Aeq(i, N+i) = - a2*delta_t;
    beq(i) = delta_t * (a1*q_dot_solar(i) + a2*(q_dot_occ(i) - q_dot_vent(i)) + a3 * T_amb(i));
end

% Inequality Constraints
lb = [T_min * ones(1,N), zeros(1,N)]';
ub = [T_max * ones(1,N), q_ac_max * ones(1,N)]';

options = optimoptions('quadprog', 'MaxIterations', 250);
[x4,fval4,exitflag4] = quadprog(H,c,[],[],Aeq,beq,lb,ub,[],options);

