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
measurements = table2array(readtable('measurements.csv', 'HeaderLines',1));

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

soma = 0;
for i = 1:2159
    s1 = A * T_b(i);
    s2 = B*[q_dot_solar(i); q_dot_occ(i) ; q_dot_ac(i) ; q_dot_vent(i); T_amb(i)];
    soma = soma + (T_b(i+1) - (s1 + s2))^2;
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

for o = 1:length(q_dot_occ)
    if q_dot_occ(o) <= 0
        lb(o) = -inf;
        ub(o) = inf;
    end
end
    
options = optimoptions('quadprog', 'MaxIterations', 1000);
[x4,fval4,exitflag4] = quadprog(H,c,[],[],Aeq,beq,lb,ub,[],options);

%% Task 4 Alternative, Maybe
N = 2160;
var  = (E2+1)/10;
H_alt = [2*var*eye(N), zeros(N); zeros(N), zeros(N)];
c_alt = [2*var*T_ref*ones(1,N), Phi'*delta_t]';

% Equality Constraints
Aeq_alt = zeros(N, 2*N);
beq_alt = zeros(N,1);
Aeq_alt(1,1) = 1;
beq_alt(1) = T_b1;
for i = 2:N
    Aeq_alt(i,i) = a3*delta_t - 1;
    Aeq_alt(i,i+1) = 1;
    Aeq_alt(i, N+i) = - a2*delta_t;
    beq_alt(i) = delta_t * (a1*q_dot_solar(i) + a2*(q_dot_occ(i) - q_dot_vent(i)) + a3 * T_amb(i));
end

[x4_alt,fval4,exitflag4] = quadprog(H_alt,c_alt,[],[],Aeq_alt,beq_alt,lb,ub,[],options);

%% Test section
aux_opt_Tbp1 = x4(1:2160);
opt_Tb = x4_alt(1:2160);

opt_qac = x4(2161:end);
opt_qac_alt = x4_alt(2161:end);

opt_Tbp1 = zeros(1, 2161);
opt_Tbp1(2:end) = aux_opt_Tbp1;
opt_Tbp1(1) = T_b1;

opt_Tb(2161) = (1-a3*delta_t)*opt_Tb(2160) + a2*delta_t*opt_qac_alt(2160) + delta_t * (a1*q_dot_solar(2160) + a2*(q_dot_occ(2160) - q_dot_vent(2160)) + a3 * T_amb(2160));

subplot(3,1,1);
plot(1:2161, opt_Tbp1, 1:2161, opt_Tb);
legend('Tbp1 in the equation','Tb in the equation')
subplot(3,1,2);
plot(1:2161, opt_Tbp1)
legend('Tbp1 in the equation')
subplot(3,1,3);
plot(1:2161, opt_Tb)
legend('Tbp in the equation')

figure;
subplot(3,1,1);
plot(1:2160, opt_qac, 1:2160, opt_qac_alt);
legend('Tbp1 in the equation','Tb in the equation')
subplot(3,1,2);
plot(1:2160, opt_qac)
legend('Tbp1 in the equation')
subplot(3,1,3);
plot(1:2160, opt_qac_alt)
legend('Tbp in the equation')

temps = [opt_Tbp1' opt_Tb];

qac = [opt_qac opt_qac_alt];