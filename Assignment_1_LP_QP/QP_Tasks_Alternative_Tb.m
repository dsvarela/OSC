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

clear D_a1 D_a2 D_a3 D_b1 D_b2 D_b3

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

Phi = Phi/delta_t;

% Auxiliary Variables
var = (E2+1)/10;
D = (A*Phi/a2 + 2*var*T_ref);


c_alt = zeros(1,N);
for i = 1:N-2
    c_alt(i) = Phi(i)/a2 - D(i+1);
end
c_alt(end-1) = Phi(N-1)/a2 - 2*var*T_ref;
c_alt(end)= Phi(N)*delta_t;
H_alt= [2*var*eye(N-1), zeros(N-1,1); zeros(1,N)];

% Lower and Upper Bound Constraints
lb_alt = [T_min * ones(1,N-1), 0]';
ub_alt = [T_max * ones(1,N-1), q_ac_max]';

for o = 2:length(q_dot_occ)
    if q_dot_occ(o) <= 0
        lb_alt(o) = -inf;
        ub_alt(o) = inf;
    end
end

% for -qack < 0
A_ineq_min = zeros(N,N);
b_ineq_min = zeros(1,N);

% for qack < 100
A_ineq_max = zeros(N,N);
b_ineq_max = zeros(1,N);
for i =  1:N-1
    if i == 1
        A_ineq_min(i,i) = 1/(delta_t*abs(a2));
        b_ineq_min(i) = q_dot_vent(i+1) - q_dot_occ(i+1) + abs(a1/a2)*q_dot_solar(i+1) - abs(a3/a2)*T_amb(i+1) + T_b1*(1+abs(a3)*delta_t)/(delta_t*abs(a2));
        
        A_ineq_max(i,i) = -1/(delta_t*abs(a2));
        b_ineq_max(i) = 100 - q_dot_vent(i+1) + q_dot_occ(i+1) - abs(a1/a2)*q_dot_solar(i+1) + abs(a3/a2)*T_amb(i+1) -T_b1*(1+abs(a3)*delta_t)/(delta_t*abs(a2));
    else
        A_ineq_min(i,i-1) = -(1+abs(a3)*delta_t)/(delta_t*abs(a2));
        A_ineq_min(i,i) = 1/(delta_t*abs(a2));
        b_ineq_min(i) = q_dot_vent(i+1) - q_dot_occ(i+1) + abs(a1/a2)*q_dot_solar(i+1) - abs(a3/a2)*T_amb(i+1);
        
        A_ineq_max(i,i-1) = (1+abs(a3)*delta_t)/(delta_t*abs(a2));
        A_ineq_max(i,i) = -1/(delta_t*abs(a2));
        b_ineq_max(i) = 100 - q_dot_vent(i+1) + q_dot_occ(i+1) - abs(a1/a2)*q_dot_solar(i+1) + abs(a3/a2)*T_amb(i+1);
    end
end
A_ineq_min(N,N) = -1;
b_ineq_min(N) = 0;
A_ineq_max(N,N) = 1;
b_ineq_max(N) = 100;

A_ineq = [A_ineq_min; A_ineq_max];
b_ineq = [b_ineq_min, b_ineq_max]';
options = optimoptions('quadprog', 'MaxIterations', 300);

[x4_alt,~] = quadprog(H_alt,c_alt,A_ineq,b_ineq,[],[],lb_alt,ub_alt,[],options);

%% Graphing of the Results
aux = delta_t*(a2*(q_dot_occ(1:end-1)- q_dot_vent(1:end-1)) + a1*q_dot_solar(1:end-1) + a3*T_amb(1:end-1));
qac4 = (x4_alt(2:end-1)-A*x4_alt(1:end-2)-aux)/a2/delta_t;

qac4 = [qac4' , x4_alt(end)]'

hours = 1:N;
subplot(2,1,1);
plot(hours, x4_alt(1:end-1));
subplot(2,1,2);
plot(hours, qac4);
%% Sum to obtain cost functions


task4sum = 0;
for i=1:2160
    task4sum = task4sum + qac4(i)*Phi(i)*delta_t +(E2+1)/10 * (x4_alt(i)-T_ref)^2;
end
