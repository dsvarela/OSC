%% Everything about this program seems fine.
clear all
close all
clc

%% Parameter Definition
% Daniel Salgado Varela - 5321263
D_a1 = 2; D_a2 = 6; D_a3 = 3;

% Rohan Chandrashekar - 5238382
D_b1 = 3; D_b2 = 8; D_b3 = 2;

%Parameter Definition
E1 = D_a1 + D_b1;
E2 = D_a2 + D_b2;
E3 = D_a3 + D_b3;

clear D_a1 D_a2 D_a3 D_b1 D_b2 D_b3

%% Task 1 -  Selecting the Air-Conditioner

pow_X = 4000; %Watts
pow_Y = 2500; %Watts
price_X = 3000; %€
price_Y = 1500; %€

%% Task 1b)
% Problem  Formulation:
%[x1,x2]=[X,Y];
% max (x1*pow_x + x2*pow_y) = min (-x1*pow_x - x2*pow_y)

cTa = [-pow_X; -pow_Y];
% X*price_X + Y*price_Y <= 24000 + 300*E1
% X + Y <= 12
Aa = [price_X, price_Y; 1, 1];
ba = [24000 + 300*E1; 12];


[x1a,fval,exitflag] = linprog(cTa,Aa,ba,[],[],[]);
X = round(x1a(1));
Y = round(x1a(2));

Max_Power = X*pow_X + Y*pow_Y;

%% Task 1c)
Maintenance_X = [200+E2 200+2*E2 200+3*E2 300+4*E2 300+5*E2 ...
    400+5*E2 500+5*E2 600+5*E2 700+5*E2 800+5*E2]';
Maintenance_Y = [50+E3 50+2*E3 100+3*E3 150+4*E3 150+5*E3 200+5*E3 ...
    250+5*E3 300+5*E3 350+5*E3 400+5*E3]';

% Problem  Formulation:
% max (x1*pow_x + x2*pow_y) = min (-x1*pow_x - x2*pow_y)
% X*total_price_X + Y*total_price_Y <= Budget
% X + Y <= 12

% Build Vectors of Total Budget and Cost for units X and Y
total_price_X =  zeros(10,1);
total_price_Y =  zeros(10,1);
Budget =  zeros(10,1);
for N = 1:10
    % Total price of X = 3000 + N year maintenance
    total_price_X(N) = 3000 + sum(Maintenance_X(1:N));
    % Total price of Y = 1500 + N year maintenance
    total_price_Y(N) = 1500 + sum(Maintenance_Y(1:N));
    Budget(N) = 24000 + 300*E1 + (N)*(4000+100*E1);
end


% % Solve Task 1c) as 10 different Linear Programming Problems.
% sol_1c = zeros(10,2);
% for N = 1:10
%     cTB = [-pow_X; -pow_Y];
%     AB = [total_price_X(N), total_price_Y(N); 1, 1];
%     bB = [Budget(N); 12];
%     [sol_1c(N,:)] = linprog(cTB,AB,bB,[],[],[0 0]);
% end

% Task 1c) - single LP
single_c = [-pow_X; -pow_Y;  %Year 1
    -pow_X; -pow_Y;  %Year 2
    -pow_X; -pow_Y;  %Year 3
    -pow_X; -pow_Y;  %Year 4
    -pow_X; -pow_Y;  %Year 5
    -pow_X; -pow_Y;  %Year 6
    -pow_X; -pow_Y;  %Year 7
    -pow_X; -pow_Y;  %Year 8
    -pow_X; -pow_Y;  %Year 9
    -pow_X; -pow_Y]; %Year 10
single_A = zeros(20,20);
single_b = zeros(20,1);
for i = 1:10
    single_A(i,2*i-1) = total_price_X(i);
    single_A(i,2*i) = total_price_Y(i);
    single_A(10+i, 2*i-1) = 1;
    single_A(10+i, 2*i) = 1;
    single_b(i) = Budget(i);
    single_b(10+i) = 12;
end

[x1c] = linprog(single_c,single_A,single_b,[],[]);

%% Post Processing of the Obtained Solutions.
x1c_organized = zeros(10,2);
for i = 1:10
    x1c_organized(i,:) = [x1c(2*i-1), x1c(2*i)];
end
for i = 1:10
    % x and y values are not integers, so they can be rounded up or down,
    % which will result in not optimal solutions (and unfeasible ones)
    
    % x and y agragate all possible combinations of rounded up or down
    % values of our optimal solutions
    x = [floor(x1c_organized(i,1)), ceil(x1c_organized(i,1))];
    y = [floor(x1c_organized(i,2)), ceil(x1c_organized(i,2))];
    units = zeros(2,2);
    price = zeros(2,2);
    max_pow = 0;
    tot_pow = 0;
    
    % we now cycle through all possible solutions pairs x,y to delete the
    % ones that are not within the constraints.
    for j=1:2
        for k = 1:2
            units(j,k) = x(j) + y(k);
            price(j,k) =total_price_X(i)*x(j) + total_price_Y(i)*y(k);
            if units(j,k) > 12 || price(j,k) > Budget(i)
                price(j,k) = 0;
                units(j,k) = 0;
            end
        end
    end
    
    % for the feasible integer solutions found, we cycle through them to
    % see which will generate the most power. That will be our opitmal
    % solution.
    for j=1:2
        for k = 1:2
            if units(j,k) ~= 0
                tot_pow = pow_X*x(j) + pow_Y*y(k);
                if tot_pow >= max_pow
                    max_pow = tot_pow;
                    opt_x = x(j);
                    opt_y = y(k);
                    left_over = Budget(i)-price(j,k);
                end
            end
        end
    end
    
    x1c_organized(i,3) = opt_x;
    x1c_organized(i,4) = opt_y;
    x1c_organized(i,5) = max_pow;
    x1c_organized(i,6) = left_over;
end

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

clear measurements

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

% Mean Square Error Computation for results discussion
mse = 0;
for i = 1:2159
    s1 = A * T_b(i);
    s2 = B*[q_dot_solar(i); q_dot_occ(i) ; q_dot_ac(i) ; q_dot_vent(i); T_amb(i)];
    mse = mse + (T_b(i+1) - (s1 + s2))^2;
end

clear mse

%% Task 4 as Described in Section 4.2.
N = 2160; % Horizon of Steps
T_b1= 22.43; % Given Temperature at step 1
q_ac_max = 100; % Maximum Power Output of Air-Conditioning
T_min = 15; % Minimum Temperature Permmited
T_max = 28; % Maximum Temperature Permmited
T_ref = 22; % Reference Temperature

Phi = Phi/delta_t; %convert to KWh to kW, so it makes units wise

var  = (E2+1)/10; % Cost factor for variations from T_ref.

% Constants Matrix for Square Variables
H = [2*var*eye(N-1), zeros(N-1,N); zeros(N, N-1), zeros(N)];

% Constants Matrix for Linear Variables
c = zeros(2*N-1,1);
for i = 1:N-1
    c(i) = -2*var*T_ref;
    c(N+i-1) = Phi(i)'*delta_t;
end
c(2*N-1) = Phi(N)'*delta_t;

% Lower and Upper Bound Constraints
lb = [T_min * ones(1,N-1), zeros(1,N)]';
ub = [T_max * ones(1,N-1), q_ac_max * ones(1,N)]';

% Cycle through q_dot_occ to find values equal to zero. Lift constraint on
% corresponding lower and upper bounds.
for i = 1:N-1
    if q_dot_occ(i+1) <= 0
        lb(i) = -inf;
        ub(i) = inf;
    end
end

% Equality Constraints, A Matrix and b Vector.
Aeq_alt = zeros(N-1, 2*N-1);
beq_alt = zeros(N-1,1);
Aeq_alt(1,1) = 1;
Aeq_alt(1,N) = - a2*delta_t;
beq_alt(1) = (1-a3*delta_t)*T_b1 + delta_t * (a1*q_dot_solar(1) + a2*(q_dot_occ(1) - q_dot_vent(1)) + a3 * T_amb(1));
for i = 2:N-1
    Aeq_alt(i,i-1) = a3*delta_t - 1;
    Aeq_alt(i,i) = 1;
    Aeq_alt(i, N+i-1) = - a2*delta_t;
    beq_alt(i) = delta_t * (a1*q_dot_solar(i) + a2*(q_dot_occ(i) - q_dot_vent(i)) + a3 * T_amb(i));
end

% Solution to the problem
options = optimoptions('quadprog', 'MaxIterations', 500);
[x4,~] = quadprog(H,c,[],[],Aeq_alt,beq_alt,lb,ub,[],options);

%% Task 4 as Described in Section 4.3.
% 
% % Auxiliary Variables
% var = (E2+1)/10;
% D = (A*Phi/a2 + 2*var*T_ref);
% 
% % Linear constants vector construction
% c_alt = zeros(1,N);
% for i = 1:N-2
%     c_alt(i) = Phi(i)/a2 - D(i+1);
% end
% c_alt(end-1) = Phi(N-1)/a2 - 2*var*T_ref;
% c_alt(end)= Phi(N)*delta_t;
% 
% % Quadratic constants vector construction
% H_alt= [2*var*eye(N-1), zeros(N-1,1); zeros(1,N)];
% 
% % Lower and Upper Bound Constraints
% lb_alt = [T_min * ones(1,N-1), 0]';
% ub_alt = [T_max * ones(1,N-1), q_ac_max]';
% 
% for o = 2:length(q_dot_occ)
%     if q_dot_occ(o) <= 0
%         lb_alt(o) = -inf;
%         ub_alt(o) = inf;
%     end
% end
% 
% % Inequalites matrices construction (done in two separate matrices that 
% % are then joined in the end)
% % for -qack < 0
% A_ineq_min = zeros(N,N);
% b_ineq_min = zeros(1,N);
% 
% % for qack < 100
% A_ineq_max = zeros(N,N);
% b_ineq_max = zeros(1,N);
% for i =  1:N-1
%     if i == 1
%         A_ineq_min(i,i) = 1/(delta_t*abs(a2));
%         b_ineq_min(i) = q_dot_vent(i+1) - q_dot_occ(i+1) + abs(a1/a2)*q_dot_solar(i+1) - abs(a3/a2)*T_amb(i+1) + T_b1*(1+abs(a3)*delta_t)/(delta_t*abs(a2));
%         
%         A_ineq_max(i,i) = -1/(delta_t*abs(a2));
%         b_ineq_max(i) = 100 - q_dot_vent(i+1) + q_dot_occ(i+1) - abs(a1/a2)*q_dot_solar(i+1) + abs(a3/a2)*T_amb(i+1) -T_b1*(1+abs(a3)*delta_t)/(delta_t*abs(a2));
%     else
%         A_ineq_min(i,i-1) = -(1+abs(a3)*delta_t)/(delta_t*abs(a2));
%         A_ineq_min(i,i) = 1/(delta_t*abs(a2));
%         b_ineq_min(i) = q_dot_vent(i+1) - q_dot_occ(i+1) + abs(a1/a2)*q_dot_solar(i+1) - abs(a3/a2)*T_amb(i+1);
%         
%         A_ineq_max(i,i-1) = (1+abs(a3)*delta_t)/(delta_t*abs(a2));
%         A_ineq_max(i,i) = -1/(delta_t*abs(a2));
%         b_ineq_max(i) = 100 - q_dot_vent(i+1) + q_dot_occ(i+1) - abs(a1/a2)*q_dot_solar(i+1) + abs(a3/a2)*T_amb(i+1);
%     end
% end
% A_ineq_min(N,N) = -1;
% b_ineq_min(N) = 0;
% A_ineq_max(N,N) = 1;
% b_ineq_max(N) = 100;
% 
% A_ineq = [A_ineq_min; A_ineq_max];
% b_ineq = [b_ineq_min, b_ineq_max]';
% 
% % Solution to the problem
% options = optimoptions('quadprog', 'MaxIterations', 300);
% [x4_alt,~] = quadprog(H_alt,c_alt,A_ineq,b_ineq,[],[],lb_alt,ub_alt,[],options);

%% Graphing of the Results
T_b_opt = [T_b1 , x4(1:N-1)']';
q_dot_ac_opt = x4(N:end);

hours = 1:N;
figure;
subplot(2,1,1);
plot(hours, T_b_opt);
title('Plot of the Optimized Values for the Building''s Temperature over a Horizon of N=2160 Steps.','Fontsize', 20);
xlabel('Time [h]', 'Fontsize', 14);
ylabel('Temperature [ºC]', 'Fontsize', 14, 'fontweight','bold');
legend({'$T_b$'}, 'Location','northeast', 'Interpreter', 'Latex', 'Fontsize', 20);
set(gca, 'Fontsize', 14)
grid on
subplot(2,1,2);
plot(hours, q_dot_ac_opt);
title('Plot of the Optimized Values for the Building''s Air-Conditioning Power over a Horizon of N=2160 Steps.', 'Fontsize', 20);
legend({'$\dot{q}_{ac}$'}, 'Location','northeast', 'Interpreter', 'Latex', 'Fontsize', 20);
xlabel('Time [h]', 'Fontsize', 14, 'fontweight','bold');
ylabel('Power [kW]', 'Fontsize', 14, 'fontweight','bold');
set(gca, 'Fontsize', 14)
grid on
%% Sum to obtain cost functions

Tbp1 = T_b_opt*(1-a3*delta_t) + delta_t*(a1*q_dot_solar(1:N) + a2*(q_dot_occ(1:N) - q_dot_vent(1:N)+ q_dot_ac_opt) + a3*T_amb(1:N));
cost_per_step = zeros(1,N);
total_cost_fn = 0;
for i=1:N
cost_per_step(i) = q_dot_ac_opt(i)*Phi(i)*delta_t +(E2+1)/10 * (T_b_opt(i)-T_ref)^2;
total_cost_fn = total_cost_fn + cost_per_step(i);
end
figure;
plot(hours, cost_per_step)
title('Plot of the Cost of Air-Conditioning for the Building over a Horizon of N=2160 Steps.', 'Fontsize', 20);
legend({'$Price$'}, 'Location','northeast', 'Interpreter', 'Latex', 'Fontsize', 20);
xlabel('Time [h]', 'Fontsize', 14, 'fontweight','bold');
ylabel('Price [€]', 'Fontsize', 14, 'fontweight','bold');
set(gca, 'Fontsize', 14)
grid on