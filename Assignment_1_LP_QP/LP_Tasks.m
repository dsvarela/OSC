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

%% Task 1 -  Selecting the Air-Conditioner

pow_X = 4000; %Watts
pow_Y = 2500; %Watts
price_X = 3000; %€
price_Y = 1500; %€

%% a)
% Note: Task 1a) was solved without considering maintenance costs.
%[x1,x2]=[X,Y];
% max (x1*pow_x + x2*pow_y) = min (-x1*pow_x - x2*pow_y)
cTa = [-pow_X; -pow_Y];
% X*price_X + Y*price_Y <= 24000 + 300*E1
% X + Y <= 12
Aa = [price_X, price_Y; 1, 1];
ba = [24000 + 300*E1; 12];


[task1a_sol,fval,exitflag,output,lambda] = linprog(cTa,Aa,ba,[],[],[]);
X = round(task1a_sol(1));
Y = round(task1a_sol(2));

Max_Power = X*pow_X + Y*pow_Y;

%% b)
Maintenance_X = [200+E2 200+2*E2 200+3*E2 300+4*E2 300+5*E2 400+5*E2 500+5*E2 600+5*E2 700+5*E2 800+5*E2];
Maintenance_Y = [50+E3 50+2*E3 100+3*E3 150+4*E3 150+5*E3 200+5*E3 250+5*E3 300+5*E3 350+5*E3 400+5*E3];

% Problem  Formulation:
% max (x1*pow_x + x2*pow_y) = min (-x1*pow_x - x2*pow_y)
% X*total_price_X + Y*total_price_Y <= Budget
% X + Y <= 12

task1b_sol = zeros(10,2);
total_price_X =  zeros(10,1);
total_price_Y =  zeros(10,1);
Budget =  zeros(10,1);
for N = 1:10
    % Total price of X = 3000 + N year maintenance
    total_price_X(N) = 3000 + sum(Maintenance_X(1:N));
    % Total price of Y = 1500 + N year maintenance
    total_price_Y(N) = 1500 + sum(Maintenance_Y(1:N));
    Budget(N) = 24000 + 300*E1 + (N)*(4000+100*E1);
    
    cTB = [-pow_X; -pow_Y];
    AB = [total_price_X(N), total_price_Y(N); 1, 1];
    bB = [Budget(N); 12];
    [task1b_sol(N,:),fval,exitflag,output,lambda] = linprog(cTB,AB,bB,[],[],[0 0]);
end

%Task 1c) - single LP

stoopid_c = [-pow_X; -pow_Y;  %Year 1
            -pow_X; -pow_Y;  %Year 2
            -pow_X; -pow_Y;  %Year 3
            -pow_X; -pow_Y;  %Year 4
            -pow_X; -pow_Y;  %Year 5
            -pow_X; -pow_Y;  %Year 6
            -pow_X; -pow_Y;  %Year 7
            -pow_X; -pow_Y;  %Year 8
            -pow_X; -pow_Y;  %Year 9
            -pow_X; -pow_Y]; %Year 10
stoopid_A = zeros(20,20);
stoopid_b = zeros(20,1);
for i = 1:10
    stoopid_A(i,2*i-1) = total_price_X(i);
    stoopid_A(i,2*i) = total_price_Y(i);
    stoopid_A(10+i, 2*i-1) = 1;
    stoopid_A(10+i, 2*i) = 1;
    stoopid_b(i) = Budget(i);
    stoopid_b(10+i) = 12;
end
lb = zeros(1,20);
[sol,fval,exitflag,output,lambda] = linprog(stoopid_c,stoopid_A,stoopid_b,[],[],lb);