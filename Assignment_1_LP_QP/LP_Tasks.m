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