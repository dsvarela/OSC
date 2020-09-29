f = 10;
a_1 = f*a1;
a_2 = f*a2;
a_3 = f*a3;

A = 1-a_3 * delta_t;
B = delta_t * [a_1 a_2 a_2 -a_2 a_3];

tbp1_est = zeros(2159,1);
sum = 0;
for i = 1:2159
    s1 = A * T_b(i);
    s2 = B*[q_dot_solar(i); q_dot_occ(i) ; q_dot_ac(i) ; q_dot_vent(i); T_amb(i)];
    tbp1_est(i) = (s1 + s2);
    sum = sum + (T_b(i+1)-tbp1_est(i))^2;
end

tbp1 = T_b(2:end);

time = 1:2159;

hold on
plot(time, abs(tbp1_est-tbp1));
legend('.9a');