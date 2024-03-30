Cm = 1.2; % N*m/A
J = 0.04; % kg*m^2
Ce = 1.2; % V*s
Re = 4; % Om
Rf = 0.9; % Om
Le = 0.01; % Henrys
Lf = 0.01; % Henrys
b = 0.1; % N*m/(rad/s)
Mmin = 1; % N*m
Mmax = 1.1; % N*m

kampl = 10; % силовой преобразователь
ktach = 0.9; % V*s тахогенератор
RF = 20*10^3; % Om сопротивление фильтра
CF = 5.1*10^-6; % F емкость фильтра
% R1 = 10*10^3; % Om
R1 = 680; % Om
R2 = 2*10^3; % Om
R3 = R2;
TF = RF*CF;
Rrot = Re+Rf;
Lrot = Le+Lf;
Trot = Lrot/Rrot;
TM = J*Rrot/(Ce*Cm);
Umax = 151.5539;
Umin = 13.7304;

