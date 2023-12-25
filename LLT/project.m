clear all
clc
close all

%% =============== Givens ===============
%Lambda = 0:0.1:1;
v_inf = 176;
Lambda = 0.5;
delta = [];
for lambda = Lambda
AR = 5+1/3;
b = 12.192;
S = b^2/AR;
Cr = 2*b/AR/(1+lambda);
Ct = lambda * Cr;
N = 50;
n = 1:N;
theta = (n'*pi / (N+1));
y = -cos(theta)*b/2;
c = Cr - 2*abs(y)*(Cr-Ct)/b;
a0 = 5.5 + (5.5-5.8)/(0-b/2)*abs(y); % must be column
AoA = deg2rad(180/pi); % rad
twist = deg2rad(0); % rad
sweep = 0; % rad
%alpha_0 = deg2rad(0); % rad
%alpha_g = AoA + 2*abs(y)*twist/b;
%alpha = alpha_g - alpha_0;
alpha = deg2rad(5.5+(5.5-3.5)/(0-b/2)*abs(y));



%% =============== Solve ===============
AIC = (4*b./a0./c.*sin(theta*n) + n.*sin(theta*n)./sin(theta));
An = linsolve(AIC, alpha);

CL = An(1) * pi * AR;
delta = [delta sum(n(2:end).*(An(2:end)/ An(1))'.^2)];
CDi = CL^2/pi/AR*(1+delta);
tau = mean(alpha)/An(1) - pi*AR/a0 - 1;
if(AR<4)
  a = a0*cos(sweep) ./ (sqrt(1 + (a0*cos(sweep)/pi/AR).^2) + a0*cos(sweep)/pi/AR);
else
  a = a0*cos(sweep) ./ (1 + a0*cos(sweep)/pi/AR);
end
e = 1/(1+delta);
w = v_inf * sum(n.*An'.*sin(theta*n)./sin(theta), 2);
Gamma = 2*b*v_inf*sum(An'.*(sin(theta*n)), 2);
alpha_ind = w/v_inf;
cl = a.*(alpha);


end % endfor


f1 = figure;
plot(y, Gamma); hold on; grid on; title("\Gamma vs. y", 'interpreter', 'latex'); xlabel("b (m.)", 'interpreter', 'latex'); ylabel("$\{Gamma}$", 'interpreter', 'latex'); xlim([y(1) y(end)]);

f2 = figure;
plot(y, alpha_ind); hold on; grid on; title("$\{alpha}_{ind}$ vs. y", 'interpreter', 'latex'); xlabel("b (m.)", 'interpreter', 'latex'); ylabel("$\{Alpha}_{ind}$ (rad)", 'interpreter', 'latex'); xlim([y(1) y(end)]);

f3 = figure;
plot(y, cl); hold on; grid on; title("Local lift coefficient vs. y", 'interpreter', 'latex'); xlabel("b (m.)", 'interpreter', 'latex'); ylabel("C_l (rad)", 'interpreter', 'latex'); xlim([y(1) y(end)]);

if(length(Lambda)>1)
  figure;
  plot(Lambda, delta); hold on; grid on; title(["$\{alpha}$ vs. $\{lambda}$ for AR = ", num2str(AR)], 'interpreter', 'latex'); xlabel("\{lambda}", "interpreter", "latex"); ylabel("$\{delta}$", "interpreter", "latex");
end

