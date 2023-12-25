clear all
clear
clc, clc, clc, clc, clc
close all

h = 0.457; % height of the wind tunnel test section
A = h*h; % Area of the wind tunnel test section
c = 0.152; % chord
b = h; % span
t = 0.12*c; % max thickness

%% Read Data
data = importdata('data.csv').data;
alpha = data(:,1); p_inf_s=data(:,2); p_inf=data(:,3); xu=data(:,40); xl=data(:,41); xr=data(:,42);
xu = xu(~isnan(xu)); xl = xl(~isnan(xl)); xr = xr(~isnan(xr)); xu=xu'; xl=xl'; xr=xr';
Pu = data(:, 4:15); % Pressure on upper surface
Pl = data(:, 16:26); % Pressure on lower surface
Pr = data(:, 27:39); % Rake Pressure Readings

ex(6) = experiment();

f1 = figure;
colors = {'b','r','y','g','k',[.73 .2 .9]};

for(i=1:6)
  k = 1+(i-1)*4;
  ex(i).alpha = mean(alpha(k:k+3));
  ex(i).p_inf_s = mean(p_inf_s(k:k+3));
  ex(i).p_inf = mean(p_inf(k:k+3));
  ex(i).pu = mean(Pu(k:k+3, :));
  ex(i).pl = mean(Pl(k:k+3, :));
  ex(i).pw = mean(Pr(k:k+3, :));

  % calculate Cp
  ex(i).cpu = (ex(i).pu-ex(i).p_inf_s) ./ (ex(i).p_inf-ex(i).p_inf_s);
  ex(i).cpl = (ex(i).pl-ex(i).p_inf_s) ./ (ex(i).p_inf-ex(i).p_inf_s);

  %% calculate Cn
  for(j=1:10) % lower
    ex(i).cn = ex(i).cn + (xl(j+1)-xl(j))*(ex(i).cpl(j+1)+ex(i).cpl(j)) / 2.0;
  end
  for(j=1:11) % upper
    ex(i).cn = ex(i).cn - (xu(j+1)-xu(j))*(ex(i).cpu(j+1)+ex(i).cpu(j)) / 2.0;
  end
  %% Cl
  ex(i).cl = ex(i).cn*cos(ex(i).alpha*pi/180);

  %% Cm
  xref = 0.25;
  for j=1:10 % lower
    ex(i).cm = ex(i).cm + (ex(i).cpl(j+1)+ex(i).cpl(j))/2*(xl(j+1)-xl(j))*(xl(j)+(xl(j+1)-xl(j))/(ex(i).cpl(j+1)+ex(i).cpl(j))*ex(i).cpl(j+1)-xref);
  end
  for j=1:11 % upper
    ex(i).cm = ex(i).cm - (ex(i).cpu(j+1)+ex(i).cpu(j))/2*(xu(j+1)-xu(j))*(xu(j)+(xu(j+1)-xu(j))/(ex(i).cpu(j+1)+ex(i).cpu(j))*ex(i).cpu(j+1)-xref);
  end

  %% Cd
  P_T = max(ex(i).pw(1), ex(i).pw(end)); % assume the freestream total pressure same as the reading at the furthest point in the rake
  u_v_inf = sqrt((ex(i).pw-ex(i).p_inf_s) ./ (P_T - ex(i).p_inf_s));
  fun = u_v_inf.*(1-u_v_inf);
  for(j=1:12)
    ex(i).cd = ex(i).cd + (fun(j+1)+fun(j))/2 * (xr(j+1) - xr(j)) / 1000;
  end
  ex(i).cd = 2/c * ex(i).cd;
  ex(i).cl_theory = 2*pi*ex(i).alpha*pi/180; % cl from thin airfoil theory

  %% Wind Tunnel Corrections
  K1 = 0.76; % K1 is a constant (for an aerofoil mounted across the entire width of a wind tunnel test section, K1 = 0.76)
  V = 0.7*t*c*b;
  epsilon_sb = K1*V/A^(3/2); % solid blockage ratio
  epsilon_wb = c/2/h*ex(i).cd; % wake blockage ratio
  sigma = pi^2/48*(c/h)^2; % σ accounts for the curvature effect
  epsilon = epsilon_sb + epsilon_wb; % total blockage ratio

  if(epsilon>0.05)
    ex(i).cl = ex(i).cl * (1-sigma-2*epsilon);
    ex(i).cd = ex(i).cd * (1-3*epsilon_sb-2*epsilon_wb);
  end
  fprintf("\nFor alpha = %d deg:\n  Cn = %d\n  Cl = %d\n  Cm = %d\n  Cd = %d\n  epsilon = %d\n\n", ex(i).alpha, ex(i).cn, ex(i).cl, ex(i).cm, ex(i).cd, epsilon);

  %% Plot Figures
  subplot(2,1,1); hold on;
  plot(xu, ex(i).cpu, 'color', colors{i}, '-x');
  subplot(2,1,2); hold on;
  plot(xl, ex(i).cpl, 'color', colors{i}, '-o');
end

subplot(2,1,1); hold on; grid on; xlabel("$x/c$", 'interpreter', 'latex'); ylabel("$C_{p,u}$", 'interpreter', 'latex'); title(strcat("$C_p$ Upper vs. $x/c$ for different $\{alpha}$"), 'interpreter', 'latex'); set(gca, 'YDir','reverse'); set(gca, "FontSize", 16);xlim([0 1]); ylim([-4.2, .4]); leg1 = [['$c_{p,u}$ $\alpha=$';'$c_{p,u}$ $\alpha=$';'$c_{p,u}$ $\alpha=$';'$c_{p,u}$ $\alpha=$';'$c_{p,u}$ $\alpha=$';'$c_{p,u}$ $\alpha=$';], [num2str(ex(1).alpha);num2str(ex(2).alpha);num2str(ex(3).alpha);num2str(ex(4).alpha);num2str(ex(5).alpha);num2str(ex(6).alpha)]];
legend(leg1, 'interpreter', 'latex');

subplot(2,1,2); hold on; grid on; xlabel("$x/c$", 'interpreter', 'latex'); ylabel("$C_{p,l}$", 'interpreter', 'latex'); title(strcat("$C_p$ Lower vs. $x/c$ for different $\{alpha}$"), 'interpreter', 'latex'); set(gca, 'YDir','reverse'); set(gca, "FontSize", 16);xlim([0 1]); ylim([-.6, 1.3]); leg2 = [['$c_{p,l}$ $\alpha=$';'$c_{p,l}$ $\alpha=$';'$c_{p,l}$ $\alpha=$';'$c_{p,l}$ $\alpha=$';'$c_{p,l}$ $\alpha=$';'$c_{p,l}$ $\alpha=$'], [num2str(ex(1).alpha);num2str(ex(2).alpha);num2str(ex(3).alpha);num2str(ex(4).alpha);num2str(ex(5).alpha);num2str(ex(6).alpha)]];
legend(leg2, 'interpreter', 'latex', 'location', 'southeast');

%% Load XFoil data
xfoil = importdata("xFoil_data").data(6:end,:);
f2 = figure;
subplot(2,2,1); xlabel("$\{alpha}$", 'interpreter', 'latex'); ylabel("$C_l$", 'interpreter', 'latex'); title("$C_l$ vs. $\{alpha}$", 'interpreter', 'latex'); set(gca, "FontSize", 16); hold on; grid on;
plot([ex(1:6).alpha]', [ex(1:6).cl]', 'r-x'); plot(xfoil(:,1), xfoil(:,2), 'b'); plot([ex(1:6).alpha]', [ex(1:6).cl_theory]', 'k-o');
legend("Wind Tunnel Data", "Xfoil Analysis", "Thin Airfoil Theory", 'Location','northwest');
xlim([0 15]);

%f3 = figure; hold on; grid on;
subplot(2,2,2); xlabel("$\{alpha}$", 'interpreter', 'latex'); ylabel("$C_d$", 'interpreter', 'latex'); title("$C_d$ vs. $\{alpha}$", 'interpreter', 'latex'); set(gca, "FontSize", 16); xlim([0 15]); ylim([-0.05 0.25]); hold on; grid on;
plot([ex(1:6).alpha]', [ex(1:6).cd]', 'r-x'); plot(xfoil(:,1), xfoil(:,3), 'b'); plot([ex(1:6).alpha]', zeros(1,6), 'k-o');
legend("Wind Tunnel Data", "Xfoil Analysis", "Thin Airfoil Theory", 'Location','northwest');

%f4 = figure; hold on; grid on;
subplot(2,2,3); xlabel("$\{alpha}$", 'interpreter', 'latex'); ylabel("$C_m$", 'interpreter', 'latex'); title("$C_m$ vs. $\{alpha}$", 'interpreter', 'latex'); set(gca, "FontSize", 16); xlim([0 15]); hold on; grid on;
plot([ex(1:6).alpha]', [ex(1:6).cm]', 'r-x'); plot(xfoil(:,1), xfoil(:,5), 'b'); plot([ex(1:6).alpha]', zeros(1,6), 'k-o');
legend("Wind Tunnel Data", "Xfoil Analysis", "Thin Airfoil Theory", 'Location','northwest');

%f5 = figure; hold on; grid on;
subplot(2,2,4); xlabel("$\{alpha}$", 'interpreter', 'latex'); ylabel("$C_l/C_d$", 'interpreter', 'latex'); title("Glide ratio $C_l/C_d$ vs. $\{alpha}$", 'interpreter', 'latex'); set(gca, "FontSize", 16); hold on; grid on;
plot([ex(1:6).alpha]', [ex(1:6).cl]'./[ex(1:6).cd]', 'r-x'); plot(xfoil(:,1), xfoil(:,2)./xfoil(:,3), 'b'); xlim([0 15]);
legend("Wind Tunnel Data", "Xfoil Analysis");

print(f1, "f1.pdf", "-r600", '-fillpage', "-landscape");
print(f2, "f2.pdf", "-r600", "-fillpage", "-landscape");
