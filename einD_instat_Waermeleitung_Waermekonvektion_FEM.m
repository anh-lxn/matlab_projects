clear; clc; close all;

%% 1) Parameter

l1 = 0.3; % m
l2 = 0.15; % m
l3 = 0.15; % m

% Innenwand
lam1 = 1.4; % W/mK
rho1 = 2300; % kg/m^3
cp1 = 880; % J/kgK

% Dämmung
lam2 = 0.04; % W/mK
rho2 = 30; % kg/m^3
cp2 = 840; % J/kgK

% Außenwand
lam3 = 0.8; % W/mK
rho3 = 1800; % kg/m^3
cp3 = 850; % J/kgK

T0 = 2000 + 273.15; % K
Tu = -273.15 + 273.15; % K
alpha = 25; % W/(m^2 K)

dt = 100.0; % s
t_end = 100000; % s
steps = round(t_end/ dt)  

%% 2) Netz

nume1 = 5;
nume2 = 10;
nume3 = 20;

x1 = linspace(0, l1, nume1 + 1);
x2 = linspace(l1, l1 + l2, nume2 + 1);
x3 = linspace(l1 + l2, l1 + l2 + l3, nume3 + 1);

% Matrix erstellen und Doppelte entfernen
X = [x1, x2(2:end), x3(2:end)]; %% alle Knoten
anz_knoten = length(X);
anz_elemente = anz_knoten - 1;

%% 3) Matrizen

M = zeros(anz_knoten, anz_knoten);   % Massenmatrix
K = zeros(anz_knoten, anz_knoten);   % Steifigeitsmatrix
f = zeros(anz_knoten, 1);            % Lastvektor
T = zeros(anz_knoten, 1);            % Temperaturvektor (Temperaturen an allen Knoten)

%% 4) Assemblierung

for element = 1:anz_elemente
    knoten = [element, element+1];
    xe1 = X(element);
    xe2 = X(element + 1);
    s = xe2 - xe1; % Länge des Elements

    xm = (xe2 + xe1) / 2;

    if xm <= l1
        lame = lam1;
        rhoe = rho1;
        cpe = cp1;
    elseif (xm > l1) && (xm <= l1 + l2)
        lame = lam2;
        rhoe = rho2;
        cpe = cp2;
    else
        lame = lam3;
        rhoe = rho3;
        cpe = cp3;
    end 
    
    Ke = lame / s * [1 -1; -1 1 ];
    Me = rhoe*cpe*s/6 * [2 1; 1 2];
    
    K(knoten, knoten) = K(knoten, knoten) + Ke;
    M(knoten, knoten) = M(knoten, knoten) + Me;
end

%% 5) Konvektion einabuen
K(end,end) = K(end,end) + alpha;
f(end, 1) = f(end, 1) + alpha * Tu;

%% 6) Anfangsbedingung
T(1:end, 1) = T0;

% Speicher
T_hist = zeros(anz_knoten, steps + 1);
T_hist(:, 1) = T;
time = linspace(0, t_end, steps + 1);

% 7) Randbedingungen
T(1, 1) = T0;

% 8) Zeitintegration und Lösung

A = M/dt + K;

for n=1:steps
    b = M/dt * T + f;
    
    A_sys = A;
    b_sys = b;


    A_sys(1,:) = 0;
    A_sys(1,1) = 1;
    b_sys(1,1) = T0;
    T_neu = A_sys \ b_sys; % = inv(A_sys) * b_sys

    T = T_neu;

    T_hist(:, n+1) = T;
end

%% =========================
%  9) PLOTS
%  =========================

% Temperaturprofile zu ausgewählten Zeiten
figure;
hold on; grid on; box on;

idx1 = 1;
idx2 = round(0.1*steps) + 1;
idx3 = round(0.3*steps) + 1;
idx4 = round(0.6*steps) + 1;
idx5 = steps + 1;

plot(X, T_hist(:,idx1)-273.15, 'LineWidth', 1.8, 'DisplayName', sprintf('t = %.1f s', time(idx1)));
plot(X, T_hist(:,idx2)-273.15, 'LineWidth', 1.8, 'DisplayName', sprintf('t = %.1f s', time(idx2)));
plot(X, T_hist(:,idx3)-273.15, 'LineWidth', 1.8, 'DisplayName', sprintf('t = %.1f s', time(idx3)));
plot(X, T_hist(:,idx4)-273.15, 'LineWidth', 1.8, 'DisplayName', sprintf('t = %.1f s', time(idx4)));
plot(X, T_hist(:,idx5)-273.15, 'LineWidth', 1.8, 'DisplayName', sprintf('t = %.1f s', time(idx5)));

xline(l1, '--k', 'Schichtgrenze 1');
xline(l1 + l2, '--k', 'Schichtgrenze 2');

xlabel('x [m]');
ylabel('Temperatur [°C]');
title('Instationäre 1D-Wärmeleitung in 3-Schicht-Wand (FEM)');
legend('Location', 'best');

%% =========================
%  10) ANIMATION
%  =========================

figure;
for n = 1:steps+1
    plot(X, T_hist(:,n)-273.15, '-o' , 'LineWidth', 2);
    grid on; box on;
    xlabel('x [m]');
    ylabel('Temperatur [°C]');
    title(sprintf('Temperaturverlauf, t = %.1f s', time(n)));
    xlim([0 l1+l2+l3]);
    ylim([min(T_hist(:))-273.15-10, max(T_hist(:))-273.15+10]);

    xline(l1, '--k');
    xline(l1 + l2, '--k');

    drawnow;
end
