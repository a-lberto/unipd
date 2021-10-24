%% Implementazione moto indisturbato
clc; clear all

v = 1; % m/s
alpha = pi/3; %alpha = -pi/2+pi*rand(1);

% Punto a piacere per il quale calcolare il potenziale e la funzione di
% corrente
x0 = 1; % x0 = -3+6*rand(1);
y0 = 1; % y0 = -3+6*rand(1);

% Valori cercati dei potenziali
phi = [-4*v:1:4*v];
psi = [-4*v:1:4*v];

% Punti ascissa per disegnare le rette
x = linspace(-3,3,7);

%% Calcolo potenziale e corrente del punto a piacere
psi0 = v*(y0*cos(alpha)-x0*sin(alpha))
phi0 = v*(x0*cos(alpha)+y0*sin(alpha))

%% Grafica
% Origine
plot(0,0,'k')
hold on

% Linee di corrente
for i=1:length(psi)
    if abs(cos(alpha))>10e-10
        g = @(psi, x) (psi/v + x*sin(alpha))/cos(alpha);
        y = g(psi(i), x);
        p1 = plot(x,y, 'b');
    else
        y = x;
        p1 = plot(zeros(length(y),1)+psi(i)/v, y, 'b');
    end
end

% Linee potenziali
for i=1:length(phi)
    if abs(sin(alpha))>10e-10
        f = @(phi, x) (phi/v - x*cos(alpha))/sin(alpha);
        y = f(phi(i), x);
        p2 = plot(x,y, 'r');
    else
        y = x;
        p2 = plot(zeros(length(y),1)+phi(i)/v, y, 'r');
    end
end

% Stampo punto e linee corrispondenti
plot(x0, y0, 'ko')

if abs(cos(alpha))>10e-10
    g = @(psi, x) (psi/v + x*sin(alpha))/cos(alpha);
    y = g(psi0, x);
    p3 = plot(x,y, 'c');
else
    y = x;
    p3 = plot(zeros(length(y),1)+psi0/v, y, 'c');
end

if abs(sin(alpha))>10e-10
    f = @(phi, x) (phi/v - x*cos(alpha))/sin(alpha);
    y = f(phi0, x);
    p4 = plot(x,y, 'g');
else
    y = x;
    p4 = plot(zeros(length(y),1)+phi0/v, y, 'g');
end

% Grafico
plot(x,zeros(length(x),1), 'k')
plot(zeros(length(x),1),x, 'k')
legend([p1 p2 p3 p4], '\Psi', '\Phi', '\Psi_0', '\Phi_0')
axis square
axis([-3 3 -3 3])
hold off