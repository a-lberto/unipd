clc; clear all

%% Implementazione doublet
% Dati iniziali
k = 1;

% Coordinate punto incognito
x0 = -1;
y0 = 1;

%% Individuo valori potenziale e corrente
% Controllo che il punto non sia nella sorgente
if x0==0 && y0==0
    print('Il punto è singolare')
else
    psi0 =-k/2/pi*y0/(x0^2+y0^2)
    phi0 = k/2/pi*x0/(x0^2+y0^2)
end

%% Grafica
% Intervallo funzione potenziale e di corrente in modo da riempire il
% grafico
scale = ceil(10*k);
psi=linspace(-k,k,8*scale);
phi=linspace(-k,k,8*scale);

% Angoli per disegnare le circonferenze
th = linspace(-pi,pi,300);

% Coefficienti
rphi = k/4/pi./phi;
rpsi = k/4/pi./psi;

% Origine
plot(0,0,'k')
hold on

for i=1:length(psi)
    % Circonferenze potenziali
    p1 = plot(rphi(i)+rphi(i)*cos(th),rphi(i)*sin(th),'r');
end

for i=1:length(psi)
    % Circonferenze corrente
    p2 = plot(rpsi(i)*cos(th),-rpsi(i)+rpsi(i)*sin(th),'b');
end

% Punto incognito
plot(x0,y0, 'ko')
rphi0 = k/4/pi/phi0;
rpsi0 = k/4/pi/psi0;

p3 = plot(rphi0+rphi0*cos(th),rphi0*sin(th),'g');
p4 = plot(rpsi0*cos(th),-rpsi0+rpsi0*sin(th),'c');

plot([-3 3],[0 0], 'k-');
plot([0 0],[-3 3], 'k-');

legend([p1 p2 p3 p4], '\Phi', '\Psi', '\Phi_0', '\Psi_0')
axis([-3 3 -3 3])
axis square
hold off
