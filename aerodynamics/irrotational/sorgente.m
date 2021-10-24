%% Implementazione moto generato da sezione filamento sorgente
clc; clear all

% Punto di origine e portata sorgente
x0 = 1; % x0 = -3+6*rand(1);
y0 = 1; % y0 = -3+6*rand(1);
Q = 3;

% Punto dove individuare valori corrente e potenziale
x1 = -2; %x1 = -3+6*rand(1);
y1 = -1; %x2 = -3+6*rand(1);

% % Esempio 1
% x0=0;y0=1;Q=1;
% x1=-2;y1=0;

% % Esempio 2
% x0=1;y0=1;Q=1;
% x1=-2;y1=0;

% % Esempio 3
% x0=0;y0=0;Q=1;
% x1=-2;y1=0;

%% Individuo valori potenziale e corrente
% Controllo che il punto non sia nella sorgente
if x0==x1 && y0==y1
    print('Il punto � singolare')
else
    if x1<x0
        if y1<y0
            psi1=Q/2/pi*(atan((y1-y0)/(x1-x0))-pi)
        else
            psi1=Q/2/pi*(atan((y1-y0)/(x1-x0))+pi)
        end
    else
        psi1=Q/2/pi*atan((y1-y0)/(x1-x0))
    end
   phi1=Q/2/pi*log(sqrt((x0-x1)^2+(y0-y1)^2))
end

% Intervallo funzione potenziale e di corrente in modo da riempire il
% grafico
psi=[-Q/2:1/4:Q/2];
phi=[-Q/3:1/6:Q/3];

%% Grafica
% Larghezza grafico e raggio per costruire le rette dentro al grafico
x=linspace(-3,3,7);
r=(max(x)+abs(x0)+abs(y0))*sqrt(2);

% Angoli per disegnare le circonferenze
th = linspace(-pi,pi,100);

% Origine
plot(0,0,'k')
hold on

for i=1:length(psi)
    % Circonferenze potenziali
    p1 = plot(x0+exp(2*pi*phi(i)/Q)*cos(th), y0+exp(2*pi*phi(i)/Q)*sin(th), 'r');
    
    % Rette corrente
    p2 = plot([x0 x0+r*cos(2*pi*psi(i)/Q)],[y0 y0+r*sin(2*pi*psi(i)/Q)],'b');
end

% Punto incognito
plot(x1,y1, 'bo')

% Linee potenziale incognite
p3 = plot(x0+exp(2*pi*phi1/Q)*cos(th), y0+exp(2*pi*phi1/Q)*sin(th), 'g');
p4 = plot([x0 x0+r*cos(2*pi*psi1/Q)],[y0 y0+r*sin(2*pi*psi1/Q)],'c');

% Grafico
plot(x,zeros(length(x),1), 'k')
plot(zeros(length(x),1),x, 'k')
legend([p1 p2 p3 p4], '\Phi', '\Psi', '\Phi_1', '\Psi_1')
axis square
axis([-3 3 -3 3])
hold off