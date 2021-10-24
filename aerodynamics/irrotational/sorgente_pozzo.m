clc; clear all;

% Implementazione pozzo e sorgente
a = 1;
Q = 1;

% Coordinate sorgente e pozzo
xS = -a;
yS = 0;
xP = a;
yP = 0;

% Coordinate punto incognito
x0 = 1;
y0 = 2;

%% Individuo valori potenziale e corrente
% Controllo che il punto non sia nella sorgente
if (x0==a || x0==-a) && y0==0
    print('Il punto è singolare')
else
    if x0<0
        if y0<0
            psi0=-Q/2/pi*atan((2*y0*a)/(x0^2+y0^2-a^2)-pi)
        else
            psi0=-Q/2/pi*atan((2*y0*a)/(x0^2+y0^2-a^2)+pi)
        end
    else
        psi0=-Q/2/pi*atan((2*y0*a)/(x0^2+y0^2-a^2))
    end
    phi0 = Q/4/pi*log((y0^2+(a+x0)^2)/(y0^2+(a-x0)^2))
end

%% Grafica
% Intervallo funzione potenziale e di corrente in modo da riempire il
% grafico
scale = ceil(10*Q);
psi=linspace(-Q/2,Q/2,2*scale);
phi=linspace(-Q/2,Q/2,2*scale);

% Angoli per disegnare le circonferenze
th = linspace(-pi,pi,300);

% Coefficienti
n = exp(-4.*pi.*phi./Q);
m = 1./tan(-2.*pi.*psi./Q);
rphi = a*2.*sqrt(n)./(n-1);
rpsi = a.*sqrt(1+m.^2);

% Origine
plot(0,0,'k')
hold on

for i=1:length(psi)
    % Circonferenze potenziali
    p1 = plot(-a*(n(i)+1)/(n(i)-1)+rphi(i)*cos(th), rphi(i)*sin(th), 'r');
end

for i=1:length(psi)
    % Circonferenze corrente
    p2 = plot(rpsi(i)*cos(th),m(i)*a+rpsi(i)*sin(th),'b');
end


% Punto incognito
plot(x0,y0, 'ko')
plot([-3 3], [0 0], 'k-')
plot([0 0], [-3 3], 'k-')
n0 = exp(-4*pi*phi0/Q);
m0 = 1/tan(-2*pi*psi0/Q);
rphi0 = a*2*sqrt(n0)/(n0-1);
rpsi0 = a*sqrt(1+m0^2);

p3 = plot(-a*(n0+1)/(n0-1)+rphi0*cos(th), rphi0*sin(th), 'g');
p4 = plot(rpsi0*cos(th),m0*a+rpsi0*sin(th),'c');

legend([p1 p2 p3 p4], '\Phi', '\Psi', '\Phi_0', '\Psi_0')
axis([-3 3 -3 3])
axis square
hold off