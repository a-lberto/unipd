clc; clear all

%% Implementazione corpo piano rankine
x=linspace(-3,3,1000);
r=3;
th=linspace(0,pi,300);

% Prendo angoli piu' vicini agli estremi dove diverge
th=pi/2+pi/2.*cos(th);
th=pi/2+pi/2.*cos(th);

Q=3;
v=1;

figure(1)
plot([0 pi],[0 0], 'k')
hold on

%% r in funzione di theta
psi =[-Q/2-1.5 -Q/2-1 -Q/2-0.5 -Q/2-0.01 -Q/2+0.01 -0.01 0.01 0.5 1 1.5];
for psii=psi
    c = [0 (-psii+max(psi))/(max(psi)-min(psi)) (-psii+max(psi))/(max(psi)-min(psi))];
    plot(th, (psii/v+Q/2/v*(1-th./pi))./sin(th), 'Color', c)
end
axis([0 pi -5 5])
hold off

phi=[-2*Q*v:1/2:2*Q*v];
psiUp = [0.0001:1/2:2*Q*v];
psiDown = [-Q/2-0.001:-1/2:-Q/2-2*Q*v];
psiIn = [-Q/2+0.0001:0.9998/2:-0.0001];

%% Piano di Rankine
% Evidenzio le tre regioni
figure(2)
plot(0,0)
hold on

% Funzione potenziale
q=x;
for i=1:length(phi)
    for j=1:length(x)
        % Controllo che l'interno della radice sia positivo o nullo
        q(j)=exp(4*pi/Q*(phi(i)-v*x(j)))-x(j)^2;
        if q(j)<0
            q(j)=0;
        end
    end
    plot(x,sqrt(q), 'r-')
    plot(x,-sqrt(q), 'r-')
end

% Corrente sopra
psi=psiUp;
for i=1:length(psi)
    y=psi(i)/v+Q/2/v*(1-th./pi);
    plot(y./tan(th),y, 'b')
end

% Corrente sotto
psi=psiDown;
for i=1:length(psi)
    y=psi(i)/v+Q/2/v*(1-th./pi);
    plot(y./tan(th),y, 'b')
end

% Corrente all'interno
psi=psiIn;
for i=1:length(psi)
    y=psi(i)/v+Q/2/v*(1-th./pi);
    plot(y./tan(th),y, 'b')
end

% Punto di ristagno
plot(-Q/2/pi/v, 0, 'go')
plot([0 0], [-3 3], 'k-')
axis([-3 3 -3 3])
axis square
hold off
