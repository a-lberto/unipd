%% Implementazione flusso attorno a un cilindro con circuitazione
clc; clear all;
set(gcf,'renderer','Painters');

v = 3;
k = 3;
R=k/2/pi/v;

% Circuitazione nulla
L = 0;

% Circuitazione moderata
L = 4*pi*R*v/2;

% Circuitazione cuspide
L = 4*pi*R*v;

% Punti di ristagno esterni al moto
L = 5*pi*R*v;

th=linspace(pi/2+0.0001,-3/2*pi-0.0001,500);

r=linspace(0,10*R,500);

[r,th]=meshgrid(r,th);

phi = v.*r.*cos(th).*(1+R^2./r.^2)-L/2/pi.*th;
psi = v.*r.*sin(th).*(1-R^2./r.^2)+L/2/pi*log(r);

X=r.*cos(th);
Y=r.*sin(th);

levels = -k*v:1/4:k*v;

contour(X,Y, phi, levels, 'r')
hold on
contour(X,Y, psi, levels, 'b')

% Punti di ristagno
if L/4/pi/R/v<=1
    % Se sono sul cilindro
    y0=-L/4/pi/v;
    x0=sqrt(R^2-y0^2);
    plot([-x0 x0],[y0 y0], 'ko')
    psi0=L/2/pi*log(R);
    contour(X, Y, psi, [psi0 psi0], 'k')
else
    % Se e' fuori dal cilindro
    x0=0;
    y0=-L/4/pi/v*(1+sqrt(1-(4*pi*v*R/L)^2));
    th0=-pi/2;
    r0=abs(y0);
    psi0=v*r0*sin(th0)*(1-R^2/r0^2)+L/2/pi*log(r0);
    plot(x0, y0, 'ko')
    contour(X, Y, psi, [psi0 psi0], 'k')
end

th=linspace(pi/2+0.0001,-3/2*pi-0.0001,500);
plot(R*cos(th), R*sin(th), 'g')

axis square
axis([-0.5 0.5 -0.5 0.5])
% pbaspect([2 1 1])
hold off