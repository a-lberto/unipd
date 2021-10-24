clc; clear all;
%% Implementazione metodo dei pannelli

% Cifre NACA 4 cifre
% Dividere prima, e seconda da terza e quarta
NACA = [0 0 12];
c=1;

vinf=1;
alpha=5/180*pi;

m = NACA(1)/100;

if NACA(2)==0
    p=1e-16;
else
    p = NACA(2)/10;
end

s = NACA(3)/100;

% Funzione linea media
lineaMedia = @(t) (m/p^2 * (2*p.*t-t.^2))           .*((t>=0)&(t<=p))+...
                  (m/(1-p)^2*(1-2*p+2*p.*t-t.^2))   .*((t>=p)&(t<=1));

% Funzione spessori
A = [0.2969 -0.1260 -0.3537 0.2843 -0.1015];
spessore = @(t) s/0.20/c.*...
    ( A(1).*sqrt(t) ...
    + A(2).*t ...
    + A(3).*t.^2 ...
    + A(4).*t.^3 ...
    + A(5).*t.^4 );

% Coordinate ascisse
x=linspace(0,c,300);
t=x/c;

% Numero pannelli
for nPannelli = 60
nPannelli
% Grafico profilo ala
% Linea media
figure(1)
plot(t,lineaMedia(t), 'b:')
hold on

% Linee spessori
plot(t,lineaMedia(t)+spessore(t), 'k')
plot(t,lineaMedia(t)-spessore(t), 'k')


th = linspace(0,-2*pi, nPannelli+1)';
% plot(0.5+0.5*cos(th),0.5*sin(th), 'ro');

% Punti per pannelli
x = 0.5+0.5*cos(th);
y = 0.5*sin(th);
n=zeros(length(x), 1);
for i=1:length(th)
    if sin(th(i))>0
        y(i) = lineaMedia(x(i))+spessore(x(i));
    else
        y(i) = lineaMedia(x(i))-spessore(x(i));
    end
    n(i)=i;
end
plot(x,y, 'bs-')


l=zeros(nPannelli,1);
sin_th=l;
cos_th=l;
xs=l;
ys=l;
clear th; th=l;
% Varie proprieta' pannelli
for i=1:nPannelli
    % Lunghezza pannello
    l(i) = sqrt((x(i+1)-x(i))^2+(y(i+1)-y(i))^2);
    
    % Seno inclinazione pannello
    sin_th(i)=(y(i+1)-y(i))/l(i);
    
    % Coseno inclinazione pannello
    cos_th(i)=(x(i+1)-x(i))/l(i);
    
    % Ascissa punto di controllo
    xs(i) = (x(i)+x(i+1))/2;
    
    % Ordinata punto di controllo
    ys(i) = (y(i)+y(i+1))/2;
    
    % Theta
    th(i)=atan2((y(i+1)-y(i)),(x(i+1)-x(i)));
end
clear sin_th cos_th

figure(1)
plot(xs, ys, 'rx')
axis square
axis([0-0.1 1+0.1 -0.5-0.1 0.5+0.1])
title('Discretizzazione')

% Etichette nodi
 if length(n)<30
    r=0.05;
    for i=1:length(th)
            text(xs(i),ys(i), num2str(n(i)));
    end
end
% Velocita' dovuta a sorgenti e vortici
beta=zeros(nPannelli);
rr=zeros(nPannelli,nPannelli+1);
logrr=zeros(nPannelli);

for i=1:nPannelli
    for j=1:nPannelli
        xx(i,j)  =  (xs(i)-x(j))*cos(th(j))+(ys(i)-y(j))*sin(th(j));
        yy(i,j)  = -(xs(i)-x(j))*sin(th(j))+(ys(i)-y(j))*cos(th(j));
        
        % r(i,j+1)^2
        rr(i,j+1) = sqrt((xs(i)-x(j+1))^2 + (ys(i)-y(j+1))^2);
        
        % r(i,j)^2
        rr(i,j) = sqrt((xs(i)-x(j))^2 + (ys(i)-y(j))^2);
        
        logrr(i,j)=log(rr(i,j+1)/rr(i,j));
        
        if i==j
            beta(i,j) = pi;
        else
            beta(i,j) = atan2(yy(i,j),(xx(i,j)-l(j)))-atan2(yy(i,j),xx(i,j));
%             beta(i,j) = acos((rr(i,j)^2+rr(i,j+1)^2-l(j)^2)/(2*rr(i,j)*rr(i,j+1)));
        end
        
%         figure(2)
%         plot(x,y,'k:')
%         hold on
%         
%         % Segmenti
%         plot([x(i) x(i+1)],[y(i) y(i+1)], 'b')
%         plot(xs(i), ys(i), 'rx')
%         plot([x(j) x(j+1)], [y(j) y(j+1)], 'r')
%         if i==j
%             plot([x(i) x(i+1)],[y(i) y(i+1)], 'g')
%         end
%         plot([x(j) xs(i)], [y(j) ys(i)], 'r:')
%         plot([x(j+1) xs(i)], [y(j+1) ys(i)], 'b:')
%         axis([0 1 -0.5 0.5])
%         axis square
%         hold off
        
    end
end

%% Soluzione del sistema
A=zeros(nPannelli+1);
b=zeros(nPannelli+1,1);

% Compilazione vettore dati b
% bij
for i=1:nPannelli
   b(i) = vinf*(sin(th(i)-alpha)); %ok
end

% Condizione di Kutta
b(end) = - vinf*(cos(th(1)-alpha)) - vinf*(cos(th(end)-alpha));

% Matrice A
for i=1:nPannelli
    for j=1:nPannelli
       % A(i,j)
       A(i,j) = 1/2/pi*(sin(th(i)-th(j))*log(rr(i,j+1)/rr(i,j)) ...
                       +cos(th(i)-th(j))*beta(i,j));       
    end
    
    % A(1, N+1)
    for j=1:nPannelli
        A(i,nPannelli+1) = A(i,nPannelli+1) + ...
           1/2/pi*(cos(th(i)-th(j))*log(rr(i,j+1)/rr(i,j)) ...
                  -sin(th(i)-th(j)*beta(i,j)));
    end
end

% A(N+1, j)
for j=1:nPannelli
    for k=[1, nPannelli]
        A(nPannelli+1,j)=A(nPannelli+1,j)+ ...
            1/2/pi*(sin(th(k)-th(j))*beta(k,j)...
                   -cos(th(k)-th(j))*log(rr(k,j+1)/rr(k,j)));
    end
end

% A(N+1, N+1)
for j=1:nPannelli
    for k=[1, nPannelli]
        A(nPannelli+1,nPannelli+1)=A(nPannelli+1,nPannelli+1)+ ...
            1/2/pi*(sin(th(k)-th(j))*log(rr(k,j+1)/rr(k,j)) ...
                  + cos(th(k)-th(j)*beta(k,j)));
    end
end
xSol=A\b;
gamma = xSol(end);
q = xSol(1:end-1);

% gamma=0.1;
% q=zeros(1,nPannelli)';
% gamma=0;
xx=ones(length(xs),1);

phi=@(csi,eta) vinf*(csi*cos(alpha)+eta*(sin(alpha)))...              % Indisturbato
              +1/2/pi*(q'*1/2*log((xs-csi*xx).^2+(ys-eta*xx).^2))...  % Sorgente
              -1/2/pi*(gamma*xx'*atan2((ys-eta*xx),(xs-csi*xx)));     % Vortice

% psi=@(csi,eta) vinf*(eta*cos(alpha)-csi*(sin(alpha))) ...
%               +1/2/pi*(q'*atan2((eta*xx-ys),(csi*xx-xs)))...
%               +1/2/pi*(gamma*xx'*1/2*log((xs-csi*xx).^2+(ys-eta*xx).^2));

% Griglia per calcolo contour
x=linspace(-0.5,1.5,200);
y=linspace(-1,1,200);
[X,Y] = meshgrid(x,y);

% % Plot del grafico
% for i=1:length(x)
%     for j=1:length(x)
%         PHI(i,j)=phi(X(i,j), Y(i,j));
%     end
% end
% contour(X,Y,PHI, linspace(-10,20, 100))
% axis square

hold off

% %% Calcolo delle velocita'
% u=zeros(nPannelli,1);
% v=u;
% for i=1:nPannelli
%     u(i)=vinf*cos(alpha);
%     v(i)=vinf*sin(alpha);
%     for j=1:nPannelli
%         u(i)=u(i)+q(j)*us(i,j)+gamma*uv(i,j);
%         v(i)=v(i)+q(j)*vs(i,j)+gamma*vv(i,j);
%     end
% end

% vt=zeros(nPannelli,1);
% for i=1:nPannelli
% %     vt(i)=vinf*(sin_th(i)*sin(alpha)+cos_th(i)*cos(alpha));
%     for j=1:nPannelli
%         vt(i)=vt(i)+...
%             q(j)*( ...
%               - sinDif(j,i)*vs(i,j) ...
%               - cosDif(j,i)*us(i,j) ...
%             ) ...
%           + gamma*( ... 
%                - sinDif(j,i)*vv(i,j) ...
%                - cosDif(j,i)*uv(i,j) ...
%             );
%     end
% end

vt=zeros(nPannelli,1);
for i=1:nPannelli
    vt(i)=vinf*cos(th(i)-alpha);
    for j=1:nPannelli
        vt(i)=vt(i)+...
            q(j)/2/pi*( ...
                sin(th(i)-th(j))*beta(i,j) ...
              - cos(th(i)-th(j))*log(rr(i,j+1)/rr(i,j)) ...
            ) ...
          + gamma/2/pi*( ... 
                sin(th(i)-th(j))*log(rr(i,j+1)/rr(i,j)) ...
               + cos(th(i)-th(j))*beta(i,j) ...
            );
    end
end

% vt=zeros(nPannelli,1);
% for i=1:nPannelli
% %     vt(i)=vinf*(sin_th(i)*sin(alpha)+cos_th(i)*cos(alpha));
%     for j=1:nPannelli
%         vt(i)=
%     end
% end

Cp=zeros(nPannelli, 1);
for i=1:nPannelli
    Cp(i)=1-(vt(i)^2)/vinf^2;
end


figure(2)
plot(xs,ys, 'k');
hold on
p1 = plot(xs(1:floor(length(xs)/2)), Cp(1:floor(length(Cp)/2)), 'bx-');
p2 = plot(xs(ceil(length(xs)/2):end), Cp(ceil(length(xs)/2):end), 'rx-');
title('Cp')
legend([p1 p2], 'Intradosso', 'Estradosso')
% axis([0-0.1 1+0.1 -1.5 1.5])
hold off
pause(1)

end



% %% Controllo se le coordinate xxij, yyij e rij sono corrette
% subplot(1,2,2)
% th = linspace(0,2*pi,100);
% for i=1:nPannelli
%     for j=1:nPannelli
%         % Profilo
%         plot(x,y,'k:')
%         hold on
%         
%         % Segmenti
%         plot([x(i) x(i+1)],[y(i) y(i+1)], 'b')
%         plot(xs(i), ys(i), 'rx')
%         plot([x(j) x(j+1)], [y(j) y(j+1)], 'r')
%         if i==j
%             plot([x(i) x(i+1)],[y(i) y(i+1)], 'g')
%         end
%         
%         % Raggi distanza
%         R=[ cos_th(i) sin_th(i),...
%            -sin_th(i) cos_th(i)];
%         xsi = xx(i,j)*cos_th(i)-yy(i,j)*sin_th(i)+x(j);
%         ysi = xx(i,j)*sin_th(i)+yy(i,j)*cos_th(i)+y(j);
%         
%         plot(x(j),y(j), 'kx');
%         plot(x(j)+rr(i,j)*cos(th),y(j)+rr(i,j)*sin(th), 'r')
%         
%         plot(x(j+1),y(j+1), 'kx');
%         plot(x(j+1)+rr(i,j+1)*cos(th),y(j+1)+rr(i,j+1)*sin(th), 'b')
%         
%         plot([x(j) xsi], [y(j) ysi], 'r:')
%         plot([x(j+1) xsi], [y(j+1) ysi], 'b:')
%         axis([0 1 -0.5 0.5])
%         axis square
%         hold off
%         pause(2/nPannelli^2)
%     end
% end