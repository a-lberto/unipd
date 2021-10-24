%% Inizializzazione
clc;
clear all;
set(gcf,'renderer','Painters');

%% Definizione punti noti
A = [10 7];
B = [99 100];
C = [89 100];
D = [10  39.5];
G = [51.5 54.5];
F = [64 59];
rG = 47.5;
rF = 19.5;
x = [A(1) B(1) C(1) D(1) F(1) G(1)];
y = [A(2) B(2) C(2) D(2) F(2) G(2)];

%% Individuazione punti curva mozzo
th = linspace(-pi/2,0, 50);
xG = G(1)+rG*cos(th);
yG = G(2)+rG*sin(th);

%% Individuazione punto di tangenza
CF = sqrt((C(2)-F(2))^2+(C(1)-F(1))^2);
FT = rF;
thCF = atan((C(2)-F(2))/(C(1)-F(1)));
thCT = asin(FT/CF);
thFT = -(pi/2-thCF-thCT);

%% Individuazione punti curva corona
th = linspace(thFT,-pi/2, 50);
xF = F(1)+rF*cos(th);
yF = F(2)+rF*sin(th);

T(1) = F(1)+rF*cos(thFT);
T(2) = F(2)+rF*sin(thFT);

%% Punti profilo
xTurbo = [A(1) G(1) xG B(1) B(1) C(1) T(1) xF D(1) A(1)];
yTurbo = [A(2) A(2) yG G(2) B(2) C(2) T(2) yF D(2) A(2)];


for ddr=100


%% Discretizzazione
ddr                            % Numero di sezioni r
ddz = ddr;                           % Numero di sezioni z
dr = (B(2)-A(2))/ddr;                % Deltar
dz = (B(1)-A(1))/ddz;                % Deltaz         


%% Inizializzazione grafico
plot(xTurbo, yTurbo, 'k-')
hold on
axis([A(1)-5 B(1)+5 A(2)-5 B(2)+5])
axis square

% % Griglia
% for i=0:ddz
%     plot([A(1)+i*dz A(1)+i*dz],[A(2) B(2)], 'c-');
% end
% 
% for i=0:ddr
%     plot([A(1) B(1)],[A(2)+i*dr A(2)+i*dr], 'c-');
% end

%% Massimi e minimi

% Trovo minimi definendo la funzione F a tratti e calcolandone i valori per
% ogni dz
Fx = linspace(A(1),B(1), ddr+1);
Fy = Fx;
for i=1:length(Fx)
    if Fx(i)>=A(1) && Fx(i)<=G(1)
        Fy(i) = A(2);
    end
    if Fx(i)>G(1) && Fx(i)<=B(1)
        Fy(i) = -sqrt(rG^2-(Fx(i)-G(1))^2)+G(2);
    end
end
% E li indico nel grafico
% plot(Fx, Fy, 'ro');

% Trovo massimi definendo la funzione G a tratti e calcolandone i valori
% per ogni dz
Gx = linspace(A(1),B(1), ddz+1);
Gy = Gx;
for i=1:length(Gx)
    if Gx(i)>=D(1) && Gx(i)<=F(1)
        Gy(i) = D(2);
    end
    if Gx(i)>F(1) && Gx(i)<=T(1)
        Gy(i) = -sqrt(rF^2-(Gx(i)-F(1))^2)+F(2);
    end
    if Gx(i)>T(1) && Gx(i)<=C(1)
        Gy(i) = (C(2)-T(2))/(C(1)-T(1))*(Gx(i)-T(1))+T(2);
    end
    if Gx(i)>C(1)
        Gy(i) = C(2);
    end
end
% E li indico nel grafico
% plot(Gx, Gy, 'bo');

% Trovo minimi definendo la funzione H a tratti e calcolandone i valori
% per ogni dr
Hy = linspace(A(2), B(2), ddz+1);
Hx = Hy;
for i=1:length(Hx)
    if Hy(i)>=A(2) && Hy(i)<=D(2)
        Hx(i) = A(1);
    end
    if Hy(i)>D(2) && Hy(i)<T(2)
        Hx(i) = sqrt(rF^2-(Hy(i)-F(2))^2)+F(1);
    end
    if Hy(i)>T(2)
        Hx(i) = (C(1)-T(1))/(C(2)-T(2))*(Hy(i)-T(2))+T(1);
    end
end
% plot(Hx, Hy, 'bs');

% Trovo minimi definendo la funzione J a tratti e calcolandone i valori
% per ogni dz
Jy = linspace(A(2), B(2), ddr+1);
Jx = Jy;
for i=1:length(Jx)
    if Jy(i)>=A(2) && Jy(i)<=G(2)
        Jx(i) = sqrt(rG^2-(Jy(i)-G(2))^2)+G(1);
    end
    if Hy(i)>G(2)
        Jx(i) = B(1);
    end
end
% plot(Jx, Jy, 'bo');




%% Individuazione punti
% Come criterio per la scelta dei punti devono essere interni al campo
% fluido oppure essere adiacenti a destra, sinistra, sopra o sotto un punto
% interno al fluido; oppure al contorno

m = 1;
n = 1;
min = B(2)*ones(1,ddz+1);
for z=A(1):dz:B(1)
    for r=A(2):dr:B(2)
        interno =  r>=Fy(m)-dr && r<=Gy(m)+dr;        % Se il punto e' interno al campo fluido
        inferiore = m>1     && r<Fy(m) && r>Fy(m-1);  % Se il punto e' appena sotto al campo fluido
        superiore = m<ddz+1 && r>Gy(m) && r<=Gy(m+1); % Se il punto e' appena sopra al campo fluido
        
        mozzo = r<=Fy(m) || z > B(1)-dz/2;              % Se il punto e' del mozzo
        corona = r>=Gy(m) && z < C(1);                % Se il punto e' della corona
        
        if  (interno||inferiore||superiore)
            n0(n) = n; % Numero nodo  
            z0(n) = z; % Coordinata x
            r0(n) = r; % Coordinata y
            n1(n) = 1;
            n2(n) = 1;
            n3(n) = 1;
            n4(n) = 1;
            lambda1(n) = 1;
            lambda2(n) = 1;
            lambda3(n) = 1;
            lambda4(n) = 1;
            
            % Coordinate Z e R generalizzate
            Z(n) = round((z0(n)-A(1))/dz+1);
            R(n) = round((r0(n)-A(2))/dr+1);
            
%             plot(z0(n), r0(n), 'kx');
%             text(z,r,num2str(n));
%             if interno
%                plot(z,r,'rx') 
%             end
%             
%             if mozzo
%                 plot(z,r,'bx')
%             end
%             
%             if corona
%                 plot(z,r,'gx')
%             end
            
            % Assegno i valori di psi ai punti
            if mozzo
                psi0(n) = 0;
                psi0old(n) = 0;
%                 plot(z0(n), r0(n), 'c+')                
            else if corona
                psi0(n) = 100;
                psi0old(n) = 100;
%                 plot(z0(n), r0(n), 'c+')
            else if interno
                psi0(n) = 50;
                psi0old(n) = -1;
            end
            end
            end
            
            
            % Il vettore min(m) tiene il valore minimo di y per ogni
            % colonna m
            if r0(n) < min(m)
                min(m) = r0(n);
            end
            n=n+1; % Incrementa il numero etichetta del punto
       end
    end
    
    % Il vettore max(m) tiene il valori massimo di y per ogni colonna m
    max(m) = r0(n-1);
    m=m+1; % Incrementa il numero di colonna di riferimento
end
punti = n-1;



%% Individuazione punti adiacenti ad ogni punto del campo fluido
% Per i punti al contorno, come adiacenti dalla parte verso l'esterno si
% considerano il punto stesso, in modo che la differenza divisa verso
% quella direzione sia nulla

toll=10e-3;
m = 1;
for n=1:punti

    % Assegno punti inferiori
    n4(1) = n0(1);    
    if n>1 && z0(n) > z0(n-1)
        n4(n) = n0(n);
        
        % Dato che questa condizione mi dice la colonna successiva
        % incremento anche m
        m = m+1;
    end
    
    % Assegno punti superiori: se va alla colonna successiva vuol dire che
    % e' il punto piu' alto
    n2(punti) = n0(punti); 
    if n<punti && z0(n) < z0(n+1)
        n2(n) = n0(n);
    end
    
    % Assegno punti a sinistra: controllo a partire dalla posizione n-1 il
    % punto precedente finche' non ne trovo uno che ha lo stesso r e z
    % minore di dz. Le maggioranze sono fatte perche' non sono uguali di
    % cifre dell'ordine di 10e-14
    for o=n-1:-1:1
        if abs(r0(o)-r0(n))<toll && abs(z0(o)-(z0(n)-dz))<toll
            n3(n)= n0(o);
            break
        else
            n3(n) = n0(n);
        end
    end
    
    % Assegno punti a destra: controllo a partire dalla posizione n+1 il
    % punto successivo finche' non ne trovo uno che ha lo stesso r e z
    % maggiore di dz
    for o=n+1:1:punti
        if abs(r0(o)-r0(n))<toll && abs(z0(o)-(z0(n)+dz))<toll
            n1(n)= n0(o);
            break
        else
            n1(n) = n0(n);
        end
    end
    
    % Assegno punti sopra: sfrutto massimi per colonna calcolati prima
    if r0(n) < max(m)
        if n2(n) == 1
            n2(n) = n+1;
        end
    end
    
    % Assegno punti sotto: sfrutto minimi per colonna calcolati prima
    if r0(n) > min(m)
        if n4(n) == 1
            n4(n) = n-1;
        end
    end
end



%% Applico il metodo delle differenze divise senza lambda
% % Per ogni punto sfrutto l'equazione delle differenze divise
% scarto = 1;
% toll=10e-10;
% N=0;
% while scarto > toll && N < 3000
%     scarto = 0;
%     for n=1:punti
%         psi0old(n)=psi0(n);
%         psi24=(psi0(n2(n))+psi0(n4(n)))/(dr^2);
%         psi13=(psi0(n1(n))+psi0(n3(n)))/(dz^2);
%         psir=1/r0(n)*(psi0(n2(n))-psi0(n4(n)))/(2*dr);
%         psi0(n)=(psi24+psi13-psir)/(2*(1/(dr^2)+1/(dz^2)));
%         
%         % Condizioni di contorno vanno mantenute
%         if psi0old(n)==100 || psi0old(n)==0
%             psi0(n)=psi0old(n);
%         end
%         
%         % Mantengo lo scarto massimo in tutto il campo fluido
%         if abs(psi0(n)-psi0old(n))>scarto
%             scarto = abs(psi0(n)-psi0old(n));
%         end
%     end
%     
%     N = N + 1;
% end

%% Calcolo i lambda di tutti i punti
for n=1:punti
    % Lambda 1
    if z0(n1(n))>Jx(R(n)) && z0(n) < Jx(R(n))
        lambda1(n) = (Jx(R(n))-z0(n))/dz;
%         plot(z0(n)+lambda1(n)*dz,r0(n), 'ro')
%         plot(z0(n),r0(n), 'r+')
    end
    
    % Lambda 2
    if r0(n2(n))>Gy(Z(n)) && r0(n)<Gy(Z(n))
        lambda2(n) = (Gy(Z(n))-r0(n))/dr;
%         plot(z0(n),r0(n)+lambda2(n)*dr, 'r^')
%         plot(z0(n),r0(n), 'r+')
    end
    
    % Lambda 3
    if z0(n3(n))<Hx(R(n)) && z0(n)>Hx(R(n))
        lambda3(n) = (z0(n)-Hx(R(n)))/dz;
%         plot(z0(n)-lambda3(n)*dz,r0(n), 'rs')
%         plot(z0(n),r0(n), 'r+')
    end
    
    % Lambda 4
    if r0(n4(n))<Fy(Z(n)) && r0(n)>Fy(Z(n))
        lambda4(n) = (r0(n)-Fy(Z(n)))/dr;
%         plot(z0(n),r0(n)-lambda4(n)*dr, 'r^')
%         plot(z0(n),r0(n), 'r+')
    end
end


%% Applico le differenze divise con i lambda
scarto = 1;
toll=10e-5;
N=0;
while scarto > toll && N < 3000
    scarto = 0;
    for n=1:punti
        psi0old(n)=psi0(n);
        p1 = psi0(n1(n));
        p2 = psi0(n2(n));
        p3 = psi0(n3(n));
        p4 = psi0(n4(n));
        l1 = lambda1(n);
        l2 = lambda2(n);
        l3 = lambda3(n);
        l4 = lambda4(n);
        
        p24 = (p2/l2+p4/l4)/(l2+l4)/(dr^2);
        p13 = (p1/l1+p3/l3)/(l1+l3)/(dz^2);
        pr = (p2-p4)/(2*r0(n)*dr*(l2+l4));
        lambdas = (l1*l3*dz^2+l2*l4*dr^2)/(l1*l2*l3*l4*dr^2*dz^2);
        
        psi0(n) = (p24+p13-pr)/lambdas;
        
        
        % Condizioni di contorno vanno mantenute
        if psi0old(n)==100 || psi0old(n)==0
            psi0(n)=psi0old(n);
        end
        
        % Mantengo lo scarto massimo in tutto il campo fluido
        if abs(psi0(n)-psi0old(n))>scarto
            scarto = abs(psi0(n)-psi0old(n));
        end
    end
    
    N = N + 1;
end

% %% Valori di psi attuali
% for n=1:punti
%     text(z0(n),r0(n), sprintf('%2.0f',psi0(n)));
% end

%% Cerco i punti per cui un valore di psi e' uguale o compreso
for psix=[5 50 95]
    clear psiXz psiXr
    m=1;
    p2=0;
    p0=0;
    i=1;
    for n=1:punti
        if n>1 && z0(n) > z0(n-1)
            m=m+1;
        end
        if z0(n)<G(1)+rF*sqrt(2) && psix>=psi0(n) && psix<=psi0(n2(n))
            p2 = psi0(n2(n));
            p0 = psi0(n);
            rx = (psix-p0)/(p2-p0)*(r0(n2(n))-r0(n))*lambda2(n)*lambda4(n2(n))+r0(n)+dr*(1-lambda4(n2(n)));
            psiXr(i) = rx;
            psiXz(i) = z0(n);
            i = i+1;
        end
        if z0(n)>=G(1)+rF*sqrt(2) && psix>=psi0(n) && psix<=psi0(n3(n))
            psi3 = psi0(n3(n));
            p0 = psi0(n);
            zx = (psix-p0)/(psi3-p0)*(z0(n3(n))-z0(n))*lambda3(n)*lambda1(n3(n))+z0(n)-dz*(1-lambda1(n3(n)));
            psiXz(i) = zx;
            psiXr(i) = r0(n);
            i = i+1;
        end
    end
%     plot(psiXz, psiXr, 'Color', [1-0.01*psix 0 0.01*psix]);
    plot(psiXz, psiXr, 'b');
    
end

%% Test per il punto n: mostra posizione punti adiacenti
% % Mi aspetto una disposizione a croce
% n = floor(rand*(punti));
% n = 112;
% plot(z0(n), r0(n), 'k+')
% plot(z0(n1(n))+dz/4, r0(n1(n)),        'ko')
% plot(z0(n2(n)),      r0(n2(n))+dr/4,   'k^')
% plot(z0(n3(n))-dz/4, r0(n3(n)),        'ks')
% plot(z0(n4(n)),      r0(n4(n))-dr/4,   'kx')
% 
% % Didascalia in basso a destra del punto ispezione
% text(B(1), A(2), num2str(n),'HorizontalAlignment', 'center');
% text(B(1)-dz/4, A(2)-dr/2, num2str(Z(n)),'HorizontalAlignment', 'center');
% text(B(1)+dz/4, A(2)-dr/2, num2str(R(n)),'HorizontalAlignment', 'center');

%% Chiusura plot
hold off;
pause(0.1)
end