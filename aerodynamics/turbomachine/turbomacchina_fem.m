%% Pulizia memoria
clear all
clc

%% Definizione punti noti
A = [10 7];
B = [99 100];
C = [89 100];
D = [10  39.5];
G = [51.5 54.5];
F = [64 59];
rG = 47.5;
rF = 19.5;

%% Punti curve mozzo e corona
% Individuazione punti curva mozzo
th = linspace(-pi/2,0, 50);
xG = G(1)+rG*cos(th);
yG = G(2)+rG*sin(th);

% Individuazione punto di tangenza
CF = sqrt((C(2)-F(2))^2+(C(1)-F(1))^2);
FT = rF;
thCF = atan((C(2)-F(2))/(C(1)-F(1)));
thCT = asin(FT/CF);
thFT = -(pi/2-thCF-thCT);

% Individuazione punti curva corona
th = linspace(thFT,-pi/2, 50);
xF = F(1)+rF*cos(th);
yF = F(2)+rF*sin(th);

T(1) = F(1)+rF*cos(thFT);
T(2) = F(2)+rF*sin(thFT);

%% Punti profilo
xTurbo = [A(1) G(1) xG B(1) B(1) C(1) T(1) xF F(1) D(1) A(1)];
yTurbo = [A(2) A(2) yG G(2) B(2) C(2) T(2) yF D(2) D(2) A(2)];

%% Riassegno variabili
xA = A(1);
yA = A(2);
xB = B(1);
yB = B(2);
xC = C(1);
yC = C(2);
xF = F(1);
yF = F(2);
xD = D(1);
yD = D(2);
xF1 = xF(1);
xF2 = T(1);
yF2 = T(2);
xG = G(1);
yG = G(2);
alpha = pi/2+thFT;

tabellafile = fopen('iter_fem.txt','w');
fprintf(tabellafile,'%6s \t %12s \t %10s \t %10s \t %10s \n\n','ndEta','ndCsi', 'nElementi', 'nNodi', 'Tempo calcolo');

for ndEta=8
tic
%% Discretizzazione
% Suddivisione trasversale
ndEta
% Suddivisione longitudinale
ndCsi = 5*ndEta

%% Curva corona
% Lunghezze tratti curva corona
s1c = xF1-xD;
s2c = rF*alpha;
s3c = sqrt((xC-xF2)^2+(yC-yF2)^2);

% Lunghezza curva corona
lcc = s1c+s2c+s3c;

% Parametrizzazioni delle ascisse e ordinate della curva della corona
pXC = @(s) ((s>=0)&(s<s1c)).*       (xD+s) ...
    +      ((s>=s1c)&(s<s1c+s2c)).* (xF+rF*sin((s-s1c)/rF)) ...
    +      (s>=s1c+s2c).*           (xF2+((xC-xF2)/s3c).*(s-(s1c+s2c)));
pYC = @(s) ((s>=0)&(s<s1c)).*       (yD) ...
    +      ((s>=s1c)&(s<s1c+s2c)).* (yF-rF*cos((s-s1c)/rF)) ...
    +      (s>=s1c+s2c).*           (yF2+((yC-yF2)/s3c).*(s-(s1c+s2c)));

% Punti della curva
sCC = linspace(0, lcc, ndCsi+1);
xCC = pXC(sCC);
yCC = pYC(sCC);

%% Curva mozzo
% Lunghezze tratti curva mozzo
s1m = xG-xA;
s2m = rG*pi/2;
s3m = yB-yG;

% Lunghezza curva mozzo
lcm = s1m+s2m+s3m;

% Parametrizzazioni delle ascisse e ordinate della curva del mozzo
pXC = @(s) ((s>=0)&(s<s1m)).*        (xA+s) ...
    +      ((s>=s1m)&(s<s1m+s2m)).*  (xG+rG*sin((s-s1m)/rG)) ... 
    +      (s>=s1m+s2m).*            (xB);
pYC = @(s) ((s>=0)&(s<s1m)).*        (yA) ...
    +      ((s>=s1m)&(s<s1m+s2m)).*  (yG-rG*cos((s-s1m)/rG)) ...
    +      (s>=s1m+s2m).*            (yG+(s-(s1m+s2m)));

% Punti della curva
sCM = linspace(0, lcm, ndCsi+1);
xCM = pXC(sCM);
yCM = pYC(sCM);

%% Inizializzazione grafico
plot(xTurbo, yTurbo, 'k-')
hold on
axis([min(xTurbo)-5 max(xTurbo)+5 min(yTurbo)-5 max(yTurbo)+5])
axis square


%% Suddivisione trasversale
z=zeros((ndCsi+1)*(ndEta+1),1);
r=zeros((ndCsi+1)*(ndEta+1),1);

i=1;
for n=1:ndCsi+1
    for m=0:ndEta
        z(i)=(xCC(n)-xCM(n))/ndEta*m+xCM(n);
        r(i)=(yCC(n)-yCM(n))/ndEta*m+yCM(n);
%         plot(z(i),r(i), 'kx')
        i=i+1;
    end
end

%% Costruzione matrice nodi, ascisse e ordinate, psi
N = zeros(ndEta+1, ndCsi+1);
zN = zeros(ndEta+1, ndCsi+1);
rN = zeros(ndEta+1, ndCsi+1);
psiN = zeros(ndEta+1, ndCsi+1);

i=1;
for n=1:ndCsi+1
   for m=1:ndEta+1
       N(m,n)=i;
       zN(m,n)=z(i);
       rN(m,n)=r(i);
       
       if m==1
           psiN(m,n)=0;
       else
           if m==ndEta+1
               psiN(m,n)=100;
           else
               psiN(m,n) = NaN;
           end
       end
       i=i+1;
   end
end

%% Costruzione matrice coefficienti
a = zeros(N(ndEta+1,ndCsi+1),N(ndEta+1,ndCsi+1));

for n=1:ndCsi+1
    for m=1:ndEta+1
        l = N(m,n);
        if n>1
            if m>1 % ------------------------In basso a sinistra
                
                % Peso N = 3
                k = 3;
                m_i = [m-1 m-1 m m];
                n_i = [n-1 n n n-1];
                z_i = diag(zN(m_i,n_i))';
                r_i = diag(rN(m_i,n_i))';
                
                h = 3;
                a(l,l) = a(l,l) + aIntGauss(z_i, r_i, h, k);
                
                h = 2;
                a(l,N(m-1,n)) = a(l, N(m-1,n)) + aIntGauss(z_i, r_i, h, k);
                
                h = 4;
                a(l,N(m,n-1)) = a(l,N(m,n-1)) + aIntGauss(z_i, r_i, h, k);
                
                h = 1;
                a(l,N(m-1,n-1)) = a(l,N(m-1,n-1)) + aIntGauss(z_i, r_i, h, k);
                
                clear z_i r_i;
            end
            if m<ndEta+1 % ------------------In alto a sinistra
                
                % Peso N = 2
                k = 2;
                m_i = [m m m+1 m+1];
                n_i = [n-1 n n n-1];
                z_i = diag(zN(m_i,n_i))';
                r_i = diag(rN(m_i,n_i))';
                
                h = 2;
                a(l,l) = a(l,l) + aIntGauss(z_i, r_i, h, k);
                
                h = 3;
                a(l,N(m+1,n)) = a(l, N(m+1,n)) + aIntGauss(z_i, r_i, h, k);
                
                h = 1;
                a(l,N(m,n-1)) = a(l,N(m,n-1)) + aIntGauss(z_i, r_i, h, k);
                
                h = 4;
                a(l,N(m+1,n-1)) = a(l, N(m+1,n-1)) + aIntGauss(z_i, r_i, h, k);
                
                clear z_i r_i;
            end
        end
        if n<ndCsi+1
            if m>1 % ------------------------In basso a destra
                
                % Peso N = 4
                k = 4;
                m_i = [m-1 m-1 m m];
                n_i = [n n+1 n+1 n];
                z_i = diag(zN(m_i,n_i))';
                r_i = diag(rN(m_i,n_i))';
                
                h = 4;
                a(l,l) = a(l,l) + aIntGauss(z_i, r_i, h, k);
                
                h = 1;
                a(l,N(m-1,n)) = a(l, N(m-1,n)) + aIntGauss(z_i, r_i, h, k);
                
                h = 3;
                a(l,N(m,n+1)) = a(l,N(m,n+1)) + aIntGauss(z_i, r_i, h, k);
                
                h = 2;
                a(l,N(m-1,n+1)) = a(l,N(m-1,n+1)) + aIntGauss(z_i, r_i, h, k);
                
                clear z_i r_i;
            end
            if m<ndEta+1 % ------------------In alto a destra
                
                % Peso N = 1        
                k = 1;
                m_i = [m m m+1 m+1]; 
                n_i = [n n+1 n+1 n];
                z_i = diag(zN(m_i,n_i))';
                r_i = diag(rN(m_i,n_i))';
                
                h = 1;
                a(l,l) = a(l,l) + aIntGauss(z_i, r_i, h, k);
                
                h = 4;
                a(l,N(m+1,n)) = a(l, N(m+1,n)) + aIntGauss(z_i, r_i, h, k);
                
                h = 2;
                a(l,N(m,n+1)) = a(l,N(m,n+1)) + aIntGauss(z_i, r_i, h, k);
                
                h = 3;
                a(l,N(m+1,n+1)) = a(l,N(m+1,n+1)) + aIntGauss(z_i, r_i, h, k); 
                
                clear z_i r_i;
            end
        end
    end
end

%% Costruzione vettore valori psi
psi = zeros((ndCsi+1)*(ndEta+1),1);
i = 1;
for n=1:ndCsi+1
    for m=1:ndEta+1
        psi(i)=psiN(m,n);
        i=i+1;
    end
end

%% Vettore soluzione
b = zeros(length(psi),1);

% Inserimento colonne note nel vettore soluzione
for i=1:length(psi)
    if (psi(i)==0) || (psi(i) ==100)
        b=b-a(:,i)*psi(i);
    end
end

%% Risoluzione sistema Ax=b
% Rimozione righe e colonne psi noti
aold = a;
bold = b;
j=1;
i=1;
while true
    if i==length(psi)+1
        break;
    else
        if psi(i) == 0 || psi(i) == 100
            if j==1
                anew=aold(2:end,2:end);
                bnew=bold(2:end);
            else
                if j==length(aold(:,1))
                    anew=aold(1:(end-1),1:(end-1));
                    bnew=bold(1:(end-1));
                else 
                    anew=aold([1:(j-1),(j+1):end],[1:(j-1),(j+1):end]);
                    bnew=bold([1:(j-1),(j+1):end]);
                end
            end
            i=i+1;
        else
            j=j+1;
            i=i+1;
        end
    end
    aold=anew;
    bold=bnew;
end
asol=aold;
bsol=bold;

tempoCalcolo = toc;

% Calcolo psi incognite
opts.SYM=true;
x=linsolve(asol,bsol, opts);

%% Completamento matrice psiN
i=1;
for n=1:ndCsi+1
    for m=1:ndEta+1
        if isnan(psiN(m,n))
            psiN(m,n) = x(i);
            i=i+1;
        end
    end
end

%% Curve di livello della psi


% Coordinate punti elemento parente
csi_i = [-1 1 1 -1];
eta_i = [-1 -1 1 1];

% Funzione di forma bidimensionale lineare
Ni = @(csi, eta, i) 1/4*(1+csi_i(i)*csi)*(1+eta_i(i)*eta);

% Ricostruzione funzione interpolata dati i 4 nodi
xCsiEta = @(csi, eta, x_i) Ni(csi, eta, 1)*x_i(1) + ...
                           Ni(csi, eta, 2)*x_i(2) + ...
                           Ni(csi, eta, 3)*x_i(3) + ...
                           Ni(csi, eta, 4)*x_i(4);

% Scelta valori curve livello
for psiX = [5:5:95]
    
    % Individuazione punti al di sotto della curva
    clear minPsi
    minPsi = zeros(ndCsi+1,1);
    for n=1:ndCsi+1
        for m=1:ndEta+1
            if psiN(m,n)>psiX
                minPsi(n) = m-1;
                break
            end
        end
    end
    
    % Individuazione punti curva con ausilio elemento parente
    clear rPsiX zPsiX
    zPsiX = zeros(1,5);
    rPsiX = zPsiX;
    
    for n=1:ndCsi
        m = minPsi(n);

        %Elemento in alto a destra
        m_i = [m m m+1 m+1]; 
        n_i = [n n+1 n+1 n];
        z_i = diag(zN(m_i,n_i))';
        r_i = diag(rN(m_i,n_i))';

        % Semplificazione riferimenti a psi
        psi1 = psiN(m_i(1), n_i(1));
        psi2 = psiN(m_i(2), n_i(2));
        psi3 = psiN(m_i(3), n_i(3));
        psi4 = psiN(m_i(4), n_i(4));

        % Suddivisione ascisse in elemento parente
        csi = linspace(-1,1,20);
        eta = zeros(1,20);
        for i=1:20
            % Ordinata noto valore ascissa e funzione psi
            % Valida per funzioni di forma lineari
            eta(i) = ( ...
                      4*psiX...
                    - (1-csi(i))*(psi1+psi4)...
                    - (1+csi(i))*(psi2+psi3)...
                  ) / (...
                      (1+csi(i))*(psi3-psi2)...
                    + (1-csi(i))*(psi4-psi1)...
                  );

            zPsiX(i) = xCsiEta(csi(i), eta(i), z_i);
            rPsiX(i) = xCsiEta(csi(i), eta(i), r_i);
        end
        plot(zPsiX, rPsiX, 'r')
    end
end
hold off



A = [ndEta, ndCsi, ndEta*ndCsi, (ndEta+1)*(ndCsi+1), toc];

fprintf(tabellafile,'%i \t %i \t %i \t %i \t %1.4f \n\n',A);
fnameeps = sprintf('%i_%i_%1.3f.eps', ndEta, (ndEta+1)*(ndCsi+1),toc);
fnamejpg = sprintf('%i_%i_%1.3f.png', ndEta, (ndEta+1)*(ndCsi+1),toc);
saveas(gcf,fnameeps,'epsc');
saveas(gcf,fnamejpg,'png');

end
fclose(tabellafile);

clc
type iter_fem.txt