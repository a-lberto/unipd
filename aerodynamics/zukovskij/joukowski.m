clc; clear all;

th = linspace(0,2*pi,100);

thx0 = 175/180*pi;

m = 0.3;
delta = 135/180*pi;
b = 2;
R = sqrt((m*sin(delta))^2+(b-m*cos(delta))^2);

i=1;
 for thx=linspace(0+thx0,2*pi+thx0,181)
    subplot(1, 2, 1)
    plot([-3 3], [0 0], 'k:')
    hold on
    plot([0 0], [-3 3], 'k:')

    % Circonferenza cerchio generatore
    plot(m*cos(delta)+R*cos(th),m*sin(delta)+R*sin(th));

    % Centro C
    plot(m.*cos(delta),m.*sin(delta), 'ko');

    % O-C
    plot([0 m.*cos(delta)],[0 m.*sin(delta)]);

    % B1-C
    plot([b m.*cos(delta)],[0 m.*sin(delta)]);
    
    % C2
    xC2 = -b*m*cos(delta)/(b-2*m*cos(delta));
    yC2 = b*m*sin(delta)/(b-2*m*cos(delta));
    plot(xC2, yC2, 'ko');
    
    % O-C2
    plot([0 xC2], [0 yC2]);
    
    % Cerchio centro C2
    rC2 = sqrt((b-xC2)^2+yC2^2);
    plot(xC2+rC2*cos(th), yC2+rC2*sin(th));
    
    % C-P
    plot([m*cos(delta) m*cos(delta)+R*cos(thx)], [m*sin(delta) m*sin(delta)+R*sin(thx)]);

    % O-P
    plot([0 m*cos(delta)+R*cos(thx)],[0 m*sin(delta)+R*sin(thx)]);

    % Punto P
    xP = R*cos(thx)+m*cos(delta);
    yP = R*sin(thx)+m*sin(delta);
    plot(xP, yP, 'bo');
    rho=sqrt(xP^2+yP^2);
    
    ni = atan2(yP,xP);
    
%     % Controllo rho e ni coincide con P
%     plot(rho*cos(ni), rho*sin(ni), 'bx');

    % Punto P2
    xP2 = b^2/rho*cos(-ni);
    yP2 = b^2/rho*sin(-ni);
    plot(xP2, yP2, 'bo');
    
    % O-P2;
    plot([0 xP2], [0 yP2]);

    % Punti B1 e B2
    plot([-b b], [0 0], 'kx');
    
    % Punto T
    xT(i) = xP+xP2;
    yT(i) = yP+yP2;
    
    
    plot([xP xT(i) xP2], [yP yT(i) yP2]);
    plot(xT(i),yT(i), 'ro');

    axis square
    axis(2.5*b*[-1 1 -1 1])
    hold off
    
    % Costruzione profilo
    
    subplot(1, 2, 2)
    plot([-3 3], [0 0], 'k:')
    hold on
    plot([0 0], [-3 3], 'k:')
    
    
    plot(xT(1:i), yT(1:i), 'b-');
    plot(xT(i),yT(i), 'ro');
    hold off
    axis square
    axis(2.5*b*[-1 1 -1 1])
    i=i+1;
    saveas(gcf,sprintf('./img/anim%03.0f.png', (thx-thx0)/pi*180))
    pause(0.001)
end