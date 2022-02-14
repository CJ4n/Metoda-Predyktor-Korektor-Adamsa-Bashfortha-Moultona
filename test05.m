% Autor: Jan Cichomski 313201
% Skrypt sprawdza co się dzieje z rozwiązaniem równania y''+3y'+3y+x^2=0,
% gdy wyłączymy korektor w metodzie predyktor-korektor
% Adamsa-Bashforha-Moultona.

clearvars
close all
x0=0;
xMax=6;
a=@(x)1;
b=@(x)3;
c=@(x)3;
d=@(x)x^2;
y0=1;
dy0=1;
% Dokładne rozwiązanie równania y''+3y'+3y+x^2=0
sol=@(x)(2*x)/3 + (13*exp(-(3*x)/2).*cos((3^(1/2)*x)/2))/9 - x.^2/3 +...
    (5*3^(1/2)*exp(-(3*x)/2).*sin((3^(1/2)*x)/2))/3 - 4/9;

Y(1,1) = x0;
Y(1,2) = y0;
Y(1,3) = dy0;
ddy = @(x,y,dy)(-dy*b(x)-y*c(x)-d(x))/a(x);

figure(1)
hold on
args=linspace(x0,xMax,1000);
plot(args,sol(args),'LineWidth',3);
legend("Dokładne rozwiązanie y''+3y'+3y+x^2=0",'Location','southwest')
ylim([-18 4])
title("Metoda Adamsa-Bashfortha-Moultona bez korektora w zależności od N")
figure(2)
hold on
plot(args,sol(args),'LineWidth',3);
legend("Dokładne rozwiązanie y''+3y'+3y+x^2=0",'Location','southwest')
ylim([-18 4])
title("Metoda Adamsa-Bashfortha-Moultona w zależności od N")

for k=3:13
    N = 2*k;
    h = (xMax-x0)/N;
    
    for i=1:3
        Y(i+1,:) = runge_kutta(ddy,h,Y(i,:));
    end
    % Metdoa Adamsa-Bashfortha-Moultona
    for i=4:N
        Y(i+1,:) = adams_bashforth(ddy,h,Y(i-3:i,:));
        Y(i+1,:) = adams_moulton(ddy,h,Y(i-2:i+1,:));
    end   
    ABMout = Y(:,2);
    % Metdoa Adamsa-Bashfortha
    for i=4:N
        Y(i+1,:) = adams_bashforth(ddy,h,Y(i-3:i,:));
        
    end
    ABout = Y(:,2);
    
    args = linspace(x0,xMax,N+1);  
    figure(1)
    plot(args,ABout,'DisplayName',"N="+int2str(N));    
    figure(2)
    plot(args,ABMout,'DisplayName',"N="+int2str(N));
end

figure(1)
movegui([300 550]);
figure(2)
movegui([900 550]);
fprintf("Jak widać, wyłączenie korektora powoduje,\nże rozwiązanie jest ")
fprintf("bardzo niestabilne w porównaniu\ndo rozwiązania z włączonym ")
fprintf("korektorem. Nawet dla\nN=20 rozwiązanie jest nadal bardzo nie ")
fprintf("dokładne \ni niestabilne. Dopiero dla N=24 rozwiązanie")
fprintf("\nbez korektora można uznać z zadowalająco dokładne.\n")
fprintf("Z kolei rozwiązanie z korektorem jest zadowalająco\n")
fprintf("dokładne już dla N=8.\n")

