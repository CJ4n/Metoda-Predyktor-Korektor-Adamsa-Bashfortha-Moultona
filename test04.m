% Autor: Jan Cichomski 313201
% Skrypt testuje stopień zbieżności funkcji ABM4_Main, poprzez obliczanie
% iloczynów błędów globalnych dla N i 2N punktów podziału. Spodziewana
% wartość iloczynu to 16.
clearvars
close all
hold on

x0=0;
xMax=5;
a=@(x)1;
b=@(x)1;
c=@(x)2;
d=@(x)exp(x);
y0=1;
dy0=1;
% Dokładne rozwiązanie równania y''+y'+2y=exp(x).
sol=@(x)(exp(-x/2).*(35*cos((7.^(1/2)*x)/2) + 15*7.^(1/2)*...
    sin((7^(1/2)*x)/2) - 7*exp((3*x)/2).*cos((7.^(1/2)*x)/2).^2 -...
    7*exp((3*x)/2).*sin((7.^(1/2)*x)/2).^2))/28;

globalError=zeros(1,7);
for k=2:7
    N=2^(k+3);
    % Przybliżenie wartości y w zadanych N punktach metodą
    % predyktor-koretkor Adamsa-Bashfortha-Moultona rzędu 4.
    Y = ABM4_Main(a,b,c,d,x0,y0,dy0,xMax,N);
        
    % Obliczenie dokładnych wartości funkcji w zadanych N punktach.
    args = linspace(x0,xMax,N+1);
    exact = sol(args);
    % Obliczenie błędu globalnego.
    globalError(k) = double(max(abs(exact-Y')));
    if(k>=3)
        plot(N,double(globalError(k-1))/double(globalError(k)),'x')
    end   
end

title("Iloczyn błędów globalnych ((i-1)-ego do i-tego) od liczby punktów")
xlabel("N")
ylabel("Iloczyn błędów globalnych ((i-1)-ego do i-tego)")
ylim([0 20])
fprintf("--------------TEST04--------------\n")
fprintf("Iloczyn błędów globalnych dla N i 2N punktów podziału:\n")
for i=3:7
    fprintf("N = %d oraz N = %d jest równy   %f \n",2^(i-1+3),2^(i+3),...
        double(globalError(i-1))/double(globalError(i)))
end
fprintf("gdzie N to liczba punktów podziału,\nw których ")
fprintf("obliczamy wartości.\n")


