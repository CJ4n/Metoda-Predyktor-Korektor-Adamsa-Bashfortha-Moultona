function Y = runge_kutta(ddy,h,Y)
% 
% Funkcja runge_kutta implementuje jeden krok metody Rungegu-Kutty rzędu 
% 4 ze wzorem kalsycznym.
% in:
%   ddy - uchwyt do równania różniczkowego postaci:
%         ddy=(x,y,dy)b(x)*dy+c(x)*y+d(x), gdize a,b,c są funkcjami od x.
%   h - długość kroku
%   Y - wektor postaci [x_i,y_i,dy_i] z początkowymi przybliżeniami 
%       wartości funkcji y_i i jej pochodnej dy_i oraz wartością x_i 
%       w jakim wartości te były obliczone.
% out:
%    Y - wektor postaci [x_i+1,y_i+1,dy_i+1], gdzie x_i+1=x_i+h, y_i+1
%        i dy_i+1 to odpowiednio wartości funkcji y i jej pochodnej
%        w punkcie x_i+1

% Funkcja pomocnicza reprezentująca pochodą wektora [x,y,y'], 
% czyli [1,y',y'']
F = @(Y_k)[1;Y_k(3);ddy(Y_k(1),Y_k(2),Y_k(3))]';

% Implementacjia jednego kroku metody Rungego-Kutty rzędu 4.

K0 = F(Y);
K1 = F(Y+h/2*K0);
K2 = F(Y+h/2*K1);
K3 = F(Y+h*K2);
Y = Y+(K0+2*K1+2*K2+K3)/6*h;

