function Y = adams_moulton(ddy,h,Y)
% 
% Funkcja adams_moulton implementuje jeden krok metody Adamsa-Moultona
% rzędu 4.
% in:
%   ddy - uchwyt do równania różniczkowego postaci:
%   ddy=(x,y,dy)b(x)*dy+c(x)*y+d(x), gdzie a, b, c są funkcjami od x.
%   h - długość kroku
%   Y - wektor postaci [x0,y,dy] z czterema początkowymi przybliżeniami
%   wartości funkcji i jej pochodnej oraz wartością x, w jakim wartości
%   te były obliczone. Y =
%   [x_i-2 y(x_i-2) y'(x_i-2)]
%   [x_i-1 y(x_i-1) y'(x_i-1)]
%   [x_i   y(x_i)   y'(x_i)  ]
%   [x_i+1 y(x_i+1) y'(x_i+1)]
% out:
%    Y - wektor postaci [x_i+1 y(x_i+1) y'(x_i+1)], gdzie x_i+1=x_i+h,
%    y(x_i+1) i  y'(x_i+1) to odpowiednio wartości funkcji y i jej 
%    pochodnej w punkcie x_i+1

% Funkcja pomocnicza reprezentująca pochodną wektora [x,y,y'], 
% czyli [1 y' y'']
F = @(Y_k)[1;Y_k(3);ddy(Y_k(1),Y_k(2),Y_k(3))]';

% Implementacji metody Adamsa-Moultona rzędu 4.
Y = Y(3,:)+h/24*(9*F(Y(4,:))+19*F(Y(3,:))-5*F(Y(2,:))+F(Y(1,:)));

