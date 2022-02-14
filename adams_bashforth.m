function Y = adams_bashforth(ddy,h,Y)
% Autor: Jan Cichomski 313201
% 
% Funkcja adams_bashforth implementuje jeden krok metody Adamsa-Bashfortha
% rzędu 4.
% in:
%   ddy - uchwyt do równania różniczkowego postaci:
%   ddy=(x,y,dy)b(x)*dy+c(x)*y+d(x), gdize a, b, c są funkcjami od x.
%   h - długość kroku
%   Y - wektor z czterema początkowymi przybliżeniami wartości
%   funkcji i jej pochodnej oraz wartością x, w jakim wartości te były
%   obliczone. Y =
%   [x_i-3 y(x_i-3) y'(x_i-3)]
%   [x_i-2 y(x_i-2) y'(x_i-2)]
%   [x_i-1 y(x_i-1) y'(x_i-1)]
%   [x_i   y(x_i)   y'(x_i)  ]
% out:
%    Y - wektor postaci [x_i+1 y(x_i+1) y'(x_i+1)], gdzie x_i+1=x_i+h,
%    y(x_i+1) i  y'(x_i+1) to odpowiednio wartości funkcji y i jej 
%    pochodnej w punkcie x_i+1

% Funkcja pomocnicza reprezentująca pochodną wektora [x,y,y'], 
% czyli [1 y' y'']
F = @(Y_k)[1;Y_k(3);ddy(Y_k(1),Y_k(2),Y_k(3))]';

% Implementacja metody Adamsa-Bashfortha, 
Y = Y(4,:)+h/24*(55*F(Y(4,:))-59*F(Y(3,:))+37*F(Y(2,:))-9*F(Y(1,:)));

