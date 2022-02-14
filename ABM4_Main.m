function y=ABM4_Main(a,b,c,d,x0,y0,dy0,xMax,N)
% Autor: Jan Cichomski 313201
%
% Funkcja ABM4_Main przybliża wartości funkcji y w N równoodległych
% punktach pomiędzy x0 a xMax, włącznie z x0 i xMax, gdzie y to dokładne 
% rozwiązanie równania różniczkowego drugiego stopnia postaci
% a(x)*y''+b(x)*y'+c(x)*y+d(x) = 0. Funkcja wykorzystuje metodę
% predyktor-korektor Adamsa-Bashfortha-Moultona rzędu 4-tego, gdzie 
% korektorem jest metoda Adamsa-Moultona rzędy 4-tego, a predyktorem 
% metoda Adamsa-Bashfortha również rzędu 4-tego. Do przybliżenia wartości
% funkcji y w pierwszych trzech kolejnych punktach wykorzystywana jest
% metoda Runegego-Kutty 4 rzędu ze wzorem klasyczny.
% in:
%   a, b, c, d - są uchwytami do funkcji od x, 'a' nie może się zerować na
%   przedziale [x0 xMax].
%   x0 - początek przedziału, dla którego przybliżamy wartości y.
%   y0 - warunek początkowy dla y, y(x0)=y0.
%   dy0 - warunek początkowy dla pochodnej funkcji y, y'(x0)=dy0.
%   xMax - koniec przedziału, dla którego przybliżamy wartości y.
%   N - liczba równoodległych punktów między x0 a xMax, w których
%       przybliżamy wartości funkcji y.
% out:
%   y - wektor przybliżonych wartości funkcji y w punktach x0,x1,...,xN,
%`      gdzie xi=x0+i*h, gdze h=(xMax-x0)/N

% Uchwyt funkcji reprezentujący równanie różniczkowe drugiego stopnia
% postaci y''=(-y'*b-y*c-d)/a, gdzie a, b, c, d są funkcjami od x
ddy = @(x,y,dy)-(dy*b(x)+y*c(x)+d(x))/a(x);

% Y będzie wektorem składającym się z wartości x_i, y(x_i) i y'(x_i)
% w N+1 równoodległych punktach:
%[x_0    ...   x_N]
%[y(x_0) ... y(x_N)]
%[y'(x_0)...y'(x_N)]
Y(1,1) = x0;
Y(1,2) = y0;
Y(1,3) = dy0;
% Krok metody.
h = (xMax-x0)/N;

% Startujemy metodę Adamsa-Bashfortha-Moultona obliczając wartości dla
% trzech pierwszych kroków używając metody Rungego-Kutty.
for i=1:3
    Y(i+1,:) = runge_kutta(ddy,h,Y(i,:));
end

% Wykorzystujemy metodę predyktor-korektor Adamsa-Bashfortha-Moultona
% do przybliżenia wartości y we wszystkich kolejnych N-3 punktach.
for i=4:N
    Y(i+1,:) = adams_bashforth(ddy,h,Y(i-3:i,:));
    Y(i+1,:) = adams_moulton(ddy,h,Y(i-2:i+1,:));
end

% Zwracana wartość funkcji y.
y = Y(:,2);


