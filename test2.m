% test01
% test02
% test03
% test04
% test05


clearvars
close all


x0=1;
xMax=5;
a=@(x)1;
b=@(x)0;
c=@(x)0;
d=@(x)-x^2 + 2*x - 5;

N = 5;
h=(xMax-x0)/(N-1);
ddy = @(x,y,dy)(-dy*b(x)-y*c(x)-d(x))/a(x);

Y0(1,1)=x0;
Y0(1,2)=1;
Y0(1,3)=1;



for i=1:3
Y0(i+1,:)=runge_kutta(ddy,h,Y0(i,:))
end
for i=4:N-1
Y0(i+1,:) =adams_bashforth(ddy,h,Y0)
end

syms x y(x)
dy=diff(y,x);
ode=diff(dy,x)==(-sym(d)-sym(c)*y-sym(b)*dy)/sym(a);
cond1=y(0)==1;
cond2=dy(0)==1;
sol=dsolve(ode,[cond1 cond2]);

% inicjalizujemy potrzebne zmienne

 % wyliczone poczatkowe wartosci funkcji


% uruchamiamy nasz progra

% wyliczona dokladna wartosc calki funkcji fun na przedziale [5, 4]
intergal =31/3;

% wyliczamy blad dokladnosci calki
Y0(4,2)-Y0(3,2)
error = abs(intergal - (Y0(4,2)-Y0(3,2)));

fprintf('Blad wynosi: %d\n', error);
 
 
 
 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USUNAĆ  double(max(abs(q-Y(:,2)')))/h^4
% %    
% %    
% %     clearvars
% % close all
% % x0=0;
% % xMax=5;
% % a=@(x)1;
% % b=@(x)1;
% % c=@(x)2;
% % d=@(x)exp(x);
% % y0=1;
% % dy0=1;
% % syms x y(x)
% % dy=diff(y,x);
% % ode=diff(dy,x)==(-sym(d)-sym(c)*y-sym(b)*dy)/sym(a);
% % cond1=y(0)==y0;
% % cond2=dy(0)==dy0;
% % sol=dsolve(ode,[cond1 cond2]);
% % sol
% % hold on
% % for k=4:8
% %     if k==5
% %         %                   continue
% %     end
% %     N=2^k;
% %     h=(xMax-x0)/N;
% %     out = ABM4_Main(a,b,c,d,y0,dy0,x0,xMax,N);
% %     
% %     
% %     args = linspace(x0,xMax,N);
% %     
% %     for i=1:N
% %         x=args(i);
% %         q(i)=subs(sol);
% %     end
% %     
% %     %    double(max(abs(q-out')))
% %     ratio(k) = double(max(abs(q-out')));
% %     sum((q-out'));
% %     double(ratio(k));
% %     % double(suma)
% %     %         plot(args,abs(q-out'),'x')
% %     double(ratio(k-1))/double(ratio(k))
% %     if(k>=5)
% %         plot(k,double(ratio(k-1))/double(ratio(k)),'x')
% %     end
% %     
% % end
% % 
% % 
% % 
% % 
% % 
% 
% 
% 
% % 
% % clearvars
% % close all
% % x0=0;
% % xMax=5;
% % a=@(x)1;
% % b=@(x)1;
% % c=@(x)2;
% % d=@(x)exp(x);
% % y0=1;
% % dy0=1;
% % syms x y(x)
% % dy=diff(y,x);
% % ode=diff(dy,x)==(-sym(d)-sym(c)*y-sym(b)*dy)/sym(a);
% % cond1=y(0)==y0;
% % cond2=dy(0)==dy0;
% % sol=dsolve(ode,[cond1 cond2]);
% % sol=@(x)(exp(-x/2)*(35*cos((7^(1/2)*x)/2) + 15*7^(1/2)*sin((7^(1/2)*x)/2)...
% %     - 7*exp((3*x)/2)*cos((7^(1/2)*x)/2)^2 -...
% %     7*exp((3*x)/2)*sin((7^(1/2)*x)/2)^2))/28;
% % hold on
% % 
% % for k=4:10
% %     N=2^k;
% % %     h=(xMax-x0)/N;    
% %     estination = ABM4_Main(a,b,c,d,y0,dy0,x0,xMax,N);    
% %     args = linspace(x0,xMax,N);
% %     
% %     for i=1:N
% % %         x=args(i);
% % %         q(i)=subs(sol);
% % q(i)=sol(args(i));
% %     end
% %     ratio(k) = double(max(abs(q-estination')));
% %     if(k>=5)
% %         plot(k,double(ratio(k-1))/double(ratio(k)),'x')
% %     end
% %     
% % end
% 
% 
% 
