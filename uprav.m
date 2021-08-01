%оценка управляемости и устойчивости
clear all;
clc;
ma = 11000;      %полная масса авто
ma1 = 4000;     %распределение масс груженого авто
ma2 = 7000;     %распределение масс груженого авто
L = 4;        %база авто
L1 = (ma2*L)/ma     %координата центра масс для груженого авто
L2 = (ma1*L)/ma     %координата центра масс для груженого авто
hc = (0.3)*L;  %высота центра масс
B1 = 2.116;     %колея передних колес
l0 = 0.7*B1;    %расстояние между осями шкворней
up = 30;        %передаточное число рулевого привода
nsh1 = 2;       %количество шин i-го моста
nsh2 = 4;        %количество шин i-го моста
g = 9.8;        %ускорение свободного падения
kuv = 120000;       %коэффициент сопротивления уводу шины   
rs = 0.491;       %свободный радиус колеса
kuv1 = kuv*nsh1;    
kuv2 = kuv*nsh2;

%зависимость между углами поворота  управляемых колес

tetan=[];
tetav = 5*pi/180:5*pi/180:45*pi/180;        %угол поворота внешних колес
tetan = atan(1./((l0/L)+(cos(tetav)./sin(tetav))))  %угол поворота внутренних колес
tetan

figure(1);
plot (tetav*180/pi, tetan*180/pi);
grid on;
legend('tetan');

%радиус поворота без учета увода колес

Rp = L.*((cos(tetan)./sin(tetan))+(cos(tetav)./sin(tetav)))./2  %радиус поворота авто без учета увода колес
Rp1 = L./(tan(tetan+tetav)./2)

figure(2);
plot (tetav*180/pi, Rp, tetav*180/pi, Rp1 );
grid on;
legend('Rp', 'Rp1');
xlabel('угол поворота внутреннего колеса','FontName','Arial Cyr')
ylabel('радиус поворота без учета увода колес м','FontName','Arial Cyr')

%поворачиваемость автомобиля с учетом увода колес

teta1 = 2*pi/180    %угол поворота колес
teta2 = 4*pi/180
teta3 = 6*pi/180

a=0;
for Vx=0:2:30
    a=a+1
    Rgz1(a) = L/teta1*(1+(ma1/kuv1-ma2/kuv2).*Vx.^2/L)      %радиус поворота с учетом увода колес
    Rgz2(a) = L/teta2*(1+(ma1/kuv1-ma2/kuv2).*Vx.^2/L)
    Rgz3(a) = L/teta3*(1+(ma1/kuv1-ma2/kuv2).*Vx.^2/L)
end
Vx=0:2:30
figure(3);
plot (Vx, Rgz1, Vx, Rgz2, Vx, Rgz3);
grid on;
legend('Rgz1', 'Rgz2', 'Rgz3');
xlabel('скорость м/с','FontName','Arial Cyr')
ylabel('радиус поворота м ','FontName','Arial Cyr')

%углы увода колес переднего и заднего мостов
Vx=0:2:30
buv11 = ma1*Vx.^2./(Rgz1*kuv1)      %угол увода колес
buv12 = ma1*Vx.^2./(Rgz2*kuv1)
buv13 = ma1*Vx.^2./(Rgz3*kuv1)
buv21 = ma2*Vx.^2./(Rgz1*kuv2)
buv22 = ma2*Vx.^2./(Rgz2*kuv2)
buv23 = ma2*Vx.^2./(Rgz3*kuv2)

figure(4);
plot (Vx, buv11*180/pi, Vx, buv12*180/pi, Vx, buv13*180/pi, Vx, buv21*180/pi, Vx, buv22*180/pi, Vx, buv23*180/pi);
grid on;
legend('buv11', 'buv12', 'buv13', 'buv21', 'buv22', 'buv23');
xlabel('скорость м/с','FontName','Arial Cyr')
ylabel('угол увода колес','FontName','Arial Cyr')

%коэффициент статической поворачиваемости
nupov = (ma2/kuv2)/(ma1/kuv1)
nupov = 0.9359     % недостаточная поворачиваемость

%характеристика траекторной управляемости автомобиля
Vxx = 60*1000/3600
Up = 30     %передаточное число рулевого привода
alphap=10:10:120        %угол поворота рулевого колеса
tetax = (alphap/Up)*pi/180
Rgzx = L./tetax*(1+(ma1/kuv1-ma2/kuv2).*Vxx.^2/L)
ay = Vxx^2./Rgzx        %боковое ускорение авто
K = 1/L*tan((alphap*pi/180)/Up)+ (ma*Vxx^2*(L1*kuv1-L2*kuv2))./(Rgzx*L^2*kuv1*kuv2) %кривизна траектории
Kx = 1/L*tan((alphap*pi/180)/Up)

figure(5);
plot (alphap, K, alphap, Kx);
grid on;
legend('K', 'Kx');
xlabel('угол поворота рулевого колеса','FontName','Arial Cyr')
ylabel('кривизна траектории','FontName','Arial Cyr')
title('поворачиваемость','FontName','Arial Cyr')

figure(6);
plot (alphap, Rgzx);
grid on;
legend('Rgzx');
xlabel('угол поворота рулевого колеса','FontName','Arial Cyr')
ylabel('радиус поворота','FontName','Arial Cyr')

figure(7);
plot (alphap, ay);
grid on;
legend('ay');
xlabel('угол поворота рулевого колеса','FontName','Arial Cyr')
ylabel('боковое ускорение авто','FontName','Arial Cyr')

alphapmin = ((0.72*L+0.2)*Up*10^(-2))*180/pi
alphapmax = ((0.72*L+2.6)*Up*10^(-2))*180/pi
alphapmin = 55.293
alphapmax = 96.546
alphapx = 56.667  %значение угла поворота рулевого колеса соответствующая ускорению а = 2 м/(с^2)
% => траекторная управляемость авто удовлетворяет нормативным требованиям

%устойчивость автомобиля
%критическая скорость по боковому опрокидыванию

tetaxx = 2*pi/180:1*pi/180:10*pi/180
Vkrop = sqrt(B1*L*g./(2*hc.*tetaxx))    %критическая скорость по боковому опрокидыванию
figure(8);
plot (tetaxx*180/pi, Vkrop);
grid on;
legend('Vkrop');
xlabel('угол поворота колес','FontName','Arial Cyr')
ylabel('скорость','FontName','Arial Cyr')

%критическая скорость по боковому скольжению
phiy1 = 0.2     %коэф сцепления
phiy2 = 0.4
phiy3 = 0.6
phiy4 = 0.8
Vkrphi1 = sqrt(phiy1*L*g./tetaxx)   %крит. скорость по боковому скольжению
Vkrphi2 = sqrt(phiy2*L*g./tetaxx)
Vkrphi3 = sqrt(phiy3*L*g./tetaxx)
Vkrphi4 = sqrt(phiy4*L*g./tetaxx)

figure(9);
plot (tetaxx*180/pi, Vkrphi1, tetaxx*180/pi, Vkrphi2, tetaxx*180/pi, Vkrphi3, tetaxx*180/pi, Vkrphi4);
grid on;
legend('Vkrphi1', 'Vkrphi2', 'Vkrphi3', 'Vkrphi4');
xlabel('угол поворота колес','FontName','Arial Cyr')
ylabel('скорость','FontName','Arial Cyr')

nupy = 0.5*B1/hc        %коэффициент поперечной устойчивости
Bkrop = atan(nupy)*180/pi      %критический угол косогора по боковому опрокидыванию
Bkrphi1 = atan(phiy1)*180/pi    %критический угол косогора по боковому скольжению
Bkrphi2 = atan(phiy2)*180/pi
Bkrphi3 = atan(phiy3)*180/pi
Bkrphi4 = atan(phiy4)*180/pi
Bkrphi1 = 11.31
Bkrphi2 = 21.801
Bkrphi3 = 30.964
Bkrphi4 = 38.66

phiyy=0.2:0.2:0.8
Bkrphix = atan(phiyy)*180/pi

figure(10);
plot (phiyy, Bkrphix);
grid on;
legend('Bkrphix');
xlabel('коэф сцепления','FontName','Arial Cyr')
ylabel('критический угол косогора по боковому скольжению','FontName','Arial Cyr')

 B = 2.5;        %ширина автомобиля:
 b = 0.275       %ширина шины
 mp = 9296.9      %подрессоренная масса авто
 L1 = (ma2*L)/ma     %координата центра масс для груженого авто
 L2 = (ma1*L)/ma     %координата центра масс для груженого авто
 hlya1 = 0.505
 hlya2 = 0.505
 E = atan((hlya1+hlya2)/L)*180/pi  %угол наклона оси крена
 hlya = (hc-(L1*hlya2-L2*hlya1)/L)*cos(E)   %плечо крена
 cx = 9.9486e+005        %коэф сум жесткости упругих элементов подвески
 nur = 1.15     %для листовой рессоры 
 Bp = B-2*b          %расстояние между упругими элементами
 clya = 0.5*cx*nur*Bp^2      %коэф угловой жесткости подвески
 Lyakp = hlya*mp*2/(clya-hlya*mp*g)*180/pi     %угол крена кузова

