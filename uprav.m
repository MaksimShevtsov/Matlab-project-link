%������ ������������� � ������������
clear all;
clc;
ma = 11000;      %������ ����� ����
ma1 = 4000;     %������������� ���� ��������� ����
ma2 = 7000;     %������������� ���� ��������� ����
L = 4;        %���� ����
L1 = (ma2*L)/ma     %���������� ������ ���� ��� ��������� ����
L2 = (ma1*L)/ma     %���������� ������ ���� ��� ��������� ����
hc = (0.3)*L;  %������ ������ ����
B1 = 2.116;     %����� �������� �����
l0 = 0.7*B1;    %���������� ����� ����� ��������
up = 30;        %������������ ����� �������� �������
nsh1 = 2;       %���������� ��� i-�� �����
nsh2 = 4;        %���������� ��� i-�� �����
g = 9.8;        %��������� ���������� �������
kuv = 120000;       %����������� ������������� ����� ����   
rs = 0.491;       %��������� ������ ������
kuv1 = kuv*nsh1;    
kuv2 = kuv*nsh2;

%����������� ����� ������ ��������  ����������� �����

tetan=[];
tetav = 5*pi/180:5*pi/180:45*pi/180;        %���� �������� ������� �����
tetan = atan(1./((l0/L)+(cos(tetav)./sin(tetav))))  %���� �������� ���������� �����
tetan

figure(1);
plot (tetav*180/pi, tetan*180/pi);
grid on;
legend('tetan');

%������ �������� ��� ����� ����� �����

Rp = L.*((cos(tetan)./sin(tetan))+(cos(tetav)./sin(tetav)))./2  %������ �������� ���� ��� ����� ����� �����
Rp1 = L./(tan(tetan+tetav)./2)

figure(2);
plot (tetav*180/pi, Rp, tetav*180/pi, Rp1 );
grid on;
legend('Rp', 'Rp1');
xlabel('���� �������� ����������� ������','FontName','Arial Cyr')
ylabel('������ �������� ��� ����� ����� ����� �','FontName','Arial Cyr')

%���������������� ���������� � ������ ����� �����

teta1 = 2*pi/180    %���� �������� �����
teta2 = 4*pi/180
teta3 = 6*pi/180

a=0;
for Vx=0:2:30
    a=a+1
    Rgz1(a) = L/teta1*(1+(ma1/kuv1-ma2/kuv2).*Vx.^2/L)      %������ �������� � ������ ����� �����
    Rgz2(a) = L/teta2*(1+(ma1/kuv1-ma2/kuv2).*Vx.^2/L)
    Rgz3(a) = L/teta3*(1+(ma1/kuv1-ma2/kuv2).*Vx.^2/L)
end
Vx=0:2:30
figure(3);
plot (Vx, Rgz1, Vx, Rgz2, Vx, Rgz3);
grid on;
legend('Rgz1', 'Rgz2', 'Rgz3');
xlabel('�������� �/�','FontName','Arial Cyr')
ylabel('������ �������� � ','FontName','Arial Cyr')

%���� ����� ����� ��������� � ������� ������
Vx=0:2:30
buv11 = ma1*Vx.^2./(Rgz1*kuv1)      %���� ����� �����
buv12 = ma1*Vx.^2./(Rgz2*kuv1)
buv13 = ma1*Vx.^2./(Rgz3*kuv1)
buv21 = ma2*Vx.^2./(Rgz1*kuv2)
buv22 = ma2*Vx.^2./(Rgz2*kuv2)
buv23 = ma2*Vx.^2./(Rgz3*kuv2)

figure(4);
plot (Vx, buv11*180/pi, Vx, buv12*180/pi, Vx, buv13*180/pi, Vx, buv21*180/pi, Vx, buv22*180/pi, Vx, buv23*180/pi);
grid on;
legend('buv11', 'buv12', 'buv13', 'buv21', 'buv22', 'buv23');
xlabel('�������� �/�','FontName','Arial Cyr')
ylabel('���� ����� �����','FontName','Arial Cyr')

%����������� ����������� ����������������
nupov = (ma2/kuv2)/(ma1/kuv1)
nupov = 0.9359     % ������������� ����������������

%�������������� ����������� ������������� ����������
Vxx = 60*1000/3600
Up = 30     %������������ ����� �������� �������
alphap=10:10:120        %���� �������� �������� ������
tetax = (alphap/Up)*pi/180
Rgzx = L./tetax*(1+(ma1/kuv1-ma2/kuv2).*Vxx.^2/L)
ay = Vxx^2./Rgzx        %������� ��������� ����
K = 1/L*tan((alphap*pi/180)/Up)+ (ma*Vxx^2*(L1*kuv1-L2*kuv2))./(Rgzx*L^2*kuv1*kuv2) %�������� ����������
Kx = 1/L*tan((alphap*pi/180)/Up)

figure(5);
plot (alphap, K, alphap, Kx);
grid on;
legend('K', 'Kx');
xlabel('���� �������� �������� ������','FontName','Arial Cyr')
ylabel('�������� ����������','FontName','Arial Cyr')
title('����������������','FontName','Arial Cyr')

figure(6);
plot (alphap, Rgzx);
grid on;
legend('Rgzx');
xlabel('���� �������� �������� ������','FontName','Arial Cyr')
ylabel('������ ��������','FontName','Arial Cyr')

figure(7);
plot (alphap, ay);
grid on;
legend('ay');
xlabel('���� �������� �������� ������','FontName','Arial Cyr')
ylabel('������� ��������� ����','FontName','Arial Cyr')

alphapmin = ((0.72*L+0.2)*Up*10^(-2))*180/pi
alphapmax = ((0.72*L+2.6)*Up*10^(-2))*180/pi
alphapmin = 55.293
alphapmax = 96.546
alphapx = 56.667  %�������� ���� �������� �������� ������ ��������������� ��������� � = 2 �/(�^2)
% => ����������� ������������� ���� ������������� ����������� �����������

%������������ ����������
%����������� �������� �� �������� �������������

tetaxx = 2*pi/180:1*pi/180:10*pi/180
Vkrop = sqrt(B1*L*g./(2*hc.*tetaxx))    %����������� �������� �� �������� �������������
figure(8);
plot (tetaxx*180/pi, Vkrop);
grid on;
legend('Vkrop');
xlabel('���� �������� �����','FontName','Arial Cyr')
ylabel('��������','FontName','Arial Cyr')

%����������� �������� �� �������� ����������
phiy1 = 0.2     %���� ���������
phiy2 = 0.4
phiy3 = 0.6
phiy4 = 0.8
Vkrphi1 = sqrt(phiy1*L*g./tetaxx)   %����. �������� �� �������� ����������
Vkrphi2 = sqrt(phiy2*L*g./tetaxx)
Vkrphi3 = sqrt(phiy3*L*g./tetaxx)
Vkrphi4 = sqrt(phiy4*L*g./tetaxx)

figure(9);
plot (tetaxx*180/pi, Vkrphi1, tetaxx*180/pi, Vkrphi2, tetaxx*180/pi, Vkrphi3, tetaxx*180/pi, Vkrphi4);
grid on;
legend('Vkrphi1', 'Vkrphi2', 'Vkrphi3', 'Vkrphi4');
xlabel('���� �������� �����','FontName','Arial Cyr')
ylabel('��������','FontName','Arial Cyr')

nupy = 0.5*B1/hc        %����������� ���������� ������������
Bkrop = atan(nupy)*180/pi      %����������� ���� �������� �� �������� �������������
Bkrphi1 = atan(phiy1)*180/pi    %����������� ���� �������� �� �������� ����������
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
xlabel('���� ���������','FontName','Arial Cyr')
ylabel('����������� ���� �������� �� �������� ����������','FontName','Arial Cyr')

 B = 2.5;        %������ ����������:
 b = 0.275       %������ ����
 mp = 9296.9      %�������������� ����� ����
 L1 = (ma2*L)/ma     %���������� ������ ���� ��� ��������� ����
 L2 = (ma1*L)/ma     %���������� ������ ���� ��� ��������� ����
 hlya1 = 0.505
 hlya2 = 0.505
 E = atan((hlya1+hlya2)/L)*180/pi  %���� ������� ��� �����
 hlya = (hc-(L1*hlya2-L2*hlya1)/L)*cos(E)   %����� �����
 cx = 9.9486e+005        %���� ��� ��������� ������� ��������� ��������
 nur = 1.15     %��� �������� ������� 
 Bp = B-2*b          %���������� ����� �������� ����������
 clya = 0.5*cx*nur*Bp^2      %���� ������� ��������� ��������
 Lyakp = hlya*mp*2/(clya-hlya*mp*g)*180/pi     %���� ����� ������

