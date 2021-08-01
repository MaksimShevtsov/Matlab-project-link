%Функция интегрирования системы ДУ

function dy = Solve_Diff_Ur(t,y)
global ma np rk0 g f0 J1 J2 J3 J5 J4 u31 u32 eta31 eta32 mu1 mu2 c1 c2 fixd fixst ke Mf1 Mf2 Mv1 Mv2 
global L1 L2 P1 P2 dy2 dy4 Mw Mf ttt fix mav dypred wp Mp a_vsh b_vsh c_vsh kp Mh Al Mf1max tf1 kwv

dy = zeros(7,1);
%--------------------------------------------------------------------------
%Функции внешних воздействий
wd = y(1);      %угловая скрость двигателя
nd = wd*30/pi;  %частота вращения двигателя
Me = Mp*(a_vsh + b_vsh*nd/np + c_vsh*(nd/np)*(nd/np));      %момент двигателя

if (wd <= wp)
   Mv1 = Me;
else
   Mv1 = Mp-kp*(wd-wp);     %регуляторная ветвь
end;

fv = f0*(1+(0.0216*rk0*y(5))^2);    %коэф. сопротивления качению
Mf = ma*g*fv*rk0*sign(y(5));        %момент, учитывающий сопротивление качению
Mw = kwv*Al*(rk0^3)*(y(5)^2);        %момент, учитывающий сопротивление воздуха
Mv2 = Mf + Mh + Mw;                 %приведенный момент сопротивления

if (t >= 0 && t <= tf1)              %момент трения сцепления
   Mf1 = Mf1max*t/tf1;
else
   Mf1 = Mf1max;
end;

fix = fixd + (fixst-fixd)*exp(-ke*rk0*abs(y(4)-y(5)));      %коэффициент сцепления колеса с дорогой
Mf2 = fix*mav*g*rk0;                %момент сцепления ведущих колес с дорогой

%--------------------------------------------------------------------------
%Формирование матрицы Якоби и вектора внешних воздействий
%матрица Якоби (7х7)
Jacob = [0  -mu1*L1*P1/(J1+J2*L1)  mu1*u31*L1*P1/(J1+J2*L1)                 0                      0  -L1*P1/(J1+J2*L1)  0;
         0  -mu1/(J1*L1+J2)        mu1*u31/(J1*L1+J2)                       0                      0  -1/(J1*L1+J2)      0;
         0  mu1*u31*eta31/J3       -(mu1*u31^2*eta31+mu2/(u32^2*eta32))/J3  mu2/(u32*eta32*J3)     0  (u31*eta31)/J3     -1/(u32*eta32*J3);
         0  0                      (mu2/u32)/(J4+J5*L2)                     -mu2/(J4+J5*L2)        0  0                  1/(J4+J5*L2);
         0  0                      (mu2/u32)*L2*P2/(J4*L2+J5)               -mu2*L2*P2/(J4*L2+J5)  0  0                  L2*P2/(J4*L2+J5);
         0  c1                     -c1*u31                                  0                      0  0                  0;
         0  0                      c2/u32                                   -c2                    0  0                  0];
%вектор внешних водейсвий (7х1)
b     = [(Mv1-Mf1*sign(y(1)-y(2))*(1-L1))/(J1+J2*L1);
         (Mv1*L1*P1+Mf1*sign(y(1)-y(2))*(1-L1))/(J1*L1+J2);
         0;
         (-Mv2*L2*P2-Mf2*sign(y(4)-y(5))*(1-L2))/(J4+J5*L2);
         (-Mv2+Mf2*sign(y(4)-y(5))*(1-L2))/(J4*L2+J5);
         0;
         0];
%--------------------------------------------------------------------------
%Формирование системы ДУ
dy(1) = Jacob(1,1)*y(1) + Jacob(1,2)*y(2) + Jacob(1,3)*y(3) + Jacob(1,4)*y(4) + Jacob(1,5)*y(5) + Jacob(1,6)*y(6) + Jacob(1,7)*y(7) + b(1);
dy(2) = Jacob(2,1)*y(1) + Jacob(2,2)*y(2) + Jacob(2,3)*y(3) + Jacob(2,4)*y(4) + Jacob(2,5)*y(5) + Jacob(2,6)*y(6) + Jacob(2,7)*y(7) + b(2);
dy(3) = Jacob(3,1)*y(1) + Jacob(3,2)*y(2) + Jacob(3,3)*y(3) + Jacob(3,4)*y(4) + Jacob(3,5)*y(5) + Jacob(3,6)*y(6) + Jacob(3,7)*y(7) + b(3);
dy(4) = Jacob(4,1)*y(1) + Jacob(4,2)*y(2) + Jacob(4,3)*y(3) + Jacob(4,4)*y(4) + Jacob(4,5)*y(5) + Jacob(4,6)*y(6) + Jacob(4,7)*y(7) + b(4);
dy(5) = Jacob(5,1)*y(1) + Jacob(5,2)*y(2) + Jacob(5,3)*y(3) + Jacob(5,4)*y(4) + Jacob(5,5)*y(5) + Jacob(5,6)*y(6) + Jacob(5,7)*y(7) + b(5);
dy(6) = Jacob(6,1)*y(1) + Jacob(6,2)*y(2) + Jacob(6,3)*y(3) + Jacob(6,4)*y(4) + Jacob(6,5)*y(5) + Jacob(6,6)*y(6) + Jacob(6,7)*y(7) + b(6);
dy(7) = Jacob(7,1)*y(1) + Jacob(7,2)*y(2) + Jacob(7,3)*y(3) + Jacob(7,4)*y(4) + Jacob(7,5)*y(5) + Jacob(7,6)*y(6) + Jacob(7,7)*y(7) + b(7);

dy2 = Jacob(2,1)*y(1) + Jacob(2,2)*y(2) + Jacob(2,3)*y(3) + Jacob(2,4)*y(4) + Jacob(2,5)*y(5) + Jacob(2,6)*y(6) + Jacob(2,7)*y(7) + b(2);
dy4 = Jacob(4,1)*y(1) + Jacob(4,2)*y(2) + Jacob(4,3)*y(3) + Jacob(4,4)*y(4) + Jacob(4,5)*y(5) + Jacob(4,6)*y(6) + Jacob(4,7)*y(7) + b(4);

dypred = dy;

ttt =t
