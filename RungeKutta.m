function [T,Y] = RungeKutta(tspan,y0,h)
global J2 J4 u31 u32 mu1 mu2 Mf1 Mf2 L1 L2 P1 P2 Md2 L1rez L2rez P1rez P2rez Pf1 Pf2 Wf1 Wf2 dy

t = tspan(1);
i = 0;
Y = zeros((tspan(2)-tspan(1))/h+1,length(y0));
K1 = zeros((tspan(2)-tspan(1))/h+1,length(y0));
K2 = zeros((tspan(2)-tspan(1))/h+1,length(y0));
K3 = zeros((tspan(2)-tspan(1))/h+1,length(y0));
K4 = zeros((tspan(2)-tspan(1))/h+1,length(y0));
dyRK = zeros(7);
Pf1 = zeros((tspan(2)-tspan(1))/h+1,length(y0));
Pf2 = zeros((tspan(2)-tspan(1))/h+1,length(y0));
Wf1 = zeros((tspan(2)-tspan(1))/h+1,length(y0));
Wf2 = zeros((tspan(2)-tspan(1))/h+1,length(y0));
dy = zeros((tspan(2)-tspan(1))/h+1,length(y0));
while t < tspan(2)
    i = i + 1;
    if i == 1
        T(1) = tspan(1);
        Vk = y0;
        Y(1,:) = y0;
    else
        deltaw1 = 0.1;
        deltaw2 = 0.001;
        Md1 = mu1*(Y(i-1,2)-Y(i-1,3)*u31);      %момент диссипативного элемента 1
        Md2 = mu2*(Y(i-1,3)/u32-Y(i-1,4));      %момент диссипативного элемента 2

        %функция вычисления дискретной функции L1
        if (abs(Y(i-1,1) - Y(i-1,2)) <= deltaw1) & (Mf1 >= abs(Y(i-1,6)+Md1+dyRK(2)*J2))
           L1 = 1;
           Y(i-1,2) = Y(i-1,1);
        else
           L1 = 0;
        end;
        
        %функция вычисления дискретной функции L2
        if (abs(Y(i-1,4) - Y(i-1,5)) <= deltaw2) & (Mf2>=abs(Y(i-1,7)+Md2-dyRK(4)*J4))
           L2 = 1;
           Y(i-1,4) = Y(i-1,5);
        else
           L2 = 0;
        end;

        P1 = 1;
        P2 = 1;
        
        %Реализация метода интегрирования Рунге-Кутта
        Vk = Y(i-1,:);
        K1 = rot90(h*Solve_Diff_Ur(t,Vk));
        K2 = rot90(h*Solve_Diff_Ur(t+0.5*h,Vk+0.5*K1));
        K3 = rot90(h*Solve_Diff_Ur(t+0.5*h,Vk+0.5*K2));
        K4 = rot90(h*Solve_Diff_Ur(t + h,Vk + K3));
        Y(i,:) = Vk + (K1 + 2*K2 + 2*K3 + K4)/6;
        T(i) = t;
        dyRK = rot90(Solve_Diff_Ur(t,Y(i,:)));
        L1rez(i) = L1;
        L2rez(i) = L2;
        P1rez(i) = P1;
        P2rez(i) = P2;
        Pf1(i) = abs(Mf1*(Y(i,2)-Y(i,1)));
        Pf2(i) = abs(Mf2*(Y(i,5)-Y(i,4)));
        Wf1(i) = Wf1(i-1) + Pf1(i)*h;
        Wf2(i) = Wf2(i-1) + Pf2(i)*h;
        dy(i,:) = K1/h;
    end;
    t = t + h;
end;
