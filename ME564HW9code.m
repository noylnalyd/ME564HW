%ME564 HW10 P1
clc
clear all
close all
T=1; 
t=0; 
Wr=[T^3/3 T^2/2; T^2/2 T];
expAT=[1 T; 0 1];
B=[0; 1];
x0=[1 1]'; x1=[-1 -1]';
%calculation for original system
for i=1:101; 
 expAt=[1 T-t; 0 1];  % computes expAt as a function of t (j)
 u(i)=(expAt*B)'*inv(Wr)*(x1 - expAT*x0);
 u2(i)=[T-t, 1]*[12/T^3 -6/T^2; -6/T^2 4/T]*(x1-[1 T; 0 1]*x0); 
 u3(i)=24*t-8; 
 xc(i,1:2) = [4*t^3-4*t^2+t+1; 12*t^2-8*t+1];
 tvec(i)=t; 
 t=t+0.01;
end
Wr_adj = -[T T^2/2; T^2/2 T^3/3];
expAT_adj = [1 0; T 1];
C=[1 0]; %=B';
t=0;
for j=1:101
    expAt = [1 t; 0 1];
    v(j) = -C*expAt*inv(Wr_adj)*(x1 - expAT_adj*x0);
    v3(j) = 30*t - 14; 
    pc(j,1:2) = [-15*t^2+14*t+2; 5*t^3-7*t^2-2*t+5];
    tvec(j)=t;
    t=t+0.01;  
end
A=[0 1; 0 0]; B=[0; 1]; C=[1 0]; D=[0]; 
sys=ss(A,B,C,D); 
[y,t,x]=lsim(sys,u,tvec,x0);  % replaces t above to be tvec. 
sys_adj = ss(-A',-C',B',D);
[z,t,p]=lsim(sys_adj,v,tvec,x1);
figure(1);
plot(t,x,t,xc,'o'); figure(gcf);
title('state over time for original system')
xlabel('t'); ylabel('state'); legend('x','xc')
figure(2);
plot(t,p,t,pc,'o');
title('state over time for adjoint system')
xlabel('t'); ylabel('state'); legend('p','pc')
figure(3); 
plot(t,u,t,u3,'o',t,v,t,v3,'o'); %,t,dot(B*u,p')+dot(x',C'*v)); 
legend('u','u3','v','v3');
title('the open-loop control inputs'); 
figure(4); plot(t,x,t,xc,'o',t,p,t,pc,'o'); 
figure(gcf); title('all the states'); 
Annotations
