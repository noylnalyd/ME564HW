%% Problem 4
clc, clear
syms t;
A = [-4/t, -2/t/t;...
        1, 0];
phi1 = [1/t/t; -1/t];
phi2 = [2/t/t/t; -1/t/t];
[eigv,eigs] = eig(A)
%% Problem 7
clc, clear

A = [-2, 0; 0 0];

[eigv,eigs] = eig(A)
%% Problem 2
clc, clear all
A = [0.0974 -0.1178 0.7876 -0.1168 0.0178;
    0.1291 -0.1174 1.2850 0.0302 0.0971;
    0.0528 0.1119 0.1325 0.7668 0.0637;
    0.0424 0.2647 0.2806 1.7644 0.1195];
[Um,Sm,Vm] = svd(A);

V1 = orth(A')
V1 = V1(:,1:2)
V2 = null(V1')
V = [V1,V2]

s2 = V^-1 * A' * A * V
sinv = zeros(2,2);
for i=1:2
    sinv(i,i) = 1/sqrt(s2(i,i));
end

U1 = A*V1*sinv
U2 = null(U1')
U = [U1,U2]

sigma = zeros(4,5);
for i=1:2
    sigma(i,i) = sqrt(s2(i,i))
end
Ar = U*sigma*V'

Anorm = norm(A-Ar,"fro")
Um;
Vm;
Sm;

%% Problem 3
clc, clear all, close all
A = [-2.6 0.25 -38 0;
    -0.075 -0.27 4.4 0;
    0.078 -0.99 -0.23 0.052;
    1.0 0.078 0 0];
B = [17 7;
    0.82 -3.2;
    0 0.046;
    0 0];

[eigv,eig] = eig(A)

xdot = @(t,x) A*x + B*[unit_impulse(t);0];

[t_i,x_i] = ode45(xdot,[0 3],zeros(4,1));

figure()
plot(t_i,x_i)
title("x for 3.a.i")

xdot = @(t,x) A*x + B*[0;unit_impulse(t)];

[t_ii,x_ii] = ode45(xdot,[0 3],zeros(4,1));

figure()
plot(t_ii,x_ii)
title("x for 3.a.ii")



Ac = [0 1;
    -7.29 -4.05];
Bc = [0 0 0 0;
    0 0 1 0];
Bu2 = B(:,2)*[10,0];


Aaug = [A, Bu2;
        Bc, Ac];
Baug = [B(:,1);0;0];
Aaug
Baug

xdot = @(t,x) Aaug*x + Baug*[unit_impulse(t)];
[t_aug,x_aug] = ode45(xdot,[0 3],zeros(6,1));
figure()
plot(t_aug,x_aug)
title("x_{aug} for 3.b")
legend(["x"+num2str((1:4)');("x_c"+num2str((1:2)'))])

%% Problem 5
clc, close all, clear all
M = .5;
m=.2;
b=.1;
I=.006;
g=9.8;
l=.3;

denom = I*(M+m)+M*m*l^2;

A = [0 1 0 0;
    0 -(I+m*l^2)*b/denom m^2*g*l^2/denom 0;
    0 0 0 1;
    0 -m*l*b/denom m*g*l*(M+m)/denom 0]

B = [0;
    (I+m*l^2)/denom;
    0;
    m*l/denom]

[eigv,eigs] = eig(A)
[Ua,Sa,Va] = svd(A)

phiAi = expm(A*.1)
phiAii = expm(A*1)

xdot = @(t,x) A*x + B*[hadamardStep(t)];
[t,x] = ode45(xdot,[0 1],zeros(4,1));
figure()
plot(t,x)
title("x for 5.e")
legend("$x$","$\dot{x}$","$\theta$","$\dot{\theta}$",'Interpreter','latex')

%% Problem 6
clc, clear all, close all
V1 = [1;-3];
V2 = [1;2];
L = [1 0; 0 -3];
V = [V1,V2];
A = V*L*V^-1
%% Problem 7
clc, clear all, close all
A = [-3 2; 0 1];
[eigv,eigs] = eig(A)
%% Problem 8
clc, close all, clear all
M = .5;
m=.2;
b=.1;
I=.006;
g=9.8;
l=.3;

denom = I*(M+m)+M*m*l^2;

A = [0 1 0 0;
    0 -(I+m*l^2)*b/denom m^2*g*l^2/denom 0;
    0 0 0 1;
    0 -m*l*b/denom m*g*l*(M+m)/denom 0];

B = [0;
    (I+m*l^2)/denom;
    0;
    m*l/denom];

[eigv,eigs] = eig(A)
syms t;
phi_z = expm(eigs*t)
phi_x = eigv*phi_z*eigv^-1;
phi_x = vpa(phi_x,5)

%% Helper functions
function [out] = unit_impulse(t)
    if(t<1 && t>0)
        out = 1.;
    else
        out = 0.;
    end
end
function [out] = hadamardStep(t)
    if(t>0)
        out = 1.;
    else
        out = 0.;
    end
end