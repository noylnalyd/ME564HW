%% Problem 1
clc, clear, close all;

As = cell(3,1);
As{1} = [2,3;3,1];
tht = pi/8;
As{2} = [cos(tht),-sin(tht);sin(tht),cos(tht)];
As{3} = [1,-1;2,-1];

figure();
hold on;
colors = ['r' 'g' 'b' 'c' 'm' 'y']
nvects = length(colors);
vects = zeros(2,nvects);
for i=1:nvects
    tht = 2*pi*(i-1)/(nvects);
    vects(:,i) = [cos(tht);sin(tht)];
end
for mat=1:3
    subplot(1,3,mat);
    hold on;
    title(['Plot of A_(' mat ')'])
    A = As{mat};
    vects(:,i)
    for i=1:nvects
        quiver(vects(1,i),vects(2,i),colors(i),'LineWidth',2);
    end
    Av = A*vects;
    for i=1:nvects
        quiver(Av(1,i),Av(2,i),[':' colors(i)],'LineWidth',2);
    end
end


%% Problem 5
clc, clear;
A = [2,1;1,-8];

% a.
a = zeros(2);
for i=0:4
    a = a + A^i / factorial(i);
end
a

b = zeros(2);
for i=0:4
    b = b + (-1)^i * A^i / factorial(i);
end
b = b^-1;
b

l1 = -3 - sqrt(36+4*17)/2;
l2 = -3 + sqrt(36+4*17)/2;

dif = exp(l2) - exp(l1);
tot = exp(l2) + exp(l1);
c1 = dif/(l2-l1);
c0 = exp(l1)-c1*l1;

c = [c0+2*c1,c1;c1,c0-8*c1]

[V,D,W] = eig(A);
for i=1:2
    D(i,i) = exp(D(i,i));
end
d = (W'^-1)*D*(W')

e = expm(A)