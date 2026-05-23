clear; clc;

z1 = 1e-3;
z2 = 1e-3;
z3 = -1e-3;

z = [z1; z2; z3];

x1 = 108.5e-3; y1 = 123e-3;
x2 = -109e-3; y2 = 123e-3;
x3 = 63.5e-3; y3 = -174e-3;

C = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1);

ax1 = (x3-x2)/C; ax2 = -(x3-x1)/C; ax3 = (x2-x1)/C;
ay1 = (y3-y2)/C; ay2 = -(y3-y1)/C; ay3 = (y2-y1)/C;
b1 = 1 + (x1*(y3-y2)+y1*(x2-x3))/C;
b2 = (-x1*(y3-y1)+y1*(x3-x1))/C;
b3 = (x1*(y2-y1)-y1*(x2-x1))/C;

M = [b1, b2, b3;
     ax1, ax2, ax3;
     ay1, ay2, ay3];

res = M * z;

%%
pos = [x1, y1, z1;
       x2, y2, z2;
       x3, y3, z3];
E = pos \ -ones(3,1);
res2 = [-1/E(3);
        -E(2)/E(3);
        E(1)/E(3)];

res
res2
err = norm(res - res2)

%% MS 水平
