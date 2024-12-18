Cap = 2.3;
tss1 = 0.0765;
tss2 = 0.0145;

eps = 0.01;

x = linspace(0.001, 3.5, 1000);
% x = linspace(-2.5, 2.5, 1000);


c = 20;
f = 1./(1 + exp(c*(x-Cap)));
ff = tss2 + (tss1 - tss2)*f;

figure
plot(x, ff)
xlabel('AE')
ylabel('tss')