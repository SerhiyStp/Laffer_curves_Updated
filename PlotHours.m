close all; clear;

lm = load('laborm.txt');
lf = load('laborf.txt');

figure('Position', [60,60,900,600])
plot(lm(:,1), lm(:,2), 'r')
hold on
plot(lf(:,1), lf(:,2), 'b--')
legend('men','women')