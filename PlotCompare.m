clear; clc;
kgrid = load('kgrid.txt');
exp_grid = load('expgrid.txt');
exp_grid2 = (0:1:44)';
v1 = load('vlast.txt');
v2 = load('vlast2.txt');
c1 = load('clast.txt');
c2 = load('clast2.txt');
k1 = load('klast.txt');
k2 = load('klast2.txt');
nm1 = load('nmlast.txt');
nm2 = load('nmlast2.txt');
nf1 = load('nflast.txt');
nf2 = load('nflast2.txt');


figure
plot(exp_grid, nm1(10,:))
hold on
plot(exp_grid2, nm2(10,:), 'r--')

% plot(kgrid, nf1(:,end))
% hold on
% plot(kgrid, nf2(:,end),'r--')
