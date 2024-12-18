AE = 1;
t1 = 1-(1-0.0765)/(1+0.0765);
t2 = 1-(1-0.0145)/(1+0.0145);
c_sigm = 14;
S = 2.323967451*AE;

Cbar = (t1 - t2)*(S + log(2)/c_sigm);

ny = 100;
yy = linspace(0.01, 3, ny)';

for i = 1 : ny
    tss(i,1) = t2*yy(i) - (t1-t2)/c_sigm*log(abs(1 + exp(-c_sigm*(yy(i)-S)))) + Cbar;
end

figure
plot(yy, tss)