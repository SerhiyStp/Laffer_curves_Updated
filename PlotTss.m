tss = load('tss_test.txt');
figure('Position', [60,60,1000,600])
subplot(1,2,1)
plot(tss(:,1), tss(:,2))

subplot(1,2,2)
plot(tss(:,1), tss(:,3))