wavewa;
wx = wave;
sx = s0;
waveca;
wy = wave;
sy = s0;

wxy = wx.*conj(wy);
pxy = (abs(wxy)).^2;
px = (abs(wx)).^2;
py = (abs(wy)).^2;
rxy = abs(wxy)./(sqrt(abs(wx)).*sqrt(abs(wy)));
phixy = atan(imag(wxy)./real(wxy));

n = length(s);

global_r = (sum(phixy')/n);
Cdelta = 0.776;
r_avg = (scales')*(ones(1,n));
r_avg = phixy ./ r_avg;
r_avg = dj*dt/Cdelta*sum(r_avg);

subplot('position',[0.1 0.75 0.65 0.2])
%subplot(2,1,1);
plot(t,sx);
hold on
plot(t,sy);
hold off
legend('WA','CA');
set(gca,'XLim',xlim(:))
xlabel('Time (day)')
ylabel('Count')
%xline(22,'--');%NY 21, CA 18, WA 22
title('a) COVID-19 Daily New Cases')

Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
subplot('position',[0.1 0.37 0.65 0.28])
imagesc(t,log2(period),phixy);
colorbar;
caxis([-pi/2,pi/2])
xlabel('Time (day)')
ylabel('Log_2(Period (day))')
title('b) Phase between WA and CA')

subplot('position',[0.77 0.37 0.2 0.28])
plot(global_r,log2(period));
xlabel('Wavelet Phase')
title('c) Global Wavelet Phase')
set(gca, 'YLim',log2([min(period),max(period)]), 'YDir','reverse');

subplot('position',[0.1 0.07 0.65 0.2])
plot(t,r_avg)
set(gca,'XLim',xlim(:))
set(gca, 'YLim',[-1.5, 1.5])
yline(0.0,'r--');
xlabel('Time (day)')
ylabel('Avg Phase')
title('d) Scale-average Wavelet Phase')
