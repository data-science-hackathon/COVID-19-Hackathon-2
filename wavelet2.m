function [wave,period,coi] = ...
	wavelet2(Y,dt,scale);
n1 = length(Y);
x(1:n1) = Y - mean(Y);
base2 = fix(log(n1)/log(2) + 0.4999);
x = [x,zeros(1,2^(base2+1)-n1)];
n = length(x);

k = [1:fix(n/2)];
k = k.*((2.*pi)/(n*dt));
k = [0., k, -k(fix((n-1)/2):-1:1)];

f = fft(x);
J1 = length(scale);

wave = zeros(J1,n);
wave = wave + i*wave;

for a1 = 1:J1
	[daughter,fourier_factor,coi,dofmin]=wave_bases('MORLET',k,scale(a1),-1);	
	wave(a1,:) = ifft(f.*daughter);
end

period = fourier_factor*scale;
coi = coi*dt*[1E-5,1:((n1+1)/2-1),fliplr((1:(n1/2-1))),1E-5];  % COI [Sec.3g]
wave = wave(:,1:n1);

return
