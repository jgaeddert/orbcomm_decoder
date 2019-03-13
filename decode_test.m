clear all;
close all;

load('data/1552071892p6.mat');

% compute decimation rate
fsym = 4800;        % symbol rate
k    = 4;           % desired samples per symbol
D    = fs/(k*fsym); % decimation rate

% approximate (relative) center frequencies
fc = [ -0.17066590625;
        0.19555240000;];
num_signals = length(fc);
id = 1; % signal id to process

% half-band decimator filter
m = 17; % filte semi-length
h = 0.5 * sinc([-2*m:2*m]*0.5) .* hamming(4*m+1)';

% down-sample
num_samples = length(samples);

% approximate matched filter (this mimics a RRC filter)
hp = sinc([-k*m:k*m]*1.09/k) .* hamming(2*k*m+1)';

% mix down by coarse frequency offset
y = samples .* exp(-1i*2*pi*fc(id)*[0:(num_samples-1)]);

% decimate by iteratively cutting sample rate by half
for i = 1:log2(D),
    y = filter(h,1,y)(1:2:end);
end;
n = length(y);

% apply matched filter and normalize level
r = filter(hp,1,y);
r = r / median(abs(r));

% recover carrier using 4th-order costas loop
v = zeros(1,n);
phi = 0;
dphi = -0.2;    % initial guess for carrier offset
alpha = 0.05;   % loop filter bandwidth
beta  = 0.5*alpha^2;
for i=1:n,
    % mix down
    v(i) = r(i) * exp(-1i*phi);

    % compute carrier phase error and adjust phase/frequency accordingly
    phase_error = imag(v(i).^4);
    phi  += alpha * phase_error;
    dphi +=  beta * phase_error;
    phi  += dphi;
end;

% recover timing
s = []; % output symbols
for i=1:n,
end;

return

% estimate residual carrier offset
%dphi_hat_0 = -arg(sum(r0(1:(n-1)) .* conj(r0(2:n))));
%r0 = r0 .* exp(-1i*dphi_hat_0*[0:(n-1)]);

% dumb spectrum estimate
nfft = 2400;
i0 = 0;
Y0 = Y1 = zeros(1,nfft);
R0 = R1 = zeros(1,nfft);
num_transforms = 0;
while i0 + nfft < n,
    Y0 = Y0 + abs(fft(y0([1:nfft]+i0).*hamming(nfft)',nfft)).^2;
    Y1 = Y1 + abs(fft(y1([1:nfft]+i0).*hamming(nfft)',nfft)).^2;
    R0 = R0 + abs(fft(r0([1:nfft]+i0).^4.*hamming(nfft)',nfft)).^2;
    R1 = R1 + abs(fft(r1([1:nfft]+i0).^4.*hamming(nfft)',nfft)).^2;
    num_transforms = num_transforms + 1;
    i0 = i0 + round(nfft/4);
end;

Y0 = 10*log10(fftshift(Y0)/num_transforms);
Y1 = 10*log10(fftshift(Y1)/num_transforms);
R0 = 10*log10(fftshift(R0)/num_transforms);
R1 = 10*log10(fftshift(R1)/num_transforms);
f = ([0:(nfft-1)]/nfft - 0.5); % * fsym * k;
plot(f,R0,'-x'); %f*1e-3,R1);
xlabel('Frequency [kHz]');
ylabel('PSD [dB]');
legend('signal 0','signal 1');
axis([-0.02 0.02 110 140]);
grid on;

