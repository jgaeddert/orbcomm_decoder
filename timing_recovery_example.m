% demonstrate timing correction using very simple interpolator
clear all;
close all;

% options
k           = 2;            % samples/symbol
m           = 7;            % filter semi-length
dt          = 0.3;          % fractional sample timing offset
num_symbols = 1200;         % number of symbols to simulate

% generate data symbols (QPSK)
s = exp(1i*2*pi*randi(4,1,num_symbols)/4 + 1i*pi/4);

% square-root matched filters (note: these mimic a RRC filter, but simpler to compute)
r  = 1.057; % bandwidth expansion ratio for approximating a RRC filter
ht = sinc(([-k*m:k*m] + dt)*r/k) .* hamming(2*k*m+1)' * r; % transmit (w/ offset)
hr = sinc(([-k*m:k*m]     )*r/k) .* hamming(2*k*m+1)' * r; % receive (no offset)

% interpolate symbols using transmit matched filter
x = filter(ht,1,reshape([s(:) zeros(num_symbols,k-1)].',1,[]));

% TODO: apply channel effects (noise carrier offset, etc.)
v = x;

% apply matched filter and normalize signal level
y = filter(hr,1,v);
y = y / median(abs(y));

% interpolate matched filter output
tau = dt;                   % timing correction to be applied
th  = [-2, -1, 0, 1, 2];    % time vector for interpolating filter
w   = sqrt(hamming(5)');    % window function for interpolating filter
hi  = sinc(th - tau) .* w;  % interpolating filter coefficients
z   = filter(hi,1,y);       % interpolate matched filter output

% plot results
syms_mf     = y(1:k:end);   % symbols for matched filter output (no correction)
syms_interp = z(3:k:end);   % symbols after MF and interpolator (2 sample delay in interp)

plot(syms_mf,    'x','Color',[1 1 1]*0.8,...
     syms_interp,'x','Color',[0 0.5 .2]);
axis([-1 1 -1 1]*1.5);
grid on;
axis square;
xlabel('in-phase');
ylabel('quadrature phase');
legend('Before timing correction','After timing correction');


