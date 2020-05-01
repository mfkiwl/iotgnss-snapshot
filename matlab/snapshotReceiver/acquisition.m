% [prn, codePhase, doppler, cno] = acquisition(s, fs, fi [, fDmax ])
%
% Acquires the satellite signals contained in in the input signal using a
% decimation and FFT approach.
%
% Parameters:
% s............. Input signal samples
% fs............ Sampling rate [Hz]
% fi............ Intermediate frequency [Hz]
% fDmax......... (optional) Maximum Doppler frequency to be searched [Hz]
%                (default: 4300)
%
% Returns:
% prn........... Found satellite numbers
% codePhase..... Corresponding code phase measurements [s]
% Tambg......... Ambiguity of the code phase measurement [s]
% doppler....... Corresponding Doppler measurements [Hz]
% cno........... Corresponding C/N0 estimates [dB-Hz]
%
function [prn, codePhase, Tambg, doppler, cno] = acquisition(s, fs, fi, fDmax)

addpath('codes');

if ~exist('fDmax', 'var')
    fDmax = 4.3e3;
end

% remove DC component (due to quantization e.g.)
s = s - mean(s);

% carrier frequeny [Hz]
f0 = 1.57542e9;

% generate the spreading codes
fc = 1.023e6; % chipping rate [chips/s]
c = codeGpsL1ca();
[numSv, L] = size(c);

% transform the spreading code to frequency domain using FFT (to speed up
% the correlation later on)
cF = conj(fft(c'));

% timing of the samples in s
ts = (0:length(s)-1) / fs;

T = ts(end);
fDstep = round(1/T);
fD = 0;

acqHits = struct('prn', cell(1,numSv), 'fD', cell(1,numSv), ...
    'codePhase', cell(1,numSv), 'peakRatio', cell(1,numSv));

sigma = NaN;

while fD < abs(fDmax)
    % carrier wipe-off and down-sampling of the input signal
    sLowRate = downSampleMix(s, fs, fi + fD, fc*(1+fD/f0), L);
    
    % execute the correlation in frequency domain using the FFT
    sLowRateF = fft(sLowRate) * ones(1, numSv);
    corrMag = abs(ifft(sLowRateF .* cF));
    
    if isnan(sigma)
        % noise floor estimate (assuming a Rayleigh distribution)
        % (but just do it once, because we can assume this value stays
        % constant over multiple iterations)
        sigma = median(corrMag(:))/sqrt(2*log(2));
    end
    
    % detection threshold (prob. false alarm = 1e-8)
    thresh = sigma * sqrt(-2*log(1e-8));
    
    [maxCorrMag, ixMax] = max(corrMag, [], 1);
    prnHits = find(maxCorrMag > thresh);
    for n=1:length(prnHits)
        prn_n = prnHits(n);
        
        % better codePhase by 3-point interpolation
        [mag,codePhase] = interpolateMax3(corrMag(:,prn_n), ixMax(prn_n));
        
        acqHitPeakRatio = mag / sigma;
        
        if isempty(acqHits(prn_n).peakRatio) || acqHitPeakRatio > acqHits(prn_n).peakRatio
            acqHits(prn_n).prn = prn_n;
            acqHits(prn_n).fD = fD;
            acqHits(prn_n).codePhase = codePhase;
            acqHits(prn_n).peakRatio = acqHitPeakRatio;
        end
    end
    
    % choose next Doppler frequency
    if fD > 0
        fD = -fD;
    else
        fD = -fD + fDstep;
    end
end

prn = find(arrayfun(@(x) ~isempty(x.prn), acqHits));
codePhase = 1023-arrayfun(@(x) x.codePhase, acqHits(prn)) / fc;
Tambg = ones(size(prn)) * L / fc;
doppler = arrayfun(@(x) x.fD, acqHits(prn));
cno = arrayfun(@(x) 10*log10(x.peakRatio^2/(2*T)), acqHits(prn));

end

% sLowRate = downSampleMix(s, fs, f, fc, L)
%
% Mixes input signal s to baseband with frequency f and at the same time
% reduces the sampling rate to the chipping rate fc. The chips are
% pre-integrated.
%
% Parameters:
% s.......... input signal samples
% fs......... sampling rate of the input signal [Hz]
% f.......... frequency to mix down to baseband [Hz]
% fc......... chipping rate [chips/s]
% L.......... length of the spreading code [chips]
%
% Returns:
% sLowRate... input signal reduced to rate fc
%
function sLowRate = downSampleMix(s, fs, f, fc, L)

sLowRate = zeros(L, 1);

carrPhase = 0;
carrPhaseIncr = 2*pi*f/fs;
codePhase = 0;
codePhaseIncr = fc / fs;

for n=1:length(s)
    ixLowRate = mod(floor(codePhase), L) + 1;
    sLowRate(ixLowRate) = sLowRate(ixLowRate) + s(n) * exp(1i*carrPhase);
    
    carrPhase = carrPhase + carrPhaseIncr;
    codePhase = codePhase + codePhaseIncr;
end

end

% [yi,xi] = interpolateMax3(y, xMax)
%
% Uses a 3-point interpolation to find the maximum of the input vector.
% This is using a polynom of order 2, i.e. we have
% [ y_n-1 ]   [ 1 -1 1 ]   [ a ]
% [  y_n  ] = [ 1  0 0 ] * [ b ]
% [ y_n+1 ]   [ 1  1 1 ]   [ c ]
%
% The maximum is at -b/(2*c).
%
% Parameters:
% y........... input vector
% xMax........ xMax location of maximum in vector y
%
% Returns:
% yi.......... interpolated maximum value
% xi.......... offset of the interpolated maximum (offset into y)
%
function [yi,xi] = interpolateMax3(y, xMax)

yi = y(xMax);
xi = xMax;

if 1 < xMax && xMax < length(y)
    a = yi;
    b = .5*(y(xMax+1) - y(xMax-1));
    c = .5*(y(xMax-1)+y(xMax+1)) - yi;
    xi = -b / (2*c);
    yi = a + xi*b + xi^2*c;
    xi = xi + xMax;
end

end