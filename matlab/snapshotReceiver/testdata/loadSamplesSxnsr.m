% s = loadSamplesSxnsr(filename, fs, T [, fmt ])
%
% Loads samples from an IfEN SX-NSR device
% (https://www.ion.org/gnss/upload/files/956_SX.NSR.Flyer_Aug2013_Letter.pdf)
%
% Parameters:
% filename... input filename to be parsed
% fs......... sampling rate in the file [Hz]
% T.......... length of the snapshot to be extracted [s]
% fmt........ (optional) format of the samples to be read (default: 'bit2')
%
% Returns:
% s.......... samples
%
function s = loadSamplesSxnsr(filename, fs, T, fmt)

if ~exist('fmt', 'var')
    fmt = 'bit2';
end

fid = fopen(filename, 'r');
if fid < 0
    error('Could not open input file "%s"', filename);
end

Nsamples = ceil(fs*T);
s = fread(fid, [1 Nsamples], fmt);

fclose(fid);
