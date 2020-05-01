% [eph, doy, alpha, beta] = parseRinexNavFile(filename, prn, tow)
%
% Parses a RINEX navigation file and returns the best matching ephemerides
% for the satellites in prn input vector.
%
% Parameters:
% filename...... RINEX navigation file (xxxxxxx.yyN)
% prn........... Satellite numbers to be extracted
% tow........... Current TOW for which the orbits are to be extracted [s]
%
% Returns:
% eph........... Vector of structs with the ephemerides
% doy........... Day of year
% alpha,beta.... Klobuchar model parameters (if not available empty)
%
function [eph, doy, alpha, beta] = parseRinexNavFile(filename, prn, tow)

if ~exist(filename, 'file')
    error('Navigation File doesnt exist!');
end

doy = 0;

fid = fopen(filename, 'r');
if fid < 0
    error('Could not open Navigation File');
end

line = fgetl(fid);
if ~strcmp(line(21),'N')
    % Doesn't seem to be a Navigation file
    error('This is not a Navigation File');
end

rinexVersion = str2double(line(1:9));

% extract the ionospheric correction parameters if they're present
if rinexVersion < 3 && nargout > 2
    [alpha,beta] = parseIono2(fid);
elseif rinexVersion >= 3
    [alpha,beta] = parseIono3(fid);
else
    alpha = [];
    beta = [];
end

findFirstEntry(fid);

N = length(prn);
eph = struct('satType', cell(1,N), 'PRN', cell(1,N), ...
    'doy', cell(1,N), 'a_f0', cell(1,N), ...
    'a_f1', cell(1,N), 'a_f2', cell(1,N), 'IODE_sf2', cell(1,N), ...
    'C_rs', cell(1,N), 'deltan', cell(1,N), 'M_0', cell(1,N), ...
    'C_uc', cell(1,N), 'e', cell(1,N), 'C_us', cell(1,N), ...
    'rootA', cell(1,N), 't_oe', cell(1,N), 'C_ic', cell(1,N), ...
    'Omega_0', cell(1,N), 'C_is', cell(1,N), 'i_0', cell(1,N), ...
    'C_rc', cell(1,N), 'omega', cell(1,N), 'omegaDot', cell(1,N), ...
    'iDot', cell(1,N), 'weekNumber', cell(1,N), 'health', cell(1,N), ...
    'T_GD', cell(1,N), 'time', cell(1,N), 't_oc', cell(1,N));

while ~feof(fid)
    p = parseEntry(fid, rinexVersion);
    if ~isempty(p) && any(prn==p.PRN) && p.time > tow
        ix = find(prn == p.PRN);
        if isempty(eph(ix).time) || p.time < eph(ix).time
            eph(ix) = p;
            doy = eph(ix).doy;
        end
    end
end

fclose(fid);

end

% [alpha,beta] = parseIono2(fid)
% parses the navigation to find the ionospheric alpha and beta parameters
% probably present. If the end of the header has been reached without
% finding them, empty matrices are returned (RINEX 2).
function [alpha,beta] = parseIono2(fid)
line = fgetl(fid);
alpha = [];
beta  = [];
while ~feof(fid) && isempty(strfind(line,'END OF HEADER')) && ...
        isempty(strfind(line,'ION ALPHA')) && ...
        isempty(strfind(line,'ION BETA'))
    line = fgetl(fid);
end

if ~isempty(strfind(line,'END OF HEADER'))
    % no iono corrections found
    return;
elseif ~isempty(strfind(line,'ION ALPHA'))
    alpha = sscanf(line(1:60),'%f %f %f %f');
    if length(alpha)==1
        alpha = cell2mat(textscan(line(1:60),'%f%f%f%f'));
    end
elseif ~isempty(strfind(line,'ION BETA'))
    beta = sscanf(line(1:60),'%f %f %f %f');
    if length(beta)==1
        beta = cell2mat(textscan(line(1:60),'%f%f%f%f'));
    end
end

line = fgetl(fid);
while ~feof(fid) && isempty(strfind(line,'END OF HEADER')) && ...
        isempty(strfind(line,'ION ALPHA')) && ...
        isempty(strfind(line,'ION BETA'))
    line = fgetl(fid);
end

if ~isempty(strfind(line,'END OF HEADER'))
    % no iono corrections found
    return;
elseif ~isempty(strfind(line,'ION ALPHA'))
    alpha = sscanf(line(1:60),'%f %f %f %f');
    if length(alpha)==1
        alpha = cell2mat(textscan(line(1:60),'%f%f%f%f'));
    end
elseif ~isempty(strfind(line,'ION BETA'))
    beta = sscanf(line(1:60),'%f %f %f %f');
    if length(beta)==1
        beta = cell2mat(textscan(line(1:60),'%f%f%f%f'));
    end
end

end

% [alpha,beta] = parseIono3(fid)
% parses the navigation to find the ionospheric alpha and beta parameters
% probably present. If the end of the header has been reached without
% finding them, empty matrices are returned (RINEX 3).
function [alpha,beta] = parseIono3(fid)
line = fgetl(fid);
alpha = [];
beta  = [];
while ~feof(fid) && isempty(strfind(line,'END OF HEADER'))
    if ~isempty(strfind(line,'IONOSPHERIC CORR'))
        if strcmp(line(1:4),'GPSA')
            alpha = [ str2double(line(5+(1:12))) str2double(line(17+(1:12))) ...
                str2double(line(29+(1:12))) str2double(line(41+(1:12))) ];
        elseif strcmp(line(1:4),'GPSB')
            beta  = [ str2double(line(5+(1:12))) str2double(line(17+(1:12))) ...
                str2double(line(29+(1:12))) str2double(line(41+(1:12))) ];
        end
    end
    line = fgetl(fid);
end
fseek(fid, 0, 'bof');
end


% findFirstEntry(fid)
% finds the first data-entry inside the file opened with pointer 'fid'
function findFirstEntry(fid)
while ~feof(fid)
    line = fgetl(fid);
    if ~isempty(strfind(line,'END OF HEADER'))
        return;
    end
end
end

% eph = parseEntry (line, fid)
% parses an entry in the Navigation file and stores all values in a struct
function eph = parseEntry(fid, rinexVersion)

if feof(fid)
    % no more data sets
    eph = [];
    return;
end
line = fgetl(fid);
while sum(abs(cell2mat(textscan(line(2:end),'%f%f%f%f%f%f%f%f%f%f'))))==0 && ~feof(fid)
    line = fgetl(fid);
end

eph = [];
if rinexVersion>=3
    s = cell2mat(textscan(line(2:23), '%f%f%f%f%f%f%f'));
    eph.satType = line(1);
else
    s = cell2mat(textscan(line(1:22), '%f%f%f%f%f%f%f'));
    eph.satType = 'G';
end
if size(s,2)~=7 || isnan(eph.satType)
    % Uncompatible data set
    eph = [];
    return;
end
eph.PRN = s(1);
eph.doy = datenum(s(2),s(3),s(4)) - datenum(s(2),1,1) + 1;

nx = 3+(rinexVersion>=3); % number of whitespaces in front of the data

line = regexprep(line,'[Dd]','e'); % str2double doesn't like D or d as separator for the exponent

% Sat. Clock
eph.a_f0 = str2double(line(nx+(20:38))); %[s]
eph.a_f1 = str2double(line(nx+(39:57))); %[s/s]
eph.a_f2 = str2double(line(nx+(58:76))); %[s/s^2]

% LINE 2
if feof(fid)
    % Incomplete data set
    eph = [];
    return;
end
line = fgetl(fid);
line = regexprep(line,'[Dd]','e');

eph.IODE_sf2 = str2double(line(nx+(1:19)));
eph.C_rs     = str2double(line(nx+(20:38)));  %[m]
eph.deltan   = str2double(line(nx+(39:57)));  %[rad/s]
eph.M_0      = str2double(line(nx+(58:76)));  %[rad]

% LINE 3
if feof(fid)
    % Incomplete data set
    eph = [];
    return;
end
line = fgetl(fid);
line = regexprep(line,'[Dd]','e');

eph.C_uc  = str2double(line(nx+(1:19)));   %[rad]
eph.e     = str2double(line(nx+(20:38)));  %[]
eph.C_us  = str2double(line(nx+(39:57)));  %[rad]
eph.rootA = str2double(line(nx+(58:76)));  %[m^.5]

% LINE 4
if feof(fid)
    % Incomplete data set
    eph = [];
    return;
end
line = fgetl(fid);
line = regexprep(line,'[Dd]','e');

eph.t_oe    = str2double(line(nx+(1:19)));   %[s (of GPS week)]
eph.C_ic    = str2double(line(nx+(20:38)));  %[rad]
eph.Omega_0 = str2double(line(nx+(39:57)));  %[rad]
eph.C_is    = str2double(line(nx+(58:76)));  %[rad]

% LINE 5
if feof(fid)
    % Incomplete data set
    eph = [];
    return;
end
line = fgetl(fid);
line = regexprep(line,'[Dd]','e');

eph.i_0      = str2double(line(nx+(1:19)));   %[rad]
eph.C_rc     = str2double(line(nx+(20:38)));  %[m]
eph.omega    = str2double(line(nx+(39:57)));  %[rad]
eph.omegaDot = str2double(line(nx+(58:76)));  %[rad/s]

% LINE 6
if feof(fid)
    % Incomplete data set
    eph = [];
    return;
end
line = fgetl(fid);
line = regexprep(line,'[Dd]','e');

eph.iDot       = str2double(line(nx+(1:19)));  %[rad/s]
eph.weekNumber = str2double(line(nx+(39:57))); %[]

% LINE 7
if feof(fid)
    % Incomplete data set
    eph = [];
    return;
end
line = fgetl(fid);
line = regexprep(line,'[Dd]','e');

eph.health = str2double(line(nx+(20:38)));
if strcmpi(eph.satType, 'G')
    eph.T_GD = str2double(line(nx+(39:57)));
else
    eph.T_GD = 0;%str2double(line(nx+(58:76)));
end

% LINE 8
if feof(fid)
    % Incomplete data set
    eph = [];
    return;
end
fgetl(fid);

% eph.time = s(1); %[s (of GPS week)]
eph.time = eph.t_oe; %[s (of GPS week)]
eph.t_oc = eph.t_oe;

end
