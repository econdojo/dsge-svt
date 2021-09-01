function [fh] = NBERbc(time, y, LS, LW, LC)
%
% NBERbc creates shaded areas for the NBER recessions for plot.
%
% - Usage - 
% [fh] = NBERbc(time, y, LS, LW, LC)
% 
% where
%   fh is a vector of the handles for series (optional),
%   time is horizontal axis,
%   y is a vector or matrix of series to be plotted,
%   LS is a cell array specifying the style of line,
%   LW is a vector specifying line-width,
%   LC is a cell array of line colors.
%
%   LS, LW, and LC must be in order of series in y.
%
%   The handles in fh are those returned by plot.
%   If you use NBERbc without output arguments, it will simply
%   add the shaded areas onto the plot.
%
% Example:
%   time = [1945:0.25:2005.75]  % quarterly series from 1945
%   y    = randn(length(time), 1)
%   
%   NBERbc(time, y, {'-'}, 2, {'b'}) will plot the series y from 1945:Q1
%   to 2005:Q4 with a blue line whose width is 2, together with
%   shaded areas of the NBER recessions.
%   
% - Note -
% By default, this program takes care of NBER recessions after
% the World War II. For earlier recessions, 'peak' and 'trough' should be modified.
%

% -------------------------------------------------------------------------
% Written by
% Munechika Katayama [mkatayama@lsu.edu]
% Department of Economics, Louisiana State University
% -------------------------------------------------------------------------
% Ver. 1.0 10/06/2006
% Ver. 1.1 12/10/2008 -- Added the peak in Dec 2007
% Ver. 1.2 09/20/2010 -- Added the trough in June 2009
%
% Future Implementation 
%    - Incorporate serial date number in Matlab


% ----------- Fill areas for the NBER recessions -----------

% -- Peaks and troughs of business cycle --
% -- Convention here is that 1945.0 is January of 1945. So, for example,
% -- 1949+9/12 means October of 1949.
peak   = [1945+1/12; 1948+10/12; 1953+6/12; 1957+7/12; 1960+3/12; 1969+11/12; 1973+10/12; 1980;      1981+6/12;  1990+6/12; 2001+2/12;  2007+11/12];
trough = [1945+9/12; 1949+9/12;  1954+4/12; 1958+3/12; 1961+1/12; 1970+10/12; 1975+2/12;  1980+6/12; 1982+10/12; 1991+2/12;2001+10/12;  2009+5/12];
nbc    = length(peak);

NBERx = reshape([peak peak trough trough]', nbc*4, 1);

NBERy = [-1000; 1000];
% This specifies the range of Y-axis. A huge range is reserved in order to
% handle any type of data. This range will be shrunk later, according to
% max and min of y1.

NBERy = [NBERy; flipud(NBERy)];
NBERy = kron(ones(nbc,1), NBERy);

% -- Fill the area --
fill(NBERx', NBERy', [0.8 0.8 0.8], 'LineStyle', 'none')

ymin = min(min(y)) - 0.2*abs(min(min(y)));
ymax = max(max(y)) + 0.2*abs(max(max(y)));

if ymin == ymax
    ymax = ymin+1;
end

set(gca, 'YLim', [ymin ymax], 'YTickMode', 'auto')
% This shrinks the range of the Y-axis.
% New Y-axis will be 10% below (above) the min (max) of y.

% ----------- Plotting Data -----------
hold on
[fh] = plot(time, y);

if nargin > 2
    for i = 1:length(fh)
        set(fh(i), 'LineStyle', LS{i} , 'LineWidth', LW(i), 'Color', LC{i});
    end
end

set(gca, 'XLim', [min(time) max(time)], 'XTickMode', 'auto', 'Layer', 'top')