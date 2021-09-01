function [ah, fh1, fh2] = NBERbcYY(time, y1, y2, LS, LW, LC)
%
% NBERbcYY creates shaded areas for the NBER recessions for plotyy.
%
% - Usage - 
% [ah, fh1, fh2] = NBERbcYY(time, y1, y2, LS, LW, LC)
% 
% where
%   ah is the axis handle,
%   fh1 and fh2 are the figure handles for two series,
%   time is horizontal axis for both series,
%   y1 (left scale) and y2 (right scale) are series to be plotted,
%   LS is a cell array specifying the style of line,
%   LW is a vector specifying line-width,
%   LC is a cell array of line colors.
%
%   LS, LW, and LC must be in order of series in y1 and then y2.
%
%   The handles (ah, fh1, and fh2) are those returned by plotyy.
%   If you use NBERbcYY without output arguments, it will simply
%   add the shaded areas onto the plotyy.
%
% Example:
%   time = [1945:0.25:2005.75]  % quarterly series from 1945
%   y1   = randn(length(time), 2)
%   y2   = rand(length(time, 1)
%   LS   = {'-'; '--'; ':'}
%   LW   = [2 1 2]
%   LC   = {'b'; 'm'; 'r'}
%     
%   The 1st series in y1 will be a blue-solid line (line-width = 2)
%   The 2nd series in y1 will be a pink-dashed line (line-width = 1)
%   y2 will be plotted using right-Y-axis with red-dotted line and
%   line-width will be 2.
%
% - Note -
% By default, this program takes care of NBER recessions after
% the World War II. For earlier recessions, 'peak' and 'trough' should be modified.
%

% -------------------------------------------------------------------------
% Written by
% Munechika Katayama [m1kataya@ucsd.edu]
% Department of Economics, University of Carifornia, San Diego
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

ymin = min(min(y1) - 0.1*abs(min(y1)));
ymax = max(max(y1) + 0.1*abs(max(y1)));

set(gca, 'YLim', [ymin ymax], 'YTickMode', 'auto')
% This shrinks the range of the left-Y-axis.
% New left-Y-axis will be 10% below (above) the min (max) of y1.
% This can be changed either by using ah(1) or by modifying directly
% this file.

% ----------- Plotting Data -----------
hold on
[ah, fh1, fh2] = plotyy(time, y1, time, y2);
fh = [fh1; fh2];
for i = 1:length(fh)
    set(fh(i), 'LineStyle', LS{i} , 'LineWidth', LW(i), 'Color', LC{i});
end

set(ah, 'YColor', 'k')      % Setting the color of left-Y-axis and right-Y-axis to be black
set(ah, 'XLim', [min(time) max(time)], 'XTickMode', 'auto', 'Layer', 'top')

hold off