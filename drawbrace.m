% Copyright (c) 2012, P�l N�verlid S�vik
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%%
function h = drawbrace(start, stop, width, varargin)
% DRAWBRACE   Draw curly brace on current figure
%    DRAWBRACE([X1, Y1], [X2, Y2]) draws a curly brace from the point [X1,
%    Y1] to the point [X2, Y2]
%
%    DRAWBRACE([X1,Y1], [X2,Y2], W) draws a brace using the brace width W.
%
%    DRAWBRACE(..., 'Param1', 'Value1', 'Param2', 'Value2') draws a brace 
%    using the LineSeries property values specified by 'Param1'/'Value1',
%    'Param2'/'Value2', ... .
%
%    H = DRAWBRACE(...) returns the LineSeries handle to the brace.
%
%    Example:
%       H = drawbrace([0 0], [1 1], 20, 'Color', 'k')
    
    % Get axis size
    pos = get(gca, 'Position');
    opos = get(gcf, 'Position');
    ylims = ylim;
    xlims = xlim;
    
    % Take logarithmic scale into account
    isxlog = strcmp(get(gca, 'XScale'), 'log');
    isylog = strcmp(get(gca, 'YScale'), 'log');
    if isxlog
        start(1) = log(start(1));
        stop(1) = log(stop(1));
        xlims = log(xlims);
    end
    if isylog
        start(2) = log(start(2));
        stop(2) = log(stop(2));
        ylims = log(ylims);
    end
    
    % Transform from axis to screen coordinates
    xscale = pos(3) * opos(3) / diff(xlims);
    yscale = pos(4) * opos(4) / diff(ylims);
    start = (start - [xlims(1) ylims(1)]) .* [xscale yscale];
    stop = (stop - [xlims(1) ylims(1)]) .* [xscale yscale];

    % Find standard width
    if nargin == 2
        width = norm(stop - start)/10; 
    end
    
    % Find brace points
    th = atan2(stop(2)-start(2), stop(1)-start(1));
    c1 = start + width*[cos(th) sin(th)];
    c2 = 0.5*(start+stop) + 2*width*[-sin(th) cos(th)] - width*[cos(th) sin(th)];
    c3 = 0.5*(start+stop) + 2*width*[-sin(th) cos(th)] + width*[cos(th) sin(th)];
    c4 = stop - width*[cos(th) sin(th)];
    
    % Assemble brace coordinates
    q = linspace(0+th, pi/2+th, 50)';
    t = flipud(q);
    part1x = width*cos(t+pi/2) + c1(1);
    part1y = width*sin(t+pi/2) + c1(2);
    part2x = width*cos(q-pi/2) + c2(1);
    part2y = width*sin(q-pi/2) + c2(2);
    part3x = width*cos(q+pi) + c3(1);
    part3y = width*sin(q+pi) + c3(2);
    part4x = width*cos(t) + c4(1);
    part4y = width*sin(t) + c4(2);
    x = [part1x; part2x; part3x; part4x];
    y = [part1y; part2y; part3y; part4y];
    
    % Transform back to axis coordinates
    x = x / xscale + xlims(1);
    y = y / yscale + ylims(1);
    if isxlog
        x = exp(x); end
    if isylog
        y = exp(y); end
    
    % Plot brace
    h = line(x, y);
    for i = 1:2:numel(varargin)
        set(h, varargin{i}, varargin{i+1}); end

end
