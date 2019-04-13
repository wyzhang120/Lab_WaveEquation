classdef WAVEPATH
    properties
        sx; sz; gx; gz; dx; tStart; tEnd;
        wp;
    end
    
    methods
        function obj=WAVEPATH(sx, sz, gx, gz, dx, tstart, tend, wp)
            obj.sx = sx;  obj.sz = sz;
            obj.gx = gx;  obj.gz = gz;
            obj.dx = dx;
            obj.tStart = tstart;
            obj.tEnd = tend;
            obj.wp = wp;
        end
        
        
        function hdl=plotWP(obj, varargin)
            if nargin == 2
                percentiles = varargin{1};
                clim=prctile(reshape(obj.wp, [numel(obj.wp) 1]), percentiles);
            else
                clim=prctile(reshape(obj.wp, [numel(obj.wp) 1]), [5, 95]);
            end
            hdl = figure;
            [nz, nx] = size(obj.wp); 
            x = (0:nx-1)*obj.dx; z = (0:nz-1)*obj.dx;
            imagesc(x,z, obj.wp);colormap(gray);caxis(clim);axis equal;axis tight;
            title(sprintf('wave path (sz=%.2f m, gz = %.2f m) \n time window (%.2f, %.2f) ms',...
                obj.sz, obj.gz, 1000*obj.tStart, 1000*obj.tEnd));
            xlabel('X (m)'); ylabel('Z (m)'); 
            hold on;
            plot(obj.sx, obj.sz, '*r', 'LineWidth', 1, 'MarkerSize', 10);
            plot(obj.gx, obj.gz, '<g', 'MarkerFaceColor', 'g', 'MarkerSize', 5);
            hold off;
            axis tight;            
        end
        
        function plotHist(obj, varargin)
            histogram(reshape(obj.wp, [numel(obj.wp) 1]));
            title('Histogram of wave path amplitude');
            xlabel('Amp'); ylabel('Counts');
        end
        
        
        function savefig(obj, hdl, varargin)
            p = inputParser;
            addRequired(p, 'obj');
            addRequired(p, 'hdl');
            addParameter(p, 'figtype', 'pdf', @ischar);
            parse(p, obj, hdl, varargin{:})
            set(hdl, 'Color', 'w');
            if ~exist('figDir', 'var')
                figDir = fullfile(pwd, 'Figures');
                if ~exist(figDir, 'dir')
                    mkdir(figDir);
                end
            end
            
            figname = fullfile(figDir, ...
                    sprintf('sz%.fm_wpath_gz%.fm_t%.f_%.fms.%s', ...
                    obj.sz, obj.gz, 1000*obj.tStart, 1000*obj.tEnd, p.Results.figtype));
            % export_fig is an external library from github; it should be
            % added to path before calling this function
            export_fig(figname, hdl);
        end
        
        
    end
end