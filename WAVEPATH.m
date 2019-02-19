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
            imagesc(x,z, obj.wp);colormap(gray);caxis(clim);
            title(sprintf('wave path (sz=%.2f m, gz = %.2f m)',...
                obj.sz, obj.gz));
            xlabel('X (m)'); ylabel('Z (m)'); 
            hold on;
            plot(obj.sx, obj.sz, '*r', 'LineWidth', 1, 'MarkerSize', 10);
            plot(obj.gx, obj.gz, '<g', 'MarkerFaceColor', 'g', 'MarkerSize', 5);
            hold off;
        end
        
        function plotHist(obj, varargin)
            histogram(reshape(obj.wp, [numel(obj.wp) 1]));
            title('Histogram of wave path amplitude');
            xlabel('Amp'); ylabel('Counts');
        end
        
        
        function savefig(obj, hdl)
            if ~exist('figDir', 'var')
                figDir = fullfile(pwd, 'Figures');
                if ~exist(figDir, 'dir')
                    mkdir(figDir);
                end
            end
            
            figname = fullfile(figDir, ...
                    sprintf('sz%.fm_wpath_gz%.fm.fig', obj.sz, obj.gz));
            savefig(hdl, figname);
        end
        
        
    end
end