classdef CSG
    properties
        sx;sz;gx;gz;
        seis; nbc; dx; nt; dt; wavelet; vel;
    end
    
    methods
        function obj = CSG(sx, sz, gx, gz)
            obj.sx = double(sx); obj.sz=double(sz);
            obj.gx = gx; obj.gz=gz;
        end
        
        function obj = getSeis(obj,vel,nbc,dx,nt,dt,fc,isFS)
            s=ricker(fc,dt);
            tic;
            fprintf('Computing shot gather at sx=%.2fm sz=%.2fm ... \n', obj.sx, obj.sz);
            obj.seis=a2d_mod_abc28(vel,nbc,dx,nt,dt,s,...
                obj.sx, obj.sz, obj.gx, obj.gz, isFS);
            obj.nbc = nbc; 
            obj.dx=dx; 
            obj.nt=nt; 
            obj.dt=dt; 
            obj.wavelet = s;
            obj.vel =vel;
            toc;
        end     
 
        function wp = getWP(obj, idxTrace, tMute)
            if tMute ~=0
                idxMute = round(tMute/obj.dt);
                trace = obj.seis(:,idxTrace); 
                trace(idxMute:end)=0;
            else
                trace = obj.seis(:,idxTrace);
            end
            tic;
            sprintf('Computing wave path (sz=%.2f m, gz = %.2f m) ... \n', ...
                obj.sz, obj.gz(idxTrace));
            wp = a2d_wavepath_abc28(trace, obj.vel,...
                obj.nbc, obj.dx, obj.nt, obj.dt, obj.wavelet,...
                obj.sx, obj.sz, obj.gx(idxTrace), obj.gz(idxTrace));
            toc;
        end
        
        function plotCSG(obj)
            figure;
            t=(0:obj.nt-1)*obj.dt;
            imagesc(obj.gz, t, obj.seis); colormap(gray);
            ylabel('Time (t)');
            xlabel('gz(m)'); title(sprintf('seismic gather src (z=%.2f m)', obj.sz));
        end
        
        function plotTrace(obj, idxTrace)
            figure;
            t=(0:obj.nt-1)*obj.dt;
            plot(t,obj.seis(:,idxTrace));xlabel('Time (t)'); 
            title(sprintf('seismic trace (sz=%.2f m, gz = %.2f m)', obj.sz, obj.gz(idxTrace)));
        end
        
        
        function plotWP(obj, wp, idxTrace, clim)
            figure;
            [nz, nx] = size(obj.vel); 
            x = (0:nx-1)*obj.dx; z = (0:nz-1)*obj.dx;
            imagesc(x,z,wp);colormap(gray);caxis(clim);
            title(sprintf('wave path (sz=%.2f m, gz = %.2f m)',...
                obj.sz, obj.gz(idxTrace)));
            xlabel('X (m)'); ylabel('Z (m)'); 
            hold on;
            plot(obj.sx, obj.sz, '*r', 'LineWidth', 1, 'MarkerSize', 10);
            plot(obj.gx(idxTrace), obj.gz(idxTrace), '<g', 'MarkerFaceColor', 'g', 'MarkerSize', 5);
            hold off;
        end

        
        function saveseis(obj)
            save(sprintf('seisSZ%.fm.mat', obj.sz), obj.seis);
        end
        
        function saveFig(obj, type, idxTrace)
            if ~exist('figDir', 'var')
                figDir = fullfile(pwd, 'Figures');
                if ~exist(figDir, 'dir')
                    mkdir(figDir);
                end
            end
            
            if strcmp(type, 'wavepath')
                figname = fullfile(figDir, ...
                    sprintf('%sz%.fm_%s_gz%.fm.fig', obj.sz, type, obj.gz(idxTrace)));
            else
                figname = fullfile(figDir, ...
                    sprintf('%sz_%s', type));
            end
            savefig(figname);
        end
        
        
    end
    
    methods(Static)
        function wpHist(wp)
            histogram(reshape(wp, [numel(wp) 1]));
            title('Histogram of wave path amplitude');
            xlabel('Amp'); ylabel('Counts');
        end
    end
end