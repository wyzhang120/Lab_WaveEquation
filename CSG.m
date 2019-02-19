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
 
        function wp = getWP(obj, idxTrace, varargin)
            trace = obj.seis(:,idxTrace);
            switch nargin
                case 2
                    tStart = 0;  tEnd = (obj.nt-1)*obj.dt;
                case 3
                    tStart = 0;  tEnd = varargin{1};
                    idxMute = round(tEnd/obj.dt);
                    trace(idxMute:end)=0;
                case 5
                    tStart = varagin(1);  tEnd = varargin{2};
                    idxMute1 = round(tStart/obj.dt);
                    idxMute2 = round(tEnd/obj.dt);
                    trace(idxMute2:end)=0;
                    trace(1:idxMute1)=0;
                otherwise
                    error('max 4 input variables for getWP');
            end
             
            tic;
            fprintf('Computing wave path (sz=%.2f m, gz = %.2f m) ... \n',...
                obj.sz, obj.gz(idxTrace));
            wp0 = a2d_wavepath_abc28(trace, obj.vel,...
                obj.nbc, obj.dx, obj.nt, obj.dt, obj.wavelet,...
                obj.sx, obj.sz, obj.gx(idxTrace), obj.gz(idxTrace));
            toc;
            wp = WAVEPATH(obj.sx, obj.sz, obj.gx(idxTrace), obj.gz(idxTrace), ...
                obj.dx, tStart, tEnd, wp0);
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
        
        
    end
    
    methods(Static)
        function wpHist(wp)
            histogram(reshape(wp, [numel(wp) 1]));
            title('Histogram of wave path amplitude');
            xlabel('Amp'); ylabel('Counts');
        end
    end
end