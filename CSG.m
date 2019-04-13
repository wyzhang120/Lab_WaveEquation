classdef CSG
    properties
        sx;sz;gx;gz;
        seis; nbc; dx; nt; dt; wavelet; vel;
        tTable; tFA;
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
            % idxTrace: int, index of trace
            % optional keyword input parameters
            % wavetype: string; one of below; default is 'direct'
            %           'direct', 'diving', 'refrac','diffr', 'refl'
            % tStart: float, starting time of the window to be migrated,
            %         default =0
            % tEnd: float, starting time of the window to be migrated,
            %         default =0
            % reflCoeff: float, array, reflectivity
            %           this array is nececary for diffraction and
            %           reflection waves
            trace = obj.seis(:,idxTrace);
            p = inputParser;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            expectedShapes = {'direct', 'diving', 'refrac','diffr', 'refl'};
            addRequired(p, 'obj');
            addRequired(p, 'idxTrace', validScalarPosNum);            
            addParameter(p, 'tStart', 0, validScalarPosNum);
            addParameter(p, 'tEnd', -1, validScalarPosNum);
            addParameter(p, 'wavetype', 'direct', @(x) any(validatestring(x,expectedShapes)));
            addParameter(p, 'reflCoeff', 0, @isnumeric)
            parse(p, obj, idxTrace, varargin{:})
            
            wavetype = p.Results.wavetype;
            refl_ss = p.Results.reflCoeff;
            tStart = p.Results.tStart;
            if tStart == 0
               idxMute1 = [];
            else
               idxMute1 = round(tStart/obj.dt);
            end
            
            tEnd = p.Results.tEnd;
            if tEnd < 0
                tEnd = (obj.nt-1)*obj.dt;
                idxMute2 = [];
            else
                idxMute2 = round(tEnd/obj.dt);
            end            
            
            if ~isempty(idxMute1)
                trace(1:idxMute1) = 0;
            end
            
            if ~isempty(idxMute2)
                trace(idxMute2:end) = 0;
            end
            
            tic;
            fprintf('Computing wave path (sz=%.2f m, gz = %.2f m) ... \n',...
                obj.sz, obj.gz(idxTrace));
            switch wavetype
                case {'direct', 'diving', 'refrac'}
                   wp0 = a2d_wavepath_abc28(trace, obj.vel,...
                            obj.nbc, obj.dx, obj.nt, obj.dt, obj.wavelet,...
                            obj.sx, obj.sz, obj.gx(idxTrace), obj.gz(idxTrace)); 
                case {'diffr', 'refl'}
                    if refl_ss == 0
                        error('Reflectivity keyword param "reflCoeff" not provided');
                    else
                    wp0=a2d_wavepath_abc28_refl(trace, obj.vel, refl_ss, ...
                        obj.nbc, obj.dx, obj.nt, obj.dt, obj.wavelet, ...
                        obj.sx, obj.sz, obj.gx(idxTrace), obj.gz(idxTrace));
                    end
                otherwise
                    error('wavetype not recognized');
            end            
            toc;
            wp = WAVEPATH(obj.sx, obj.sz, obj.gx(idxTrace), obj.gz(idxTrace), ...
                obj.dx, tStart, tEnd, wp0);
        end
        
        function obj = getFirstArrival(obj)
            if isempty(obj.vel)
                error('velocity model not available; run "getSeis" first');
            end
            
            if ~isempty(obj.nbc)
                pad = obj.nbc;
            else
                pad = 20;
            end
            vppad = padarray(obj.vel, [pad pad], 'replicate','both' );
            [padNz, padNx]=size(vppad);
            srcx1 = pad + round(obj.sz/obj.dx) + 1; 
            srcx2 = pad + round(obj.sx/obj.dx) + 1; 
            slnPad=1./vppad;
            fprintf('Computing time table sx=%.2fm sz=%.2fm... \n', obj.sx, obj.sz);
            tic;
            ttPad = tt(slnPad, padNz, padNx, obj.dx, srcx1, srcx2);
            toc;
            obj.tTable = ttPad(pad + 1 : padNz - pad, pad + 1 : padNx - pad);
            idxGz = round(obj.gz/obj.dx) + 1;
            obj.tFA = obj.tTable(idxGz, padNx - 2*pad);
        end
        
        function plotTimeTable(obj)
            hdl = figure;
            [nz, nx] = size(obj.vel);
            h = obj.dx;
            imagesc(h*((0:nx-1)), h*((0:nz-1)), obj.tTable);
            colormap(hdl, jet);
            axis equal;axis tight;colorbar;
            xlabel('x[m]'); ylabel('z[m]');
            title('Traveltime table [ms]');
            hold on;
            plot(obj.sx, obj.sz, '*r', 'LineWidth', 1, 'MarkerSize', 10);
            hold off;
        end
  
        function hdl = plotCSG(obj, varargin)
            % optional keyword input parameters
            % FA: string; 'y': mark first arrival time; 'n': default, no mark
            p = inputParser;
            addRequired(p, 'obj');
            addParameter(p, 'FA', 'n', @ischar);
            addParameter(p, 'ylim', 'auto');
            addParameter(p, 'idxTrace', 0);
            parse(p, obj, varargin{:}) 
            idxTrace = p.Results.idxTrace;
            hdl = figure;
            t=(0:obj.nt-1)*obj.dt;
            imagesc(obj.gz, t, obj.seis); colormap(gray);
            ylabel('Time (t)'); ylim(p.Results.ylim);
            xlabel('gz(m)'); title(sprintf('shot gather (sz=%.2f m)', obj.sz));
            if strcmp(p.Results.FA, 'y')
                if isempty(obj.tFA)
                    error('first arrival not computed yet; run "getFirstArrival" before plotting');
                end
                hold on;                
                plot(obj.gz, obj.tFA, 'r');
                hold off;
            end
            
            if idxTrace
                ntmp = 100;
                xtmp = obj.gz(idxTrace) * ones(ntmp);
                ylimTmp = ylim;
                ytmp = linspace(ylimTmp(1), ylimTmp(2), ntmp);
                hold on;
                plot(xtmp, ytmp, '--c');
                hold off;
            end
        end
        
        function hdl = plotTrace(obj, idxTrace, varargin)
            % optional keyword input parameters
            % FA: string; 'y': mark first arrival time; 'n': default, no mark
            p = inputParser;
            addRequired(p, 'obj');
            addParameter(p, 'FA', 'n', @ischar);
            addParameter(p, 'xlim', 'auto');
            addParameter(p, 'ylim', 'auto');
            addParameter(p, 'LineWidth', 1);
            addParameter(p, 'tVertical', 'n', @ischar );
            parse(p, obj, varargin{:}) 
            hdl = figure;
            t=(0:obj.nt-1)*obj.dt;
            plot(t,obj.seis(:,idxTrace), 'LineWidth', p.Results.LineWidth);xlabel('Time (s)'); 
            xlim(p.Results.xlim); ylim(p.Results.ylim);
            title(sprintf('seismic trace (sz=%.2f m, gz = %.2f m)', obj.sz, obj.gz(idxTrace)));
            if strcmp(p.Results.FA, 'y')
                if isempty(obj.tFA)
                    error('first arrival not computed yet; run "getFirstArrival" before plotting');
                end
                hold on;                
                plot(obj.tFA(idxTrace), 0, '*r', 'LineWidth', 1, 'MarkerSize', 10);
                hold off;
            end
            
            if strcmp(p.Results.tVertical, 'y')
                view([90 -90]);
                set(gca, 'xdir', 'reverse');
            end
        end

        
        function saveseis(obj, hdl, type, varargin)
            p = inputParser;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            addRequired(p, 'obj');
            addRequired(p, 'hdl');
            addRequired(p, 'type', @ischar);
            addParameter(p, 'idxTrace', 0, validScalarPosNum);
            addParameter(p, 'figtype', 'pdf', @ischar)
            parse(p, obj, hdl, type, varargin{:})
            type = p.Results.type;
            idxTrace = p.Results.idxTrace;

            if ~exist('figDir', 'var')
                figDir = fullfile(pwd, 'Figures');
                if ~exist(figDir, 'dir')
                    mkdir(figDir);
                end
            end
            switch type
                case 'gather'
                    figname_= 'gather';
                case 'trace'
                    if idxTrace == 0
                        error('keyword input idxTrace needed');
                    else
                        figname_= sprintf('trace_gz%.fm', obj.gz(idxTrace));
                    end
                otherwise
                    figname_ = '';
            end
            figname = fullfile(figDir, ...
                    sprintf('sz%.fm_%s.%s', obj.sz, figname_, p.Results.figtype));
            saveas(hdl, figname);
        end
        
        
    end
    
end