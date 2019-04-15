clear;close all;
load('Figures\wavepath53m.mat');
% trace chosen for wavepath computation
itrace = idxTrace;
trace = csg0.seis(:, itrace);
% central frequency of Ricker
fc = mod1.fc; fmax = 1./mod1.dt; nt = double(mod1.nt);
trF = fft(trace)/length(trace);
A0 = max(abs(trF));
%% low pass filter
filtOrder = 100; delayFilt = filtOrder/2;
fcut = [30, 50, 100, 150, 200, 250, 300, 400, 500];
ccode = ['b', 'g', 'y', 'm', 'c', 'r', '.b', '.g', '.y', '.m'];
nfilt = length(fcut);
trFilt = zeros(nt-delayFilt, nfilt);
for i=1:nfilt
    lpFilt = designfilt('lowpassfir','FilterOrder', filtOrder, ...
         'SampleRate', fmax, 'CutoffFrequency', fcut(i));
    temp = filter(lpFilt, trace);
    trFilt(:, i) = temp(delayFilt+1:end);
end
trFiltFFT = fft(trFilt)/length(trFilt(:,1));
bpFilt = designfilt('bandpassfir','FilterOrder', filtOrder, ...
         'SampleRate', fmax, 'CutoffFrequency1', 30, 'CutoffFrequency2', 300);
trBp = filter(bpFilt, trace);
trBp(1:delayFilt) = [];
trBpFFT = fft(trBp) / length(trBp);
%% plotting
t = mod1.dt * (0:nt-1);
wav = zeros(1, nt);
nwav = length(csg0.wavelet);
wav(1:nwav) = csg0.wavelet;
fwav = (0:nt-1) * mod1.df;
wavF = fft(wav)/length(wav);
figure(1);
subplot(1,2,1);
plot(t, wav);
xlim([0, 0.01]);
xlabel('t(s)'); ylabel('Amp');
subplot(1,2,2);
plot(fwav, abs(wavF));
xlim([0, 1000]);
xlabel('freq(Hz)'); ylabel('Amp');
title('Ricker 300 Hz');

    
plotFir(t, fwav, sprintf('Trace %d', itrace), trace, trF, ...
        trBp, trBpFFT, delayFilt, ccode(1), [30, 300]);
for i=1:nfilt
    plotFir(t, fwav, sprintf('Trace %d', itrace), trace, trF, ...
        trFilt(:,i), trFiltFFT(:,i), delayFilt, ccode(i), fcut(i))
end
lpFirOverlay(t, fwav, nfilt, itrace, trace, trFilt, trF, trFiltFFT, delayFilt, ccode)
%%
function plotFir(t, fwav, ititle, trace, trF, trFilt, trFiltFFT, delayFilt, colorcode, fcut)
    figure;
    subplot(1,3,1);
    plot(t(1:end-delayFilt), trace(1:end-delayFilt), 'k'); hold on; 
    plot(t(1:end-delayFilt), trFilt, colorcode);
    hold off;
    xlabel('t(s)'); ylabel('Amp');
    xlim([0.05, 0.1]);
    title(ititle);
    subplot(1,3,2);
    plot(fwav, abs(trF), 'k'); hold on; 
    plot(fwav(1:end-delayFilt), abs(trFiltFFT), colorcode); 
    hold off;
    xlim([0, 1000]); xticks(0:100:1000);
    xlabel('freq(Hz)'); ylabel('Amp');
    title('Amp(f)');
    subplot(1,3,3);
    plot(fwav, unwrap(angle(trF)), 'k'); hold on; 
    plot(fwav(1:end-delayFilt), unwrap(angle(trFiltFFT)), colorcode);
    hold off;
    xlim([0, 1000]);xticks(0:100:1000);
    xlabel('freq(Hz)'); ylabel('Phase');
    title('Phase(f)');
    if length(fcut) == 2
        lgd2 = sprintf('%g - %g Hz', fcut(1), fcut(2));
    elseif length(fcut) == 1
        lgd2 = sprintf('%g Hz',fcut);
    end
    legend('original', lgd2);
end
%% 
function lpFirOverlay(t, fwav, nfilt, itrace, trace, trFilt, trF, trFiltFFT, delayFilt, ccode)
    figure;
    subplot(1,3,1);
    plot(t(1:end-delayFilt), trace(1:end-delayFilt), 'k'); hold on; 
    for i=1:nfilt
        plot(t(1:end-delayFilt), trFilt(:, i), ccode(i)); 
    end
    hold off;
    xlabel('t(s)'); ylabel('Amp');
    xlim([0.05, 0.1]);
    title(sprintf('Trace %d', itrace));
    subplot(1,3,2);
    plot(fwav, abs(trF), 'k'); hold on; 
    for i=1:nfilt
    plot(fwav(1:end-delayFilt), abs(trFiltFFT(:, i)), ccode(i)); 
    end
    hold off;
    xlim([0, 1000]);xticks(0:100:1000);
    xlabel('freq(Hz)'); ylabel('Amp');
    title('Amp(f)');
    subplot(1,3,3);
    plot(fwav, unwrap(angle(trF)), 'k'); hold on; 
    for i=1:nfilt
    plot(fwav(1:end-delayFilt), unwrap(angle(trFiltFFT(:,i))), ccode(i)); 
    end
    hold off;
    xlim([0, 1000]);xticks(0:100:1000);
    xlabel('freq(Hz)'); ylabel('Phase');
    title('Phase(f)');
    legend('original', 'filtered');
end