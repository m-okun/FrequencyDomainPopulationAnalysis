data = load('populationSpikeData.mat');

raster = full(data.raster'); % units by time (1ms bins)
raster = raster(1:floor((size(raster, 1)/1000))*1000, :); % chop of the last fraction of a second
sr = 500; % sampling rate (Hz), must be several fold higher than the top frequency we are interested in (here 100 Hz).
raster = squeeze(sum(reshape(raster, 1000/sr, [], size(raster, 2)))); % convert from 1 kHz to sampling rate

pr = sum(raster, 2); % population rate
%pr = pr - mean(pr); % if we want to treat it as continuous signal

% Note that many of the clusters did not pass quality criteria, and were not analised as single, well-isolated units.
% Also, unit #1 is MUA from all spikes that were not separated into one of the other clusters.
% data.sua is true for clusters that were considered as well-isolated. The example units below are some of those.
exampleUnits = [1071 202 72 247];


%% Compute PSD and coherency with populaiton rate for example units

global CHRONUXGPU  % on our PCs with GPU this produces ~3-fold speed up. On PCs without GPU would make no difference
CHRONUXGPU = true; % relies on github.com/m-okun/chronuxGPU

mfr = mean(data.raster, 2) * 1000; % mean firing rate of each unit
mfr = full(mfr);

psd = [];
coh = [];
coh_conf = [];
phase = [];
phase_confU = [];
phase_confL = [];
coh_rateadjusted = [];
coh_confrateadjusted = [];

for pos = 1:numel(exampleUnits)
  u = find(data.units == exampleUnits(pos));
  if abs(mean(pr) < 1e-9) % i.e. if we want to treat population rate as continuous signal
    [freq, c, pspectd] = widebandCoherency(pr - raster(:,u), raster(:,u), 1/sr, ...
      'typespk1', 'c', ...
      'tapers', 5, ...
      'freqfactor', 1.333, ...
      'maxFreq', 100, ...
      'debug', pos == 1);
  else
    [freq, c, pspectd] = widebandCoherency(pr - raster(:,u), raster(:,u), 1/sr, ...
      'tapers', 5, ...
      'freqfactor', 1.333, ...
      'maxFreq', 100, ...
      'debug', pos == 1);
  end
  
  psd(end+1, :) = pspectd.power2;
  coh(end+1, :) = c.coh;
  coh_conf(end+1, :) = c.cohConf;
  phase(end+1,:) = c.phase;
  phase_confU(end+1,:) = c.phaseCu;
  phase_confL(end+1,:) = c.phaseCl;
  
  % compute rate adjusted coherence (i.e. coherence the unit would have had if its firing rate was 1 spk/s):
  kappa = (1 + (mfr(u) - 1)*mfr(u) ./ psd(end,:));
  kappa(kappa < 0) = NaN; % the values are not supposed to be negative, it sometimes happens for neurons with low firing rate,
  % presumably because of the effect refractory period has on power spectrum
  kappa = kappa.^-0.5;
  coh_rateadjusted(end+1,:) = coh(end,:) .* kappa;
  coh_confrateadjusted(end+1,:) = coh_conf(end,:) .* kappa;
  
  fprintf('.')
end
fprintf('\n')


%% Plot the results

LineWidth = 3;
scaling_power = 1/2; % scale the y-axis to show low values more prominently
opt.smooth = true;

for pos = 1:numel(exampleUnits)
  subplot(4, numel(exampleUnits), pos) % PSD plot
  if opt.smooth
    semilogx(freq, smooth(psd(pos,:)), 'k')
  else
    semilogx(freq, psd(pos,:), 'k')
  end
  ylim([0 max(psd(pos,:))]*1.2)
  box off
  
  subplot(4, numel(exampleUnits), pos+numel(exampleUnits)) % coherence plot
  if opt.smooth    
    lowerconf = max(0, smooth(coh_rateadjusted(pos, :) - coh_confrateadjusted(pos,:))).^scaling_power;
    upperconf = min(1, smooth(coh_rateadjusted(pos, :) + coh_confrateadjusted(pos,:))).^scaling_power;
  else
    lowerconf = max(0, coh_rateadjusted(pos, :) - coh_confrateadjusted(pos,:)).^scaling_power;
    upperconf = min(1, coh_rateadjusted(pos, :) + coh_confrateadjusted(pos,:)).^scaling_power;    
  end
  fill([freq, freq(end:-1:1)], [upperconf; lowerconf(end:-1:1)]', ...
    'k', 'FaceAlpha', 0.25, 'LineStyle', 'none')
  if opt.smooth    
    hold on, plot(freq, smooth(coh_rateadjusted(pos, :)).^scaling_power, 'k--', 'LineWidth', LineWidth-1)
  else
    hold on, plot(freq, coh_rateadjusted(pos, :).^scaling_power, 'k--', 'LineWidth', LineWidth-1)
  end
  set(gca, 'YTick', [0 0.04 (0.1:0.1:1)].^scaling_power)
  set(gca, 'YTickLabel', num2str([0 0.04 (0.1:0.1:1)]'))
  set(gca, 'XScale', 'log')
  box off
  
  subplot(4, numel(exampleUnits), [pos+2*numel(exampleUnits) pos+3*numel(exampleUnits)]); hold on; % phase plot
  lowerconf = bestUnwrap(phase_confL(pos,:), opt);
  upperconf = bestUnwrap(phase_confU(pos,:), opt);
  phi = phase(pos,:);
  phi(isnan(lowerconf)) = NaN;
  phi = bestUnwrap(phi, opt);
  f = freq;
  f(isnan(lowerconf)) = NaN;
  m = +Inf; % phase (to scale the plot)
  
  % break into intervals where phase is statistically significant (i.e. confidence interval is not NaN)
  lowerconf = chunkArrayWithNaNs(torow(lowerconf));
  upperconf = chunkArrayWithNaNs(torow(upperconf));
  phi = chunkArrayWithNaNs(torow(phi)); 
  f = chunkArrayWithNaNs(f);  
  for i = 1:numel(lowerconf)
    if all(lowerconf{i} + 2*pi < upperconf{i})
      lowerconf{i} = lowerconf{i} + 2*pi;
    end
    if all(lowerconf{i} > phi{i})
      phi{i} = phi{i} + 2*pi;
    end
    if all(upperconf{i} < phi{i})
      upperconf{i} = upperconf{i} + 2*pi;
    end
    fill([f{i}, f{i}(end:-1:1)], [upperconf{i}, lowerconf{i}(end:-1:1)], ...
      'k', 'FaceAlpha', 0.25, 'LineStyle', 'none');
    plot(f{i}, phi{i}, 'k', 'LineWidth', LineWidth);
    m = min(m, min(lowerconf{i}));
  end  
  
  set(gca, 'XScale', 'log', 'YTick', -2*pi:pi/2:2*pi, ...
    'YTickLabel', {'0', '\pi/2', '\pi', '-\pi/2', '0', '\pi/2', '\pi', '-\pi/2', '0'})
  for a = -2*pi:2*pi:2*pi
    hold on, plot([1e-2 1e2], a*[1 1], '--', 'Color', 0.5*[1 1 1])
  end
  ylim([m m+2*pi])
  box off
  
end

subplot(4, numel(exampleUnits), 1)
ylabel('PSD')
subplot(4, numel(exampleUnits), 1+numel(exampleUnits))
ylabel('Coherence (rate adjusted)')
subplot(4, numel(exampleUnits), [2*numel(exampleUnits)  3*numel(exampleUnits)]+1)
xlabel('Frequency (Hz)')
ylabel('Phase (rad)')
