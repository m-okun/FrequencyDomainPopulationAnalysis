% performs coherency/power analysis on point or continuous signals
% Inputs: spk1 - a column vector representing a non-constant signal in time (continuous or spike counts)
%         spk2 - a matrix where each column represents a non-constant signal in time, s.t. size(spk2, 1) = numel(spk1)
%         dt - 1/(sampling rate) for spk1, spk2
%         Additional optional inputs (to be provided as name-value pairs)
%           typespk1 - 'pb' (default) or 'c', describes spk1
%           typespk2 - 'pb' (default) or 'c', describes spk2
%           segfactor - 10 by default (duration of segments is at least this many times the 1/(highest frequency) we estimate using it, must be >= 2)
%           freqfactor - 2 by default (has to be > 1, duration of segment is at least
%                        segfactor/freqfacor times the 1/(lowest frequency) we estimate using it)
%           tapers - 3 by default (one number)
%           precomputetapers - (false by default), if true tapers are precomputed in advance (could potentially save running time)
%           debug - if true will print information on sizes of segments used
%           maxFreq - by default it is 0.5/dt
%           decimate - (true by default), if false will not decimate the signals (for low frequencies) to reduce runtime
%           monotoneFreq - false by default, if true the results won't have the same frequency twice (possible at transition across different segment durations,
%                              and can be used as a control for the effect of segment duration)
% Outputs: freq - frequencies
%          c - coherency structure with fields:
%             c.coh  - coherence between spk1 and each column of spk2
%             c.phase- phase of each column of spk2 wrt spk1, e.g. a small positive (linear) phase means spk2 precedes spk1
%             c.cohConf - 95% confidence interval for coherence
%             c.phaseCu, c.phaseCl - upper and lower confidence intervals for phase
%          pspectd - power spectral density structure with fields:
%             pspectd.power1, pspectd.conf1 - PSD and confidence for spk1 signal
%             pspectd.power2, pspectd.conf2 - PSD and confidence for spk2 signals
%
function [freq, c, pspectd] = widebandCoherency(spk1, spk2, dt, varargin)

assert(size(spk1, 2) == 1, 'widebandCoherency: spk1 must be a column vector')
assert(max(spk1) - min(spk1) > 0, 'widebandCoherency: spk1 is constant?!')
assert(numel(spk1) == size(spk2, 1), 'widebandCoherency: spk1 and spk2 do not correspond to each other in size')

inp = inputParser;
inp.addParameter('typespk1', 'pb', @(x) strcmp(x, 'pb') || strcmp(x, 'c'));
inp.addParameter('typespk2', 'pb', @(x) strcmp(x, 'pb') || strcmp(x, 'c'));
inp.addParameter('segfactor', 10, @(x) isnumeric(x) && isscalar(x) && (x>0) );
inp.addParameter('freqfactor', 2, @(x) isnumeric(x) && isscalar(x) && (x>0) );
inp.addParameter('tapers', 3, @(x) isnumeric(x) && isscalar(x) && (x>0) );
inp.addParameter('maxFreq', 0.5/dt, @isnumeric);
inp.addParameter('precomputetapers', false, @islogical);
inp.addParameter('debug', false, @islogical);
inp.addParameter('decimate', true, @islogical);
inp.addParameter('monotoneFreq', false, @islogical);

inp.parse(varargin{:});
assert(inp.Results.segfactor >= 2, 'widebandCoherency: segfactor of less than 2 not allowed')

% There are two issues:
% * To have an estimate for some particular frequency f, we need at least
% a few segments of duration 1/f. So the length of available data provides an
% automatic constraint on how low we can go.
%
% * On the other hand for (high) frequency f, no point using intervals longer than (say) 20/f.

params.Fs = 1/dt;
params.err = [1 0.05];

startFreq = inp.Results.maxFreq;

c.coh = [];
c.phase = [];
c.cohConf = [];
c.phaseCu = [];
c.phaseCl = [];
for spk2col = 1:size(spk2, 2)
  freq{spk2col} =       [];
  coh{spk2col} =        [];
  phase{spk2col} =      [];
  conf{spk2col} =       [];
  phaseCu{spk2col} =    [];
  phaseCl{spk2col} =    [];
end

if nargout > 2 % need to compute pspectd
  pspectd.f = []; % concatenate to the output variable
  pspectd.power1 = [];
  pspectd.power1_conf = [];
  pspectd.power2 = [];
  pspectd.power2_conf = [];
end

while true
  L = inp.Results.segfactor/startFreq; % segment duration for this iteration
  if inp.Results.segfactor*L > numel(spk1)*dt
    if inp.Results.debug
      fprintf('Cannot use segments of %2.1f sec. (total data duration of %2.1f sec < %1.1f * %2.1f) ==> done\n', ...
        L, numel(spk1)*dt, inp.Results.segfactor, L);
    end
    break % We cannot have enough repeats with segments of this size
  end
  
  params.fpass = [startFreq/inp.Results.freqfactor startFreq];
  if params.Fs > 500*params.fpass(1) && inp.Results.decimate
    % we're now estimating frequencies that are way lower than signals' sampling rate, 
    % so we will decimate the signals, to make the computations more efficient
    spk2 = spk2(1:floor(numel(spk1)/10)*10, :);
    spk1 = spk1(1:floor(numel(spk1)/10)*10);
    
    if strcmpi(inp.Results.typespk1, 'pb')
      spk1 = sum(reshape(spk1, 10, []))';
    else
      spk1 = decimate(spk1, 10);
      spk1 = spk1 - mean(spk1);
    end
    if strcmpi(inp.Results.typespk2, 'pb') && size(spk2, 2) > 1
      spk2 = squeeze(sum(reshape(spk2, 10, size(spk2, 1)/10, size(spk2, 2))));
    elseif strcmpi(inp.Results.typespk2, 'pb') %spk2 is just one column
      spk2 = sum(reshape(spk2, 10, []))';
    elseif strcmpi(inp.Results.typespk2, 'c')
      tmp = zeros(size(spk2, 1)/10, size(spk2, 2));
      for i = 1:size(spk2, 2)
        tmp(:,i) = decimate(spk2(:,i), 10);
        tmp(:,i) = tmp(:,i) - mean(tmp(:,i));
      end
      spk2 = tmp;
    else
      error('widebandCoherency: bad input')
    end
    if inp.Results.debug
      fprintf('Signal decimation (from %2.2f Hz to %2.2f Hz)\n', params.Fs, params.Fs/10)
    end
    params.Fs = params.Fs/10;
    dt = 1/params.Fs;
  end % signal decimation
  
  if inp.Results.precomputetapers
    % precompute tapers (can save time) for all subsequent computations
    N = size(createdatamatc(spk1, 0:L:numel(spk1)/params.Fs-L, params.Fs, [0 L]), 1);
    % N is the length of segments to which mtspectrumsegc will break spk1
    params.tapers = dpsschk(inp.Results.tapers*[1 1], N, params.Fs);
    clear N
  else
    params.tapers = [inp.Results.tapers/L L inp.Results.tapers]; % vector [W T p] where W is the bandwidth, 
                                                                 % T is the duration of the data and p is an integer 
                                                                 % such that 2TW-p tapers are used.
  end
  
  if nargout > 2 % power spectrum computation   
    switch inp.Results.typespk1
      case 'pb'
        [psd1, f, ~, ~, ~, ~, conf1] = mtspectrumsegpb(spk1, L, params, true, true);
      case 'c'
        [psd1, f, ~, ~, conf1] = mtspectrumsegc(spk1, L, params, true);
      otherwise
        error('widebandCoherency: unknown type of spk1')
    end
    assert(all(conf1(2,:)' >= psd1 | isnan(conf1(2,:)'))) % top range of confidence interval cannot be below the result itself
                                                          % if most segments have no spikes confidence can be NaN
    psd1_conf = conf1(1,:)' - psd1; % i.e. after adding it to psd1 we will be getting back the top range of confidence interval

    pspectd.f = [pspectd.f; f']; % concatenate results for the present scale to the output variable
    pspectd.power1 = [pspectd.power1 psd1'];
    pspectd.power1_conf = [ pspectd.power1_conf psd1_conf'];
    
    for spk2col = 1:size(spk2, 2)
      switch inp.Results.typespk2
        case 'pb'
          [psd2, ~, ~, ~, ~, ~, conf2] = mtspectrumsegpb(spk2(:, spk2col), L, params, true, true);
        case 'c'
          [psd2, ~, ~, ~, conf2] = mtspectrumsegc(spk2(:, spk2col), L, params, true);
        otherwise
          error('widebandCoherency: unknown type of spk2')
      end
      assert(all(conf2(2,:)' >= psd2 | isnan(conf2(2,:)'))) % if most segments have no spikes confidence can be NaN
      psd2_conf = conf2(1,:)' - psd2; 
      
      % concatenate to the output variable:
      pspectd.power2(spk2col, numel(pspectd.f)-numel(f)+1:numel(pspectd.f)) = psd2;
      pspectd.power2_conf(spk2col, numel(pspectd.f)-numel(f)+1:numel(pspectd.f)) = psd2_conf;      
    end    
  end % pspectd computation  
    
  for spk2col = 1:size(spk2, 2) % coherence computation
    assert(max(spk2(:, spk2col)) - min(spk2(:, spk2col)) > 0, 'widebandCoherency: spk2 is constant?!')
    if strcmpi(inp.Results.typespk1, 'pb') && strcmpi(inp.Results.typespk2, 'pb')
      [~,phi,~,~,~,f,zerosp] = coherencysegpb(spk1, spk2(:, spk2col), L, params, false);
      [C, ~, ~, ~, ~, ~, ~, confC] = coherencysegpb(spk1, spk2(:, spk2col), L, params, true); % recompute with averaging, to get confidence for coherence      
    elseif strcmpi(inp.Results.typespk1, 'c') && strcmpi(inp.Results.typespk2, 'c')
      assert(abs(mean(spk1)) < 1e-3 && abs(mean(spk2(:, spk2col))) < 1e-3, 'widebandCoherency: continuous data must be mean subtracted')
      [~,phi,~,~,~,f] = coherencysegc(spk1, spk2(:, spk2col), L, params, false);
      [C, ~, ~, ~, ~, ~, ~, confC] = coherencysegc(spk1, spk2(:, spk2col), L, params, true);
    elseif strcmpi(inp.Results.typespk1, 'c') && strcmpi(inp.Results.typespk2, 'pb')
      assert(abs(mean(spk1)) < 1e-3, 'widebandCoherency: continuous data must be mean subtracted')
      [~,phi,~,~,~,f,zerosp] = coherencysegcpb(spk1, spk2(:, spk2col), L, params, false);
      [C, ~, ~, ~, ~, ~, ~, confC] = coherencysegcpb(spk1, spk2(:, spk2col), L, params, true);
    elseif strcmpi(inp.Results.typespk1, 'pb') && strcmpi(inp.Results.typespk2, 'c')
      assert(abs(mean(spk2)) < 1e-3, ['widebandCoherency: continuous data must be mean subtracted' num2str(abs(mean(spk2)))])
      [~,phi,~,~,~,f,zerosp] = coherencysegcpb(spk2(:, spk2col), spk1, L, params, false);
      phi = -phi;
      [C, ~, ~, ~, ~, ~, ~, confC] = coherencysegcpb(spk2(:, spk2col), spk1, L, params, true);
    else
      error('widebandCoherency: unknown type of spk1 & spk2')
    end
    
    if numel(confC) == 1 % confC is just a number which we turn into a vector corresponding to frequencies
      confC = confC*ones(size(f));
    else % for coherency between two continuous signals confC is actually a vector (for each frequency)
      assert(numel(confC) == numel(f))
    end    
    
    % In some intervals there may be 0 spikes, and so C & phi are NaN there (this is clearly irrelevant for continuous signals)
    warning('off') % In this case circ_mean will also issue a warning which we suppress 
    if exist('zerosp', 'var')
      [phi, phiCu, phiCl] = circ_mean(phi(:, ~zerosp), [], 2);
    else % when both signals are continuous
      [phi, phiCu, phiCl] = circ_mean(phi, [], 2);
    end
    warning('on') 
    % In addition, rare cases (where the mean is on just 2 values?) were encountered where confidence intervals were complex
    phiCu(abs(phiCu - real(phiCu)) > 0) = NaN;
    phiCl(abs(phiCl - real(phiCl)) > 0) = NaN;
    
    if inp.Results.debug && spk2col == 1
      fprintf('Using segments of %2.3f sec. for frequencies in the %2.3f - %2.3f Hz range (%d tapers provide a resolution of %2.3f Hz)\n', ...
        L, params.fpass, inp.Results.tapers, params.tapers(1));
    end
    
    freq{spk2col} =       [freq{spk2col} torow(f)];
    coh{spk2col} =        [coh{spk2col} torow(C)];
    phase{spk2col} =      [phase{spk2col} torow(phi)];
    conf{spk2col} =       [conf{spk2col} torow(confC)];
    phaseCu{spk2col} =    [phaseCu{spk2col} torow(phiCu)];
    phaseCl{spk2col} =    [phaseCl{spk2col} torow(phiCl)];
  end % loop on columns of spk2
  startFreq = startFreq/inp.Results.freqfactor;
end % looping on different segment sizes

freq = cell2mat(freq'); freq = freq(1, :);
if nargout > 2
  assert(max(abs(pspectd.f - freq')) < 1e-9) % frequencies for which power was computed should be identical
end

[~,i] = sort(freq);
if inp.Results.monotoneFreq
  I = [true diff(freq(i)) > 0];
  i = i(I);
end
freq = freq(i);
c.coh =     cell2mat(coh');     c.coh = c.coh(:,i);
c.phase =   cell2mat(phase');   c.phase = c.phase(:,i);
c.phaseCu = cell2mat(phaseCu'); c.phaseCu = c.phaseCu(:,i);
c.phaseCl = cell2mat(phaseCl'); c.phaseCl = c.phaseCl(:,i);
c.cohConf    = cell2mat(conf'); c.cohConf = c.cohConf(:,i);

% It can happen that the function finds a phase (with quite a small confidence interval!)
% where there is none by construction (and the coherence was correctly identified as 0).
% Therefore we do the following:
c.phaseCu(c.coh - c.cohConf <= 0) = NaN;
c.phaseCl(c.coh - c.cohConf <= 0) = NaN;

if nargout > 2
  pspectd.f = freq;
  pspectd.power1 = pspectd.power1(i);
  pspectd.power1_conf = pspectd.power1_conf(i);
  pspectd.power2 = pspectd.power2(:,i);
  pspectd.power2_conf = pspectd.power2_conf(:,i);  
end