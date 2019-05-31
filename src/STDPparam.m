% Jan 2015
%
% timothee.masquelier@alum.mit.edu
%
% This code was used in: Masquelier T, Portelli G and Kornprobst P (2016). Microsaccades enable efficient synchrony-based coding in the retina: a simulation study. Scientific Reports. 
%
% Parameters for the STDP-based learning.


global PARAM

PARAM.goOn = 0;

if PARAM.goOn
    PARAM.nRun = 1;
%     PARAM.stdp_a_pos = 2^-4; % avoid > 2^-4
% %     PARAM.stdp_a_neg = -.9* PARAM.stdp_a_pos; % at least .85
%     PARAM.stdp_a_neg = -.875* PARAM.stdp_a_pos; % at least .85
% %     PARAM.threshold = Inf;
% %     neuron.nextFiring = Inf;
% %     PARAM.fixedFiringMode = false;
else % new computation
    clear all
    global PARAM
    PARAM.goOn = false;
    
    PARAM.dump=false;
    PARAM.beSmart=true; % save time by not computing the potential when it is estimated that the threshold cannot be reached. Rigorous only when no inhibition, and no PSP Kernel. Otherwise use at your own risks...
 
    PARAM.fixedFiringMode = false;
    PARAM.fixedFiringLatency = 10e-3;
    PARAM.fixedFiringPeriod = 150e-3;

    n=0; % spike train number
    PARAM.nRun = 1; % number of times spike train is propagated
    
    % Random generators
    PARAM.randomState = 0;
%     d = dir('../mat/rand*.mat');
%     if isempty(d)
%         PARAM.randomState=0;
%     else
%         last = d(end).name;
%         PARAM.randomState = str2num(last(5:7))+1;
%     end
%     % just to warn other threads that this random is done
%     tmp=0;
%     save(['../mat/rand' sprintf('%03d',PARAM.randomState) '.mat'],'tmp');

    PARAM.useSavedWeight = true;

    PARAM.jitter = 0*30e-3; % add jitter to input spike trains
    
    %********
    %* STDP *
    %********
    PARAM.stdp_alpha = 0.0; % STDP function: x^apha*exp(-x)
    PARAM.stdp_t_pos = 3e-3;%17e-3; %(source: Bi & Poo 2001)
    PARAM.stdp_t_neg = 3e-3;%34e-3; %(source: Bi & Poo 2001)
    
    if exist(['../mat/conv.' sprintf('%03d',PARAM.randomState) '.mat'],'file')
        load(['../mat/conv.' sprintf('%03d',PARAM.randomState) '.mat']);
        n_batch = 1+length(conv)/101;
    else
        n_batch = 1;
    end
    n_variable_rate = 5;
    rate_begin = 2^-8; % avoid > 2^-8 (Tim, Dec 2014)
    rate_end = 2^-8; % avoid > 2^-8 (Tim, Dec 2014)
    PARAM.stdp_a_pos = rate_begin^( 1 + (min(n_batch,n_variable_rate)-1)/(n_variable_rate-1) * (log(rate_end)/log(rate_begin)-1)) ...
        * PARAM.stdp_alpha^-PARAM.stdp_alpha*exp(PARAM.stdp_alpha); % avoid > 1e-2 (Tim, Jan 2011)
    PARAM.stdp_a_neg = - 1.0 * PARAM.stdp_a_pos;%- 1.001 * PARAM.stdp_a_pos; % .55 wih w_out = -.05*PARAM.stdp_a_pos
    PARAM.stdp_cut = 7;
    PARAM.minWeight = 0.0;
    PARAM.w_out = .0015 * PARAM.stdp_a_pos; % Homeostatic term. All synaptic weights are decreased by this amount at each postsynaptic spike (see Gilson, Masquelier, and Hugues. PLoS Comput Biol 2011.)
    
    %***************
    %* EPSP Kernel *
    %***************
    PARAM.tm = 20e-3; % membrane time constant (typically 10-30ms)
    PARAM.ts = 2e-3; % synapse time constant
    PARAM.epspCut = 5;% specifies after how many ms we neglect the epsp
    PARAM.tmpResolution = .05e-3;
    
    PARAM.memristor = false;
    PARAM.t_op = .0005; 
    PARAM.t_on = .0005;        
    if PARAM.memristor     % Bernabe Linares: EPSP suitable for STDP with memristors:
        tn = 3e-3;
        tp = 40e-3;
        t1 = 5e-3;
        t2 = 75e-3;
        An = 1;
        Ap = An/4;
        current = [ -An/(exp(t1/tn)-1)*(exp([PARAM.tmpResolution:PARAM.tmpResolution:t1]/tn)-1) Ap/(exp(t2/tp)-1)*(exp((t1+t2-[t1+PARAM.tmpResolution:PARAM.tmpResolution:t1+t2])/tp)-1) zeros(1,(PARAM.epspCut*PARAM.tm)/PARAM.tmpResolution) ];
    %     figure; plot(PARAM.tmpResolution*(1:length(current)),current)    
        f = exp(-[0:PARAM.tmpResolution:(t1+PARAM.epspCut*PARAM.tm)]/PARAM.tm); % filter it (~LIF)
        PARAM.epspKernel=filter(f,sum(f),current);
        PARAM.refractoryPeriod = t1+t2;     
    else   
        % Double exp (Gerstner 2002)
        PARAM.epspKernel = pspKernel(0:PARAM.tmpResolution:PARAM.epspCut*PARAM.tm,PARAM.ts,PARAM.tm);
        PARAM.refractoryPeriod = 1e-3;
        
        % Rapidly adapting thr (Fontaine et al PCB 2014)
        tau_thr = 5e-3;
        exp_filter = exp(-(0:PARAM.tmpResolution:PARAM.epspCut*tau_thr)/tau_thr);
        PARAM.epspKernel = PARAM.epspKernel-filter(exp_filter,sum(exp_filter),PARAM.epspKernel);
    end
    
    [m idx] = max(PARAM.epspKernel);
    PARAM.epspKernel = PARAM.epspKernel/m;
    %     figure; plot(PARAM.tmpResolution*(1:length(PARAM.epspKernel)),PARAM.epspKernel)
    PARAM.epspMaxTime = (idx-1)*PARAM.tmpResolution;
    
    % post synaptic spike kernel
    PARAM.usePssKernel = false;
    PARAM.pssKernel = -0.5*exp(-[0:PARAM.tmpResolution:PARAM.epspCut*200e-3]/(200e-3));
%     % time constant: tm
%     PARAM.pssKernel =   0*pspKernel(0:PARAM.tmpResolution:PARAM.epspCut*PARAM.tm,PARAM.ts/10,PARAM.tm/10) ...
%                     -   3*pspKernel(0:PARAM.tmpResolution:PARAM.epspCut*PARAM.tm,PARAM.ts,PARAM.tm) ...
%                     +   2*exp(-[0:PARAM.tmpResolution:PARAM.epspCut*PARAM.tm]/PARAM.tm);
%     % time constant: tm/2
%     PARAM.pssKernel =   -3*pspKernel(0:PARAM.tmpResolution:PARAM.epspCut*PARAM.tm/2,PARAM.ts/2,PARAM.tm/2) ...
%                     +   2*exp(-[0:PARAM.tmpResolution:PARAM.epspCut*PARAM.tm/2]/(PARAM.tm/2));
%     PARAM.pssKernel =   0*PARAM.epspKernel;

    % inhibitory postsynaptic potential (positive by convention, scaled so that max is 1)
    PARAM.ipspKernel = pspKernel(0:PARAM.tmpResolution:PARAM.epspCut*PARAM.tm,5e-3,PARAM.tm); % Wang 2002: tau_GABA = 5 ms
    PARAM.ipspKernel = PARAM.ipspKernel / (max(PARAM.ipspKernel));
    PARAM.inhibStrength = 2; % inhibition strength (in fraction of threshold)
    
%     % Simple exp
%     PARAM.epspKernel = exp(-[0:PARAM.tmpResolution:PARAM.epspCut*PARAM.tm]/PARAM.tm);
%     PARAM.epspMaxTime = 0;
%     disp(['Neglecting EPSP when below ' num2str(PARAM.epspKernel(end))])

    % figure;plot(0:PARAM.tmpResolution:PARAM.epspCut*PARAM.tm,PARAM.epspKernel);
    % return


%     %***************
%     %* Spike Train *
%     %***************
%     PARAM.nPattern = 3;
    PARAM.nAfferent = 2*80^2;
%     PARAM.oscillations = false;
%     PARAM.nCopyPasteAfferent = round( .5 * PARAM.nAfferent );
%     PARAM.dt = 1e-3;
%     PARAM.maxFiringRate = 90;
    PARAM.spontaneousActivity = 0;
%     PARAM.copyPasteDuration = 50e-3;
%     PARAM.jitter=1e-3;
%     PARAM.spikeDeletion=.0;
%     PARAM.maxTimeWithoutSpike = PARAM.copyPasteDuration;
%     PARAM.patternFreq = 1/3;
%     for idx = 1:PARAM.nPattern
%         PARAM.posCopyPaste{idx} = [];
%     end
%     PARAM.T = PARAM.nPattern*(500/PARAM.patternFreq)*PARAM.copyPasteDuration;
%     if PARAM.patternFreq>0
%         rand('state',PARAM.randomState);
%         skip = false;
%         for p = 1:round( PARAM.T / PARAM.copyPasteDuration )
%             if skip
%                 skip = false;
%             else
%                 if rand<1/(1/PARAM.patternFreq-1)
%                     idx = ceil(rand * PARAM.nPattern);
%                     PARAM.posCopyPaste{idx} = [PARAM.posCopyPaste{idx} p];
%                     skip = true; % skip next
%                 end
%             end
%         end
%     end
% %     PARAM.nCopyPaste = length(PARAM.posCopyPaste);
    
    %**********
    %* Neuron *
    %**********
    % The threshold corresponds roughly to the number of coincindent spikes we want to detect
    % Should be as big as possible (avoids FA), avoiding to get "stuck" in low
    % spike density zone of the copy-pasted pattern.
    % Appropriate value depends on number of afferent. Eg around 125 for 500
    % afferents.
    % Then initial weights should be tuned so as to reach threshold.
    
    PARAM.nNeuron = 3;
    %PARAM.threshold = 350;%450 for 40x40x2 RGC, 200 for 40x40x2 LGN
    PARAM.threshold = 60;%450 for 40x40x2 RGC, 200 for 40x40x2 LGN
    
    
%    PARAM.nuThr = round( 2.0 * PARAM.maxFiringRate/2*PARAM.nAfferent*PARAM.tm); % trace parameter for thr computation use inf not to use
    PARAM.nuThr = Inf;
%     if PARAM.goOn && PARAM.threshold==Inf
%         warning('"Going on" in infinite threshold mode has no sense. Setting goOn = false');
%         PARAM.goOn = false;
%     end
    %****************
    %* Neural codes *
    %****************
%     PARAM.realValuedPattern = false;
%     PARAM.codingMethod = 4; % 1 for poisson, 2 for LIF, 3 for intensity to phase, 4 for LIF with oscillatory drive
%     PARAM.gammaFreq = 50; % freq of oscillatory drive
%     PARAM.oscillMagnitude = 1;% magnitude of oscillatory drive
%     PARAM.oscillPhase = rand*2*pi;% phase of oscillatory drive
%     PARAM.resetPeriod = Inf;%100e-3;
%     PARAM.resetStd = 25e-3;
%     if PARAM.resetPeriod<Inf
%         if PARAM.resetStd==0
%             PARAM.resetTimes = PARAM.resetPeriod * [1:ceil(PARAM.T/PARAM.resetPeriod)];
%         else
%             PARAM.resetTimes = cumsum( PARAM.resetStd * randn(1,round(1.1*PARAM.T/PARAM.resetPeriod)) + PARAM.resetPeriod );
%         end
%     else
%         PARAM.resetTimes = [];
%     end
%     % LIF model for afferents
%     PARAM.R = 1;
%     PARAM.Vthr = 1.72;
%     PARAM.Vreset = 0;
%     PARAM.Vrest = 0;
% 	PARAM.Imin = 0;
% 	PARAM.Imax = 0.8;
%     
%     PARAM.interStimuliInterval = 1 / 16;
%     PARAM.interPatternInterval = 1 / 2;

    
    %*************
    %* Reporting *
    %*************
%     PARAM.plottingPeriod = [ [0 1]  [5 6] PARAM.nRun*PARAM.T+[-1 0] ];
%     PARAM.plottingPeriod = [ [-1 -1] [-1 -1] [-1 -1] ];
end
