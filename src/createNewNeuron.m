function neuron = createNewNeuron(PARAM,N)
%         N = round(1.75*PARAM.epspCut*PARAM.tm*length(spikeList)/spikeList(end));
neuron.epspAmplitude = zeros(1,N);
neuron.epspTime = zeros(1,N);
neuron.epspAfferent = uint16(zeros(1,N));
neuron.nEpsp = 0;
if PARAM.fixedFiringMode
    neuron.nextFiring = PARAM.fixedFiringLatency;
else
    neuron.nextFiring = Inf;
end
neuron.firingTime = zeros(1,1000000);
neuron.nFiring = 0;
neuron.alreadyDepressed = false(1,PARAM.nAfferent);
%         neuron.nInefficient = 0;
neuron.maxPotential = 0;
%         neuron.currentPotential = NaN;
neuron.trPot = 0;%PARAM.initialTr;


%         neuron(nn).initialPot = 0;
%         neuron(nn).lastInhib = NaN;

neuron.ipspTime = zeros(1,1000);
neuron.nIpsp = 0;

    if PARAM.fixedFiringMode
        neuron.weight = .5*ones(1,PARAM.nAfferent);
    else
        if PARAM.threshold == Inf % must be a post-fixedFiringMode computatio
            load('../mat/weight.mat');
            neuron.weight = weight;
    %         neuron.weight = double([neuron.weight>.5]);
        else
            %mw = .5;
            %neuron.weight = 1 - (2*(1-mw)) * rand(1,PARAM.nAfferent); % random
            neuron.weight = 0 + 1 * rand(1,PARAM.nAfferent);
            %neuron.weight = mw*ones(1,PARAM.nAfferent); % all equal
            %neuron.weight = rand(1,PARAM.nAfferent)>1-mw; % binary
            neuron.weight = max(0,min(1,neuron.weight));
        end
    end
