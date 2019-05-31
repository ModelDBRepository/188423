% Jan 2015
%
% timothee.masquelier@alum.mit.edu
%
% This code was used in: Masquelier T, Portelli G and Kornprobst P (2016). Microsaccades enable efficient synchrony-based coding in the retina: a simulation study. Scientific Reports. 
%
% This code performs the STDP-based learning on the RGC spikes.
% The scripts successively loads all the ../data/afferent.rand###.###.###.mat files, and each time launches the STDPContinuous.c mex file (see comments in its header).
% All the parameters are gathered in param.m

if ~exist('PARAM','var')
    global PARAM
end

STDPparam

%rand('state',PARAM.randomState);

timedLog(['RANDOM STATE = ' int2str(PARAM.randomState) ]);

% if ~PARAM.goOn
%     % generate spike train
%     %     rand('state',PARAM.randomState);
%     %     randn('state',PARAM.randomState);
%     if exist(['../mat/afferent.rand' sprintf('%03d',PARAM.randomState) '.mat'],'file')
%         load(['../mat/afferent.rand' sprintf('%03d',PARAM.randomState) '.mat'])
%         %         if ~exist('spikeList','var') % 'false' file
%         %             [spikeList afferentList] = generateSpikeTrain;
%         %             save(['../mat/afferent.rand' sprintf('%03d',PARAM.randomState) '.mat'],'spikeList','afferentList')
%         %         end
%     else
%         timedLogLn('Generating spike train...');
%         if PARAM.realValuedPattern
%             [spikeList afferentList patternPeriod, values, times] = generateSpikeTrainWithRealValuedPattern2;
%             save(['../mat/afferent.rand' sprintf('%03d',PARAM.randomState) '.mat'],'spikeList','afferentList','patternPeriod','values','times')
%         else
%             [spikeList afferentList] = generateSpikeTrain;
%             save('-v7.3',['../mat/afferent.rand' sprintf('%03d',PARAM.randomState) '.mat'],'spikeList','afferentList')
%         end
%     end
% end


filePath = '../data/';
% filePath = '../../../data/';
fileList = dir([filePath 'afferent.rand' sprintf('%03d',PARAM.randomState) '.*.*.mat']);
disp([int2str(length(fileList)) ' files found']);

for f=1:length(fileList)
    %disp(['Loading ' filePath fileList(f).name ])
    load([filePath fileList(f).name])
    
    if isempty(spikeList)
        continue
    end

    % init neuron
    N = round(2^6.5*length(PARAM.epspKernel)*PARAM.tmpResolution*length(spikeList)/spikeList(end));
    for nn=1:PARAM.nNeuron %neuron loop
        neuron(nn) = createNewNeuron(PARAM,N);
    end %neuron loop

    if PARAM.useSavedWeight
        
%         disp(['Used weights saved in ../mat/weight.rand' sprintf('%03d',PARAM.randomState)  '.mat' ] )
%         tmp=load(['../mat/weight.rand' sprintf('%03d',PARAM.randomState)  '.mat' ]);
%         neuron.weight = tmp.weight;
%         clear tmp
        
        dirlist = dir([filePath 'weight.t=*s.txt']);
        % dirlist = dir(['../mat/sav.weight.24.txt']);
        % dirlist = dir(['../mat/sav.weight.t=0000080.479s.txt']);
        if ~isempty(dirlist)
            disp(['Used weights saved in ' filePath dirlist(end).name ])
            weight = load([filePath dirlist(end).name]);
            for nn=1:PARAM.nNeuron
                neuron(nn).weight = weight(nn,:);
            end
        end
        delete([filePath 'weight*.txt'])
        
    end
    
%     % tmp
%     disp(['Use weights based on rate.mat'])
%     load rate
%     neuron.weight = double(rate>15);
    % truncate
     %afferentList(spikeList>1040) =[];
     %spikeList(spikeList>1040) =[];
     
    
%      disp('Throwing away first 1s...')
%     afferentList(spikeList<468.7500) =[];
%     spikeList(spikeList<468.7500) =[];

    % add jitter
    if PARAM.jitter>0
        disp(['Adding a ' int2str(1000*PARAM.jitter) ' ms jitter to input spike trains'])
        spikeList = spikeList + PARAM.jitter * randn(size(spikeList));
        [spikeList idx] = sort(spikeList);
        afferentList = afferentList(idx);
        afferentList = afferentList(spikeList>0);
        spikeList = spikeList(spikeList>0);
    end
    
    % add spontaneous activity
    if PARAM.spontaneousActivity>0
        disp(['Adding ' num2str(PARAM.spontaneousActivity) ' Hz of Poisson spontaneous activity'])
        poissonSpikeList = cumsum(exprnd(1/PARAM.spontaneousActivity/PARAM.nAfferent,round(1.01*spikeList(end)*PARAM.spontaneousActivity*PARAM.nAfferent),1));
        poissonSpikeList(find(poissonSpikeList>spikeList(end),1,'first'):end) = []; % remove extra
        poissonAfferentList = uint16( floor(PARAM.nAfferent*rand(size(poissonSpikeList))) );
        spikeList = [ spikeList ;  poissonSpikeList];
        afferentList = [ afferentList ;  poissonAfferentList];
        clear poissonSpikeList
        clear poissonAfferentList
        [spikeList idx] = sort(spikeList);
        afferentList = afferentList(idx);        
    end    
    
    % PARAM.threshold = Inf;
    % stop = 0;
    for r=1:PARAM.nRun
        tic

        n=n+1;

        timedLogLn(['Run ' int2str(n) ' (' int2str(length(spikeList)) ' iterations ~ ' int2str(2e-9*length(spikeList)*length(neuron)*(ceil( PARAM.epspMaxTime / PARAM.tmpResolution )+1)) ' min )'])

        %     t=(n-1)*PARAM.T; % simulation time
%         if n>1
%             spikeList = spikeList+PARAM.T;% shift spike list
% 
%             maxLastFiring = 0;
%             for nn=1:length(neuron)
%                 if neuron(nn).nFiring>0 && neuron(nn).firingTime(neuron(nn).nFiring)>maxLastFiring
%                     maxLastFiring = neuron(nn).firingTime(neuron(nn).nFiring);
%                 end
%             end
%             spikeList = max(spikeList,maxLastFiring); % this is to avoid inserting epsp before last firing
%             %         for n=1:length(neuron) %neuron loop
%             %             neuron(n).nextFiring = Inf;
%             %         end
%             for nn=1:length(neuron) % make sure epsp are chronologic
%                 neuron(nn).nEpsp=0;
%                 neuron(nn).maxPotential=0;
%                 if PARAM.fixedFiringMode
%                     neuron(nn).nextFiring = (n-1)*PARAM.T + PARAM.fixedFiringLatency;
%                 else
%                     neuron(nn).nextFiring = Inf;
%                 end
%             end
%         end

        if PARAM.dump


            %         % beginning
            %         nSpike = round(1/3/PARAM.T*length(spikeList));
            %         delete('dump.txt');
            %         neuron=STDPContinuous(neuron,spikeList(1:nSpike),afferentList(1:nSpike)-1,true,false,PARAM);
            %         copyfile('dump.txt','dump.beginning.txt')
            %         break % no need to go further

            % end
            if r==PARAM.nRun
                nSpike = round((PARAM.T-1)/PARAM.T*length(spikeList));
                neuron=STDPContinuous(neuron,spikeList(1:nSpike),afferentList(1:nSpike)-1,false,PARAM.beSmart,PARAM);
                delete('dump.txt');
                neuron=STDPContinuous(neuron,spikeList(nSpike+1:end),afferentList(nSpike+1:end)-1,true,false,PARAM);
                copyfile('dump.txt','dump.end.txt')
            else
                neuron=STDPContinuous(neuron,spikeList,afferentList-1,false,PARAM.beSmart,PARAM); % C indexes start at 0
            end

        else
%             neuron=STDPContinuous(neuron,spikeList',afferentList',false,PARAM.beSmart,PARAM); % C indexes start at 0
            neuron=STDPContinuous(neuron,reshape(spikeList,[1 length(spikeList)]),reshape(afferentList,[1 length(afferentList)]),false,PARAM.beSmart,PARAM); % C indexes start at 0
        end

        %     % add new 'virgin' neurons
        %     if r<PARAM.nRun-1
        %         for nn=1:PARAM.nNeuron %neuron loop
        %             neuron = [ neuron createNewNeuron(PARAM,N) ];
        %         end
        %     end

        disp(' ');
        toc

        if sum([neuron.nFiring]) == 0
            warning('Neurons do not fire')
            break;
        end
        
        disp(['Avg firing rate =' num2str(mean([neuron(:).nFiring])/(spikeList(end)-spikeList(1))) ])
        disp(['Max firing rate =' num2str(max([neuron(:).nFiring])/(spikeList(end)-spikeList(1))) ])

    end % run loop

    if max([neuron.nFiring])>length(neuron(1).firingTime)
        warning('Increase firingTime array size')
    end

    % % save
    % if PARAM.fixedFiringMode
    %     weight = neuron.weight;
    %     save('../mat/weight.mat','weight');
    % end

    % save all
    c = clock;
    timeTag = [sprintf('%02.0f',c(2)) '.' sprintf('%02.0f',c(3)) '.' sprintf('%02.0f',c(4)) '.' sprintf('%02.0f',c(5)) '.' sprintf('%02.0f',c(6)) ];

    % perf;
    % for nn=1:length(neuron)
    %     timeTag = [timeTag '.n' sprintf('%02.0f',nn) '.HR' sprintf('%03.0f',100*neuron(nn).HR) '.FA' sprintf('%03.2f',neuron(nn).FA)  ];
    % end

    % if stop
    %     timeTag = [timeTag '.S'];
    % end
    % timeTag = [timeTag '.nF' int2str(neuron.nFiring)];

    % % dispMultiPattern
    % multiPattern
    % writeLatencies


    % % tags useful for multi patterns
    % hr = -sort(-max(HR,[],2));
    % timeTag = [timeTag '.pat'];
    % for pat=1:PARAM.nPattern
    %     timeTag = [timeTag '.' sprintf('%03.0f',100*hr(pat)) ];
    % end
    % hr = -sort(-max(HR,[],1));
    % timeTag = [timeTag '.neur'];
    % for neur=1:length(neuron)
    %     timeTag = [timeTag '.' sprintf('%03.0f',100*hr(neur)) ];
    % end

    % tags useful for mono pattern
    % timeTag = [timeTag '.HR.' sprintf('%03.0f',100*HR) '.FA.' sprintf('%f',FA) ];
    
    if exist([filePath 'conv.' sprintf('%03d',PARAM.randomState) '.mat'],'file')
        load([filePath 'conv.' sprintf('%03d',PARAM.randomState) '.mat']);
    else
        conv = [];
    end
    conv = [ conv convergence(filePath) ];
    save([filePath 'conv.' sprintf('%03d',PARAM.randomState) '.mat'],'conv');

end

%clear spikeList
%clear afferentList
clear patternPeriod
clear values
clear times

% disp(['Saving results in ' filePath 'matlab.rand' sprintf('%03d',PARAM.randomState) '.' timeTag '.mat']);
% save([filePath 'matlab.rand' sprintf('%03d',PARAM.randomState) '.' timeTag '.mat']);
% disp(timeTag)
% save

if PARAM.fixedFiringMode
    weight = neuron.weight;
    save([filePath 'weight.mat'],'weight');
end

if PARAM.dump
    showPot
end

disp([ int2str(mean(sum(reshape([neuron.weight],length(neuron(1).weight),length(neuron))>0.5))) ' selected weights on avg'])
% dispMultiPattern

% displayResult
% perf


