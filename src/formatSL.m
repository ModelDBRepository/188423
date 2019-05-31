% Jan 2015
%
% timothee.masquelier@alum.mit.edu
%
% This code was used in: Masquelier T, Portelli G and Kornprobst P (2016). Microsaccades enable efficient synchrony-based coding in the retina: a simulation study. Scientific Reports. 
% 
% This script loads the spikes.spk text file produced by Virtual Retina,
% and stores them in files ../data/afferent.rand###.###.###.mat.

randState = 0; %seed for random generator (if any) and ref number for the computation
%n = uint16(80); % crop size
%delta = uint16(n);% crop shift

% %--------------------------------------
timedLogLn('formatSL')

%timedLog('Loading spikes...')

% load('../mat/spikes.spk')
% spikeList = spikes;
% clear spikes;
%
%spikeList = [];
%count = 1;
for i=0:7
    load(['../data/fr_' int2str(i) '/spikes.spk'])
    
    transient_period_to_remove = 0;
    spikes(spikes(:,2)<transient_period_to_remove,:) = [];
    spikes(:,2) = spikes(:,2)-transient_period_to_remove;
    
    if i==0
        %offset = ceil(spikes(end,2)/5e-3)*5e-3;
        offset = 0*15000/8;
        %offset = 1000/8-transient_period_to_remove;
    end
    afferentList = uint16(spikes(:,1));
    spikeList = spikes(:,2)+i*offset;  
    clear spikes;
    
%     afferentList(spikeList>=1000)=[];
%     spikeList(spikeList>=1000)=[];
    
    fileName = ['afferent.rand' sprintf('%03d',randState) '.' sprintf('%03d',0) '.' sprintf('%03d',i) '.mat'];
    disp(['Saving in ../data/' fileName])
    save(['../data/' fileName],'spikeList','afferentList')

end

