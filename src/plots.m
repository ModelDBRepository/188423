weight = zeros(PARAM.nAfferent,PARAM.nNeuron);
for n=1:PARAM.nNeuron
    weight(:,n) = neuron(n).weight';
end
disp([num2str(mean(sum(weight>0.5))) ' synapses selected on avg' ])


% run ../../plotRF
figure
colors = get(gcf,'DefaultAxesColorOrder');
subplot(3,1,1)
for n=1:PARAM.nNeuron
    if neuron(n).nFiring>0
        [b, x] = hist(neuron(n).firingTime(1:neuron(n).nFiring),100);
        plot(x,b/(x(2)-x(1)),'Color',colors(mod(n-1,7)+1,:))
        hold on
    end
end
%axis([0 1 .5 PARAM.nNeuron+.5])
%axis 'auto x'

subplot(3,1,2)
maxTime = 0;
for n=1:PARAM.nNeuron
    if neuron(n).nFiring>0
        maxTime = max(maxTime,neuron(n).firingTime(neuron(n).nFiring));
        plot(neuron(n).firingTime(1:neuron(n).nFiring),n*ones(1,neuron(n).nFiring),'.')
        hold on
    end
end
axis([maxTime-10 maxTime .5 PARAM.nNeuron+.5])

subplot(3,1,3)
for n=1:PARAM.nNeuron
    plot(neuron(n).firingTime(1:neuron(n).nFiring),n*ones(1,neuron(n).nFiring),'.')
    hold on
end
axis([0 10 .5 PARAM.nNeuron+.5])

%disp(['Avg firing rate =' num2str(mean([neuron(:).nFiring])/maxTime) ])

figure('Name','RF','Color','white')
clear result
% sigmaC = .15/.1; % cat
% sigmaS = .8/.1; % cat
interp = 1;
sigmaC = interp * .05/.05; % primate
sigmaS = interp * .15/.05; % primate
sz = round(8*sigmaS);
dog = DoG(sz,sigmaC,sigmaS);
N = sqrt(size(weight,1)/2); % RF size

for n=1:PARAM.nNeuron
%     subplot(6,6,n)
    RF = zeros(interp*N);
    for polarite=[-1 1]
        w = reshape(neuron(n).weight((1:N^2)+(polarite+1)/2*N^2),N,N);
        interp_r = zeros(interp*size(w,1),size(w,2));
        interp_r(1:interp:end,:) = w;
        interp_rc = zeros(interp*size(w));
        interp_rc(:,1:interp:end) = interp_r;
        RF = RF + conv2(interp_rc,+polarite*dog,'same');
    end
    c = mod((n-1),6)+1;
    l = floor((n-1)/6)+1;
    result((l-1)*interp*(N+5)+1:(l-1)*interp*(N+5)+interp*N,(c-1)*interp*(N+5)+1:(c-1)*interp*(N+5)+interp*N) =  RF;%RF(end:-1:1,:);
%     imagesc(RF(end:-1:1,:));
%     image((RF(end:-1:1,:)/2.5+.5)*64);
end
imagesc(result);
colormap gray
axis off          % Remove axis ticks and numbers
axis image        % Set aspect ratio to obtain square pixels

co = get(gcf,'DefaultAxesColorOrder');
for n=1:PARAM.nNeuron
    rectangle('Position',[(n-1)*interp*(N+5)+1,1, size(RF)-ones(1,2)],'LineWidth', 5,'EdgeColor',co(n,:))
end

