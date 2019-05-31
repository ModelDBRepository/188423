% Jan 2015
%
% timothee.masquelier@alum.mit.edu
%
% This code was used in: Masquelier T, Portelli G and Kornprobst P (2016). Microsaccades enable efficient synchrony-based coding in the retina: a simulation study. Scientific Reports. 
%
% This code implements a variation of the model by Engbert et al PNAS 2011
% The main difference is that the drift is a random walk (i.e. U is used for microsaccades, but not to chose the direction of the drift, which is completely random)
% Saves the trajectory in ../data/trajectory.mat, which contains a nx3 natrix
% (i,j,flag), where the flag is 1 for microsaccade take off, -1 for
% microsaccade landing, and NaN for the drift.


clear all
L = 501; % grid size
resolution =20;
i0 = (L-1)/2; % neutral
j0 = (L-1)/2; % neutral
epsilon = 1e-2*resolution^-2; % relaxing factor
lambda = 4; % weight of the potential
hc = 87; % threshold

N_traj = 3e6; % we record the last N_traj steps (the first steps are discarded)
N = 50/epsilon + N_traj; % total number of steps
%N = 1e5 + N_traj; % total number of steps
trajectory = NaN * ones(N_traj,3);

% init
J = repmat(1:L,[L 1]);
I = J';
U = lambda  * resolution^-2 * ( (I-i0-1).^2 + (J-j0-1).^2 );
% figure
% imagesc(U);
% colorbar
% return

H = zeros(size(U));

% up down left right
moves = [[-1 0];[1 0];[0 -1];[0 1]];
%local = zeros(1,size(moves,1));

i=i0;
j=j0;
lastMS = NaN;
isi = [];
orientation = [];
magnitude = [];
count = zeros(L,L);
k=1;

% Use a kernel if you want multiple nodes to sink.
sigma = .1*resolution;
%sigma2 = 2;
x = -round(4*sigma):round(4*sigma);
x = repmat(x,[length(x) 1]);
y = x';
kernel = exp(-(x.^2+y.^2)/(2*sigma^2));
%kernel =  (size(x,1)-1)/2-max(abs(x),abs(y));
%kernel = 1 - tanh( (x.^2+y.^2).^.5 - 4 );
%kernel = 1/sigma*exp(-(x.^2+y.^2)/(2*sigma^2)) - 1/sigma2*exp(-(x.^2+y.^2)/(2*sigma2^2));
%kernel = ones(11,11);
%kernel = kernel((size(kernel,1)+1)/2,:);
%kernel = kernel/sum(kernel(:));

% figure
% imagesc(kernel)
% colorbar
% return

% %figure
% minMag = 10;
% x = -minMag:minMag;
% x = repmat(x,[length(x) 1]);
% y = x';
% kernelInf = ( (x.^2+y.^2).^.5<=10 ) * Inf;

refr = 2; % refractory period for microsaccades

while k<=N;
%     plot(j,L+1-i,'o')
%     hold on
%     pause

    if mod(k,round(N/10))==0
        timedLog(['Iteration #' int2str(k) ])
        %figure;imagesc(H);colorbar
    end
    
%     % up down left right
%     for m=1:size(moves,1)
%         local(m) = (1+0.1*randn) * (  U(i+moves(m,1),j+moves(m,2))+H(i+moves(m,1),j+moves(m,2)) );
%     end
%     candidate = find(local==min(local));
%     choice = candidate( ceil(rand*length(candidate)) );

    % random walk
    choice = floor(4*rand)+1;
    
%     % persistence
%     if exist('choice','var') &&  trajectory(mod(k-2,N_traj)+1,3)~=-1 % saccade landing
%         local(choice) = 2.0*local(choice);
%     end
% 
%     cs = [0 cumsum(local)];
%     choice = find(rand*cs(end)>cs,1,'last');
    
    
    i = i + moves(choice,1);
    j = j + moves(choice,2);
    trajectory(mod(k-1,N_traj)+1,:) = [i-(i0+1),j-(j0+1),NaN];
    k = k+1;


    if exist('kernel','var') % nodes in a neighbourhood sink
        H(i-(size(kernel,1)-1)/2:i+(size(kernel,1)-1)/2,j-(size(kernel,2)-1)/2:j+(size(kernel,2)-1)/2) = ...
        H(i-(size(kernel,1)-1)/2:i+(size(kernel,1)-1)/2,j-(size(kernel,2)-1)/2:j+(size(kernel,2)-1)/2) + kernel;
%         if choice>2
%             H(i-(length(kernel)-1)/2:i+(length(kernel)-1)/2,j) = H(i-(length(kernel)-1)/2:i+(length(kernel)-1)/2,j)+kernel';
%         else
%             H(i,j-(length(kernel)-1)/2:j+(length(kernel)-1)/2) = H(i,j-(length(kernel)-1)/2:j+(length(kernel)-1)/2)+kernel;
%         end
        
    else
        H(i,j) = H(i,j)+1; % local node sinks   
        H(i,j) = H(i,j)/(1-epsilon); % just to cancel the global relaxing step
    end
    
    H = (1-epsilon)*H; % global relaxing step
    
    if H(i,j) >= hc  && (isnan(lastMS) || k-lastMS >= refr )% microsaccade
        
        U1 = 1.25e-3 * lambda * resolution^-4 * ( (I-i).^2 .* (J-j).^2 );
        globalPot = H + U + U1;

        if exist('kernelInf','var')
            globalPot(i-(size(kernelInf,1)-1)/2:i+(size(kernelInf,1)-1)/2,j-(size(kernelInf,2)-1)/2:j+(size(kernelInf,2)-1)/2) = ...
            globalPot(i-(size(kernelInf,1)-1)/2:i+(size(kernelInf,1)-1)/2,j-(size(kernelInf,2)-1)/2:j+(size(kernelInf,2)-1)/2) + kernelInf;
        end

        
        candidate = find(globalPot==min(globalPot(:)));
        choice = candidate( ceil(rand*length(candidate)) );
        [theta,rho] = cart2pol(J(choice)-j,i-I(choice));
        if rho>1 % real saccade
            trajectory(mod(k-2,N_traj)+1,3) = 1; % flag to mark take off
            if k>N-N_traj
                orientation(end+1) = theta;
                magnitude(end+1) = rho;
                if ~isnan(lastMS)
                    isi(end+1) = k-lastMS;
                end
                lastMS = k;
            end
            
            i = I(choice);
            j = J(choice);
            trajectory(mod(k-1,N_traj)+1,:) = [i-(i0+1),j-(j0+1),-1]; % last digit = flag to mark landing
        else
            i = I(choice);
            j = J(choice);
            trajectory(mod(k-1,N_traj)+1,:) = [i-(i0+1),j-(j0+1),NaN]; % false saccade
        end
        k = k+1;
        
        if exist('kernel','var') % nodes in a neighbourhood sink
            H(i-(size(kernel,1)-1)/2:i+(size(kernel,1)-1)/2,j-(size(kernel,2)-1)/2:j+(size(kernel,2)-1)/2)=H(i-(size(kernel,1)-1)/2:i+(size(kernel,1)-1)/2,j-(size(kernel,2)-1)/2:j+(size(kernel,2)-1)/2)+kernel;
        else
            H(i,j) = H(i,j)+1; % local node sinks   
            H(i,j) = H(i,j)/(1-epsilon); % just to cancel the global relaxing step
        end
        
        H = (1-epsilon)*H; % global relaxing step
               
%         plot(trajectory(mod((k-N_traj+1:k)-1,N_traj)+1,2),L+1-trajectory(mod((k-N_traj+1:k)-1,N_traj)+1,1))
%         axis([0 L+1 0 L+1])
%         drawnow
%         
%         % if ~isnan(lastMS)
    end
    %if length(isi)>1e3
%      if true%k>=N-N_traj
%         % count(i,j) = count(i,j)+1;
%         imagesc(H);
%         hold on
%         plot(trajectory(mod((k-1-N_traj+1:(k-1))-1,N_traj)+1,2)+(j0+1),trajectory(mod((k-1-N_traj+1:(k-1))-1,N_traj)+1,1)+(i0+1),'k')
%         plot(trajectory(mod(k-1-1,N_traj)+1,2)+(j0+1),trajectory(mod(k-1-1,N_traj)+1,1)+(i0+1),'ok')
%         %axis([0 L+1 0 L+1])
%         colorbar
%         pause
%      end

end

% figure
% imagesc(count)
% colorbar

% reorder
idx = mod(k-1+(-N_traj:-1),N_traj)+1;
trajectory = trajectory(idx,:);

%_____________________________________________________________________________________________________________________________________
%saving
%save
save('../data/trajectory.mat','trajectory')

disp([int2str(length(isi)+1) ' saccades']);
disp(['mean(isi)=' int2str(mean(isi))])
disp(['min(isi)=' int2str(min(isi))])
disp(['mean(magnitude)=' int2str(mean(magnitude))])
disp(['min(magnitude)=' int2str(min(magnitude))])

disp(['Proportion of <8 = ' num2str(sum(isi<8)/length(isi))])
cut = find( find(trajectory(:,3)==1)>size(trajectory,1)/2 , 1 , 'first' )-1;
disp(['Proportion of <8 (1st half) = ' num2str(sum(isi(1:cut)<8)/cut)])
disp(['Proportion of <8 (2nd half) = ' num2str(sum(isi(cut+1:end)<8)/(length(isi)-cut))])


% disp(['Proportion of consecutive ones (1st half) = ' num2str(sum(isi(1:round(length(isi)/2))==refr)/length(isi)*2)])
% disp(['Proportion of consecutive ones (2nd half) = ' num2str(sum(isi(round(length(isi)/2)+1:end)==refr)/length(isi)*2)])




