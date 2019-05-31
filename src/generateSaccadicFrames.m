% Jan 2015
%
% timothee.masquelier@alum.mit.edu
%
% This code was used in: Masquelier T, Portelli G and Kornprobst P (2016). Microsaccades enable efficient synchrony-based coding in the retina: a simulation study. Scientific Reports. 
%
% Generates frame sequences from ../data/trajectory.mat and a single image.
%
% By definiton, microsaccades in trajectory.mat last one time step.
% The time step in Virtual Retina is 5ms (i.e. 200 frames/s).
% We thus interpolate the microsaccades trajectory x5 such that
% microsaccade duration is 25ms (in the biological range).
%
% The drift is a random walk with time step 5ms. The size of the grid is
% chosen so that the resulting Brownian motion has a diffusion constant of 40 arcmin^2 / s, see Kuang et al Current Biol 2012
%
% We save the new interpolated trajectory in ../data/interpolated_trajectory.mat
%
% n_pool should be <=  number of core in your machine (multithread)
n_pool = 4;

load ../data/trajectory.mat

im = imread('../img/face.jpg'); % change path if neccessary
N_interpolate = 5; % number of frames in the MS sections are multiplicated by this number

D = 40;  % arcmin^2 / s, see Kuang et al Current Biol 2012
dt = 5e-3; % inter-frame interval
pixels_per_degree = 25;
scale = sqrt(2*dt*D)/60*pixels_per_degree;

trajectory(:,1:2) = scale*trajectory(:,1:2); 

count = 1;
meanLum = mean(im(:));


interpolated_trajectory = NaN*ones(15000/dt,3);

for t=1:size(trajectory,1)-1
    if mod(count,round(size(interpolated_trajectory,1)/10))==0
        timedLog(['Interpolating... ' num2str(count/size(interpolated_trajectory,1)) ])
    end

    if trajectory(t,3)==1 % saccade take off
        n_interpolate = N_interpolate;
    else % interpolate
        n_interpolate = 1;
    end
    for n=1:n_interpolate
        shift = (n-1)/n_interpolate*(trajectory(t+1,1:2)-trajectory(t,1:2))+trajectory(t,1:2);
                
        interpolated_trajectory(count,1:2) = shift;
        if n>1 
            interpolated_trajectory(count,3) = NaN;
        else
            interpolated_trajectory(count,3) = trajectory(t,3);
        end
            
        count = count+1;
%         fr(:) = meanLum;
        if count>size(interpolated_trajectory,1)
            break
        end
    end
    if count>size(interpolated_trajectory,1)
        break
    end
end

interpolated_trajectory(count:end,:) = [];

save('../data/interpolated_trajectory.mat','interpolated_trajectory')
%return

%clear trajectory

batch = ceil(size(interpolated_trajectory,1)/n_pool);

matlabpool(n_pool)

parfor b=1:n_pool
    xform = eye(3);
    for t= (b-1)*batch+1 : min(size(interpolated_trajectory,1),b*batch)
        if mod(t-(b-1)*batch,round(batch/10))==0
            timedLog(['Batch ' int2str(b) ' - Frame #' int2str(t) ])
        end
        xform(3,1:2) = interpolated_trajectory(t,[2 1]);        
        fr = imtransform(im, maketform('affine',xform), 'XData',[1 size(im,2)],'YData',[1 size(im,1)],'FillValues',meanLum);

        imwrite(fr,['../img/frame/fr_' sprintf('%07d',t-1) '.jpg' ],'jpg'); % change path if neccessary
    end
end

matlabpool close

% ---------------------------
% Plotting (can be commented)

figure
subplot(1,2,1)
plot(trajectory(:,1),trajectory(:,2),'+-c')
hold on
plot(trajectory(trajectory(:,3)==1,1),trajectory(trajectory(:,3)==1,2),'dk')
plot(trajectory(trajectory(:,3)==-1,1),trajectory(trajectory(:,3)==-1,2),'or')

subplot(1,2,2)
plot(interpolated_trajectory(:,1),interpolated_trajectory(:,2),'+-c')
hold on
plot(interpolated_trajectory(interpolated_trajectory(:,3)==1,1),interpolated_trajectory(interpolated_trajectory(:,3)==1,2),'dk')
plot(interpolated_trajectory(interpolated_trajectory(:,3)==-1,1),interpolated_trajectory(interpolated_trajectory(:,3)==-1,2),'or')
