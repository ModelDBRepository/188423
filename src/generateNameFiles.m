% Jan 2015
%
% timothee.masquelier@alum.mit.edu
%
% This code was used in: Masquelier T, Portelli G and Kornprobst P (2016). Microsaccades enable efficient synchrony-based coding in the retina: a simulation study. Scientific Reports. 
% 
% Generates n_file files with list of frames to be processed by Virtual
% Retina. This is useful to launch several threads of Virtual Retina (if multiple cores are available).

load ../data/interpolated_trajectory.mat

n_file = 8;
n_frame = length(interpolated_trajectory);
batch = ceil(n_frame/n_file);


for n=1:n_file
    n
    fid = fopen(['../img/frame/file_name_' int2str(n-1) '.txt'],'w');
    for i= (n-1)*batch+1 : min(n_frame,n*batch)
        fprintf(fid,'../img/frame/fr_%07d.jpg\n',i-1);    
    end
    fclose(fid);
end
