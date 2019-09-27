% 270919
% Rasmus Petersen
%
% Transpose (rotate by 90 deg) images in video
%
% Cd to the directory where the video files are located

ff = dir('*.avi')

for i = 1:numel(ff)
    disp(['Processing ' ff(i).name ':'])
    vObj = VideoReader([ff(i).name]);
    disp('...Reading...')
    v = read(vObj,[1 inf]);
    nframes = size(v,4);  
    vnew = permute(v,[2 1 3 4]);
    fname = ['trans_' ff(i).name];
    vObjnew = VideoWriter(fname);
    open(vObjnew)
    disp('...Writing...')
    writeVideo(vObjnew,vnew);
    close(vObjnew)    
end

