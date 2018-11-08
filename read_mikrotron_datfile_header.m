function [ data ] = read_mikrotron_datfile_header( fid )

% rsp 080713 (exactly as in dat2mat)

data.offset = fread(fid,1,'uint32');
data.header = fread(fid,1,'uint32');

data.header_sig = fread(fid,20,'char=>char');
data.record_start = fread(fid,30,'char=>char');
data.camera_name = fread(fid,100,'char=>char');

data.header_sig;
data.record_start;
data.camera_name;

data.camera_man = fread(fid,100,'char=>char');
data.camera_model = fread(fid,100,'char=>char');
data.camera_firmware = fread(fid,100,'char=>char');
data.camera_serial = fread(fid,100,'char=>char');
data.usercomment = fread(fid,1024,'char=>char');

data.hack = fread(fid,2,'char=>char');

data.camera_count = fread(fid,1,'uint32');
data.xoffset = fread(fid,1,'uint32');
data.yoffset = fread(fid,1,'uint32');
data.width = fread(fid,1,'uint32');
data.height = fread(fid,1,'uint32');
data.imagesize = fread(fid,1,'uint32');
data.framerate = fread(fid,1,'uint32') ;     % fps
data.exposuretime = fread(fid,1,'uint32');   % muS
data.dataformat = fread(fid,1,'uint32');


data.bayer = fread(fid,3,'double');
data.gamma = fread(fid,3,'double');
fseek(fid,1672,-1);
data.nframes = fread(fid,1,'uint64');
data.startframe = fread(fid,1,'uint64');
data.triggerframe = fread(fid,1,'uint64');
data.triggertick = fread(fid,1,'uint64');
data.internal = fread(fid,1,'uint64');
data.internal = fread(fid,1,'uint32');
data.imageblitz = fread(fid,4,'uint32');
data.irig = fread(fid,1,'uint32');
data.tickcountfreq = fread(fid,1,'uint64');


end

