function [angles, curvatures]=variables_from_tr4(tr4file,s,Nwhisker)
%%%%%%%%%%%%%%%%%%
% Aim: Extract kinematic variables and curvature from tr4files
% for 2D tracking: Azimuth angle and curvature
% for 3D tracking: Azimuth angle, elevation angle, horizontal curvature,
% vertical curvature, 3D curvature
% 
% Inputs:
% tr4file: Name of the tr4file as 'tr4filename.tr4' or 'tr4filename.tr4_2D'
% s: Evaluate bezier curves at s=s* 0<=s*<=1 . s=0 indicates base of the whisker
% s=1 the tip. 
% Nwhisker:  vector with Number of the whiskers being tracked. If tracking single
% whisker, Nwhisker=1.
%
% Outputs:
%
% for 2D tracking
% angles= Matrix containing horizontal angle for all whiskers. Rows represent
% time and columns different whiskers
% curvatures=Matrix containing  horizontal cuvature for all whiskers. Rows represent
% time and columns different whiskers

%for 3D tracking
% angles= Matrix containing angles for all whiskers. Rows represent
% time, columns different whiskers and third dimension represent (1)
% Horizontal angle (2) Vertical angle.
% curvatures=Matrix containing cuvatures for all whiskers. Rows represent
% time, columns different whiskers and third dimension represent (1)
% Horizontal angle (2) Vertical angle (3) 3D curvature


%% For 2D tracking

if strcmp(tr4file(end-5:end),'tr4_2D')
display('2d tracking')
load(tr4file,'-mat','whisker')
handles.whisker = whisker;
clear whisker

else
display('3d tracking')

load(tr4file,'-mat','whisker','calib')
handles.whisker = whisker;
handles.calib = calib;
clear whisker calib
end


for w = Nwhisker

    idx = find(handles.whisker(w).tracked);
    theta = zeros(3,length(handles.whisker(w).tracked));
    curv_hc = zeros(3,length(handles.whisker(w).tracked));  % head-centred
    curv3 = zeros(1,length(handles.whisker(w).tracked));
    for i = 1:length(idx)
        fr = idx(i);
        r = squeeze(handles.whisker(w).r3all(fr,:,:));
        theta(:,fr) = base_angle3(r,s);
        curv_hc(1,fr) = curvature(r([2 3],:),s);
        curv_hc(2,fr) = curvature(r([1 3],:),s);
        curv_hc(3,fr) = curvature(r([1 2],:),s);
        curv3(fr) = curvature3(r,s);
   end
    clear i fr
    
    azimuth(w,idx) = theta(3,idx);
    elevation(w,idx) = theta(1,idx);
    kappa_h(w,idx)=curv_hc(3,idx)';
    kappa_v(w,idx)=curv_hc(1,idx)';
    kappa_3d(w,idx)=curv3(idx);

    
end

if strcmp(tr4file(end-5:end),'tr4_2D')
   angles=azimuth';
   curvatures=kappa_h'; 
else
angles(:,:,1)=azimuth';
angles(:,:,2)=elevation';

curvatures(:,:,1)=kappa_h';
curvatures(:,:,2)=kappa_v';
curvatures(:,:,3)=kappa_3d';
clear w azimuth elevation

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function theta = base_angle3(r,t)
% angle of tangent vector to Bezier defined by control points {r} at point
% t
% theta(1) is in the x-y plane
% theta(2) is in the y-z plane
% theta(3) is in the x-z plane (not sure meaningful)

dBdt = bezierdtval(r,t);    % tangent vector
theta = zeros(3,length(t));
theta(3) = atan2(-dBdt(2,:),dBdt(1,:))*(180/pi);        % azimuth - whisker pointing caudal is 0'.
theta(1) = 180-atan2(-dBdt(2,:),dBdt(3,:))*(180/pi);    % elevation - whisker pointing ventral is 0'.
theta(2) = atan2(-dBdt(3,:),dBdt(1,:))*(180/pi);
clear dBdt
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kappa = curvature(r,t)

dBdt = bezierdtval(r,t);
d2Bdt2 = bezierdt2val(r,t);
kappa = (dBdt(1,:).*d2Bdt2(2,:)-dBdt(2,:).*d2Bdt2(1,:)) ./ (dBdt(1,:).^2+dBdt(2,:).^2).^(3/2);
clear dBdt d2Bt2
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kappa = curvature3(r,t)
% extension of the kappa formula to 3 dimensions (wikipedia curvature page)
dBdt = bezierdtval(r,t);
d2Bdt2 = bezierdt2val(r,t);

if numel(t)>1
    error('Generalise the code!')
end

kappa = norm(cross(dBdt,d2Bdt2))/norm(dBdt).^3;
% equivalent but more cumbersome formula:
% kappa = sqrt((d2Bdt2(3,:).*dBdt(2,:)-d2Bdt2(2,:).*dBdt(3,:)).^2+(d2Bdt2(1,:).*dBdt(3,:)-d2Bdt2(3,:).*dBdt(1,:)).^2+(d2Bdt2(2,:).*dBdt(1,:)-d2Bdt2(1,:).*dBdt(2,:)).^2)./ ...
%     (dBdt(1,:).^2+dBdt(2,:).^2+dBdt(3,:).^2).^(3/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

function dBdt = bezierdtval(bez,t)
% Evaluate 1st deriv of bezier curve wrt t, with parameters bez(:,1),bez(:,2),... at points t
% size(t) = [1,N], where N is number of points at which to evaluate function.  Typically, t = [0,1]
% size(bez) = [2,order+1]
order = size(bez,2)-1;
switch order
    case 2
        p0 = bez(:,1);
        p1 = bez(:,2);
        p2 = bez(:,3);
        %         dBdt = -p0*2*(1-t) + 2*p1*(1-2*t) + 2*p2*t;
        dBdt = 2*(p0-2*p1+p2)*t + 2*(-p0+p1)*ones(size(t));
    otherwise
        error('Write more code!')
end
clear p0 p1 p2
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dB2dt2 = bezierdt2val(bez,t)
% Evaluate 2nd deriv of bezier curve wrt t, with parameters bez(:,1),bez(:,2),... at points t
% size(t) = [1,N], where N is number of points at which to evaluate function.  Typically, t = [0,1]
% size(bez) = [2,order+1]
order = size(bez,2)-1;
switch order
    case 2
        p0 = bez(:,1);
        p1 = bez(:,2);
        p2 = bez(:,3);
        dB2dt2 = 2*(p0-2*p1+p2) * ones(size(t));
    otherwise
        error('Write more code!')
end
clear p0 p1 p2
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
