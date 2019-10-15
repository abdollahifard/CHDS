function [simul]=myDS(Is,ti,path_sim,params)
%% Description:
% DS search is used,
% At each simulation node if the most compatible data event of the ti is
% not consistent in some neighboring points, those points are removed and
% added to simulation path again,
% If a hard data is among the inconsistent neighbors, no neighbors are
% removed and instead the synthesis of current node is postponed. 


%beta = params.beta;
search_radius=params.search_radius;% 20 %maximum extension of the data events
n = params.n;% 20;                 %maximum number of points in the data event
%alpha = params.alpha; %1.5         % maximum allowed simulation nodes is alpha*numel(path_sim)
%Manhattan distance, 2= continuous with euclidean distance)

ind_hd = find(~isnan(Is));
simul = Is;
ti_size = [size(ti,1) size(ti,2) size(ti,3)]; %to make sure size_ti is a vector of size 3
simul_size = [size(simul,1),size(simul,2),size(simul,3)];

%thr = inf;%%%%%%%%%%%%%%%%%%%%%%???
N = numel(path_sim);
cnt=0;
%looping simulation nodes
t1=toc;
while (numel(path_sim)>0)
    cnt=cnt+1;
    simnod = path_sim(1);
    path_sim(1)=[];
    %finding where we are in the simulation
    [ysim,xsim,zsim] = ind2sub(simul_size,simnod);
    ys=max(ysim-search_radius,1);
    ye=min(ysim+search_radius,simul_size(1));
    xs=max(xsim-search_radius,1);
    xe=min(xsim+search_radius,simul_size(2));
    zs=max(zsim-search_radius,1);
    ze=min(zsim+search_radius,simul_size(3));
    dev0 = simul(ys:ye,xs:xe,zs:ze);
    [iii] = find(~isnan(dev0));
    [y,x,z]=ind2sub(size(dev0),iii);
    y=y+ys-1;
    x=x+xs-1;
    z=z+zs-1;
    lag_mag = sqrt((y-ysim).^2+(x-xsim).^2+(z-zsim).^2);
    in_radius = lag_mag<=search_radius;
    y=y(in_radius);
    x=x(in_radius);
    z=z(in_radius);
    lag_mag=lag_mag(in_radius);
    if numel(x)==0
        y_match=ceil(size(ti,1)*rand);
        x_match=ceil(size(ti,2)*rand);
        z_match=ceil(size(ti,3)*rand);
        simul(ysim,xsim,zsim) = ti(y_match,x_match,z_match);
        continue
    elseif numel(x)>n
        [~,ind]=sort(lag_mag);
        lag_mag = lag_mag(ind(1:n));
        y = y(ind(1:n));
        x = x(ind(1:n));
        z = z(ind(1:n));
    end
    ymm=minmax([y;ysim]'); 
    xmm=minmax([x;xsim]');    
    zmm=minmax([z;zsim]');
    dev = nan(diff(ymm)+1,diff(xmm)+1,diff(zmm)+1);
    ind_dev = sub2ind(size(dev),y-ymm(1)+1,x-xmm(1)+1,z-zmm(1)+1);
    try
    inds_dev =  sub2ind(simul_size,y,x,z);
    catch
        keyboard
    end
    dev(ind_dev) = simul(inds_dev);
    mask = ~isnan(dev);
    dev(~mask)=0;
    
    J=ssd(ti,dev,mask);
    
	if mod(cnt,1000)==0
        t2=toc;
        fprintf('Fine Sim: iter %d, no. dev nodes %d,t_av of last 1000 iters %d\n',...
		cnt,numel(ind_dev),(t2-t1)/1000);
        t1=t2;
	end
    [~,match_ind] = min(J(:));
    [y_match,x_match,z_match] = ind2sub(size(J),match_ind);
    y_match = y_match+ysim-ymm(1);
    x_match = x_match+xsim-xmm(1);
    z_match = z_match+zsim-zmm(1);
    simul(ysim,xsim,zsim) = ti(y_match,x_match,z_match);
    
end

