function [simul,cnt]=f_CHDS(Is,ti,path_sim,params)
JJ=[];


%% Description:
% DS search is used,
% At each simulation node if the most compatible data event of the ti is
% not consistent in some neighboring points, those points are removed and
% added to simulation path again,
% If a hard data is among the inconsistent neighbors, no neighbors are
% removed and instead the synthesis of current node is postponed. 
if numel(unique(ti(:)))>20
    is_continuous = true;
    dr_ti = max(ti(:))-min(ti(:));% ti dynamic range
    t_std = dr_ti/15;
else
    is_continuous = false;
    t_std = 0;
end

alpha = params.alpha;
search_radius=params.search_radius;% 20 %maximum extension of the data events
n = params.n;% 20;                 %maximum number of points in the data event
beta = params.beta; %1.5         % maximum allowed simulation nodes is beta*numel(path_sim)
m = params.m;
%Manhattan distance, 2= continuous with euclidean distance)

ind_hd = find(~isnan(Is));
simul = Is;
simul_fs = zeros(size(Is));%sum of values
simul_fs2 = zeros(size(Is));% sum of squared values
simul_fn = zeros(size(Is));% num of values

%ti_size = [size(ti,1) size(ti,2) size(ti,3)]; %to make sure size_ti is a vector of size 3
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
    [ysim,xsim,zsim] = ind2sub(simul_size,simnod);
    % extracting a box of size (2r+1)^3 (if valid) containing candidate
    %     nodes of dev (r=search_radius):
    ys=max(ysim-search_radius,1);% start
    ye=min(ysim+search_radius,simul_size(1)); % end
    xs=max(xsim-search_radius,1);
    xe=min(xsim+search_radius,simul_size(2));
    zs=max(zsim-search_radius,1);
    ze=min(zsim+search_radius,simul_size(3));
    dev0 = simul(ys:ye,xs:xe,zs:ze);
    % y,x,z will show the coordinates of nodes of final data-event in the SG
    indnn = find(~isnan(dev0));
    [y,x,z] = ind2sub(size(dev0),indnn); 
    y=y+ys-1;% coordinates of candidate dev points in SG
    x=x+xs-1;
    z=z+zs-1;
    % removing points outside the search radius:
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
        fill_neighborhood;
        continue
    elseif numel(x)>n
        % keep n closest nodes:
        [~,ind]=sort(lag_mag);
        %lag_mag = lag_mag(ind(1:n));
        y = y(ind(1:n));
        x = x(ind(1:n));
        z = z(ind(1:n));
    end
    ymm=minmax([y;ysim]'); 
    xmm=minmax([x;xsim]');    
    zmm=minmax([z;zsim]');
    % forming final data event (considering the n and r constraints)
    dev = nan(diff(ymm)+1,diff(xmm)+1,diff(zmm)+1);
    dev_h = dev;
    ind_dev = sub2ind(size(dev),y-ymm(1)+1,x-xmm(1)+1,z-zmm(1)+1);
    ind_dev_sg =  sub2ind(simul_size,y,x,z);
    dev(ind_dev) = simul(ind_dev_sg);
    dev_h(ind_dev) = Is(ind_dev_sg);% dev containing only hard data
    
    
    
    % Weighting strategy: 
    %    1- giving higher weights to hard data
    %    2- if cnt>alpha*N then give higher wieghts to points which their
    %      values is different from the mean of the dev (rare values)
    
    weights = zeros(size(dev));
    bind_dev = ~isnan(dev);% binary index of not-nan values in simul
    bind_dev_h = ~isnan(dev_h);% binary index of not-nan values in Is
    weights(bind_dev)=1;
    weights(bind_dev_h)=params.hd_w;
    dev_vals = dev(bind_dev); 
    %dev_vals_h = dev(bind_dev); % may contain nans
    if cnt>alpha*N
        mean_dev = mean(dev_vals);
        sig_dev = std(dev_vals);
        if sig_dev == 0
            w1 = ones(size(dev_vals));
        else
            w1 = abs(dev_vals-mean_dev)+1/3*sig_dev;
        end
        weights(bind_dev) = w1.*weights(bind_dev);
    end
    %weights = weights/sum(weights(:));
    
    
    dev(~bind_dev)=-100;
    %tic
    
    J=ssd(ti,dev,weights);
	if mod(cnt-1,1000)==0
        t2=toc;
        fprintf('Coarse Sim: iter %d, filled %d out of %d, no. dev nodes %d,t_av of last 1000 iters %d\n',...
		cnt-1, sum(sum(sum(~isnan(simul)))),N,sum(sum(sum(bind_dev))),(t2-t1)/1000);
        t1=t2;
	end
    %[~,match_ind] = min(J(:));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % IT IS VERY IMPORTANT NOT TO CHOOSE THE SAME MATCH PARTICULARLY WHEN
    %    CNT>alpha*n
    % pick one of the best matches randomly:
    mm = min(J(:));
    match_ind = find(J==mm);
    if (numel(match_ind)<3)&&(cnt>alpha*N)
        [~ , match_ind]=sort(J(:));
        match_ind=match_ind(1:3);        
    end   
    match_ind = match_ind(ceil(rand*numel(match_ind)));
    [y_match,x_match,z_match] = ind2sub(size(J),match_ind);
    match = ti(y_match:y_match+size(dev,1)-1,x_match:x_match+size(dev,2)-1 ...
        ,z_match:z_match+size(dev,3)-1);
    
    
    y_match = y_match+ysim-ymm(1);
    x_match = x_match+xsim-xmm(1);
    z_match = z_match+zsim-zmm(1);
    
    % remove inconsistent points from the simulation grid iff there exists
    % no hard data among them, otherwise dont remove any simulated node and
    % postpone the simulation of current node by putting it in a random
    % position in the path_sim
    if is_continuous
        [indr] = find((abs(match-dev)>dr_ti/5)&bind_dev);
    else
        [indr] = find((match~=dev)&bind_dev);
    end
    [yr,xr,zr]= ind2sub(size(dev),indr);
    yr=ymm(1)+yr-1;
    xr=xmm(1)+xr-1;
    zr=zmm(1)+zr-1;
    ind_to_remove = sub2ind(simul_size,yr,xr,zr);
%     if ~isempty(ind_to_remove)
%         keyboard
%     end
    if cnt<beta*N
        if all(~ismember(ind_to_remove,ind_hd))
            simul(ind_to_remove)=nan;
            for kk = 1:numel(ind_to_remove)
                rk = ceil(rand*numel(path_sim));
                if ~isempty(path_sim)
                    path_sim = [path_sim(1:rk-1);ind_to_remove(kk);path_sim(rk:end)];
                else
                    path_sim=ind_to_remove;
                end
            end
            simul(simnod) = ti(y_match,x_match,z_match);
            fill_neighborhood;
        else
            rk = ceil(rand*numel(path_sim));
            if ~isempty(path_sim)
                path_sim = [path_sim(1:rk-1);simnod;path_sim(rk:end)];
            else
                path_sim=simnod;
            end
        end
    else
        simul(simnod) = ti(y_match,x_match,z_match);
        fill_neighborhood;
    end
    if params.disp
        if mod(cnt,5)==0
            cnt
            if is_continuous
                imshow(uint8(simul));drawnow
            else
%                 imshow(uint8((simul+1)*70));drawnow
                II=simul*128+127;
                II(isnan(II))=0;
                imshow(uint8(II));drawnow
                JJ(:,:,end+1)=uint8(II);
                
            end
        end
    end    
end


simul_std = sqrt(simul_fs2./(simul_fn-1)-simul_fs.^2./((simul_fn-1).*(simul_fn)));
B=(simul_std<=t_std)&(isnan(simul));
simul(B) = simul_fs(B)./simul_fn(B); 




function fill_neighborhood
    rr = floor(m);
    yr1 = min([ysim-1,y_match-1,rr]);
    xr1 = min([xsim-1,x_match-1,rr]);
    yr2 = min([simul_size(1)-ysim,size(ti,1)-y_match,rr]);
    xr2 = min([simul_size(2)-xsim,size(ti,2)-x_match,rr]);
%     if size(ti,3)==1
%         zr1=0;zr2=0;
%     else
        zr1 = min([zsim-1,z_match-1,rr]);
        zr2 = min([simul_size(3)-zsim,size(ti,3)-z_match,rr]);
%     end
    simul_fs(ysim-yr1:ysim+yr2,xsim-xr1:xsim+xr2,zsim-zr1:zsim+zr2)=...
        simul_fs(ysim-yr1:ysim+yr2,xsim-xr1:xsim+xr2,zsim-zr1:zsim+zr2)+...
        ti(y_match-yr1:y_match+yr2,x_match-xr1:x_match+xr2,z_match-zr1:z_match+zr2);
    simul_fs2(ysim-yr1:ysim+yr2,xsim-xr1:xsim+xr2,zsim-zr1:zsim+zr2)=...
        simul_fs2(ysim-yr1:ysim+yr2,xsim-xr1:xsim+xr2,zsim-zr1:zsim+zr2)+...
        (ti(y_match-yr1:y_match+yr2,x_match-xr1:x_match+xr2,z_match-zr1:z_match+zr2)).^2;
    simul_fn(ysim-yr1:ysim+yr2,xsim-xr1:xsim+xr2,zsim-zr1:zsim+zr2)=...
        simul_fn(ysim-yr1:ysim+yr2,xsim-xr1:xsim+xr2,zsim-zr1:zsim+zr2)+1;

end
end



