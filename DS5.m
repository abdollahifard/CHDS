function [simul,cnt]=DS5(Is,ti,path_sim1,path_sim2,params)
%% Description:
% DS search is used,
% At each simulation node if the most compatible data event of the ti is
% not consistent in some neighboring points, those points are removed and
% added to simulation path again,
% If a hard data is among the inconsistent neighbors, no neighbors are
% removed and instead the synthesis of current node is postponed.

search_radius=params.search_radius;% 20 %maximum extension of the data events
n = params.n;% 20;                 %maximum number of points in the data event
f = params.f;%0.5;                 %maximum fraction of scanned training image
t = params.t; %0.05;               %distance threshold (between 0 and 1)
distance_type= params.distance_type;%1;%type of variable (0=categorical, 1= continuous with
alpha = params.alpha; %2         % maximum allowed simulation nodes is alpha*numel(path_sim)
beta = params.beta;   %1.2
%Manhattan distance, 2= continuous with euclidean distance)

ind_hd = find(~isnan(Is));
simul = Is;
%ti=load(tifile);
%data = load(datafile);
ti_size = [size(ti,1) size(ti,2) size(ti,3)]; %to make sure size_ti is a vector of size 3
simul_size = [size(simul,1),size(simul,2),size(simul,3)];
sizeyxti = ti_size(1)*ti_size(2);
sizeyxsim = simul_size(1)*simul_size(2);
%tic

%defining path in simulation
%path_sim = randperm(simul_size(1)*simul_size(2)*simul_size(3));

bestmin1 = zeros(size(path_sim1,2),1);
nbtries1 = zeros(size(path_sim1,2),1);
bestmin2 = zeros(size(path_sim2,2),1);
nbtries2 = zeros(size(path_sim2,2),1);
thr = inf;
N1 = numel(path_sim1);
N2 = numel(path_sim2);
N=N1+N2;
cnt=0;
%looping simulation nodes
while (numel(path_sim1)>0)||(numel(path_sim2)>0)
    if numel(path_sim1)>0
        simnod = path_sim1(1);
        path_sim1(1)=[];
        is_path1 = 1;
    else
        simnod = path_sim2(1);
        path_sim2(1)=[];
        is_path1 = 0;
    end
    %finding where we are in the simulation
    [ysim,xsim,zsim] = findcoord(simnod,simul_size(1),simul_size(2));
    
    %finding distance with all previously simulated points
    all_sim_pts = find(isfinite(simul)==1);
    [yt,xt,zt] = findcoord(all_sim_pts,simul_size(1),simul_size(2));
    d = sqrt((ysim-yt).^2+(xsim-xt).^2+(zsim-zt).^2);
    
    %sorting
    [d,s] = sort(d);
    all_sim_pts = all_sim_pts(s);
    
    %checking how many points are within outer search radius
    nb_in_search_radius = sum(d<=search_radius);
    
    %taking the closest n or less
    if nb_in_search_radius < n
        informed_nodes_counter = nb_in_search_radius;
    else
        informed_nodes_counter = n;
    end
    
    if informed_nodes_counter == 0
        vect_to_informed_nodesy = 0;
        vect_to_informed_nodesx = 0;
        vect_to_informed_nodesz = 0;
        data_event_sim = 0;
        data_event_Is = 0;%????????????????????????????????
    else
        %finding vectors to each point and data event
        %for i = 1:informed_nodes_counter
        [yt,xt,zt] = findcoord(all_sim_pts(1:informed_nodes_counter),simul_size(1),simul_size(2));
        vect_to_informed_nodesy = yt-ysim;
        vect_to_informed_nodesx = xt-xsim;
        vect_to_informed_nodesz = zt-zsim;
        
        %indice of the data point in the simul
        id_sim = (sizeyxsim).*(zt-1)+simul_size(1).*(xt-1)+yt;
        data_event_sim = simul(id_sim);
        data_event_Is = Is(id_sim);
        %end
    end
    
    %defining dark zone related to the extension of the data event
    yshiftmax = ti_size(1)-max(vect_to_informed_nodesy);
    yshiftmin = -min(vect_to_informed_nodesy);
    xshiftmax = ti_size(2)-max(vect_to_informed_nodesx);
    xshiftmin = -min(vect_to_informed_nodesx);
    zshiftmax = ti_size(3)-max(vect_to_informed_nodesz);
    zshiftmin = -min(vect_to_informed_nodesz);
    
    %taking care not to go out of the training image
    if yshiftmax > ti_size(1)
        yshiftmax = ti_size(1);
    end
    if yshiftmin < 0
        yshiftmin = 0;
    end
    if xshiftmax > ti_size(2)
        xshiftmax = ti_size(2);
    end
    if xshiftmin < 0
        xshiftmin = 0;
    end
    if zshiftmax > ti_size(3)
        zshiftmax = ti_size(3);
    end
    if zshiftmin < 0
        zshiftmin = 0;
    end
    
    %defining a random path
    window = [(yshiftmax-yshiftmin) (xshiftmax-xshiftmin) (zshiftmax-zshiftmin)];
    window_size = window(1)*window(2)*window(3);
    window_sizeyx = window(1)*window(2);
    path_ti = randperm(window_size);
    
    %scanning ti
    mindist = inf;  %initial best distance is set to inf. Updated with every best distance encountered
    nb_of_tries = ceil(window_size*f); %number of tries in the ti
    weights = (1+(params.hd_w-1)*(~isnan(data_event_Is)));
    if cnt>beta*N
        mean_dev = mean(data_event_sim);
        sig_dev = std(data_event_sim);
        if sig_dev == 0
            w1 = ones(size(data_event_sim));
        else
            w1 = abs(data_event_sim-mean_dev)+1/3*sig_dev;
        end
        weights = w1.*weights;
    end
    weights = weights/sum(weights);
    
    for i = 1:nb_of_tries
        
        %finding coordinates of the point in the random path and shifting
        zti = ceil(path_ti(i)./window_sizeyx);
        id2d = (path_ti(i)-(window_sizeyx.*(zti-1)));
        xti = ceil(id2d./window(1));
        yti = id2d-((xti-1).*window(1));
        
        yti = yti+yshiftmin;
        xti = xti+xshiftmin;
        zti = zti+zshiftmin;
        
        %vector of indices of informed points in the ti
        id = 1+(sizeyxti.*((zti+vect_to_informed_nodesz)-1)+ti_size(1).*((xti+vect_to_informed_nodesx)-1)+(yti+vect_to_informed_nodesy-1));
        data_event_ti = ti(id);
        
        %evaluating the distance between the data event in the simulation and the one in the ti.
        if distance_type == 0              %for categorial variable
            distance = sum(weights.*(data_event_sim~=data_event_ti));  %using sum instead of mean: takes less time
        elseif distance_type == 1          %for continuous variable with euclidean distance
            distance = sum(((data_event_sim-data_event_ti).^2).*weights);
        elseif distance_type == 2          %continuous variable, weighting with inverse square distance
            distance = sum(((data_event_sim-data_event_ti).^2).*weights);
        end
        
        %checking if the distance is under the minimum distance found for thes node in the simulation
        if distance < mindist
            mindist = distance;
            bestpoint = sizeyxti.*(zti-1)+ti_size(1).*(xti-1)+yti;
        end
        %if distance under t, the best point is accepted.
        if mindist <= thr*(i/nb_of_tries)
            break
        end
    end
    
    if is_path1
        nbtries1(simnod) = i;
        bestmin1(simnod) = mindist;
        thr = mean(bestmin1(1:simnod))*t;
    else
        nbtries2(simnod) = i;
        bestmin2(simnod) = mindist;
        thr = mean(bestmin2(1:simnod))*t;
    end
    
    %if enough of the ti has been scanned, take the best point so far (minimum distance)
    %     if mindist<=t%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [yti,xti,zti] = findcoord(bestpoint,size(ti,1),size(ti,2));
    x = xti+vect_to_informed_nodesx;
    y = yti+vect_to_informed_nodesy;
    
    id = size(ti,1).*(x-1)+y;
    data_event_ti = ti(id);
    vx = vect_to_informed_nodesx (data_event_ti~=data_event_sim);
    vy = vect_to_informed_nodesy (data_event_ti~=data_event_sim);
    ind_to_remove = sub2ind(size(simul),ysim+vy, xsim+vx);
    %     if (cnt>4225)&&(~isempty(ind_to_remove))
    %         keyboard
    %     end
    % remove inconsistent points from the simulation grid iff there exists
    % no hard data among them, otherwise dont remove any simulated node and
    % postpone the simulation of current node by putting it in a random
    % position in the path_sim
    if (cnt>alpha*N)
        simul(simnod) = ti(bestpoint);
    else
        if is_path1
            if all(~ismember(ind_to_remove,ind_hd))
                simul(ind_to_remove)=nan;
                for kk = 1:numel(ind_to_remove)
                    rk = ceil(rand*numel(path_sim1));
                    if ~isempty(path_sim1)
                        path_sim1 = [path_sim1(1:rk-1);ind_to_remove(kk);path_sim1(rk:end)];
                    else
                        path_sim1=ind_to_remove;
                    end
                end
                simul(simnod) = ti(bestpoint);
            else
                rk = ceil(rand*numel(path_sim1));
                if ~isempty(path_sim1)
                    path_sim1 = [path_sim1(1:rk-1);simnod;path_sim1(rk:end)];
                else
                    path_sim1=ind;
                end
            end
        else
            if all(~ismember(ind_to_remove,ind_hd))
                simul(ind_to_remove)=nan;
                for kk = 1:numel(ind_to_remove)
                    rk = ceil(rand*numel(path_sim2));
                    if ~isempty(path_sim2)
                        path_sim2 = [path_sim2(1:rk-1);ind_to_remove(kk);path_sim2(rk:end)];
                    else
                        path_sim2=ind_to_remove;
                    end
                end
                simul(simnod) = ti(bestpoint);
            else
                rk = ceil(rand*numel(path_sim2));
                if ~isempty(path_sim2)
                    path_sim2 = [path_sim2(1:rk-1);simnod;path_sim2(rk:end)];
                else
                    path_sim2=simnod;
                end
            end
        end
    end
    
    
    cnt=cnt+1;
    if params.disp
        if mod(cnt,100)==0
            cnt
            imshow(simul+.5);drawnow
        end
    end
    %     end
    
    
end

%evalquality = mean(bestmin);
%evaltries = mean(nbtries);
simul=reshape(simul,simul_size(1),simul_size(2),simul_size(3));
%toc


