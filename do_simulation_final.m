function [im] = do_simulation_final(Is,TI,params)
%% Description:
% DS search is used,
% The order of scanning of sim_nodes are determined either by considering their
%    distance from hard data by RANDOMLY sampling from a non-uniform distribution or by
%    random permutations of all nodes,
%
% At each simulation node if the most compatible data event of the ti is
%    not consistent in some neighboring points, those points are removed and
%    added to simulation path again,
% If a hard data is among the inconsistent neighbors, no neighbors are
%    removed and instead the synthesis of current node is postponed.
%% Parameters:
% params.n : maximum number of neighboring points considered in a data-event
% params.search_radius : maximum radius of data-event around simulation node
% params.m : in each m*m window in the SG only one point is simulated using DS

%% other tunable parameters (will be set to defaults if not set)
% params.simul_type: simulation type: 1- Fast CHDS(default) 2- CHDS
%   3- ordinary DS (with ssd search)
% params.hd_w : weights for hard data in match finding (>1)
% params.beta : after simulating beta*N points the removing procedure
%   stops, the remaining sim_nodes are filled by ordinary DS without
%   removing inconsistent points resulting in a fast completion of DS. (N
%   is the total no. of sim-nodes to be filled by DS).
%   By setting
% params.alpha : after simulating alpha.N points the definition of distance
%   is weighted by w=abs(data_event_sim-mean_dev)+1/3*sig_dev;
% params.apply_priority : if one the points are scanned by random sampling
%   from a non-uniform pdf giving more chance to points neat to hard data
% params.hr_type: hierarchy type of sim_nodes:
%   0 : no hierarchy,
%   1 : an interlaced regular subset including 0,25 of sim_nodes are
%   synthesized first, then we continue with the remaining nodes,
%   2 : similar to 1, but the first set of synthesized nodes are considered
%   as hard conditioning data for the remainder of procedure.
% params.disp : shows the procedure of completion if set to one.

% Setting default values:
if ~isfield(params,'simul_type')
    params.simul_type=1;
end
if ~isfield(params,'hd_w')
    params.hd_w = 10;
end
if ~isfield(params,'hr_type')
    params.hr_type = 0;
end
if ~isfield(params,'apply_priority')
    params.apply_priority = 1;
end
if ~isfield(params,'beta')
    params.beta = 2;
end
if ~isfield(params,'alpha')
    params.alpha = 1.2;
end
if ~isfield(params,'disp')
    params.disp = 0;
end

Is=double(Is);
TI=double(TI);

% ind1=find(Is==1);
% ind2=find(Is==0);
% [x1,y1]=ind2sub(size(Is),ind1);
% [x2,y2]=ind2sub(size(Is),ind2);



tic
if size(TI,3)==1
    m = params.m; % simulaiton path consists of one sample in each m*m window
    % padding Is to a multiple of m+ padding a strip of width m at each side
    size_Is_initial = size(Is);
    Is = padarray(Is,ceil(size(Is)/m)*m-size(Is),nan,'post');
    %Is = padarray(Is,[m,m],nan);
    
    %% defining the simulation path:
    
    if (params.simul_type ==1)|(params.simul_type ==2) % Fast CHDS or CHDS
        % forming the path:
        x_path = ceil(m/2):m:size(Is,1);
        y_path = ceil(m/2):m:size(Is,2);
        [x_path,y_path] = meshgrid(x_path,y_path);
        B_1st_path = logical(zeros(size(x_path)));
        B_1st_path ([2:2:end],[2:2:end])= 1;
        Is_col=im2col(Is',[m,m],'distinct');
        if size(Is_col,1)~=1
            B_to_remove = any(~isnan(Is_col));% hard data to be removed from path
        else
            B_to_remove = (~isnan(Is_col));
        end
        B_to_remove = reshape(B_to_remove, size(x_path));
        
        x_path1 = x_path (B_1st_path&(~B_to_remove));
        y_path1 = y_path (B_1st_path&(~B_to_remove));
        x_path2 = x_path ((~B_1st_path)&(~B_to_remove));
        y_path2 = y_path ((~B_1st_path)&(~B_to_remove));
        x_path1 = x_path1(:);
        y_path1 = y_path1(:);
        x_path2 = x_path2(:);
        y_path2 = y_path2(:);
        
        % adding random translations to regularly sampled points:
        % x_path1=x_path1+ceil(rand(size(x_path1))*m)-ceil(m/2);%%%%%%%%%%%%%%%%%%%%
        % y_path1=y_path1+ceil(rand(size(y_path1))*m)-ceil(m/2);
        % x_path2=x_path2+ceil(rand(size(x_path2))*m)-ceil(m/2);%%%%%%%%%%%%%%%%%%%%
        % y_path2=y_path2+ceil(rand(size(y_path2))*m)-ceil(m/2);
        c = mod(m,2)*.5;
        x_path1=x_path1+ceil((rand(size(x_path1))-.5)*m/2-c);%%%%%%%%%%%%%%%%%%%%
        y_path1=y_path1+ceil((rand(size(x_path1))-.5)*m/2-c);
        x_path2=x_path2+ceil((rand(size(x_path2))-.5)*m/2-c);%%%%%%%%%%%%%%%%%%%%
        y_path2=y_path2+ceil((rand(size(x_path2))-.5)*m/2-c);
        if (params.hr_type == 0)
            sim_points = sub2ind(size(Is),[x_path1;x_path2],[y_path1;y_path2]);
            if (params.apply_priority == 0)
                path_sim = sim_points(randperm(numel(sim_points)));
            else
                path_sim = non_uniform_randperm(Is, sim_points, m);
            end
        else
            sim_points1=sub2ind(size(Is),x_path1,y_path1);
            sim_points2=sub2ind(size(Is),x_path2,y_path2);
            if params.apply_priority == 0
                path_sim1 = sim_points1(randperm(numel(sim_points1)));
                path_sim2 = sim_points2(randperm(numel(sim_points2)));
            else
                path_sim1 = non_uniform_randperm(Is, sim_points1, m);
                path_sim2 = non_uniform_randperm(Is, sim_points2, m);
            end
        end
        
        
        %% Densification using Direct Sampling:
        if params.simul_type==2% CHDS
            if params.hr_type == 0
                [simul,cnt]=CHDS(Is,TI,path_sim,params);
            elseif params.hr_type == 1
                [simul,cnt]=DS5(Is,TI,path_sim1,path_sim2,params);
            else
                [simul,cnt1]=CHDS(Is,TI,path_sim1,params);
                [simul,cnt2]=CHDS(simul,TI,path_sim2,params);
                cnt=cnt1+cnt2;
            end
        elseif params.simul_type==1% Fast CHDS
            %if params.hr_type == 0
                [simul,cnt]=f_CHDS(Is,TI,path_sim,params);
            %elseif params.hr_type == 1
                %[simul,cnt]=DS5(Is,TI,path_sim1,path_sim2,params);
            %else
            %    [simul,cnt1]=CHDS(Is,TI,path_sim1,params);
            %    [simul,cnt2]=CHDS(simul,TI,path_sim2,params);
            %    cnt=cnt1+cnt2;
            %end
        end
        fprintf('no. of iterations of coarse simulation phase is %d\n',cnt)
        
        %im_ds=simul;
        path_sim=find(isnan(simul));
        ind=randperm(numel(path_sim));
        path_sim=path_sim(ind);
        im = myDS(simul,TI,path_sim,params);

    elseif params.simul_type == 3 % ordinary DS with ssd search
        path_sim=find(isnan(Is));
        ind=randperm(numel(path_sim));
        path_sim=path_sim(ind);
        [im]=myDS(Is,TI,path_sim,params);
        %im_ds=[];im_idw=[];
    end
    %im=im(m+1:end-m,m+1:end-m);
    im=im(1:size_Is_initial(1),1:size_Is_initial(2));
    
    
    
    fprintf('total time is %f seconds \n',toc)
    %if params.disp
    figure
    imshow(im)
    %     hold on;
    %     plot(y1,x1,'ro','markersize',3)
    %     plot(y2,x2,'bo','markersize',3)
    %end
    
    
else
    m = params.m; % simulaiton path consists of one sample in each m*m window
    % padding Is to a multiple of m+ padding a strip of width m at each side
    size_Is_initial = size(Is);
    Is = padarray(Is,ceil(size(Is)/m)*m-size(Is),nan,'post');
    %Is = padarray(Is,[m,m,m],nan);
    
    %% defining the simulation path:
    
    if (params.simul_type ==1)||(params.simul_type ==2) % Fast CHDS or CHDS
        % forming the path:
        x_path = ceil(m/2):m:size(Is,1);
        y_path = ceil(m/2):m:size(Is,2);
        z_path = ceil(m/2):m:size(Is,3);
        [x_path,y_path,z_path] = meshgrid(x_path,y_path,z_path);
        B_1st_path = false(size(x_path));
        B_1st_path ([2:2:end],[2:2:end],[2:2:end])= 1;
        B_to_remove=false(size(x_path));
        for ii=1:size(Is,1)/m
            for jj=1:size(Is,2)/m
                for kk=1:size(Is,3)/m
                    
                    temp=Is((ii-1)*m+1:ii*m,(jj-1)*m+1:jj*m,(kk-1)*m+1:kk*m);
                    if any(~isnan(temp))
                        B_to_remove(ii,jj,kk)=true;
                    end
                end
            end
        end
        
        x_path1 = x_path (B_1st_path&(~B_to_remove));
        y_path1 = y_path (B_1st_path&(~B_to_remove));
        z_path1 = z_path (B_1st_path&(~B_to_remove));
        x_path2 = x_path ((~B_1st_path)&(~B_to_remove));
        y_path2 = y_path ((~B_1st_path)&(~B_to_remove));
        z_path2 = z_path ((~B_1st_path)&(~B_to_remove));
        x_path1 = x_path1(:);
        y_path1 = y_path1(:);
        z_path1 = z_path1(:);
        x_path2 = x_path2(:);
        y_path2 = y_path2(:);
        z_path2 = z_path2(:);
        % adding random translations to regularly sampled points:
        % x_path1=x_path1+ceil(rand(size(x_path1))*m)-ceil(m/2);%%%%%%%%%%%%%%%%%%%%
        % y_path1=y_path1+ceil(rand(size(y_path1))*m)-ceil(m/2);
        % x_path2=x_path2+ceil(rand(size(x_path2))*m)-ceil(m/2);%%%%%%%%%%%%%%%%%%%%
        % y_path2=y_path2+ceil(rand(size(y_path2))*m)-ceil(m/2);
        c = mod(m,2)*.5;
        x_path1=x_path1+ceil((rand(size(x_path1))-.5)*m/2-c);%%%%%%%%%%%%%%%%%%%%
        y_path1=y_path1+ceil((rand(size(x_path1))-.5)*m/2-c);
        z_path1=z_path1+ceil((rand(size(x_path1))-.5)*m/2-c);
        x_path2=x_path2+ceil((rand(size(x_path2))-.5)*m/2-c);%%%%%%%%%%%%%%%%%%%%
        y_path2=y_path2+ceil((rand(size(x_path2))-.5)*m/2-c);
        z_path2=z_path2+ceil((rand(size(x_path2))-.5)*m/2-c);
        if (params.hr_type == 0)
            sim_points = sub2ind(size(Is),[x_path1;x_path2],[y_path1;y_path2],[z_path1;z_path2]);
            if (params.apply_priority == 0)
                path_sim = sim_points(randperm(numel(sim_points)));
            else
                path_sim = non_uniform_randperm(Is, sim_points, m);
            end
        else
            sim_points1=sub2ind(size(Is),x_path1,y_path1,z_path1);
            sim_points2=sub2ind(size(Is),x_path2,y_path2,z_path2);
            if params.apply_priority == 0
                path_sim1 = sim_points1(randperm(numel(sim_points1)));
                path_sim2 = sim_points2(randperm(numel(sim_points2)));
            else
                path_sim1 = non_uniform_randperm(Is, sim_points1, m);
                path_sim2 = non_uniform_randperm(Is, sim_points2, m);
            end
        end
        
        
        %% Densification using Direct Sampling:
        if params.simul_type==2% CHDS
            if params.hr_type == 0
                [simul,cnt]=CHDS(Is,TI,path_sim,params);
            elseif params.hr_type == 1
                [simul,cnt]=DS5(Is,TI,path_sim1,path_sim2,params);
            else
                [simul,cnt1]=CHDS(Is,TI,path_sim1,params);
                [simul,cnt2]=CHDS(simul,TI,path_sim2,params);
                cnt=cnt1+cnt2;
            end
        elseif params.simul_type==1% Fast CHDS
            %if params.hr_type == 0
                [simul,cnt]=f_CHDS(Is,TI,path_sim,params);
            %elseif params.hr_type == 1
                %[simul,cnt]=DS5(Is,TI,path_sim1,path_sim2,params);
            %else
            %    [simul,cnt1]=CHDS(Is,TI,path_sim1,params);
            %    [simul,cnt2]=CHDS(simul,TI,path_sim2,params);
            %    cnt=cnt1+cnt2;
            %end
        end
        fprintf('no. of iterations of coarse simulation phase is %d\n',cnt)
        
        %im_ds=simul;
        
        path_sim=find(isnan(simul));
        ind=randperm(numel(path_sim));
        path_sim=path_sim(ind);
        im = myDS(simul,TI,path_sim,params);
        %im_idw=[];
    elseif params.simul_type == 3 % ordinary DS with ssd search
        path_sim=find(isnan(Is));
        ind=randperm(numel(path_sim));
        path_sim=path_sim(ind);
        [im]=myDS(Is,TI,path_sim,params);
        %im_ds=[];im_idw=[];
    end
    %im=im(m+1:end-m,m+1:end-m,m+1:end-m);
    im=im(1:size_Is_initial(1),1:size_Is_initial(2),1:size_Is_initial(3));
   
end