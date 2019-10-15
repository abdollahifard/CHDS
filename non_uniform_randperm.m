function path_sim = non_uniform_randperm(Is, sim_points,m)
% sim_points are linear indices of points to be permutaed,
% Is includes hard conditioning data, points close to hard data are given
%    more chance to be selected earlier in the simulation path.
% points within a circle with radius 2*m centered at any hard data are given 10
%   times more probability

[x_hd,y_hd] = find(~isnan(Is));
[x_sim , y_sim] = ind2sub(size(Is),sim_points);
if isempty(x_hd)
    path_sim = sim_points(randperm(numel(sim_points)));
else
    [X_hd,X_sim] = meshgrid(x_hd,x_sim);
    [Y_hd,Y_sim] = meshgrid(y_hd,y_sim);
    d = (X_hd - X_sim).^2 + (Y_hd - Y_sim).^2;
    p = zeros(size(d));
    p (d<4*m) = 30;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Priority = sum(p,2)+1;
    path_sim = zeros(size(sim_points));
    for i = 1 : numel(sim_points)
        CP = cumsum (Priority);
        rr = rand*CP(end);
        ind = find (CP>rr , 1, 'first');
        path_sim(i) = sim_points(ind);
        sim_points (ind) = [];
        Priority(ind) = [];
    end
end

