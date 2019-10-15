function myViewGridCategorical(I)
hold on;axis off
% ind=find(I==0);
% [y,x,z]=ind2sub(size(I),ind);
% plot3(y ,x ,z ,'.','color',[1,1,0]);
Q = unique(I);
n=numel(Q);


ind=find(I==Q(2));
[y,x,z]=ind2sub(size(I),ind);
plot3(y ,x ,z ,'.','color',[0,0,1]);

if n>=3
    ind=find(I==Q(3));
    [y,x,z]=ind2sub(size(I),ind);
    plot3(y ,x ,z ,'.','color',[0,1,0]);
end

if n>=4
    ind=find(I==Q(4));
    [y,x,z]=ind2sub(size(I),ind);
    plot3(y ,x ,z ,'.','color',[1,0,0]);
end

if n>=5
    ind=find(I==Q(5));
    [y,x,z]=ind2sub(size(I),ind);
    plot3(y ,x ,z ,'.','color',[0,1,1]);
end

if n>=6
    ind=find(I==Q(6));
    [y,x,z]=ind2sub(size(I),ind);
    plot3(y ,x ,z ,'.','color',[1,0.5,0]);
end
if n>=7
    ind=find(I==Q(7));
    [y,x,z]=ind2sub(size(I),ind);
    plot3(y ,x ,z ,'.','color',[1,0,1]);
end
if n>=8
    ind=find(I==Q(8));
    [y,x,z]=ind2sub(size(I),ind);
    plot3(y ,x ,z ,'.','color',[0,0,0]);
end
if n>=9
    ind=find(I==Q(9));
    [y,x,z]=ind2sub(size(I),ind);
    plot3(y ,x ,z ,'.','color',[0.5,0.5,0.5]);
end
if n>=10
    ind=find(I==Q(10));
    [y,x,z]=ind2sub(size(I),ind);
    plot3(y ,x ,z ,'.','color',[0.25,0.25,0.25]);
end
