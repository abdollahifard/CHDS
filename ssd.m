function z = ssd(I, Pattern, mask)
% mask determines location of valid values
% Pattern should not be nanor inf  in invalid points,
if sum(sum(sum(mask~=0)))/numel(mask)>1/15
    if size(I,3)==1
        a2 = filter2(mask, I.^2, 'valid');
        b2 = sum(sum((Pattern.^2).*mask));
        ab = filter2(Pattern.*mask, I, 'valid').*2;
        z = sqrt((a2 - ab + b2)/sum(sum(mask)));
    else
        mask = mask(end:-1:1,end:-1:1,end:-1:1);
        Pattern = Pattern(end:-1:1,end:-1:1,end:-1:1);
        s=size(mask);s(3)=size(mask,3);
        mask=double(mask);
        a2 = convn(I.^2, mask, 'valid');
        b2 = sum(sum(sum((Pattern.^2).*mask)));
        ab = convn(I, Pattern.*mask, 'valid').*2;
        z = sqrt((a2 - ab + b2)/sum(sum(sum(mask))));
    end
else
    if size(I,3)==1
        a2 = filter2(mask, I.^2, 'valid');
        b2 = sum(sum((Pattern.^2).*mask));
        ab = filter2(Pattern.*mask, I, 'valid').*2;
        z = sqrt((a2 - ab + b2)/sum(sum(mask)));
        
        
        
        s=size(mask);
        [dx,dy]=find(mask~=0);
        %dx=dx-1;dy=dy-1;dz=dz-1;
        a2=zeros(size(I));
        ab=zeros(size(I));
        b2 = sum(sum(sum((Pattern.^2).*mask)));
        for i=1:numel(dx)
            J=circshift(I,[-(dx(i)-1),-(dy(i)-1)]);
            a2=a2+J.^2*mask(dx(i),dy(i));
            ab=ab+2*Pattern(dx(i),dy(i))*mask(dx(i),dy(i))*J;
        end
        z = sqrt((a2 - ab + b2)/sum(sum(mask)));
        z = z(1:end-s(1)+1,1:end-s(2)+1);
    else
        s=size(mask);s(3)=size(mask,3);
        [ind]=find(mask~=0);
        [dx,dy,dz]=ind2sub(size(mask),ind);
        %dx=dx-1;dy=dy-1;dz=dz-1;
        a2=zeros(size(I));
        ab=zeros(size(I));
        b2 = sum(sum(sum((Pattern.^2).*mask)));
        for i=1:numel(dx)
            J=circshift(I,[-(dx(i)-1),-(dy(i)-1),-(dz(i)-1)]);
            a2=a2+J.^2*mask(dx(i),dy(i),dz(i));
            ab=ab+2*Pattern(dx(i),dy(i),dz(i))*mask(dx(i),dy(i),dz(i))*J;
        end
        z = sqrt((a2 - ab + b2)/sum(sum(sum(mask))));
        z = z(1:end-s(1)+1,1:end-s(2)+1,1:end-s(3)+1);
    end
    
end

