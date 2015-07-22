function  [drawMSD]= drawMSD( MSD, n )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
MSDbyt=cell(n,1);
drawMSD=zeros(n,3);
for i=1:n
    MSDbyt{i}=[];
    for j=1:size(MSD,2)
        if size(MSD{j},1)>=i
            MSDbyt{i}=[MSDbyt{i};MSD{j}(i)];
        end
    end
    drawMSD(i,1)=mean(MSDbyt{i});
    drawMSD(i,2)=std(MSDbyt{i})/sqrt(size(MSDbyt{i},1));
    drawMSD(i,3)=size(MSDbyt{i},1);
   
end          
end

