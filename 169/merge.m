tra=[];
n=0;
for i=1:20
    filename=sprintf('169-%i.mat.',i);
    load(filename);
    addtra=zeros(size(results,1),4);
    addtra(:,4)=results(:,4)+n;
    addtra(:,1:3)=results(:,1:3);
    tra=[tra;addtra];
    n=tra(end,4);
end