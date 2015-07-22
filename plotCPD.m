function  plotCPD( CPD168, CPD169, sort )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
plot(CPD168{sort,1},CPD168{sort,2},'.r')
hold on
plot(CPD169{sort,1},CPD169{sort,2},'.g')
legend('168KO','169WT','Location','best')
if sort==1
    title('all','FontName','Times New Roman','FontWeight','Bold','FontSize',16)
end
if sort==2
    title('random','FontName','Times New Roman','FontWeight','Bold','FontSize',16)
end
if sort==3
    title('confined','FontName','Times New Roman','FontWeight','Bold','FontSize',16)
end
if sort==4
    title('directed','FontName','Times New Roman','FontWeight','Bold','FontSize',16)
end
xlabel('diffusion rate (¦Ìm^2/s)','FontName','Times New Roman','FontSize',14)
ylabel('Cumulative probability distribution','FontName','Times New Roman','FontSize',14,'Rotation',90)
axis([-6 0 0 1])
set(gca, 'XTick', [-6:1:0])
set(gca,'XTickLabel',{'10^-^6','10^-^5','10^-^4','10^-^3','10^-^2','10^-^1','10^0'})
hold off
end

