function [ CPD ] = drawCPD( diffusion, tom)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tom=[sum(tom,2),tom];
diffall=diffusion.*tom(:,1);
diffrand=diffusion.*tom(:,2);
diffconf=diffusion.*tom(:,3);
diffdire=diffusion.*tom(:,4);
diffall=diffall(find(diffall>0));
diffrand=diffrand(find(diffrand>0));
diffconf=diffconf(find(diffconf>0));
diffdire=diffdire(find(diffdire>0));
diffall=sort(diffall);
diffrand=sort(diffrand);
diffconf=sort(diffconf);
diffdire=sort(diffdire);
diffall=log10(diffall);
diffrand=log10(diffrand);
diffconf=log10(diffconf);
diffdire=log10(diffdire);
diffallY=[1:size(diffall,1)]'/size(diffall,1);
diffrandY=[1:size(diffrand,1)]'/size(diffrand,1);
diffconfY=[1:size(diffconf,1)]'/size(diffconf,1);
diffdireY=[1:size(diffdire,1)]'/size(diffdire,1);
CPD=cell(4,2);
CPD{1,1}=diffall;
CPD{2,1}=diffrand;
CPD{3,1}=diffconf;
CPD{4,1}=diffdire;
CPD{1,2}=diffallY;
CPD{2,2}=diffrandY;
CPD{3,2}=diffconfY;
CPD{4,2}=diffdireY;



end

