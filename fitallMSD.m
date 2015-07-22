function [cleanfitresults, tom] = fitallMSD(inputfile, dim, fitend, randomstd, complex, plotoutput)
%-------------------------------------------------------------------------
%
% EXAMPLE:
%       [cleanfitresults, tom] = fitallMSD('Fab2IC_combined_results.mat','2Dxy', 10, 0.01, 0.6, 0);
% PURPOSE:
%       use fit2DMSD.m to go through all MSDs of tracks found in inputfile
%       or to cycle through them one by one via keyboard
%
% INPUTS:
%       inputfile = contains the tracks
%       keyboard  = 0 (no interaction, no plot), 1 (one MSD plot after the
%                   other by hitting a key)
%       tlag      = time lag between frames
%       fitstart  = first step to use for fit (usually 1)
%       fitend = length of MSD to fit to (should be > 5)
%       randomstd  = a standard deviation of all R^2 values for
%                    random, caged and directed is calculated. If the std
%                    is below the value 'randomstd', the msd will be tagged
%                    as 'random' movement.
%       complex    = if none of the Rsq is "better" than this value
%                    call the type of movement "complex"
%
% OUTPUTS:
%       results = structure with quite a lot of info
%
% CALLS:
%       CSVwithheadread.m
%       fit2DMSD.m
%       MSDtype.m
%
% HSTORY:
%       20080605 TOB
% TODO:
%       nothing
%-------------------------------------------------------------------------

%% load data
% [MSDsheadrow, MSDs] = CSVwithheadread(inputfile);
load(inputfile);

if plotoutput == 0
    h = waitbar(0,['Processing track: ', num2str(0),'/',num2str(MSDs(end,1))],...
        'Name',['Processing file: ', inputfile]);
end

for i = 1:MSDs(end,1)

    %% select one track
    MSD = MSDs((MSDs(:,1)==i),:);
    if isempty(MSD) == 0;
        %% fit
        if plotoutput == 0
            waitbar(i/MSDs(end,1),h,['Processing track: ', num2str(i),'/',num2str(MSDs(end,1))]);
            singleresults = fitoneMSD(MSD, dim, fitend, randomstd, complex, 0);
            fitresults{i,1} = singleresults;

        elseif plotoutput == 1
            singleresults = fitoneMSD(MSD, dim, fitend, randomstd, complex, 1);
            fitresults{i,1} = singleresults;
            w = waitforbuttonpress;
        else
        end
    end

end

if plotoutput == 0
    close(h);
end

%% cleanup
%z = 0;
%for i = 1:size(fitresults,1)
%    if isempty(fitresults{i,1}) == 0
%        if strcmp(fitresults{i,1}.movement,'track too short') == 0
%            z = z + 1;
%            cleanfitresults{z,1} = fitresults{i,1};
%        end
%    end
%end
            cleanfitresults = fitresults;



%% evaluate & plot
random = zeros(size(cleanfitresults,1),1); confined = zeros(size(cleanfitresults,1),1); directed = zeros(size(cleanfitresults,1),1); complex = zeros(size(cleanfitresults,1),1);

for i = 1:size(cleanfitresults,1);
    if strcmp(cleanfitresults{i,1}.movement, 'random'); 
        random(i) =  1; 
    end
    if strcmp(cleanfitresults{i,1}.movement, 'caged');                   % changed from 'confined' to 'caged' - FYL 
        confined(i) = 1; 
    end
    if strcmp(cleanfitresults{i,1}.movement, 'directed'); 
        directed(i) =  1; 
    end
    if strcmp(cleanfitresults{i,1}.movement, 'complex'); 
        complex(i) =  1; 
    end
end

tom            = [random,confined,directed,complex];
sumtom=sum(tom);
tom_frac       = (sumtom./sum(sumtom))';

bar(tom_frac,'group');
set(gca,'XTickLabel',{['random  ('  , num2str(sum(random)),   ')'], ...
                      ['confined  (', num2str(sum(confined)), ')'], ...
                      ['directed  (', num2str(sum(directed)), ')'], ...
                      ['complex  (' , num2str(sum(complex)),  ')']}); 


