%% Sort Tracks according to Diffusion Mode
% INPUT: 
%           Multiple .mat results from processim.m
%           Parameters
% FUNCTION: 
%           Calculates MSD of all the tracks in the .mat results
%           Sorts the tracks according to their diffusion mode:
%           Directed, Random, Confined, and Complex
% CALLS:
%           fitallMSD, fitoneMSD, msd_directed_2Dfun, msd_random_2Dfun,
%           msd_caged_2Dfun, MSDalldim
% OUTPUT:
%           Saves a results file containing the output from fitallMSD, and
%           the number of tracks sorted into each mode of diffusion
% HISTORY:
%           20130303 FYL combined the functions of wanli_msd_analysis.m,
%           fitallMSD.m, and fitoneMSD.m
%           20130307 FYL fixed bug in parameter input

%% Load files
clear all;

[filename pathname] = uigetfile('*.mat' , 'Please select the result file(s) to analyze.' , 'MultiSelect', 'on');
n_file = size(filename, 2);

% All tracks from the loaded files would be put into [tracks]
tracks = [];

for i = 1 : n_file
    fullpath = [pathname cell2mat(filename(i))];
    temp = load(fullpath);
    tracks = [tracks ; temp.tracks];
end

% [tracks] is a 4-column matrix, col 1 & col 2 are xy coordinates, col 3 is
% the frame number (not reset), col 4 is the track number (not renumbered)

%% Renumber tracks (from wanli_msd.analysis.m)
% Placing the renumbered indices in col 5
z = 1;
tracks(1 , 5) = 1;
for i = 1 : size(tracks , 1) - 1
    if tracks(i + 1 , 4) == tracks(i , 4)
        tracks(i + 1 , 5) = z;
    else
        z = z + 1;
        tracks(i + 1 , 5) = z;
    end
end

%% Reset the frame number for each track (from wanli_msd.analysis.m)
% Placing the reset frame numbers in col 6
for i = 1 : tracks(end , 5)
    % Select one track
    track = tracks((tracks(: , 5) == i) , :);
    
    % Reset the frame number (subtracted by first frame number in the
    % track)
    track(: , 6) = track(: , 3) - track(1 , 3) + 1;
    
    % Put vales back next to selected track in col 6
    tracks((tracks(: , 5) == i) , 6) = track(: , 6); 
end

%% Calculate SDs (from wanli_msd.analysis.m)
% Placing the SDs in col 7
for trn = 1 : tracks(end , 5)
    % Select one track
    track = tracks((tracks(: , 5) == trn) , :);

    % Pre-allocate some variables
    sd2Dxy  = zeros(size(track , 1) - 1 , 1);

    % calcuate sd
    for j = 1 : size(track , 1)-1
        Xdisp = (track(j + 1 , 1) - track(j , 1)) .^ 2;
        Ydisp = (track(j + 1 , 2) - track(j , 2)) .^ 2;
        sd2Dxy(j , 1) = Xdisp + Ydisp;
    end

    % Put SD back into tracks file
    tracks((tracks(: , 5) == trn) , 7) = [NaN ; sd2Dxy];
end

%% Calculated MSDs (from wanli_msd.analysis.m)

% Pre-allocate some variables
MSDs = [];

for i = 1 : 10
    SDs{1 , i}  = [];
end

for trn = 1 : tracks(end,5)

    % Select one track and  reshaped to [frame# Xpos Ypos Zpos(0)]
    track     = tracks((tracks(: , 5) == trn) , :);
    track2Dxy = [track(:,6), track(:,1), track(:,2), zeros(size(track,1),1)];

    % Calculate msd
    [msd2Dxy , SDarray3D , SDarray2D , SDarray1D] = MSDalldim(track2Dxy , 3);
    
    %             track#                             frame#          nSD              MDS         dMSD
    msdoutput  = [trn * ones(size(msd2Dxy , 1) , 1), msd2Dxy(: , 1), msd2Dxy(: , 4) , msd2Dxy(:,2 : 3)];
    MSDs = [MSDs; msdoutput];

    % Output SDarrays for 1-10 steps
    for i = 1 : 10
        if size(SDarray2D , 2) > i - 1
            SDs{1 , i} = [SDs{1 , i} ; SDarray2D{1 , i}];
            SDs{1 , i} = sort(SDs{1 , i});
        end
    end

end

for i = 1 : 10
    SDs{1 , i}(: , 2) = linspace(0 , 1 , size(SDs{1 , i} , 1));
end

% Saves file to temp.mat
save('temp.mat');

%% Fit to diffusion models (calling fitallMSD.m)
% Ask whether you want to see each track's MSD plotted
plotoutput = questdlg('Would you like to see each track''s MSD plotted?','MSD Plots','Yes','No', 'No');
plotoutput = strcmp(plotoutput , 'Yes');

defaultparam = questdlg('Would you like to use default parameters for sorting?','Sorting Parameter','Yes','No', 'Yes');

%       Info on fitallMSD.m:
%       fitallMSD(inputfile, dim, fitend, randomstd, complex, plotoutput)
%       fitend     = length of MSD to fit to (should be > 5)
%       randomstd  = a standard deviation of all R^2 values for
%                    random, caged and directed is calculated. If the std
%                    is below the value 'randomstd', the msd will be tagged
%                    as 'random' movement.
%       complex    = if none of the Rsq is "better" than this value
%                    call the type of movement "complex"

if strcmp(defaultparam , 'Yes') == 1
    [cleanfitresults, tom] = fitallMSD('temp.mat','2Dxy', 10, 0.01, 0.6, plotoutput);
    fitting_params = {10 , 0.01 , 0.6};
else
    prompt = {'Enter length of frames to fit to (should be >5):' , 'Enter cutoff value for std of all R^2 values:' , 'Enter cutoff value for R^2 of complex movement:'};
    dlg_title = 'Custom Parameters';
    num_lines = 1;
    def = {'10' , '0.01' , '0.6'};
    fitting_params = inputdlg(prompt , dlg_title , num_lines , def);
    [cleanfitresults, tom] = fitallMSD('temp.mat' , '2Dxy', str2num(cell2mat(fitting_params(1))) , str2num(cell2mat(fitting_params(2))) , str2num(cell2mat(fitting_params(3))) , plotoutput);
end

% tom = [random confined directed complex];

%% Obtain Mean MSD for each group
random = 0; confined = 0; directed = 0; complex = 0;
MSDrandom = cell(1); MSDconfined = cell(1); MSDdirected = cell(1); MSDcomplex = cell(1);

for i = 1:size(cleanfitresults,1);
    if strcmp(cleanfitresults{i,1}.movement, 'random'); 
        random = random + 1; 
        MSDrandom{random} = cleanfitresults{i,1}.originalMSD;
    end
    if strcmp(cleanfitresults{i,1}.movement, 'caged');
        confined = confined + 1; 
        MSDconfined{confined} = cleanfitresults{i,1}.originalMSD;
    end
    if strcmp(cleanfitresults{i,1}.movement, 'directed'); 
        directed = directed + 1; 
        MSDdirected{directed} = cleanfitresults{i,1}.originalMSD;
    end
    if strcmp(cleanfitresults{i,1}.movement, 'complex'); 
        complex = complex + 1; 
        MSDcomplex{complex} = cleanfitresults{i,1}.originalMSD;
    end
end

%% Save to file
variables = {'tracks' , 'filename' , 'MSDs' , 'fitting_params' , 'cleanfitresults' , 'tom' , 'MSDrandom' , 'MSDconfined' , 'MSDdirected' , 'MSDcomplex'};
uisave(variables);
    

