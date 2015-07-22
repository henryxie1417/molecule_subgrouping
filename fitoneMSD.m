function results = fitoneMSD(onemsd, dim, fitend, randomstd, complex, plotoutput)
%-------------------------------------------------------------------------
%
% EXAMPLE:
%       results = fitoneMSD('NKLFasLegfp_fixed_allMSDs.csv', '3D', 0, 0.015, 0.33, 0);
%
% PURPOSE:
%       fit three models for diffusion (random, caged, directed) to
%       a msd and calculate some "goodness of fit" statistics
% INPUTS:
%       onemsd     = one track of te MSDs file
%       dim        = dimension to fit: 3D, 2Dxy, 2Dxz, or 2Dyz
%       fitend     = end position of msd to fit, chose 0 to fit whole msd
%       randomstd  = a standard deviation of all R^2 values for
%                    random, caged and directed is calculated. If the std
%                    is below the value 'randomstd', the msd will be tagged
%                    as 'random' movement. How to define the randomstd? we
%                    set randomstd as 0.015.
%       complex    = if none of the Rsq is "better" than this value
%                    call the type of movement "complex", we set it as
%                    0.33.
%       plotoutput = 0 = no plot, 1 = plot
%
% OUTPUTS:
%       results   = structure with quite a lot of info
%
% CALLS:
%       nothing
%
% HISTORY:
%       20080603 TOB
%       20080821 Dongfang corrected the 3D_caged fitting (line 152). The
%       value of results.caged_Dcage = x(2) and results.caged_Dves
%       =x(3).
%
% TODO:
%       - put parameter guesses in inputs?
%       - add weighting acording to dMSDs
%
%-------------------------------------------------------------------------

%% check size of trajectory
if size(onemsd,1) < fitend
    results.movement = 'track too short';
else

    %% shape input data
    if fitend == 0; 
        xdata = onemsd(1:end,2);
    else
        xdata = onemsd(1:fitend,2);
    end
    results.originaltlag = onemsd(:,2);
    
    results.dim = dim;
    if strcmp(dim, '3D') == 1
        if fitend == 0
            ydata = onemsd(1:end,6);
        else
            ydata = onemsd(1:fitend,6);
        end
        results.originalMSD = onemsd(:,6);
        
    elseif strcmp(dim, '2Dxy') == 1
        if fitend == 0
            ydata = onemsd(1:end,4);
        else
            ydata = onemsd(1:fitend,4);
        end
        results.originalMSD = onemsd(:,4);
        
    elseif strcmp(dim, '2Dxz') == 1
        if fitend == 0
            ydata = onemsd(1:end,10);
        else
            ydata = onemsd(1:fitend,10);
        end
        results.originalMSD = onemsd(:,10);
        
    elseif strcmp(dim, '2Dyz') == 1
        if fitend == 0
            ydata = onemsd(1:end,12);
        else
            ydata = onemsd(1:fitend,12);
        end
        results.originalMSD = onemsd(:,12);
        
    else
        error(['Dimension ', dim, ' not defined. Choose 3D, 2Dxy, 2Dxz, or 2Dyz.']);
    end


%% fit to random model
    if strcmp(dim, '2Dxy') == 1 || strcmp(dim, '2Dxz') == 1 || strcmp(dim, '2Dyz') == 1
        param_min  = [0,    0    ];
        param_init = [0.001,0.001];
        param_max  = [inf,  inf  ];
        options     = optimset('Maxfunevals',10000,'TolFun',0.000000001, 'TolX',0.00000001,'Display','off');
        [x,resnorm,residuals,exitflag,output] = lsqcurvefit(@msd_random_2Dfun,param_init,xdata,ydata, param_min, param_max, options);

        results.random_fit      = msd_random_2Dfun(x,xdata);
        results.random_offset   = x(1);
        results.random_D        = x(2);
        results.random_exitflag = exitflag;
        results.random_output   = output;
    else
        param_min  = [0,    0    ];
        param_init = [0.001,0.001];
        param_max  = [inf,  inf  ];
        options     = optimset('Maxfunevals',10000,'TolFun',0.000000001, 'TolX',0.00000001,'Display','off');
        [x,resnorm,residuals,exitflag,output] = lsqcurvefit(@msd_random_3Dfun,param_init,xdata,ydata, param_min, param_max, options);

        results.random_fit      = msd_random_3Dfun(x,xdata);
        results.random_offset   = x(1);
        results.random_D        = x(2);
        results.random_exitflag = exitflag;
        results.random_output   = output;
    end

    % P-value R^2 and chi^2
    [H,results.random_P,KSSTAT] = kstest2(ydata,results.random_fit);

    SSE                       = norm(residuals)^2;           % sum of squared errors
    SST                       = norm(ydata - mean(ydata))^2; % total sum of squares;
    results.random_Rsq        = 1 - SSE ./ SST;              % see regstats for even more statistics

    results.random_chi  = sum(((ydata - results.random_fit).^2)./ydata);
    % see regstats and Residual Analysis for even more statistics


%% fit to caged model
    if strcmp(dim, '2Dxy') == 1 || strcmp(dim, '2Dxz') == 1 || strcmp(dim, '2Dyz') == 1
        param_min  = [0,0,0];
        param_init = [0.001,0.001,0.3];
        param_max  = [inf,inf,inf];
        options     = optimset('Maxfunevals',10000,'TolFun',0.000000001, 'TolX',0.00000001,'Display','off');
        [x,resnorm,residuals,exitflag,output] = lsqcurvefit(@msd_caged_2Dfun,param_init,xdata,ydata, param_min, param_max, options);

        results.caged_fit      = msd_caged_2Dfun(x,xdata);
        results.caged_offset   = x(1);
        results.caged_D        = x(2);
        results.caged_L        = x(3);
        results.caged_exitflag = exitflag;
        results.caged_output   = output;
    else
        param_min  = [0,    0,    0,    0  ];
        param_init = [0.001,0.001,0.001,0.3];
        param_max  = [inf,  inf,  inf,  inf];
        options     = optimset('Maxfunevals',10000,'TolFun',0.000000001, 'TolX',0.00000001,'Display','off');
        [x,resnorm,residuals,exitflag,output] = lsqcurvefit(@msd_caged_3Dfun,param_init,xdata,ydata, param_min, param_max, options);

        results.caged_fit      = msd_caged_3Dfun(x,xdata);
        results.caged_offset   = x(1);
        %results.caged_Dves     = x(2); results.caged_D = x(2);
        results.caged_Dcage    = x(2);
        results.caged_Dves     = x(3);
        results.caged_R        = x(4);
        results.caged_exitflag = exitflag;
        results.caged_output   = output;
    end

    % P-value R^2 and chi^2
    [H,results.caged_P,KSSTAT] = kstest2(ydata,results.caged_fit);

    SSE                       = norm(residuals)^2;           % sum of squared errors
    SST                       = norm(ydata - mean(ydata))^2; % total sum of squares;
    results.caged_Rsq   = 1 - SSE ./ SST;              % see regstats for even more statistics

    results.caged_chi   = sum(((ydata - results.caged_fit).^2)./ydata);


%% fit to directed model
    if strcmp(dim, '2Dxy') == 1 || strcmp(dim, '2Dxz') == 1 || strcmp(dim, '2Dyz') == 1
        param_min  = [0,0,0];
        param_init = [0.001,0.001,0.01];
        param_max  = [inf,inf,inf];
        options     = optimset('Maxfunevals',10000,'TolFun',0.000000001, 'TolX',0.00000001,'Display','off');
        [x,resnorm,residuals,exitflag,output] = lsqcurvefit(@msd_directed_2Dfun,param_init,xdata,ydata, param_min, param_max, options);

        results.directed_fit      = msd_directed_2Dfun(x,xdata);
        results.directed_offset   = x(1);
        results.directed_D        = x(2);
        results.directed_v        = x(3);
        results.directed_exitflag = exitflag;
        results.directed_output   = output;
    else
        param_min  = [0,0,0];
        param_init = [0.001,0.001,0.01];
        param_max  = [inf,inf,inf];
        options     = optimset('Maxfunevals',10000,'TolFun',0.000000001, 'TolX',0.00000001,'Display','off');
        [x,resnorm,residuals,exitflag,output] = lsqcurvefit(@msd_directed_3Dfun,param_init,xdata,ydata, param_min, param_max, options);

        results.directed_fit      = msd_directed_3Dfun(x,xdata);
        results.directed_offset   = x(1);
        results.directed_D        = x(2);
        results.directed_v        = x(3);
        results.directed_exitflag = exitflag;
        results.directed_output   = output;
    end


    % P-value R^2 and chi^2
    [H,results.directed_P,KSSTAT] = kstest2(ydata,results.directed_fit);

    SSE                       = norm(residuals)^2;           % sum of squared errors
    SST                       = norm(ydata - mean(ydata))^2; % total sum of squares;
    results.directed_Rsq   = 1 - SSE ./ SST;              % see regstats for even more statistics

    results.directed_chi  = sum(  ((ydata - results.directed_fit).^2)./ydata);

%% type of movement (my version - under developement :-)

    % round Rsq values to two decimals
    random_Rsq   = round(results.random_Rsq.*100)/100;
    caged_Rsq    = round(results.caged_Rsq.*100)/100;
    directed_Rsq = round(results.directed_Rsq.*100)/100;

    % if standard deviation for all Rsq is below the value randomstd
    % decide to call the type of movement 'random'
    % Note: the smaller this is set, the less tracks will be 'random'
    if max(ydata) < 0.2 % if the max MSD is below 0.2, then assign to immboile movement. 
        results.movement = 'immobile';
    elseif random_Rsq < complex && caged_Rsq < complex && directed_Rsq < complex
        results.movement = 'complex';
    elseif std([random_Rsq caged_Rsq directed_Rsq]) < randomstd
        results.movement = 'random';
    elseif caged_Rsq > random_Rsq && caged_Rsq > directed_Rsq
        results.movement = 'caged';
    elseif directed_Rsq > random_Rsq && directed_Rsq > caged_Rsq
        results.movement = 'directed';
    else
        results.movement = 'something is wrong';
    end


    
%% type of movement (Dongfangs original version)

    %     if ((results.caged_L)^2 < 0) && (results.directed_v > 0);
    %         if results.random_P > results.directed_P
    %             results.movement = 'random';
    %         end
    %         if results.directed_P > results.random_P
    %             results.movement = 'directed';
    %         end
    %         if results.directed_P == results.random_P
    %             if results.random_chi < results.directed_chi
    %                 results.movement = 'random';
    %             end
    %             if results.directed_chi < results.random_chi
    %                 results.movement = 'directed';
    %             end
    %         end
    %
    %     elseif ((results.directed_v)^2 < 0) && (results.caged_L > 0);
    %         if results.random_P > results.caged_P
    %             results.movement = 'random';
    %         end
    %         if results.caged_P > results.random_P
    %             results.movement = 'caged';
    %         end
    %         if results.caged_P == results.random_P
    %             if results.random_chi < results.caged_chi
    %                 results.movement = 'random';
    %             end
    %             if results.caged_chi < results.random_chi
    %                 results.movement = 'caged';
    %             end
    %         end
    %
    %
    %     elseif ((results.directed_v)^2 < 0) && ((results.caged_L)^2 < 0)
    %         results.movement = 'random';
    %
    %     elseif results.caged_L > 0 && results.directed_v > 0
    %         if results.caged_P > results.directed_P && results.caged_P > results.random_P
    %             results.movement = 'caged';
    %         end
    %         if results.directed_P > results.caged_P && results.directed_P > results.random_P
    %             results.movement = 'directed';
    %         end
    %         if results.random_P > results.directed_P && results.random_P > results.caged_P
    %             results.movement = 'random';
    %         end
    %
    %         if (results.random_P == results.directed_P) &&...
    %                 (results.random_P == results.caged_P) &&...
    %                 (results.directed_P == results.caged_P)
    %             if results.caged_chi < results.random_chi && results.caged_chi < results.directed_chi
    %                 results.movement = 'caged';
    %             end
    %             if results.random_chi < results.caged_chi && results.random_chi < results.directed_chi
    %                 results.movement = 'random';
    %             end
    %             if results.directed_chi < results.random_chi && results.directed_chi < results.caged_chi
    %                 results.movement = 'directed';
    %             end
    %         end
    %     end

    % in case no assignemenyt was done for 'movement'
    if isfield(results, 'movement') == 0
        results.movement = 'something is wrong here';
    end


    %% plot
    if plotoutput == 1
        figure(1);
        clf;

        subplot(1,3,[1:2])
        plot(results.originaltlag, results.originalMSD,'ko-'); hold on;
        plot(xdata, results.random_fit,'r-');   hold on;
        plot(xdata, results.caged_fit,'g:'); hold on;
        plot(xdata, results.directed_fit,'b-.'); hold on;
        axis([0 25*results.originaltlag(1) 0 50])
        title(['track =', num2str(onemsd(1,2)), '   movement = ', results.movement]);
        legend('data','random','caged','directed','Location','best');

        subplot(1,3,3)
        info = {...
            sprintf('D')...
            ,sprintf('%5.5f \t random',results.random_D)...
            ,sprintf('%5.5f \t caged',results.caged_D)...
            ,sprintf('%5.5f \t directed',results.directed_D)...
            ,sprintf(' ')...
            ,sprintf('offset')...
            ,sprintf('%0.5f \t random',results.random_offset)...
            ,sprintf('%0.5f \t caged',results.caged_offset)...
            ,sprintf('%0.5f \t directed' ,results.directed_offset)...
            ,sprintf(' ')...
            ,sprintf('R^2')...
            ,sprintf('%0.5f \t random',results.random_Rsq)...
            ,sprintf('%0.5f \t caged',results.caged_Rsq)...
            ,sprintf('%0.5f \t directed',results.directed_Rsq)...
            ,sprintf(' ')...
            ,sprintf('chi^2')...
            ,sprintf('%0.5f \t random',results.random_chi)...
            ,sprintf('%0.5f \t caged',results.caged_chi)...
            ,sprintf('%0.5f \t directed',results.directed_chi)...
            ,sprintf(' ')...
            ,sprintf('P')...
            ,sprintf('%0.5f \t random',results.random_P)...
            ,sprintf('%0.5f \t caged',results.caged_P)...
            ,sprintf('%0.5f \t directed',results.directed_P)...
            ,sprintf(' ')...
            ,sprintf('exitflag')...
            ,sprintf('%0.5f \t random',results.random_exitflag)...
            ,sprintf('%0.5f \t caged',results.caged_exitflag)...
            ,sprintf('%0.5f \t directed',results.directed_exitflag)};
        annotation('textbox',get(gca,'Position'),'String',info,'LineStyle','none');
        %text(0.1,0.95,str1,'HorizontalAlignment','Left','VerticalAlignment','Top');
        
%         subplot(1,3,3)
%          plot3(track1(1:17,1),track1(1:17,2),track1(1:17,3));
%         
        
        
        axis off;
    end
end


