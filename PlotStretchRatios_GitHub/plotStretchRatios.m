clc; clear; close all;

samples = {'A1', 'A2', 'B1', 'B2'};

for i = 1:length(samples)

    sampleType = samples{i};
    nExperiments = 3;

    % Set reference lengts in milimeters
     switch sampleType
        case 'A1'; x0 = 134; y0 = 134;
        case 'A2'; x0 = 333; y0 = 333;
        case 'B1'; x0 = 147; y0 = 147;
        case 'B2'; x0 = 216; y0 = 216;
    end

    dataFolder = 'exp';   % folder where your .txt files are stored

    plotFunction(sampleType, nExperiments, x0, y0, dataFolder);

end

function plotFunction(sampleType, nExperiments, x0, y0, dataFolder)
% plot_stretch_ratios
% -------------------
% sampleType     : string, e.g. 'A', 'B', 'C'
% nExperiments   : number of experimental datasets
% x0, y0         : original length in x and y
% dataFolder     : folder containing data files ('' if same dir)
%
% File naming convention:
%   Experimental: <sampleType>_exp_<i>.txt   (i = 1,2,...,nExperiments)
%   Theoretical : <sampleType>_thr.txt

    if nargin < 5
        dataFolder = '';
    end

    % Ensure folder ends with slash if not empty
    if ~isempty(dataFolder)
        if dataFolder(end) ~= filesep
            dataFolder = [dataFolder filesep];
        end
    end

    figure; hold on;axis equal

    % --- Plot experimental datasets (scatter) ---
    for i = 1:nExperiments
        
        fname = sprintf('%s%s_exp%d.txt', dataFolder, sampleType, i);
        if ~isfile(fname)
            warning('File not found: %s', fname);
            continue
        end
        
        data = load(fname);  % expects: [x  y]

        % For sample A and B
        if ~strcmp(sampleType, 'B2') 
            x = (data(:,1) + data(:,2)) / 2;
            y = (data(:,3) + data(:,4)) / 2;
        else
            % Sample C has two columns only
            x = data(:,1);
            y = data(:,2);    
        end

        lambda_x = x ./ x0;
        lambda_y = y ./ y0;
        
        scatter(lambda_x, lambda_y, 180, 'filled', 'DisplayName', sprintf('Experiment %d', i));

        % Point label
        for j=1:length(x)
            % text(lambda_x(j), lambda_y(j), num2str(j));
        end
    end

    % --- Plot theoretical curve (line) ---
    fname_thr = sprintf('%s%s_thr.txt', dataFolder, sampleType);
    if isfile(fname_thr)

        thr = load(fname_thr); % expects [x  y]
        lambda_x_t = thr(:,1);
        lambda_y_t = thr(:,2);

        % lambda_x_t = x_t ./ x0
        % lambda_y_t = y_t ./ y0

        plot(lambda_x_t, lambda_y_t, '-', 'LineWidth', 2, ...
            'DisplayName', 'Theoretical');
    else
        warning('Theoretical file not found: %s', fname_thr);
    end

    % xlabel('\lambda_x = x / x_0');
    % ylabel('\lambda_y = y / y_0');
    % xlim([0.611 1.08]);
    % ylim([0.707 1]);
    xlim([0.66 1.05])
    ylim([0.61 1.05])

    title(sprintf('Stretch ratio plot for sample %s', sampleType), 'FontSize', 18);
    legend('Location','best','FontSize',24);
    set(gca, 'FontSize', 48, 'LabelFontSizeMultiplier', 1, ...
        'TitleFontSizeMultiplier', 1, 'LineWidth', 1);
    hold off;
    % title(sprintf('Stretch ratio plot for sample %s', sampleType));
    % legend('Location','best');
    % hold off;

end
