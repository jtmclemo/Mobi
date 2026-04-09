function Mobi_smoke_test()

    dataFolder = fullfile(pwd, 'Final_CSVs');

    if ~exist(dataFolder, 'dir')
        parentDataFolder = fullfile(pwd, '..', 'Final_CSVs');
        if exist(parentDataFolder, 'dir')
            dataFolder = parentDataFolder;
        else
            error('Expected demo data folder not found: %s', dataFolder);
        end
    end

    files = dir(fullfile(dataFolder, '*.csv'));
    names = string({files.name});
    names = names(~startsWith(names, ["Mobi_", "bchs4397_", "group_"]));

    if numel(names) < 2
        error('Need at least two candidate CSV files in %s.', dataFolder);
    end

    selectedNames = names(1:2);
    selectedFiles = fullfile(dataFolder, selectedNames);

    X_all = [];
    ids = selectedNames(:);

    for k = 1:numel(selectedFiles)
        [X, ~, ~, ~] = Mobi_parsing(selectedFiles(k), 'Identity', strings(0,1));

        if size(X,1) ~= 1
            error('Smoke-test file produced %d feature rows instead of 1: %s', ...
                size(X,1), selectedFiles(k));
        end

        X_all = [X_all; X]; %#ok<AGROW>
    end

    tdaCfg = struct();
    tdaCfg.makePlots = false;
    tdaCfg.percentile = 75;
    tdaCfg.thresholdRule = "percentile";
    tdaCfg.distanceMetric = "euclidean";

    tdaResult = Mobi_tda(X_all, ids, tdaCfg);

    optCfg = struct();
    optCfg.elRankThreshold = 32.7;
    optCfg.familyMinOverlap = 7;
    optCfg.redundancyMinOverlap = 7;

    optResult = Mobi_optimization(selectedFiles(:), 3, [], optCfg);

    fprintf('Mobi smoke test passed.\n');
    fprintf('  Files checked: %d\n', numel(selectedFiles));
    fprintf('  Feature columns: %d\n', size(X_all,2));
    fprintf('  TDA clusters: %d\n', numel(unique(tdaResult.clusters)));
    fprintf('  Selected peptides: %d\n', height(optResult.selectedPeptides));
end
