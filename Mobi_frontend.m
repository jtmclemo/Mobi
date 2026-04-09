function Mobi_frontend()

    cfg = Mobi_config();
    clc;

    fprintf('============================================================\n');
    fprintf('%s (%s, %s)\n', cfg.general.programName, cfg.general.version, cfg.general.dateLabel);
    fprintf('============================================================\n\n');

    % HLA source selection
    fprintf('Choose HLA source for parsing/scoring:\n');
    fprintf('  1. Global/common\n');
    fprintf('  2. Region-based\n');
    fprintf('  3. Ancestry-group\n');
    fprintf('  4. OptiType file\n\n');

    hlaSourceChoice = input('Enter choice [1-4]: ', 's');
    hlaSourceChoice = strtrim(hlaSourceChoice);

    useReference   = false;
    useOptiType    = false;
    referenceHLAs  = strings(0,1);
    referenceInfo  = struct();
    optitypeTbl    = table();
    optitypeLookup = [];

    switch hlaSourceChoice
        case '1'
            useReference = true;
            referenceInfo = choose_global_reference_panel(cfg);
            referenceHLAs = referenceInfo.HLAs;

            fprintf('\nUsing built-in Mobi reference.\n');
            fprintf('Mode: %s\n', referenceInfo.ModeLabel);
            fprintf('Number of HLA alleles selected: %d\n\n', numel(referenceHLAs));

        case '2'
            useReference = true;
            referenceInfo = choose_reference_region_panel(cfg);
            referenceHLAs = referenceInfo.HLAs;

            fprintf('\nUsing built-in Mobi reference.\n');
            fprintf('Mode: %s\n', referenceInfo.ModeLabel);
            fprintf('Selected region: %s\n', referenceInfo.FilterLabel);
            fprintf('Number of HLA alleles selected: %d\n\n', numel(referenceHLAs));

        case '3'
            useReference = true;
            referenceInfo = choose_reference_ancestry_panel(cfg);
            referenceHLAs = referenceInfo.HLAs;

            fprintf('\nUsing built-in Mobi reference.\n');
            fprintf('Mode: %s\n', referenceInfo.ModeLabel);
            fprintf('Selected ancestry group: %s\n', referenceInfo.FilterLabel);
            fprintf('Number of HLA alleles selected: %d\n\n', numel(referenceHLAs));

        case '4'
            useOptiType = true;

            fprintf('Please select OptiType file.\n');
            [optiFile, optiPath] = uigetfile( ...
                {'*.csv;*.xlsx;*.xls', 'OptiType files (*.csv, *.xlsx, *.xls)'}, ...
                'Select OptiType file');

            if isequal(optiFile, 0)
                error('No OptiType file selected.');
            end

            optiFullPath = fullfile(optiPath, optiFile);
            fprintf('Using OptiType file: %s\n\n', optiFullPath);

            [~,~,ext] = fileparts(optiFullPath);

            if any(strcmpi(ext, {'.xlsx', '.xls'}))
                validate_optitype_spreadsheet_state(optiFullPath);
            end

            if strcmpi(ext, '.csv') || any(strcmpi(ext, {'.xlsx','.xls'}))
                optitypeTbl = readtable(optiFullPath, ...
                    'TextType', 'string', ...
                    'ReadVariableNames', false);

                if width(optitypeTbl) < 7
                    error('OptiType file must have at least 7 columns: filename plus 6 HLA columns.');
                end

                optitypeTbl = optitypeTbl(:, 1:7);
                optitypeTbl.Properties.VariableNames = ...
                    ["FileName","HLA1","HLA2","HLA3","HLA4","HLA5","HLA6"];
            else
                error('Unsupported OptiType file format: %s', ext);
            end

            fprintf('OptiType variable names detected:\n');
            print_string_list(string(optitypeTbl.Properties.VariableNames), '  ');
            optitypeLookup = build_optitype_lookup(optitypeTbl);

        otherwise
            error('Invalid HLA source choice. Please enter 1, 2, 3, or 4.');
    end

    % Select NetMHCpan folder
    fprintf('Please select NetMHCpan folder.\n');
    folder = uigetdir(pwd, 'Select folder containing NetMHCpan CSV files');

    if isequal(folder, 0)
        error('No NetMHCpan folder selected.');
    end

    fprintf('\nUsing NetMHCpan folder: %s\n\n', folder);

    % Start run log
    timestampTag = datestr(now, cfg.frontend.timestampFormat);

    if cfg.frontend.useDiaryLog
        logPath = fullfile(folder, ...
            sprintf('%s%s%s', cfg.frontend.logPrefix, timestampTag, cfg.frontend.logExt));
        logFID = fopen(logPath, 'w');
        if logFID < 0
            error('Could not open log file for writing: %s', logPath);
        end
    else
        logPath = "";
        logFID = -1;
    end

    cleanupObj = onCleanup(@() safely_close_log(logFID));

    % Run header
    log_only(logFID, '==================== RUN HEADER ====================\n');
    log_only(logFID, 'Program: %s\n', cfg.general.programName);
    log_only(logFID, 'Version: %s\n', cfg.general.version);
    log_only(logFID, 'Date label: %s\n', cfg.general.dateLabel);
    log_only(logFID, 'Working folder: %s\n', folder);

    if useReference
        fprintf('Selected HLA source mode: %s\n', referenceInfo.ModeLabel);
        fprintf('Selected HLA filter: %s\n', referenceInfo.FilterLabel);
        fprintf('Selected HLA count: %d\n\n', numel(referenceHLAs));
        log_only(logFID, 'HLA source mode: %s\n', referenceInfo.ModeLabel);
        log_only(logFID, 'HLA filter label: %s\n', referenceInfo.FilterLabel);
        log_only(logFID, 'HLA count used: %d\n', numel(referenceHLAs));
    elseif useOptiType
        fprintf('Selected HLA source mode: OptiType\n');
        fprintf('Selected HLA filter: Patient-specific\n\n');
        log_only(logFID, 'HLA source mode: OptiType\n');
        log_only(logFID, 'HLA filter label: Patient-specific\n');
    end

    if strlength(logPath) > 0
        fprintf('Full transcript log file: %s\n\n', logPath);
        log_only(logFID, 'Txt log file: %s\n', logPath);
    end
    log_only(logFID, '====================================================\n\n');

    % Find candidate NetMHCpan CSV files
    files = dir(fullfile(folder, '*.csv'));

    if isempty(files)
        error('No CSV files found in selected folder.');
    end

    allNames = string({files.name})';
    excludePrefixes = string(cfg.frontend.excludePrefixes(:));

    excludeMask = false(size(allNames));
    for p = 1:numel(excludePrefixes)
        excludeMask = excludeMask | startsWith(allNames, excludePrefixes(p));
    end

    files = files(~excludeMask);

    if isempty(files)
        error(['No valid NetMHCpan CSV files remained after excluding ', ...
               'program-generated output CSVs.']);
    end

    log_only(logFID, 'Found %d candidate NetMHCpan CSV files after filtering.\n\n', length(files));

    % Initialize storage
    X_all          = [];
    ids_all        = strings(0,1);
    summaryTbl_all = table();
    feature_names  = strings(0,1);
    sampleFiles    = strings(0,1);
    id_col         = cfg.parsing.idColumn;

    skippedFiles = strings(0,1);
    skipReasons  = strings(0,1);

    % Parse NetMHCpan files
    for k = 1:length(files)

        filename = fullfile(folder, files(k).name);
        print_both(logFID, 'Processing file %d of %d: %s\n', k, length(files), files(k).name);

        try
            sampleFileName = string(files(k).name);

            if useReference
                typedHLAs      = referenceHLAs;
                hlaSourceLabel = string(referenceInfo.ModeLabel);
                hlaFilterLabel = string(referenceInfo.FilterLabel);
            else
                typedHLAs = get_optitype_hlas(optitypeLookup, sampleFileName);

                if isempty(typedHLAs)
                    log_only(logFID, 'WARNING: No OptiType row found for %s. Skipping.\n', files(k).name);
                    skippedFiles(end+1,1) = sampleFileName; %#ok<AGROW>
                    skipReasons(end+1,1)  = "No matching OptiType row"; %#ok<AGROW>
                    continue;
                end

                hlaSourceLabel = "OptiType";
                hlaFilterLabel = "Patient-specific";
            end

            [X, feature_names, ~, summaryTbl] = Mobi_parsing(filename, id_col, typedHLAs);

            if size(X,1) ~= 1
                log_only(logFID, 'WARNING: File %s produced %d rows instead of 1. Skipping.\n', ...
                    files(k).name, size(X,1));

                skippedFiles(end+1,1) = sampleFileName; %#ok<AGROW>
                skipReasons(end+1,1)  = "Parser did not return exactly one row"; %#ok<AGROW>
                continue;
            end

            X_all = [X_all; X]; %#ok<AGROW>
            ids_all(end+1,1) = sampleFileName; %#ok<AGROW>
            sampleFiles(end+1,1) = string(filename); %#ok<AGROW>

            summaryTbl.FileName  = repmat(sampleFileName, height(summaryTbl), 1);
            summaryTbl.HLASource = repmat(hlaSourceLabel, height(summaryTbl), 1);
            summaryTbl.HLAFilter = repmat(hlaFilterLabel, height(summaryTbl), 1);
            summaryTbl.HLAsUsed  = repmat(strjoin(typedHLAs, ";"), height(summaryTbl), 1);

            summaryTbl_all = union_and_append_tables(summaryTbl_all, summaryTbl);

        catch ME
            log_only(logFID, 'WARNING: Skipping file %s because of error:\n%s\n', files(k).name, ME.message);
            skippedFiles(end+1,1) = string(files(k).name); %#ok<AGROW>
            skipReasons(end+1,1)  = string(ME.message); %#ok<AGROW>
        end
    end

    log_only(logFID, '\nAccepted files: %d\n', numel(ids_all));
    log_only(logFID, 'Skipped files:  %d\n\n', numel(skippedFiles));

    if ~isempty(skippedFiles)
        log_only(logFID, 'Skipped file summary:\n');
        log_only(logFID, '%s', evalc('print_two_column_pairs(''File'', skippedFiles, ''Reason'', skipReasons);'));
        log_only(logFID, '\n');
    end

    if size(X_all,1) < 2
        error('Need at least 2 valid files to run TDA.');
    end

    log_only(logFID, 'Features used in TDA:\n');
    log_only(logFID, '%s', evalc('print_string_list(feature_names(:), ''  '');'));
    log_only(logFID, '\n');

    % Build baseline TDA / optimization configs from cfg
    tdaCfg = struct();
    tdaCfg.makePlots       = cfg.tda.makePlots;
    tdaCfg.distanceMetric  = cfg.tda.distanceMetric;
    tdaCfg.weights13       = cfg.tda.defaultWeights13;
    tdaCfg.thresholdRule   = cfg.tda.thresholdRule;
    tdaCfg.percentile      = cfg.tda.defaultPercentile;
    tdaCfg.mergeSingletons = true;
    tdaCfg.maxOptimizationGroupSize = Inf;
    tdaCfg.minOptimizationGroupSize = 1;

    optCfg = struct();
    optCfg.maxPeptides = cfg.optimization.maxPeptides;
    optCfg.elRankThreshold = cfg.optimization.elRankThreshold;
    optCfg.familyMinOverlap = cfg.optimization.familyMinOverlap;
    optCfg.redundancyMinOverlap = cfg.optimization.redundancyMinOverlap;
    optCfg.normalizeWeights = cfg.optimization.normalizeWeights;
    optCfg.weights = cfg.optimization.weights;

    % Run TDA
    log_only(logFID, '\nRunning TDA on combined dataset...\n');
    tTDA = tic;
    tdaResult = Mobi_tda(X_all, ids_all, tdaCfg);
    mainTDATime = toc(tTDA);
    log_only(logFID, 'Main TDA runtime: %.2f seconds\n', mainTDATime);

    if isfield(tdaResult, 'chosenPercentile')
        log_only(logFID, 'Chosen epsilon percentile: %.2f\n', tdaResult.chosenPercentile);
    elseif isfield(tdaCfg, 'percentile')
        log_only(logFID, 'Configured epsilon percentile: %.2f\n', tdaCfg.percentile);
    end

    fileNumbers = (1:length(ids_all))';

    print_both(logFID, '\nMerge events (in order of epsilon):\n');
    edges = tdaResult.edges;
    nSamples = length(ids_all);

    parent = 1:nSamples;
    members = cell(nSamples,1);
    for s = 1:nSamples
        members{s} = fileNumbers(s);
    end

    for k = 1:size(edges,1)
        i = edges(k,1);
        j = edges(k,2);
        d = edges(k,3);

        ri = find_root(parent, i);
        rj = find_root(parent, j);

        if ri ~= rj
            leftGroup   = sort(members{ri});
            rightGroup  = sort(members{rj});
            mergedGroup = sort([leftGroup; rightGroup]);

            print_both(logFID, 'epsilon = %.3f : %s <--> %s  =>  %s\n', ...
                d, ...
                group_to_string(leftGroup), ...
                group_to_string(rightGroup), ...
                group_to_string(mergedGroup));

            parent(rj) = ri;
            members{ri} = mergedGroup;
            members{rj} = [];
        end
    end

    print_both(logFID, 'Chosen clustering threshold: %.4f\n\n', tdaResult.threshold);

    Tclusters = table(fileNumbers, ids_all, tdaResult.clusters, ...
        'VariableNames', {'FileNumber','File','Cluster'});
    print_both(logFID, 'Cluster assignments:\n');
    clusterText = evalc('print_cluster_table(Tclusters);');
    print_text_both(logFID, clusterText);
    print_both(logFID, '\n');

    rawGroupLabels  = tdaResult.clusters;
    rawUniqueGroups = unique(rawGroupLabels);

    log_only(logFID, 'Cluster size summary:\n');
    for g = 1:length(rawUniqueGroups)
        thisGroup = rawUniqueGroups(g);
        idx = (rawGroupLabels == thisGroup);
        nums = find(idx)';
        log_only(logFID, '  Group %d: %d files (%s)\n', ...
            thisGroup, sum(idx), format_file_list(nums));
    end
    log_only(logFID, '\n');

    if isfield(tdaResult, 'optimizationGroups')
        practicalGroupLabels  = tdaResult.optimizationGroups;
        practicalUniqueGroups = unique(practicalGroupLabels);

        log_only(logFID, 'Practical optimization-group size summary:\n');
        for g = 1:length(practicalUniqueGroups)
            thisGroup = practicalUniqueGroups(g);
            idx = (practicalGroupLabels == thisGroup);
            nums = find(idx)';
            log_only(logFID, '  Practical Group %d: %d files (%s)\n', ...
                thisGroup, sum(idx), format_file_list(nums));
        end
        log_practical_grouping_audit(logFID, tdaResult, false);
        log_only(logFID, '\n');
    else
        practicalGroupLabels  = rawGroupLabels;
        practicalUniqueGroups = rawUniqueGroups;
    end

    % Prompt for optimization
    doOptimization = false;
    mainOptimizationTime = 0;
    allResults = struct();
    robustness = [];

    while true
        userInput = input('Would you like to continue with neopeptide optimization? (y/n): ', 's');

        if strcmpi(userInput, 'y') || strcmpi(userInput, 'yes')
            doOptimization = true;
            break;
        elseif strcmpi(userInput, 'n') || strcmpi(userInput, 'no')
            doOptimization = false;
            break;
        else
            fprintf('Please enter y or n.\n');
        end
    end

    if doOptimization

        log_only(logFID, '\nProceeding with neopeptide optimization...\n\n');

        maxPeptidesInput = input( ...
            sprintf('Maximum number of peptides to select [default %d]: ', ...
            cfg.optimization.maxPeptides), 's');

        if isempty(strtrim(maxPeptidesInput))
            maxPeptides = cfg.optimization.maxPeptides;
        else
            maxPeptides = str2double(maxPeptidesInput);
        end

        if isnan(maxPeptides) || maxPeptides <= 0 || mod(maxPeptides,1) ~= 0
            error('Maximum number of peptides must be a positive integer.');
        end

        minGroupSizeInput = input( ...
            'Minimum group size for practical subgroup optimization [default 3]: ', 's');

        if isempty(strtrim(minGroupSizeInput))
            minGroupSizeForPracticalOptimization = 3;
        else
            minGroupSizeForPracticalOptimization = str2double(minGroupSizeInput);
        end

        if isnan(minGroupSizeForPracticalOptimization) || ...
                minGroupSizeForPracticalOptimization <= 0 || ...
                mod(minGroupSizeForPracticalOptimization,1) ~= 0
            error('Minimum practical group size must be a positive integer.');
        end

        % Re-run TDA once with practical grouping constraints now known
        tdaCfg.maxOptimizationGroupSize = maxPeptides;
        tdaCfg.minOptimizationGroupSize = minGroupSizeForPracticalOptimization;
        tdaCfg.makePlots = false;
        tdaResult = Mobi_tda(X_all, ids_all, tdaCfg);

        % Refresh raw and practical labels after practical regrouping
        rawGroupLabels  = tdaResult.clusters;
        rawUniqueGroups = unique(rawGroupLabels);

        if isfield(tdaResult, 'optimizationGroups')
            practicalGroupLabels  = tdaResult.optimizationGroups;
            practicalUniqueGroups = unique(practicalGroupLabels);
        else
            practicalGroupLabels  = rawGroupLabels;
            practicalUniqueGroups = rawUniqueGroups;
        end

        defaultWeights = cfg.optimization.weights;

        fprintf('Enter custom optimization weights or press Enter to keep defaults.\n');
        log_only(logFID, 'Enter custom optimization weights or press Enter to keep defaults.\n');

        wPatient    = input(sprintf('Weight for patient coverage [default %.2f]: ', ...
                       defaultWeights.wPatientCoverage), 's');
        wHLA        = input(sprintf('Weight for HLA coverage [default %.2f]: ', ...
                       defaultWeights.wHLACoverage), 's');
        wBinding    = input(sprintf('Weight for binding [default %.2f]: ', ...
                       defaultWeights.wBinding), 's');
        wPrevalence = input(sprintf('Weight for prevalence [default %.2f]: ', ...
                       defaultWeights.wPrevalence), 's');
        wFamily     = input(sprintf('Weight for family novelty [default %.2f]: ', ...
                       defaultWeights.wFamilyNovelty), 's');
        wRedundancy = input(sprintf('Weight for redundancy penalty [default %.2f]: ', ...
                       defaultWeights.wRedundancy), 's');

        WEIGHTS = defaultWeights;

        if ~isempty(strtrim(wPatient)),    WEIGHTS.wPatientCoverage = str2double(wPatient); end
        if ~isempty(strtrim(wHLA)),        WEIGHTS.wHLACoverage     = str2double(wHLA); end
        if ~isempty(strtrim(wBinding)),    WEIGHTS.wBinding         = str2double(wBinding); end
        if ~isempty(strtrim(wPrevalence)), WEIGHTS.wPrevalence      = str2double(wPrevalence); end
        if ~isempty(strtrim(wFamily)),     WEIGHTS.wFamilyNovelty   = str2double(wFamily); end
        if ~isempty(strtrim(wRedundancy)), WEIGHTS.wRedundancy      = str2double(wRedundancy); end

        weightValues = [ ...
            WEIGHTS.wPatientCoverage, ...
            WEIGHTS.wHLACoverage, ...
            WEIGHTS.wBinding, ...
            WEIGHTS.wPrevalence, ...
            WEIGHTS.wFamilyNovelty, ...
            WEIGHTS.wRedundancy];

        if any(isnan(weightValues)) || any(weightValues < 0)
            error('All optimization weights must be numeric and nonnegative.');
        end

        if sum(weightValues) == 0
            error('At least one optimization weight must be positive.');
        end

        if cfg.optimization.normalizeWeights
            weightValues = weightValues / sum(weightValues);
            WEIGHTS.wPatientCoverage = weightValues(1);
            WEIGHTS.wHLACoverage     = weightValues(2);
            WEIGHTS.wBinding         = weightValues(3);
            WEIGHTS.wPrevalence      = weightValues(4);
            WEIGHTS.wFamilyNovelty   = weightValues(5);
            WEIGHTS.wRedundancy      = weightValues(6);
        end

        optCfg.maxPeptides = maxPeptides;
        optCfg.weights = WEIGHTS;
        optCfg.minGroupSizeForPracticalOptimization = minGroupSizeForPracticalOptimization;
        optCfg.referenceEvaluation = build_reference_evaluation_config(useReference, referenceInfo);

        print_both(logFID, '\nUsing normalized optimization weights:\n');
        print_both(logFID, '  Patient coverage:   %.4f\n', WEIGHTS.wPatientCoverage);
        print_both(logFID, '  HLA coverage:       %.4f\n', WEIGHTS.wHLACoverage);
        print_both(logFID, '  Binding:            %.4f\n', WEIGHTS.wBinding);
        print_both(logFID, '  Prevalence:         %.4f\n', WEIGHTS.wPrevalence);
        print_both(logFID, '  Family novelty:     %.4f\n', WEIGHTS.wFamilyNovelty);
        print_both(logFID, '  Redundancy penalty: %.4f\n', WEIGHTS.wRedundancy);
        print_both(logFID, 'Practical max group size: %d\n', maxPeptides);
        print_both(logFID, 'Minimum practical group size: %d\n\n', minGroupSizeForPracticalOptimization);

        log_only(logFID, 'Updated practical optimization-group size summary:\n');
        for g = 1:length(practicalUniqueGroups)
            thisGroup = practicalUniqueGroups(g);
            idx = (practicalGroupLabels == thisGroup);
            nums = find(idx)';
            log_only(logFID, '  Practical Group %d: %d files (%s)\n', ...
                thisGroup, sum(idx), format_file_list(nums));
        end
        log_practical_grouping_audit(logFID, tdaResult, true);
        log_only(logFID, '\n');

        tOpt = tic;
        baseInputForSuite = struct();
        baseInputForSuite.sampleFiles = sampleFiles;
        baseInputForSuite.sampleNames = ids_all;

        optimizationSuite = Mobi_run_optimization_suite( ...
            baseInputForSuite, optCfg, tdaResult, "all");
        allResults = optimizationSuite.allResults;

        log_optimization_suite(logFID, optimizationSuite);
        print_practical_optimization_suite(logFID, optimizationSuite);

        mainOptimizationTime = toc(tOpt);
        log_only(logFID, '\nMain optimization runtime: %.2f seconds\n', mainOptimizationTime);

        log_only(logFID, '\nGlobal + honest + practical optimizations complete.\n\n');
    else
        log_only(logFID, '\nNeopeptide optimization skipped for the main run. Continuing to robustness prompt.\n\n');
    end

    % Optional robustness check
    robustnessMode = prompt_robustness_mode(cfg, mainTDATime, mainOptimizationTime);

    if robustnessMode ~= 0

        log_only(logFID, '\nPreparing robustness analysis...\n');

        baseInput = struct();
        baseInput.X = X_all;
        baseInput.sampleNames = ids_all;
        baseInput.sampleFiles = sampleFiles;
        baseInput.groupLabels = rawGroupLabels;
        baseInput.uniqueGroups = rawUniqueGroups;
        baseInput.folder = folder;
        baseInput.mainTDAResult = tdaResult;
        baseInput.didOptimization = doOptimization;
        baseInput.optCfg = optCfg;

        runFns = struct();
        runFns.runTDA = @(baseInput_, tdaCfg_) run_tda_wrapper(baseInput_, tdaCfg_);
        runFns.runOptimization = @(baseInput_, optCfg_, tdaResult_) run_optimization_wrapper(baseInput_, optCfg_, tdaResult_);

        robustnessText = evalc('robustness = Mobi_robustness(runFns, baseInput, cfg, robustnessMode);');
        log_only(logFID, '%s', robustnessText);

    else
        log_only(logFID, '\nRobustness check skipped.\n');
    end

    if cfg.frontend.enableStructuredCSVExport
        outputDir = export_structured_outputs( ...
            folder, timestampTag, cfg, summaryTbl_all, skippedFiles, skipReasons, ...
            ids_all, fileNumbers, rawGroupLabels, practicalGroupLabels, ...
            Tclusters, tdaResult, allResults, robustness, logPath);
        log_only(logFID, '\nStructured output folder: %s\n', outputDir);
    end

    print_both(logFID, '\nRun complete.\n');

end


% STRUCTURED EXPORTS
function outputDir = export_structured_outputs( ...
    folder, timestampTag, cfg, summaryTbl_all, skippedFiles, skipReasons, ...
    ids_all, fileNumbers, rawGroupLabels, practicalGroupLabels, ...
    Tclusters, tdaResult, allResults, robustness, logPath)

    outputDir = fullfile(folder, sprintf('%s%s', cfg.frontend.outputPrefix, timestampTag));

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    fprintf('\nExporting structured output files...\n');

    runInfo = table( ...
        string(cfg.general.programName), ...
        string(cfg.general.version), ...
        string(cfg.general.dateLabel), ...
        string(folder), ...
        string(logPath), ...
        string(datestr(now, 'yyyy-mm-dd HH:MM:SS')), ...
        'VariableNames', {'Program','Version','DateLabel','InputFolder','LogPath','ExportedAt'});
    safe_writetable(runInfo, fullfile(outputDir, 'run_info.csv'));

    if istable(summaryTbl_all) && ~isempty(summaryTbl_all)
        safe_writetable(summaryTbl_all, fullfile(outputDir, 'parsed_feature_summary.csv'));
    end

    skippedTbl = table(string(skippedFiles(:)), string(skipReasons(:)), ...
        'VariableNames', {'File','Reason'});
    safe_writetable(skippedTbl, fullfile(outputDir, 'skipped_files.csv'));

    safe_writetable(Tclusters, fullfile(outputDir, 'tda_raw_clusters.csv'));

    if ~isempty(ids_all)
        practicalTbl = table( ...
            fileNumbers(:), ...
            string(ids_all(:)), ...
            rawGroupLabels(:), ...
            practicalGroupLabels(:), ...
            'VariableNames', {'FileNumber','File','RawCluster','PracticalGroup'});
        safe_writetable(practicalTbl, fullfile(outputDir, 'tda_practical_groups.csv'));
    end

    if isfield(tdaResult, 'barcode') && ~isempty(tdaResult.barcode)
        barcodeTbl = table(tdaResult.barcode(:,1), tdaResult.barcode(:,2), ...
            'VariableNames', {'Birth','Death'});
        safe_writetable(barcodeTbl, fullfile(outputDir, 'tda_barcode.csv'));
    end

    tdaSummary = table( ...
        string(tdaResult.thresholdRule), ...
        tdaResult.threshold, ...
        tdaResult.chosenPercentile, ...
        string(tdaResult.distanceMetric), ...
        mat2str(tdaResult.weights), ...
        'VariableNames', {'ThresholdRule','Threshold','ChosenPercentile','DistanceMetric','Weights'});
    safe_writetable(tdaSummary, fullfile(outputDir, 'tda_summary.csv'));

    [selectedTbl, historyTbl, featureTbl] = collect_optimization_export_tables(allResults);
    referenceMetricTbl = collect_reference_metric_export_table(allResults);

    if istable(selectedTbl) && ~isempty(selectedTbl)
        safe_writetable(selectedTbl, fullfile(outputDir, 'optimization_selected_peptides.csv'));
    end

    if istable(historyTbl) && ~isempty(historyTbl)
        safe_writetable(historyTbl, fullfile(outputDir, 'optimization_selection_history.csv'));
    end

    if istable(featureTbl) && ~isempty(featureTbl)
        safe_writetable(featureTbl, fullfile(outputDir, 'optimization_peptide_features.csv'));
    end

    if istable(referenceMetricTbl) && ~isempty(referenceMetricTbl)
        safe_writetable(referenceMetricTbl, fullfile(outputDir, 'optimization_reference_metrics.csv'));
    end

    if isstruct(robustness) && isfield(robustness, 'summaryTable') && istable(robustness.summaryTable)
        safe_writetable(robustness.summaryTable, fullfile(outputDir, 'robustness_summary.csv'));
    end

    save(fullfile(outputDir, 'Mobi_run_workspace.mat'), ...
        'cfg', 'summaryTbl_all', 'skippedFiles', 'skipReasons', ...
        'ids_all', 'fileNumbers', 'rawGroupLabels', 'practicalGroupLabels', ...
        'tdaResult', 'allResults', 'robustness', '-v7.3');
end


function [selectedTbl, historyTbl, featureTbl] = collect_optimization_export_tables(allResults)

    selectedTbl = table();
    historyTbl = table();
    featureTbl = table();

    if isempty(allResults) || ~isstruct(allResults) || ~isfield(allResults, 'result')
        return;
    end

    for i = 1:numel(allResults)
        if isempty(allResults(i).result) || ~isstruct(allResults(i).result)
            continue;
        end

        label = string(allResults(i).label);
        mode = string(allResults(i).mode);
        group = allResults(i).group;
        fileList = strjoin(string(allResults(i).selectedFileNames(:)), ";");

        res = allResults(i).result;

        if isfield(res, 'selectedPeptides') && istable(res.selectedPeptides) && ~isempty(res.selectedPeptides)
            T = add_export_context(res.selectedPeptides, mode, label, group, fileList);
            selectedTbl = union_and_append_tables(selectedTbl, T);
        end

        if isfield(res, 'selectionHistory') && istable(res.selectionHistory) && ~isempty(res.selectionHistory)
            T = add_export_context(res.selectionHistory, mode, label, group, fileList);
            historyTbl = union_and_append_tables(historyTbl, T);
        end

        if isfield(res, 'peptideFeatureTable') && istable(res.peptideFeatureTable) && ~isempty(res.peptideFeatureTable)
            T = add_export_context(res.peptideFeatureTable, mode, label, group, fileList);
            featureTbl = union_and_append_tables(featureTbl, T);
        end
    end
end


function metricsTbl = collect_reference_metric_export_table(allResults)

    metricsTbl = table();

    if isempty(allResults) || ~isstruct(allResults) || ~isfield(allResults, 'result')
        return;
    end

    for i = 1:numel(allResults)
        if isempty(allResults(i).result) || ~isstruct(allResults(i).result) || ...
                ~isfield(allResults(i).result, 'referenceMetrics')
            continue;
        end

        metrics = allResults(i).result.referenceMetrics;
        if ~isstruct(metrics) || ~isfield(metrics, 'didEvaluate') || ~metrics.didEvaluate
            continue;
        end

        T = table( ...
            string(allResults(i).mode), ...
            string(allResults(i).label), ...
            allResults(i).group, ...
            numel(allResults(i).files), ...
            metrics.numSelectedPeptides, ...
            metrics.numHitHLAs, ...
            metrics.numReferencePopulations, ...
            metrics.coveragePercent, ...
            metrics.expectedHits, ...
            metrics.conditionalMeanHits, ...
            metrics.pc90, ...
            'VariableNames', { ...
                'OptimizationMode', ...
                'OptimizationLabel', ...
                'Group', ...
                'NumFiles', ...
                'NumSelectedPeptides', ...
                'NumHitHLAs', ...
                'NumReferencePopulations', ...
                'CoveragePercent', ...
                'ExpectedHits', ...
                'ConditionalMeanHits', ...
                'PC90'});

        metricsTbl = [metricsTbl; T]; %#ok<AGROW>
    end
end


function T = add_export_context(T, mode, label, group, fileList)
    T.OptimizationMode = repmat(mode, height(T), 1);
    T.OptimizationLabel = repmat(label, height(T), 1);
    T.Group = repmat(group, height(T), 1);
    T.InputFiles = repmat(fileList, height(T), 1);
end


function safe_writetable(T, outPath)
    try
        writetable(T, outPath);
        fprintf('  Wrote %s\n', outPath);
    catch ME
        warning('Could not write %s: %s', outPath, ME.message);
    end
end


function refEval = build_reference_evaluation_config(useReference, referenceInfo)
    refEval = struct();
    refEval.enabled = false;
    refEval.referenceTable = table();
    refEval.modeLabel = "";
    refEval.filterLabel = "";

    if ~useReference || ~isstruct(referenceInfo) || ~isfield(referenceInfo, 'ReferenceTable')
        return;
    end

    refEval.enabled = true;
    refEval.referenceTable = referenceInfo.ReferenceTable;
    refEval.populationWeightMode = "sample_size";
    refEval.hitCountMode = "multiplicity";

    if isfield(referenceInfo, 'ModeLabel')
        refEval.modeLabel = string(referenceInfo.ModeLabel);
    end

    if isfield(referenceInfo, 'FilterLabel')
        refEval.filterLabel = string(referenceInfo.FilterLabel);
    end
end


function log_reference_metrics(logFID, result)
    if ~isstruct(result) || ~isfield(result, 'referenceMetrics')
        return;
    end

    metrics = result.referenceMetrics;
    if ~isstruct(metrics) || ~isfield(metrics, 'didEvaluate') || ~metrics.didEvaluate
        return;
    end

    log_only(logFID, ...
        ['Reference population effectiveness: coverage %.2f%% | expected hits %.2f | ', ...
         'conditional mean hits %.2f | PC90 %.2f\n'], ...
        metrics.coveragePercent, metrics.expectedHits, metrics.conditionalMeanHits, metrics.pc90);
end


function print_reference_metrics(logFID, result)
    if ~isstruct(result) || ~isfield(result, 'referenceMetrics')
        return;
    end

    metrics = result.referenceMetrics;
    if ~isstruct(metrics) || ~isfield(metrics, 'didEvaluate') || ~metrics.didEvaluate
        return;
    end

    print_both(logFID, ...
        ['Reference population effectiveness: coverage %.2f%% | expected hits %.2f | ', ...
         'conditional mean hits %.2f | PC90 %.2f\n'], ...
        metrics.coveragePercent, metrics.expectedHits, metrics.conditionalMeanHits, metrics.pc90);
end


% OPTIMIZATION SUITE PRINTING
function log_optimization_suite(logFID, suite)

    log_only(logFID, '\n============================================================\n');
    log_only(logFID, 'GLOBAL COHORT OPTIMIZATION\n');
    log_only(logFID, '============================================================\n');
    log_group_result(logFID, suite.global);

    log_only(logFID, '\n============================================================\n');
    log_only(logFID, 'HONEST PER-CLUSTER OPTIMIZATION (RAW TDA GROUPS)\n');
    log_only(logFID, '============================================================\n');

    for i = 1:numel(suite.honestGroups)
        log_group_result(logFID, suite.honestGroups(i));
    end

    log_only(logFID, '\nOptimization suite peptide-set semantics:\n');
    log_only(logFID, '  globalSelectedPeptides: %d peptides\n', numel(suite.globalSelectedPeptides));
    log_only(logFID, '  honestSelectedPeptideUnion: %d peptides\n', numel(suite.honestSelectedPeptideUnion));
    log_only(logFID, '  practicalSelectedPeptideUnion: %d peptides\n', numel(suite.practicalSelectedPeptideUnion));
    log_only(logFID, '  selectedPeptideUnionAcrossModes: %d peptides\n', numel(suite.selectedPeptideUnionAcrossModes));
end


function log_practical_grouping_audit(logFID, tdaResult, printCompactToConsole)

    if ~isfield(tdaResult, 'practicalGroupingAudit') || isempty(tdaResult.practicalGroupingAudit)
        log_only(logFID, 'Practical grouping audit: unavailable.\n');
        return;
    end

    audit = tdaResult.practicalGroupingAudit;
    summary = audit.summary;

    if printCompactToConsole
        print_both(logFID, 'Practical grouping audit summary:\n');
        print_both(logFID, '  Small-group merges: %d\n', summary.NumMerges);
        print_both(logFID, '  Oversized-group splits: %d\n', summary.NumSplits);
        print_both(logFID, '  Repair passes: %d\n', summary.NumRepairSteps);
        print_both(logFID, '  Files changed practical group relative to initial raw grouping: %d\n\n', ...
            summary.NumFilesChanged);
    else
        log_only(logFID, 'Practical grouping audit summary:\n');
        log_only(logFID, '  Small-group merges: %d\n', summary.NumMerges);
        log_only(logFID, '  Oversized-group splits: %d\n', summary.NumSplits);
        log_only(logFID, '  Repair passes: %d\n', summary.NumRepairSteps);
        log_only(logFID, '  Files changed practical group relative to initial raw grouping: %d\n', ...
            summary.NumFilesChanged);
    end

    log_only(logFID, '\nPractical grouping audit events:\n');
    events = audit.events;
    for i = 1:numel(events)
        evidence = "not applicable";
        if isfinite(events(i).DistanceEvidence)
            evidence = sprintf('%.6g', events(i).DistanceEvidence);
        end

        log_only(logFID, '  Step %d | %s | %s\n', ...
            events(i).StepNumber, events(i).Stage, events(i).ActionType);
        log_only(logFID, '    Reason: %s\n', events(i).Reason);
        log_only(logFID, '    Affected files: %s\n', format_file_list(events(i).AffectedFileIndices));
        log_only(logFID, '    Affected names: %s\n', strjoin(events(i).AffectedFileNames(:)', ', '));
        log_only(logFID, '    Raw groups: %s\n', format_numeric_list(events(i).SourceRawGroups));
        log_only(logFID, '    Practical before: %s\n', format_numeric_list(events(i).SourcePracticalGroupsBefore));
        log_only(logFID, '    Practical after: %s\n', format_numeric_list(events(i).TargetPracticalGroupsAfter));
        log_only(logFID, '    Previous group sizes: %s\n', format_numeric_list(events(i).PreviousGroupSizes));
        log_only(logFID, '    New group sizes: %s\n', format_numeric_list(events(i).NewGroupSizes));
        log_only(logFID, '    Distance/split evidence: %s\n', evidence);
    end

    log_only(logFID, '\nRaw-to-practical grouping mapping:\n');
    mapping = audit.rawToPracticalMapping;
    for i = 1:height(mapping)
        log_only(logFID, ...
            '  File %d | %s | raw %d | initial practical %d | final practical %d | moved %d | %s\n', ...
            mapping.FileNumber(i), mapping.FileName(i), mapping.RawCluster(i), ...
            mapping.InitialPracticalGroup(i), mapping.PracticalGroup(i), ...
            mapping.MovedFromInitialPracticalGroup(i), mapping.AuditReasonCode(i));
    end
end


function print_practical_optimization_suite(logFID, suite)

    print_both(logFID, '\n============================================================\n');
    print_both(logFID, 'PRACTICAL PER-CLUSTER OPTIMIZATION\n');
    print_both(logFID, '============================================================\n');

    for i = 1:numel(suite.practicalGroups)
        groupResult = suite.practicalGroups(i);

        if ~groupResult.didRun
            print_both(logFID, '\nSkipping %s (%s).\n', groupResult.label, groupResult.skipReason);
            continue;
        end

        print_both(logFID, '\n=== %s ===\n', groupResult.label);
        print_both(logFID, 'Practical group size: %d\n', numel(groupResult.files));
        print_both(logFID, 'Files: %s\n', format_name_list(groupResult.selectedFileNames));

        if isfield(groupResult.result, 'selectedPeptides') && ...
                istable(groupResult.result.selectedPeptides) && ...
                ~isempty(groupResult.result.selectedPeptides)
            print_both(logFID, 'Selected peptides:\n');
            print_text_both(logFID, evalc('print_generic_table(groupResult.result.selectedPeptides);'));
        else
            print_both(logFID, 'No peptides selected for %s.\n', groupResult.label);
        end

        print_reference_metrics(logFID, groupResult.result);
    end
end


function log_group_result(logFID, groupResult)

    if ~isstruct(groupResult) || ~isfield(groupResult, 'label') || groupResult.label == ""
        return;
    end

    if ~groupResult.didRun
        log_only(logFID, '\nSkipping %s (%s).\n', groupResult.label, groupResult.skipReason);
        return;
    end

    log_only(logFID, '\n=== %s ===\n', groupResult.label);
    log_only(logFID, 'Group size: %d\n', numel(groupResult.files));
    log_only(logFID, 'Files: %s\n', format_name_list(groupResult.selectedFileNames));

    if isfield(groupResult.result, 'selectedPeptides') && ...
            istable(groupResult.result.selectedPeptides) && ...
            ~isempty(groupResult.result.selectedPeptides)
        log_only(logFID, 'Selected peptides:\n');
        log_only(logFID, '%s', evalc('print_generic_table(groupResult.result.selectedPeptides);'));
    else
        log_only(logFID, 'No peptides selected for %s.\n', groupResult.label);
    end

    log_reference_metrics(logFID, groupResult.result);
end


% ROBUSTNESS WRAPPERS / PROMPTS
function mode = prompt_robustness_mode(cfg, mainTDATime, mainOptimizationTime)

[fastRuns, fastSec] = estimate_robustness_time(cfg, 1, mainTDATime, mainOptimizationTime);
[standardRuns, standardSec] = estimate_robustness_time(cfg, 2, mainTDATime, mainOptimizationTime);
[fullRuns, fullSec] = estimate_robustness_time(cfg, 3, mainTDATime, mainOptimizationTime);

fprintf('\nOptional robustness check:\n');
fprintf('  0. None\n');
fprintf('  1. Fast (~%d runs, estimated %s)\n', fastRuns, format_duration(fastSec));
fprintf('  2. Standard (~%d runs, estimated %s)\n', standardRuns, format_duration(standardSec));
fprintf('  3. Full (~%d runs, estimated %s)\n', fullRuns, format_duration(fullSec));
  
    while true
        s = input('Choose robustness level [0-3]: ', 's');
        s = strtrim(s);

        if isempty(s)
            mode = 0;
            return;
        end

        mode = str2double(s);
        if ismember(mode, [0 1 2 3])
            return;
        end

        fprintf('Please enter 0, 1, 2, or 3.\n');
    end
end


function [nRuns, estSeconds] = estimate_robustness_time(cfg, robustnessMode, mainTDATime, mainOptimizationTime)

    baseThr = cfg.optimization.elRankThreshold;
    baseMax = cfg.optimization.maxPeptides;

    switch robustnessMode
        case 0
            nRuns = 0;
            estSeconds = 0;

        case 1
            nRuns = 3;

        case 2
            nRuns = 10;

        case 3
            basePct = cfg.tda.defaultPercentile;
            percentiles = unique(max(1, min(99, [basePct-10, basePct, basePct+10])));
            rankThresholds = unique(max(0.1, min(100, [baseThr, 40])));
            maxPeptides = unique(max(1, round([baseMax, baseMax+2])));
            nWeightSets = 7;

            nRuns = numel(percentiles)*numel(rankThresholds)*numel(maxPeptides)*nWeightSets;

        otherwise
            error('Invalid robustness mode.');
    end

    estSeconds = nRuns * (mainTDATime + mainOptimizationTime);
end


function s = format_duration(sec)

    if sec < 60
        s = sprintf('%.0f sec', sec);
    elseif sec < 3600
        s = sprintf('%.1f min', sec/60);
    else
        s = sprintf('%.2f hr', sec/3600);
    end
end


function tdaResult = run_tda_wrapper(baseInput, tdaCfg)

    try
        tdaResult = Mobi_tda(baseInput.X, baseInput.sampleNames, tdaCfg);

        if ~isfield(tdaResult, 'chosenPercentile') && isfield(tdaCfg, 'percentile')
            tdaResult.chosenPercentile = tdaCfg.percentile;
        end

    catch ME
        error('run_tda_wrapper failed: %s', ME.message);
    end
end


function optResult = run_optimization_wrapper(baseInput, optCfg, tdaResult)

    suite = Mobi_run_optimization_suite(baseInput, optCfg, tdaResult, "robustness");

    optResult = struct();
    optResult.suite = suite;
    optResult.selectedPeptideUnionAcrossModes = suite.selectedPeptideUnionAcrossModes;
    optResult.practicalSelectedPeptideUnion = suite.practicalSelectedPeptideUnion;
    optResult.globalSelectedPeptides = suite.globalSelectedPeptides;
    optResult.honestSelectedPeptideUnion = suite.honestSelectedPeptideUnion;
    optResult.signaturePeptideMode = "selectedPeptideUnionAcrossModes";
end


% REFERENCE PANEL HELPERS
function referenceInfo = choose_global_reference_panel(cfg)

    topN = prompt_topN(cfg);
    T = Mobi_reference();
    validate_reference_vars(T, ["Locus","Allele","AlleleFrequency","SampleSize"]);

    loci = cfg.reference.defaultLoci;
    panel = build_common_hla_panel(T, loci, topN);

    fprintf('\nBuilt-in global/common panel summary:\n');
    print_generic_table(panel);

    referenceInfo = struct();
    referenceInfo.ModeLabel   = "Global/common Mobi reference";
    referenceInfo.FilterLabel = "All reference populations";
    referenceInfo.Panel       = panel;
    referenceInfo.HLAs        = string(panel.Allele);
    referenceInfo.ReferenceTable = T;
end


function referenceInfo = choose_reference_region_panel(cfg)

    regionOptions = string(cfg.reference.regionOptions(:));

    fprintf('\nBuilt-in reference regions:\n');
    for i = 1:numel(regionOptions)
        fprintf('  %2d. %s\n', i, regionOptions(i));
    end
    fprintf('\n');

    regionChoice = input('Choose region from built-in reference: ', 's');
    regionChoice = strtrim(regionChoice);
    idx = str2double(regionChoice);

    if isnan(idx) || idx < 1 || idx > numel(regionOptions) || mod(idx,1) ~= 0
        error('Invalid region choice.');
    end

    regionLabel = regionOptions(idx);
    topN = prompt_topN(cfg);

    T = Mobi_reference();
    validate_reference_vars(T, ["Locus","Allele","AlleleFrequency","SampleSize"]);

    regionCol = find_first_existing_var(T, ["Region","GeographicRegion","GeographicalRegion"]);

    if regionLabel ~= "All Regions"
        if regionCol == ""
            error(['Mobi_reference() does not contain a region column. ', ...
                   'Add a column named Region, GeographicRegion, or GeographicalRegion.']);
        end
        T = T(strcmpi(strtrim(string(T.(regionCol))), regionLabel), :);
    end

    if isempty(T)
        error('No reference rows remained after filtering to region: %s', regionLabel);
    end

    loci = cfg.reference.defaultLoci;
    panel = build_common_hla_panel(T, loci, topN);

    fprintf('\nBuilt-in region panel summary:\n');
    print_generic_table(panel);

    referenceInfo = struct();
    referenceInfo.ModeLabel   = "Region-filtered Mobi reference";
    referenceInfo.FilterLabel = regionLabel;
    referenceInfo.Panel       = panel;
    referenceInfo.HLAs        = string(panel.Allele);
    referenceInfo.ReferenceTable = T;
end


function referenceInfo = choose_reference_ancestry_panel(cfg)

    T = Mobi_reference();
    validate_reference_vars(T, ["Locus","Allele","AlleleFrequency","SampleSize","AncestryGroup"]);

    ancestryCol = "AncestryGroup";
    ancestryOptions = unique(strtrim(string(T.(ancestryCol))));
    ancestryOptions = ancestryOptions(ancestryOptions ~= "");
    ancestryOptions = sort(ancestryOptions);

    if ~any(strcmpi(ancestryOptions, "All Ancestry Groups"))
        ancestryOptions = [ancestryOptions; "All Ancestry Groups"];
    end

    fprintf('\nBuilt-in ancestry groups:\n');
    for i = 1:numel(ancestryOptions)
        fprintf('  %2d. %s\n', i, ancestryOptions(i));
    end
    fprintf('\n');

    ancestryChoice = input('Choose ancestry group from built-in reference: ', 's');
    ancestryChoice = strtrim(ancestryChoice);
    idx = str2double(ancestryChoice);

    if isnan(idx) || idx < 1 || idx > numel(ancestryOptions) || mod(idx,1) ~= 0
        error('Invalid ancestry group choice.');
    end

    ancestryLabel = ancestryOptions(idx);
    topN = prompt_topN(cfg);

    if ancestryLabel ~= "All Ancestry Groups"
        T = T(strcmpi(strtrim(string(T.(ancestryCol))), ancestryLabel), :);
    end

    if isempty(T)
        error('No reference rows remained after filtering to ancestry group: %s', ancestryLabel);
    end

    loci = cfg.reference.defaultLoci;
    panel = build_common_hla_panel(T, loci, topN);

    fprintf('\nBuilt-in ancestry panel summary:\n');
    print_generic_table(panel);

    referenceInfo = struct();
    referenceInfo.ModeLabel   = "Ancestry-filtered Mobi reference";
    referenceInfo.FilterLabel = ancestryLabel;
    referenceInfo.Panel       = panel;
    referenceInfo.HLAs        = string(panel.Allele);
    referenceInfo.ReferenceTable = T;
end


function topN = prompt_topN(cfg)
    prompt = sprintf('Top alleles per locus to use [default %d]: ', cfg.reference.defaultTopN);
    topNInput = input(prompt, 's');

    if isempty(strtrim(topNInput))
        topN = cfg.reference.defaultTopN;
    else
        topN = str2double(topNInput);
    end

    if isnan(topN) || topN <= 0 || mod(topN,1) ~= 0
        error('Top alleles per locus must be a positive integer.');
    end
end


function validate_reference_vars(T, requiredVars)
    for k = 1:numel(requiredVars)
        if ~any(strcmp(string(T.Properties.VariableNames), requiredVars(k)))
            error('Mobi_reference() must contain the variable "%s".', requiredVars(k));
        end
    end
end


function varName = find_first_existing_var(T, candidates)
    varName = "";
    for k = 1:numel(candidates)
        if any(strcmp(string(T.Properties.VariableNames), candidates(k)))
            varName = candidates(k);
            return;
        end
    end
end


function panel = build_common_hla_panel(T, loci, topN)

    T = T(~ismissing(T.Allele), :);
    T = T(~ismissing(T.Locus), :);

    if ~isnumeric(T.AlleleFrequency)
        T.AlleleFrequency = str2double(string(T.AlleleFrequency));
    end

    if ~isnumeric(T.SampleSize)
        T.SampleSize = str2double(string(T.SampleSize));
    end

    T = T(~isnan(T.AlleleFrequency), :);

    out = table();

    for L = 1:numel(loci)
        locus = string(loci(L));
        S = T(strcmpi(strtrim(string(T.Locus)), locus), :);

        if isempty(S)
            continue;
        end

        alleles = unique(string(S.Allele), 'stable');
        avgFreq  = nan(numel(alleles),1);
        nStudies = zeros(numel(alleles),1);
        effN     = zeros(numel(alleles),1);

        for a = 1:numel(alleles)
            idx = string(S.Allele) == alleles(a);
            F = S.AlleleFrequency(idx);

            if any(~isnan(S.SampleSize(idx)))
                W = S.SampleSize(idx);
                good = ~isnan(F) & ~isnan(W) & W > 0;
                if any(good)
                    avgFreq(a) = sum(F(good) .* W(good)) / sum(W(good));
                    effN(a) = sum(W(good));
                else
                    avgFreq(a) = mean(F, 'omitnan');
                    effN(a) = sum(~isnan(F));
                end
            else
                avgFreq(a) = mean(F, 'omitnan');
                effN(a) = sum(~isnan(F));
            end

            nStudies(a) = sum(~isnan(F));
        end

        Tloc = table( ...
            repmat(locus, numel(alleles), 1), ...
            alleles(:), ...
            avgFreq(:), ...
            nStudies(:), ...
            effN(:), ...
            'VariableNames', {'Locus','Allele','MeanAlleleFrequency','NumRows','EffectiveWeight'} );

        Tloc = sortrows(Tloc, {'MeanAlleleFrequency','EffectiveWeight'}, {'descend','descend'});
        Tloc = Tloc(1:min(topN, height(Tloc)), :);

        out = [out; Tloc]; %#ok<AGROW>
    end

    panel = out;
end


function optitypeLookup = build_optitype_lookup(optitypeTbl)

    optitypeLookup = containers.Map('KeyType', 'char', 'ValueType', 'any');

    if isempty(optitypeTbl)
        return;
    end

    for rowIdx = 1:height(optitypeTbl)
        row = optitypeTbl(rowIdx, :);

        typedHLAs = [ ...
            string(row.HLA1); ...
            string(row.HLA2); ...
            string(row.HLA3); ...
            string(row.HLA4); ...
            string(row.HLA5); ...
            string(row.HLA6)];

        typedHLAs = unique(Mobi_utils.normalize_hla_list(typedHLAs), 'stable');
        if isempty(typedHLAs)
            continue;
        end

        exactKey = normalize_optitype_file_key(string(row.FileName), false);
        stemKey  = normalize_optitype_file_key(string(row.FileName), true);

        if strlength(exactKey) > 0 && ~isKey(optitypeLookup, char(exactKey))
            optitypeLookup(char(exactKey)) = typedHLAs;
        end

        if strlength(stemKey) > 0 && ~isKey(optitypeLookup, char(stemKey))
            optitypeLookup(char(stemKey)) = typedHLAs;
        end
    end
end


function typedHLAs = get_optitype_hlas(optitypeLookup, sampleFileName)

    typedHLAs = strings(0,1);

    if isempty(optitypeLookup)
        return;
    end

    exactKey = normalize_optitype_file_key(sampleFileName, false);
    if strlength(exactKey) > 0 && isKey(optitypeLookup, char(exactKey))
        typedHLAs = optitypeLookup(char(exactKey));
        return;
    end

    stemKey = normalize_optitype_file_key(sampleFileName, true);
    if strlength(stemKey) > 0 && isKey(optitypeLookup, char(stemKey))
        typedHLAs = optitypeLookup(char(stemKey));
    end
end


function key = normalize_optitype_file_key(fileName, dropExtension)

    key = strtrim(string(fileName));
    key = erase(key, """");
    key = erase(key, "'");
    key = regexprep(key, '\s+', '');
    key = lower(key);

    if dropExtension && strlength(key) > 0
        [~, stem, ~] = fileparts(char(key));
        key = string(stem);
    end
end


function validate_optitype_spreadsheet_state(optiFullPath)

    [folder, name, ext] = fileparts(char(optiFullPath));
    lockPath = fullfile(folder, ['~$', name, ext]);

    if exist(lockPath, 'file')
        error(['The selected OptiType spreadsheet appears to still be open or was not closed cleanly: ', ...
            '%s\nClose Excel completely, delete the "~$" lock file if it remains, or export the sheet as CSV and select that instead.'], ...
            lockPath);
    end
end


% TABLE / PRINT HELPERS
function T = union_and_append_tables(T, A)

    if isempty(T)
        T = A;
        return;
    end

    if isempty(A)
        return;
    end

    varsT = string(T.Properties.VariableNames);
    varsA = string(A.Properties.VariableNames);
    allVars = unique([varsT, varsA], 'stable');

    T = add_missing_vars_as_string(T, allVars);
    A = add_missing_vars_as_string(A, allVars);

    T = T(:, cellstr(allVars));
    A = A(:, cellstr(allVars));

    T = [T; A];
end


function T = add_missing_vars_as_string(T, allVars)
    varsT = string(T.Properties.VariableNames);
    missing = allVars(~ismember(allVars, varsT));

    for i = 1:numel(missing)
        T.(missing(i)) = repmat("", height(T), 1);
    end
end


function print_string_list(vals, prefix)
    vals = string(vals(:));
    vals = vals(vals ~= "");
    for i = 1:numel(vals)
        fprintf('%s%s\n', prefix, vals(i));
    end
end


function print_two_column_pairs(header1, col1, header2, col2)
    col1 = string(col1(:));
    col2 = string(col2(:));

    width1 = max(strlength([string(header1); col1]));
    fmtHeader = sprintf('%%-%ds | %%s\\n', width1);
    fmtRow    = sprintf('%%-%ds | %%s\\n', width1);

    fprintf(fmtHeader, header1, header2);
    fprintf('%s-+-%s\n', repmat('-',1,width1), repmat('-',1,max(strlength([string(header2); col2]))));

    for i = 1:numel(col1)
        fprintf(fmtRow, col1(i), col2(i));
    end
end


function print_cluster_table(T)
    files = string(T.File);
    nums  = T.FileNumber;
    clus  = T.Cluster;

    widthNum  = max(strlength(string([0; nums])));
    widthFile = max(strlength(["File"; files]));
    widthClus = max(strlength(string([0; clus])));

    fmtHeader = sprintf('%%-%ds  %%-%ds  %%-%ds\\n', widthNum, widthFile, widthClus);
    fmtRow    = sprintf('%%-%dd  %%-%ds  %%-%dd\\n', widthNum, widthFile, widthClus);

    fprintf(fmtHeader, char("FileNumber"), char("File"), char("Cluster"));
    fprintf('%s  %s  %s\n', repmat('-',1,widthNum), repmat('-',1,widthFile), repmat('-',1,widthClus));

    for i = 1:height(T)
        fprintf(fmtRow, nums(i), char(files(i)), clus(i));
    end
end


function print_numeric_matrix(A, headers)
    if isempty(A)
        fprintf('  [empty]\n');
        return;
    end

    ncol = size(A,2);

    if nargin < 2 || isempty(headers)
        headers = strings(1,ncol);
        for j = 1:ncol
            headers(j) = "Col" + j;
        end
    else
        headers = string(headers(:))';
    end

    strA = strings(size(A));
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            if isinf(A(i,j))
                strA(i,j) = "Inf";
            elseif isnan(A(i,j))
                strA(i,j) = "NaN";
            else
                strA(i,j) = sprintf('%.6g', A(i,j));
            end
        end
    end

    widths = zeros(1,ncol);
    for j = 1:ncol
        widths(j) = max(strlength([headers(j); strA(:,j)]));
    end

    for j = 1:ncol
        fprintf('%-*s  ', widths(j), headers(j));
    end
    fprintf('\n');

    for j = 1:ncol
        fprintf('%s  ', repmat('-',1,widths(j)));
    end
    fprintf('\n');

    for i = 1:size(strA,1)
        for j = 1:ncol
            fprintf('%-*s  ', widths(j), strA(i,j));
        end
        fprintf('\n');
    end
end


function print_generic_table(T)
    if isempty(T) || height(T) == 0
        fprintf('  [empty table]\n');
        return;
    end

    varNames = string(T.Properties.VariableNames);
    ncol = numel(varNames);
    nrow = height(T);

    cellText = strings(nrow, ncol);

    for j = 1:ncol
        col = T.(varNames(j));

        if iscell(col)
            for i = 1:nrow
                if isstring(col{i}) || ischar(col{i})
                    cellText(i,j) = string(col{i});
                elseif isnumeric(col{i}) || islogical(col{i})
                    cellText(i,j) = mat2str(col{i});
                else
                    cellText(i,j) = "<cell>";
                end
            end
        elseif isstring(col) || ischar(col)
            cellText(:,j) = string(col);
        elseif isnumeric(col) || islogical(col)
            for i = 1:nrow
                if isscalar(col(i))
                    cellText(i,j) = sprintf('%.6g', col(i));
                else
                    cellText(i,j) = mat2str(col(i));
                end
            end
        elseif iscategorical(col)
            cellText(:,j) = string(col);
        else
            for i = 1:nrow
                try
                    cellText(i,j) = string(col(i));
                catch
                    cellText(i,j) = "<value>";
                end
            end
        end
    end

    widths = zeros(1, ncol);
    for j = 1:ncol
        widths(j) = max(strlength([varNames(j); cellText(:,j)]));
    end

    for j = 1:ncol
        fprintf('%-*s  ', widths(j), varNames(j));
    end
    fprintf('\n');

    for j = 1:ncol
        fprintf('%s  ', repmat('-',1,widths(j)));
    end
    fprintf('\n');

    for i = 1:nrow
        for j = 1:ncol
            fprintf('%-*s  ', widths(j), cellText(i,j));
        end
        fprintf('\n');
    end
end


function r = find_root(parent, x)
    r = x;
    while parent(r) ~= r
        r = parent(r);
    end
end


function s = group_to_string(v)
    v = v(:)';
    if isempty(v)
        s = "{}";
        return;
    end
    s = "{" + strjoin(string(v), ",") + "}";
end


function s = format_file_list(v)
    v = v(:)';
    if isempty(v)
        s = 'files none';
        return;
    end
    s = sprintf('files %s', strjoin(string(v), ','));
end


function s = format_name_list(vals)
    vals = string(vals(:));
    vals = vals(vals ~= "");
    if isempty(vals)
        s = 'files none';
        return;
    end
    s = sprintf('files %s', strjoin(vals, ','));
end


function s = format_numeric_list(vals)
    vals = vals(:)';
    if isempty(vals)
        s = '{}';
    else
        s = "{" + strjoin(string(vals), ",") + "}";
    end
end


function print_both(logFID, varargin)
    fprintf(varargin{:});
    log_only(logFID, varargin{:});
end


function print_text_both(logFID, txt)
    fprintf('%s', txt);
    log_only(logFID, '%s', txt);
end


function log_only(logFID, varargin)
    if logFID > 0
        fprintf(logFID, varargin{:});
    end
end


function safely_close_log(logFID)
    if logFID > 0
        fclose(logFID);
    end
end
