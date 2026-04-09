function robustness = Mobi_robustness(runFns, baseInput, cfg, robustnessMode)

    arguments
        runFns struct
        baseInput struct
        cfg struct
        robustnessMode
    end

    % Validate function handles
    requiredFns = ["runTDA", "runOptimization"];
    for k = 1:numel(requiredFns)
        if ~isfield(runFns, requiredFns(k)) || ~isa(runFns.(requiredFns(k)), 'function_handle')
            error('Mobi_run_robustness:MissingFunctionHandle', ...
                'runFns.%s must be provided as a function handle.', requiredFns(k));
        end
    end

    % Normalize robustness mode
    mode = i_normalizeMode(robustnessMode);

    fprintf('\n============================================================\n');
    fprintf('OPTIONAL ROBUSTNESS CHECK\n');
    fprintf('Mode: %s\n', mode.label);
    fprintf('============================================================\n');

    if mode.value == 0
        robustness = struct();
        robustness.mode = mode;
        robustness.didRun = false;
        robustness.message = 'Robustness mode set to none.';
        fprintf('Robustness check skipped.\n');
        return;
    end

    % Build baseline configs directly from cfg
    baselineTDA = struct();
    baselineTDA.makePlots       = false;  
    baselineTDA.distanceMetric  = cfg.tda.distanceMetric;
    baselineTDA.weights13       = cfg.tda.defaultWeights13;
    baselineTDA.thresholdRule   = cfg.tda.thresholdRule;
    baselineTDA.percentile      = cfg.tda.defaultPercentile;

    baselineOpt = struct();
    baselineOpt.maxPeptides       = cfg.optimization.maxPeptides;
    baselineOpt.elRankThreshold   = cfg.optimization.elRankThreshold;
    baselineOpt.familyMinOverlap  = cfg.optimization.familyMinOverlap;
    baselineOpt.redundancyMinOverlap = cfg.optimization.redundancyMinOverlap;
    baselineOpt.normalizeWeights  = cfg.optimization.normalizeWeights;
    baselineOpt.weights           = cfg.optimization.weights;

    if isfield(baseInput, 'optCfg') && isstruct(baseInput.optCfg)
        baselineOpt = i_overlayStructFields(baselineOpt, baseInput.optCfg);
    end

    % Build sweep grid from cfg + mode
    sweep = i_buildSweepFromConfig(cfg, mode, baselineTDA, baselineOpt);

    fprintf('Baseline percentile: %.2f\n', baselineTDA.percentile);
    fprintf('Baseline EL rank threshold: %.2f\n', baselineOpt.elRankThreshold);
    fprintf('Baseline max peptides: %d\n', baselineOpt.maxPeptides);
    fprintf('Total robustness runs (including baseline): %d\n', sweep.nRuns);

    % Baseline run
    fprintf('\nRunning baseline analysis...\n');
    tStartAll = tic;

    t0 = tic;
    baselineTDAResult = runFns.runTDA(baseInput, baselineTDA);
    baselineOptResult = runFns.runOptimization(baseInput, baselineOpt, baselineTDAResult);
    baselineElapsed = toc(t0);

    baselineSignature = i_buildRunSignature(baselineTDAResult, baselineOptResult);

    % Prepare storage
    runs(sweep.nRuns) = struct( ...
        'runID', [], ...
        'label', "", ...
        'tdaCfg', struct(), ...
        'optCfg', struct(), ...
        'tdaResult', struct(), ...
        'optResult', struct(), ...
        'summary', struct(), ...
        'elapsedSec', []);

    rowCounter = 1;

    % Store baseline row
    runs(rowCounter).runID      = rowCounter;
    runs(rowCounter).label      = "baseline";
    runs(rowCounter).tdaCfg     = baselineTDA;
    runs(rowCounter).optCfg     = baselineOpt;
    runs(rowCounter).tdaResult  = baselineTDAResult;
    runs(rowCounter).optResult  = baselineOptResult;
    runs(rowCounter).summary    = i_compareToBaseline( ...
        baselineTDAResult, baselineOptResult, baselineSignature, baselineSignature);
    runs(rowCounter).elapsedSec = baselineElapsed;
    rowCounter = rowCounter + 1;

    % Sweeps
    for iRun = 1:numel(sweep.runConfigs)

        runCfg = sweep.runConfigs(iRun);

        thisTDA = baselineTDA;
        thisOpt = baselineOpt;

        thisTDA.percentile        = runCfg.percentile;
        thisOpt.elRankThreshold   = runCfg.elRankThreshold;
        thisOpt.maxPeptides       = runCfg.maxPeptides;
        thisOpt.weights           = runCfg.weights;

        fprintf('\n[%d/%d] Robustness run: %s\n', rowCounter, sweep.nRuns, runCfg.label);
        fprintf('  percentile      = %.2f\n', thisTDA.percentile);
        fprintf('  EL rank cutoff  = %.2f\n', thisOpt.elRankThreshold);
        fprintf('  max peptides    = %d\n', thisOpt.maxPeptides);
        fprintf('  weights label   = %s\n', i_weightLabel(thisOpt.weights));

        tRun = tic;
        thisTDAResult = runFns.runTDA(baseInput, thisTDA);
        thisOptResult = runFns.runOptimization(baseInput, thisOpt, thisTDAResult);
        thisElapsed = toc(tRun);

        thisSignature = i_buildRunSignature(thisTDAResult, thisOptResult);
        thisSummary = i_compareToBaseline( ...
            thisTDAResult, thisOptResult, thisSignature, baselineSignature);

        runs(rowCounter).runID      = rowCounter;
        runs(rowCounter).label      = string(runCfg.label);
        runs(rowCounter).tdaCfg     = thisTDA;
        runs(rowCounter).optCfg     = thisOpt;
        runs(rowCounter).tdaResult  = thisTDAResult;
        runs(rowCounter).optResult  = thisOptResult;
        runs(rowCounter).summary    = thisSummary;
        runs(rowCounter).elapsedSec = thisElapsed;

        rowCounter = rowCounter + 1;
    end

    % Trim unused preallocated rows
    runs = runs(1:rowCounter-1);

    % Build summary table
    summaryTable = i_buildSummaryTable(runs);

    % Overall metrics
    stableClusterCount = sum(summaryTable.SameNumGroups);
    stablePeptideCount = sum(summaryTable.PeptideSetExactMatch);
    highPeptideOverlapCount = sum(summaryTable.PeptideJaccard >= 0.80);
    exactAssignmentCount = sum(summaryTable.AssignmentExactMatch);

    totalRuns = height(summaryTable);
    totalElapsed = toc(tStartAll);

    robustness = struct();
    robustness.mode = mode;
    robustness.didRun = true;
    robustness.baseline = struct( ...
        'tdaCfg', baselineTDA, ...
        'optCfg', baselineOpt, ...
        'tdaResult', baselineTDAResult, ...
        'optResult', baselineOptResult, ...
        'signature', baselineSignature);
    robustness.sweep = sweep;
    robustness.runs = runs;
    robustness.summaryTable = summaryTable;
    robustness.overall = struct();
    robustness.overall.totalRuns = totalRuns;
    robustness.overall.totalElapsedSec = totalElapsed;
    robustness.overall.sameNumGroupsCount = stableClusterCount;
    robustness.overall.sameNumGroupsFrac = stableClusterCount / totalRuns;
    robustness.overall.exactAssignmentCount = exactAssignmentCount;
    robustness.overall.exactAssignmentFrac = exactAssignmentCount / totalRuns;
    robustness.overall.exactPeptideSetCount = stablePeptideCount;
    robustness.overall.exactPeptideSetFrac = stablePeptideCount / totalRuns;
    robustness.overall.highPeptideOverlapCount = highPeptideOverlapCount;
    robustness.overall.highPeptideOverlapFrac = highPeptideOverlapCount / totalRuns;
    robustness.signaturePeptideMode = "selectedPeptideUnionAcrossModes";

    % Print compact report
    fprintf('\n============================================================\n');
    fprintf('ROBUSTNESS SUMMARY (%s)\n', upper(mode.label));
    fprintf('============================================================\n');
    fprintf('Total runs: %d\n', robustness.overall.totalRuns);
    fprintf('Same number of groups as baseline: %d/%d (%.1f%%)\n', ...
        stableClusterCount, totalRuns, 100*robustness.overall.sameNumGroupsFrac);
    fprintf('Exact cluster assignment match:    %d/%d (%.1f%%)\n', ...
        exactAssignmentCount, totalRuns, 100*robustness.overall.exactAssignmentFrac);
    fprintf('Exact peptide set match:           %d/%d (%.1f%%)\n', ...
        stablePeptideCount, totalRuns, 100*robustness.overall.exactPeptideSetFrac);
    fprintf('Peptide Jaccard >= 0.80:           %d/%d (%.1f%%)\n', ...
        highPeptideOverlapCount, totalRuns, 100*robustness.overall.highPeptideOverlapFrac);
    fprintf('Elapsed time: %.2f seconds\n', totalElapsed);

end


% MODE NORMALIZATION
function mode = i_normalizeMode(robustnessMode)

    if isnumeric(robustnessMode)
        switch robustnessMode
            case 0
                mode.value = 0; mode.label = "none";
            case 1
                mode.value = 1; mode.label = "fast";
            case 2
                mode.value = 2; mode.label = "standard";
            case 3
                mode.value = 3; mode.label = "full";
            otherwise
                error('Mobi_run_robustness:InvalidMode', ...
                    'Numeric robustnessMode must be 0, 1, 2, or 3.');
        end
        return;
    end

    if isstring(robustnessMode) || ischar(robustnessMode)
        s = lower(strtrim(string(robustnessMode)));
        switch s
            case "none"
                mode.value = 0; mode.label = "none";
            case "fast"
                mode.value = 1; mode.label = "fast";
            case "standard"
                mode.value = 2; mode.label = "standard";
            case "full"
                mode.value = 3; mode.label = "full";
            otherwise
                error('Mobi_run_robustness:InvalidMode', ...
                    'String robustnessMode must be "none", "fast", "standard", or "full".');
        end
        return;
    end

    error('Mobi_run_robustness:InvalidModeType', ...
        'robustnessMode must be numeric or string.');
end


function base = i_overlayStructFields(base, overrides)
    names = fieldnames(overrides);
    for k = 1:numel(names)
        base.(names{k}) = overrides.(names{k});
    end
end


% BUILD SWEEP FROM CONFIG
function sweep = i_buildSweepFromConfig(~, mode, baselineTDA, baselineOpt)

    basePct = baselineTDA.percentile;
    baseThr = baselineOpt.elRankThreshold;
    baseMax = baselineOpt.maxPeptides;
    baseW   = baselineOpt.weights;
    sweep = struct();
    sweep.modeLabel = mode.label;

    switch mode.value
        case 1
            % Fast = TDA threshold sensitivity only.
            sweep.percentiles = unique(i_clipNumeric([basePct-10, basePct+10], 1, 99));
            sweep.rankThresholds = baseThr;
            sweep.maxPeptideOptions = baseMax;
            sweep.weightSets = {baseW};

        case 2
            % Standard = one-factor-at-a-time sensitivity.
            sweep.percentiles = unique(i_clipNumeric([basePct-10, basePct, basePct+10], 1, 99));
            sweep.rankThresholds = unique(i_clipNumeric([baseThr-7.7, baseThr, baseThr+7.3], 0.1, 100));
            sweep.maxPeptideOptions = unique(max(1, round([baseMax-2, baseMax, baseMax+2])));
            sweep.weightSets = i_makeWeightSetCollection(baseW, "partial");

        case 3
            % Full = exhaustive cartesian sweep.
            sweep.percentiles = unique(i_clipNumeric([basePct-10, basePct, basePct+10], 1, 99));
            sweep.rankThresholds = unique(i_clipNumeric([baseThr, 40], 0.1, 100));
            sweep.maxPeptideOptions = unique(max(1, round([baseMax, baseMax+2])));
            sweep.weightSets = i_makeWeightSetCollection(baseW, "full");
        
        otherwise
            error('Sweep builder called with unsupported mode.');
    end

    sweep.runConfigs = i_buildRunConfigList(sweep, baselineTDA, baselineOpt);
    sweep.nRuns = 1 + numel(sweep.runConfigs);

end


function runConfigs = i_buildRunConfigList(sweep, baselineTDA, baselineOpt)

    runConfigs = struct( ...
        'label', {}, ...
        'percentile', {}, ...
        'elRankThreshold', {}, ...
        'maxPeptides', {}, ...
        'weights', {});

    basePct = baselineTDA.percentile;
    baseThr = baselineOpt.elRankThreshold;
    baseMax = baselineOpt.maxPeptides;
    baseW = baselineOpt.weights;

    switch sweep.modeLabel
        case "fast"
            for i = 1:numel(sweep.percentiles)
                pctVal = sweep.percentiles(i);
                if pctVal ~= basePct
                    runConfigs(end+1) = i_runConfig( ... %#ok<AGROW>
                        sprintf('percentile %.2f', pctVal), pctVal, baseThr, baseMax, baseW);
                end
            end

        case "standard"
            for i = 1:numel(sweep.percentiles)
                pctVal = sweep.percentiles(i);
                if pctVal ~= basePct
                    runConfigs(end+1) = i_runConfig( ... %#ok<AGROW>
                        sprintf('percentile %.2f', pctVal), pctVal, baseThr, baseMax, baseW);
                end
            end

            for i = 1:numel(sweep.rankThresholds)
                thrVal = sweep.rankThresholds(i);
                if thrVal ~= baseThr
                    runConfigs(end+1) = i_runConfig( ... %#ok<AGROW>
                        sprintf('EL rank %.2f', thrVal), basePct, thrVal, baseMax, baseW);
                end
            end

            for i = 1:numel(sweep.maxPeptideOptions)
                maxVal = sweep.maxPeptideOptions(i);
                if maxVal ~= baseMax
                    runConfigs(end+1) = i_runConfig( ... %#ok<AGROW>
                        sprintf('max peptides %d', maxVal), basePct, baseThr, maxVal, baseW);
                end
            end

            for i = 1:numel(sweep.weightSets)
                wStruct = sweep.weightSets{i};
                if ~i_sameWeights(wStruct, baseW)
                    runConfigs(end+1) = i_runConfig( ... %#ok<AGROW>
                        sprintf('weights %s', i_weightLabel(wStruct)), basePct, baseThr, baseMax, wStruct);
                end
            end

        case "full"
            for iPct = 1:numel(sweep.percentiles)
                for iThr = 1:numel(sweep.rankThresholds)
                    for iMax = 1:numel(sweep.maxPeptideOptions)
                        for iW = 1:numel(sweep.weightSets)
                            pctVal = sweep.percentiles(iPct);
                            thrVal = sweep.rankThresholds(iThr);
                            maxVal = sweep.maxPeptideOptions(iMax);
                            wStruct = sweep.weightSets{iW};

                            if ~i_isSameConfig(pctVal, thrVal, maxVal, wStruct, baselineTDA, baselineOpt)
                                runConfigs(end+1) = i_runConfig( ... %#ok<AGROW>
                                    "full cartesian", pctVal, thrVal, maxVal, wStruct);
                            end
                        end
                    end
                end
            end
    end
end


function cfg = i_runConfig(label, pctVal, thrVal, maxVal, wStruct)
    cfg = struct();
    cfg.label = string(label);
    cfg.percentile = pctVal;
    cfg.elRankThreshold = thrVal;
    cfg.maxPeptides = maxVal;
    cfg.weights = wStruct;
end


% BUILD WEIGHT PERTURBATIONS
function weightSets = i_makeWeightSetCollection(baseW, level)

    % Convert to vector in fixed order for perturbation
    names = i_weightNames();
    baseVec = zeros(1, numel(names));
    for k = 1:numel(names)
        baseVec(k) = baseW.(names{k});
    end

    W = {};
    W{end+1} = i_vectorToWeightStruct(baseVec, names);

    switch lower(level)
        case "partial"
            % Slightly emphasize patient coverage
            v1 = baseVec;
            v1(1) = v1(1) + 0.08;
            v1(3) = max(0, v1(3) - 0.04);
            v1(6) = max(0, v1(6) - 0.04);
            W{end+1} = i_vectorToWeightStruct(i_normalizeIfNeeded(v1), names);

            % Slightly emphasize HLA coverage
            v2 = baseVec;
            v2(2) = v2(2) + 0.08;
            v2(5) = max(0, v2(5) - 0.04);
            v2(6) = max(0, v2(6) - 0.04);
            W{end+1} = i_vectorToWeightStruct(i_normalizeIfNeeded(v2), names);

            % Slightly emphasize family novelty / reduce redundancy
            v3 = baseVec;
            v3(5) = v3(5) + 0.06;
            v3(6) = max(0, v3(6) - 0.06);
            W{end+1} = i_vectorToWeightStruct(i_normalizeIfNeeded(v3), names);

        case "full"
            % Add the partial set first
            W = i_makeWeightSetCollection(baseW, "partial");

            % More binding-heavy
            v4 = baseVec;
            v4(3) = v4(3) + 0.10;
            v4(1) = max(0, v4(1) - 0.05);
            v4(4) = max(0, v4(4) - 0.05);
            W{end+1} = i_vectorToWeightStruct(i_normalizeIfNeeded(v4), names);

            % More prevalence-heavy
            v5 = baseVec;
            v5(4) = v5(4) + 0.10;
            v5(5) = max(0, v5(5) - 0.05);
            v5(6) = max(0, v5(6) - 0.05);
            W{end+1} = i_vectorToWeightStruct(i_normalizeIfNeeded(v5), names);

            % More redundancy penalty
            v6 = baseVec;
            v6(6) = v6(6) + 0.10;
            v6(1) = max(0, v6(1) - 0.05);
            v6(5) = max(0, v6(5) - 0.05);
            W{end+1} = i_vectorToWeightStruct(i_normalizeIfNeeded(v6), names);

        otherwise
            error('Unknown weight perturbation level.');
    end

    weightSets = i_removeDuplicateWeightStructs(W);
end


function names = i_weightNames()
    names = { ...
        'wPatientCoverage', ...
        'wHLACoverage', ...
        'wBinding', ...
        'wPrevalence', ...
        'wFamilyNovelty', ...
        'wRedundancy'};
end


function s = i_vectorToWeightStruct(v, names)
    s = struct();
    for k = 1:numel(names)
        s.(names{k}) = v(k);
    end
end


function v = i_normalizeIfNeeded(v)
    v(v < 0) = 0;
    s = sum(v);
    if s <= 0
        error('Weight vector became all zeros during robustness perturbation.');
    end
    v = v / s;
end


function W = i_removeDuplicateWeightStructs(W)
    keep = true(size(W));
    for i = 1:numel(W)
        vi = i_weightStructToVector(W{i});
        for j = i+1:numel(W)
            vj = i_weightStructToVector(W{j});
            if max(abs(vi - vj)) < 1e-12
                keep(j) = false;
            end
        end
    end
    W = W(keep);
end


function v = i_weightStructToVector(s)
    names = i_weightNames();
    v = zeros(1, numel(names));
    for k = 1:numel(names)
        v(k) = s.(names{k});
    end
end


function tf = i_sameWeights(a, b)
    tf = max(abs(i_weightStructToVector(a) - i_weightStructToVector(b))) < 1e-12;
end


function label = i_weightLabel(w)
    v = i_weightStructToVector(w);
    label = sprintf('[%.2f %.2f %.2f %.2f %.2f %.2f]', v);
end


% BASELINE COMPARISON
function signature = i_buildRunSignature(tdaResult, optResult)

    signature = struct();

    signature.groupAssignments = i_extractAssignments(tdaResult);
    signature.numGroups        = i_extractNumGroups(tdaResult, signature.groupAssignments);
    signature.groupSizes       = i_extractGroupSizes(tdaResult, signature.groupAssignments);
    signature.peptides         = i_extractPeptides(optResult);

end


function summary = i_compareToBaseline(tdaResult, optResult, thisSig, baseSig)

    summary = struct();

    summary.numGroups = i_extractNumGroups(tdaResult, thisSig.groupAssignments);
    summary.groupSizes = i_extractGroupSizes(tdaResult, thisSig.groupAssignments);

    summary.sameNumGroups = isequaln(thisSig.numGroups, baseSig.numGroups);

    if isempty(thisSig.groupAssignments) || isempty(baseSig.groupAssignments)
        summary.assignmentExactMatch = false;
        summary.assignmentAgreement = NaN;
    else
        summary.assignmentExactMatch = isequaln(thisSig.groupAssignments(:), baseSig.groupAssignments(:));
        summary.assignmentAgreement = mean(thisSig.groupAssignments(:) == baseSig.groupAssignments(:), 'omitnan');
    end

    [jaccardPeptides, exactPeptideMatch, nPeptides] = i_comparePeptideSets(thisSig.peptides, baseSig.peptides);
    summary.peptideJaccard = jaccardPeptides;
    summary.peptideSetExactMatch = exactPeptideMatch;
    summary.numSelectedPeptides = nPeptides;

    % Pull optional score if available
    summary.coverageScore = i_extractScalarIfPresent(optResult, ...
        {'coverageScore','overallCoverage','finalScore','objectiveValue','bestScore'});

    % High-level stability flag
    summary.stableOverall = ...
        summary.sameNumGroups && ...
        (summary.peptideJaccard >= 0.80 || summary.peptideSetExactMatch);

end


% EXTRACTORS
function assignments = i_extractAssignments(tdaResult)

    assignments = [];

    candidateFields = { ...
        'optimizationGroups', ...
        'optimizationGroupLabels', ...
        'clusterAssignments', ...
        'groupAssignments', ...
        'clusters', ...
        'groupLabels', ...
        'rawClusters', ...
        'assignments', ...
        'labels', ...
        'idx'};

    for k = 1:numel(candidateFields)
        f = candidateFields{k};
        if isstruct(tdaResult) && isfield(tdaResult, f)
            x = tdaResult.(f);
            if isnumeric(x) || islogical(x)
                assignments = x(:);
                return;
            end
        end
    end

    % If explicit assignments not present, try reconstructing from groups/cells
    groupFieldNames = {'groups','clusters','connectedComponents'};
    for k = 1:numel(groupFieldNames)
        f = groupFieldNames{k};
        if isstruct(tdaResult) && isfield(tdaResult, f)
            G = tdaResult.(f);
            if iscell(G)
                % If groups are cell arrays of sample indices, reconstruct labels
                maxIdx = 0;
                for i = 1:numel(G)
                    if isnumeric(G{i}) && ~isempty(G{i})
                        maxIdx = max(maxIdx, max(G{i}));
                    end
                end
                if maxIdx > 0
                    assignments = NaN(maxIdx,1);
                    for i = 1:numel(G)
                        if isnumeric(G{i})
                            assignments(G{i}) = i;
                        end
                    end
                    return;
                end
            end
        end
    end
end


function numGroups = i_extractNumGroups(tdaResult, assignments)

    numGroups = NaN;

    candidateFields = {'numGroups','nGroups','numClusters','nClusters'};
    for k = 1:numel(candidateFields)
        f = candidateFields{k};
        if isstruct(tdaResult) && isfield(tdaResult, f) && isnumeric(tdaResult.(f)) && isscalar(tdaResult.(f))
            numGroups = tdaResult.(f);
            return;
        end
    end

    if ~isempty(assignments)
        valid = assignments(~isnan(assignments));
        if ~isempty(valid)
            numGroups = numel(unique(valid));
            return;
        end
    end
end


function groupSizes = i_extractGroupSizes(tdaResult, assignments)

    groupSizes = [];

    candidateFields = {'optimizationGroupSizes','groupSizes','clusterSizes','componentSizes','rawComponentSizes'};
    for k = 1:numel(candidateFields)
        f = candidateFields{k};
        if isstruct(tdaResult) && isfield(tdaResult, f) && isnumeric(tdaResult.(f))
            groupSizes = sort(tdaResult.(f)(:).', 'descend');
            return;
        end
    end

    if ~isempty(assignments)
        valid = assignments(~isnan(assignments));
        if ~isempty(valid)
            u = unique(valid);
            sz = zeros(size(u));
            for i = 1:numel(u)
                sz(i) = sum(valid == u(i));
            end
            groupSizes = sort(sz(:).', 'descend');
        end
    end
end


function peptides = i_extractPeptides(optResult)

    peptides = strings(0,1);

    % Direct string/cell fields
    directFields = { ...
        'selectedPeptideUnionAcrossModes', ...
        'practicalSelectedPeptideUnion', ...
        'globalSelectedPeptides', ...
        'honestSelectedPeptideUnion', ...
        'selectedPeptides', ...
        'bestPeptides', ...
        'chosenPeptides'};
    for k = 1:numel(directFields)
        f = directFields{k};
        if isstruct(optResult) && isfield(optResult, f)
            x = optResult.(f);
            peptides = i_toStringColumn(x);
            peptides = unique(peptides(peptides ~= ""));
            if ~isempty(peptides)
                return;
            end
        end
    end

    % Table fields with likely peptide column
    tableFields = {'selectedPeptides','selectedTable','finalTable','bestTable','chosenTable'};
    for k = 1:numel(tableFields)
        f = tableFields{k};
        if isstruct(optResult) && isfield(optResult, f) && istable(optResult.(f))
            T = optResult.(f);
            peptideVars = {'Peptide','peptide','SelectedPeptide'};
            for j = 1:numel(peptideVars)
                if any(strcmp(T.Properties.VariableNames, peptideVars{j}))
                    peptides = i_toStringColumn(T.(peptideVars{j}));
                    peptides = unique(peptides(peptides ~= ""));
                    if ~isempty(peptides)
                        return;
                    end
                end
            end
        end
    end
end


function out = i_extractScalarIfPresent(S, fieldNames)
    out = NaN;
    if ~isstruct(S)
        return;
    end
    for k = 1:numel(fieldNames)
        f = fieldNames{k};
        if isfield(S, f) && isnumeric(S.(f)) && isscalar(S.(f))
            out = S.(f);
            return;
        end
    end
end


function [jaccardVal, exactMatch, nThis] = i_comparePeptideSets(thisPeptides, basePeptides)

    thisPeptides = unique(i_toStringColumn(thisPeptides));
    basePeptides = unique(i_toStringColumn(basePeptides));

    thisPeptides = thisPeptides(thisPeptides ~= "");
    basePeptides = basePeptides(basePeptides ~= "");

    nThis = numel(thisPeptides);

    if isempty(thisPeptides) && isempty(basePeptides)
        jaccardVal = 1;
        exactMatch = true;
        return;
    end

    inter = numel(intersect(thisPeptides, basePeptides));
    uni   = numel(union(thisPeptides, basePeptides));

    if uni == 0
        jaccardVal = NaN;
    else
        jaccardVal = inter / uni;
    end

    exactMatch = isequal(thisPeptides(:), basePeptides(:));
end


function s = i_toStringColumn(x)
    if isstring(x)
        s = x(:);
    elseif ischar(x)
        s = string({x}).';
    elseif iscellstr(x)
        s = string(x(:));
    elseif iscell(x)
        s = string(x(:));
    else
        s = strings(0,1);
    end
end


% SUMMARY TABLE
function T = i_buildSummaryTable(runs)

    n = numel(runs);

    RunID = zeros(n,1);
    Label = strings(n,1);
    Percentile = NaN(n,1);
    ELRankThreshold = NaN(n,1);
    MaxPeptides = NaN(n,1);
    WeightLabel = strings(n,1);
    NumGroups = NaN(n,1);
    GroupSizes = strings(n,1);
    SameNumGroups = false(n,1);
    AssignmentExactMatch = false(n,1);
    AssignmentAgreement = NaN(n,1);
    NumSelectedPeptides = NaN(n,1);
    PeptideJaccard = NaN(n,1);
    PeptideSetExactMatch = false(n,1);
    CoverageScore = NaN(n,1);
    StableOverall = false(n,1);
    ElapsedSec = NaN(n,1);

    for i = 1:n
        RunID(i) = runs(i).runID;
        Label(i) = string(runs(i).label);
        Percentile(i) = runs(i).tdaCfg.percentile;
        ELRankThreshold(i) = runs(i).optCfg.elRankThreshold;
        MaxPeptides(i) = runs(i).optCfg.maxPeptides;
        WeightLabel(i) = string(i_weightLabel(runs(i).optCfg.weights));
        NumGroups(i) = runs(i).summary.numGroups;
        GroupSizes(i) = string(mat2str(runs(i).summary.groupSizes));
        SameNumGroups(i) = logical(runs(i).summary.sameNumGroups);
        AssignmentExactMatch(i) = logical(runs(i).summary.assignmentExactMatch);
        AssignmentAgreement(i) = runs(i).summary.assignmentAgreement;
        NumSelectedPeptides(i) = runs(i).summary.numSelectedPeptides;
        PeptideJaccard(i) = runs(i).summary.peptideJaccard;
        PeptideSetExactMatch(i) = logical(runs(i).summary.peptideSetExactMatch);
        CoverageScore(i) = runs(i).summary.coverageScore;
        StableOverall(i) = logical(runs(i).summary.stableOverall);
        ElapsedSec(i) = runs(i).elapsedSec;
    end

    T = table( ...
        RunID, Label, Percentile, ELRankThreshold, MaxPeptides, WeightLabel, ...
        NumGroups, GroupSizes, SameNumGroups, AssignmentExactMatch, ...
        AssignmentAgreement, NumSelectedPeptides, PeptideJaccard, ...
        PeptideSetExactMatch, CoverageScore, StableOverall, ElapsedSec);
end


% UTILITIES
function tf = i_isSameConfig(pctVal, thrVal, maxVal, wStruct, baselineTDA, baselineOpt)
    tf = abs(pctVal - baselineTDA.percentile) < 1e-12 && ...
         abs(thrVal - baselineOpt.elRankThreshold) < 1e-12 && ...
         isequal(maxVal, baselineOpt.maxPeptides) && ...
         max(abs(i_weightStructToVector(wStruct) - i_weightStructToVector(baselineOpt.weights))) < 1e-12;
end


function x = i_clipNumeric(x, lo, hi)
    x = max(lo, min(hi, x));
end
