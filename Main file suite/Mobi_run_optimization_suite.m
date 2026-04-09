function suite = Mobi_run_optimization_suite(baseInput, optCfg, tdaResult, runMode)
% Mobi_run_optimization_suite
% Canonical execution path for Mobi peptide optimization modes.

    if nargin < 4 || isempty(runMode)
        runMode = "all";
    end

    runMode = lower(strtrim(string(runMode)));

    sampleFiles = string(baseInput.sampleFiles(:));

    if isfield(baseInput, 'sampleNames') && ~isempty(baseInput.sampleNames)
        sampleNames = string(baseInput.sampleNames(:));
    else
        sampleNames = sampleFiles;
    end

    if numel(sampleFiles) ~= numel(sampleNames)
        error('baseInput.sampleFiles and baseInput.sampleNames must have the same length.');
    end

    if ~isfield(optCfg, 'maxPeptides') || isempty(optCfg.maxPeptides)
        cfg = Mobi_config();
        optCfg.maxPeptides = cfg.optimization.maxPeptides;
    end

    if ~isfield(optCfg, 'weights') || isempty(optCfg.weights)
        cfg = Mobi_config();
        optCfg.weights = cfg.optimization.weights;
    end

    minGroupSize = 1;
    if isfield(optCfg, 'minGroupSizeForPracticalOptimization') && ...
            ~isempty(optCfg.minGroupSizeForPracticalOptimization)
        minGroupSize = optCfg.minGroupSizeForPracticalOptimization;
    end

    runGlobal = ismember(runMode, ["all", "main", "global_honest_practical", "global", "robustness"]);
    runHonest = ismember(runMode, ["all", "main", "global_honest_practical", "honest"]);
    runPractical = ismember(runMode, ["all", "main", "global_honest_practical", "practical", "robustness"]);

    suite = struct();
    suite.didRun = false;
    suite.runMode = runMode;
    suite.optCfg = optCfg;
    suite.global = empty_group_result();
    suite.honestGroups = empty_group_result_array();
    suite.practicalGroups = empty_group_result_array();
    suite.allResults = empty_all_results();

    resultCounter = 0;

    if runGlobal
        groupFiles = cellstr(sampleFiles);
        res = Mobi_optimization(groupFiles, optCfg.maxPeptides, optCfg.weights, optCfg);
        res = attach_reference_metrics_if_available(res, optCfg);
        suite.global = make_group_result("global", "Global cohort", 0, groupFiles, sampleNames, res);
        resultCounter = resultCounter + 1;
        suite.allResults(resultCounter) = group_result_to_all_result(suite.global);
    end

    if runHonest
        rawLabels = get_raw_group_labels(tdaResult);
        rawGroups = unique(rawLabels);

        for g = 1:numel(rawGroups)
            groupLabel = rawGroups(g);
            idx = rawLabels == groupLabel;
            groupFiles = cellstr(sampleFiles(idx));
            groupNames = sampleNames(idx);
            res = Mobi_optimization(groupFiles, optCfg.maxPeptides, optCfg.weights, optCfg);
            res = attach_reference_metrics_if_available(res, optCfg);
            groupResult = make_group_result( ...
                "honest", sprintf('Honest raw Group %d', groupLabel), groupLabel, ...
                groupFiles, groupNames, res);

            suite.honestGroups(end+1) = groupResult;
            resultCounter = resultCounter + 1;
            suite.allResults(resultCounter) = group_result_to_all_result(groupResult);
        end
    end

    if runPractical
        practicalLabels = get_practical_group_labels(tdaResult);
        practicalGroups = unique(practicalLabels);

        for g = 1:numel(practicalGroups)
            groupLabel = practicalGroups(g);
            idx = practicalLabels == groupLabel;

            if sum(idx) < minGroupSize
                groupResult = make_group_result( ...
                    "practical", sprintf('Practical Group %d', groupLabel), groupLabel, ...
                    cellstr(sampleFiles(idx)), sampleNames(idx), struct());
                groupResult.didRun = false;
                groupResult.skipReason = sprintf('Group size %d < minimum %d', sum(idx), minGroupSize);
            else
                groupFiles = cellstr(sampleFiles(idx));
                groupNames = sampleNames(idx);
                res = Mobi_optimization(groupFiles, optCfg.maxPeptides, optCfg.weights, optCfg);
                res = attach_reference_metrics_if_available(res, optCfg);
                groupResult = make_group_result( ...
                    "practical", sprintf('Practical Group %d', groupLabel), groupLabel, ...
                    groupFiles, groupNames, res);
            end

            suite.practicalGroups(end+1) = groupResult;

            if groupResult.didRun
                resultCounter = resultCounter + 1;
                suite.allResults(resultCounter) = group_result_to_all_result(groupResult);
            end
        end
    end

    suite.globalSelectedPeptides = extract_peptides_from_group_result(suite.global);
    suite.honestSelectedPeptidesByGroup = build_peptides_by_group_table(suite.honestGroups);
    suite.practicalSelectedPeptidesByGroup = build_peptides_by_group_table(suite.practicalGroups);
    suite.honestSelectedPeptideUnion = extract_union_from_groups(suite.honestGroups);
    suite.practicalSelectedPeptideUnion = extract_union_from_groups(suite.practicalGroups);
    suite.selectedPeptideUnionAcrossModes = unique([ ...
        suite.globalSelectedPeptides(:); ...
        suite.honestSelectedPeptideUnion(:); ...
        suite.practicalSelectedPeptideUnion(:)], 'stable');
    suite.allResultsTable = build_all_results_table(suite.allResults);
    suite.summary = build_suite_summary(suite);
    suite.didRun = runGlobal || runHonest || runPractical;
    Mobi_validate.optimization_suite_output(suite);
end


function labels = get_raw_group_labels(tdaResult)
    if isfield(tdaResult, 'rawClusters')
        labels = tdaResult.rawClusters(:);
    elseif isfield(tdaResult, 'groupLabels')
        labels = tdaResult.groupLabels(:);
    elseif isfield(tdaResult, 'clusters')
        labels = tdaResult.clusters(:);
    else
        error('tdaResult must contain rawClusters, groupLabels, or clusters.');
    end
end


function labels = get_practical_group_labels(tdaResult)
    if isfield(tdaResult, 'optimizationGroups')
        labels = tdaResult.optimizationGroups(:);
    else
        labels = get_raw_group_labels(tdaResult);
    end
end


function out = empty_group_result()
    out = struct( ...
        'mode', "", ...
        'label', "", ...
        'group', NaN, ...
        'files', {{}}, ...
        'selectedFileNames', strings(0,1), ...
        'result', struct(), ...
        'didRun', false, ...
        'skipReason', "");
end


function out = empty_group_result_array()
    out = repmat(empty_group_result(), 0, 1);
end


function out = empty_all_results()
    out = repmat(struct( ...
        'mode', "", ...
        'label', "", ...
        'group', NaN, ...
        'files', {{}}, ...
        'selectedFileNames', strings(0,1), ...
        'result', struct()), 0, 1);
end


function out = make_group_result(mode, label, group, files, selectedFileNames, result)
    out = empty_group_result();
    out.mode = string(mode);
    out.label = string(label);
    out.group = group;
    out.files = files;
    out.selectedFileNames = string(selectedFileNames(:));
    out.result = result;
    out.didRun = true;
    out.skipReason = "";
end


function out = group_result_to_all_result(groupResult)
    out = struct();
    out.mode = groupResult.mode;
    out.label = groupResult.label;
    out.group = groupResult.group;
    out.files = groupResult.files;
    out.selectedFileNames = groupResult.selectedFileNames;
    out.result = groupResult.result;
end


function peptides = extract_peptides_from_group_result(groupResult)
    peptides = strings(0,1);

    if ~isstruct(groupResult) || ~isfield(groupResult, 'didRun') || ~groupResult.didRun
        return;
    end

    if isfield(groupResult.result, 'selectedPeptides') && ...
            istable(groupResult.result.selectedPeptides) && ...
            any(strcmp(groupResult.result.selectedPeptides.Properties.VariableNames, 'Peptide'))
        peptides = unique(string(groupResult.result.selectedPeptides.Peptide), 'stable');
        peptides = peptides(peptides ~= "");
    end
end


function peptides = extract_union_from_groups(groups)
    peptides = strings(0,1);
    for i = 1:numel(groups)
        peptides = [peptides; extract_peptides_from_group_result(groups(i))]; %#ok<AGROW>
    end
    peptides = unique(peptides(peptides ~= ""), 'stable');
end


function T = build_peptides_by_group_table(groups)
    Mode = strings(0,1);
    Group = zeros(0,1);
    Label = strings(0,1);
    Peptide = strings(0,1);

    for i = 1:numel(groups)
        peptides = extract_peptides_from_group_result(groups(i));
        for p = 1:numel(peptides)
            Mode(end+1,1) = groups(i).mode; %#ok<AGROW>
            Group(end+1,1) = groups(i).group; %#ok<AGROW>
            Label(end+1,1) = groups(i).label; %#ok<AGROW>
            Peptide(end+1,1) = peptides(p); %#ok<AGROW>
        end
    end

    T = table(Mode, Group, Label, Peptide);
end


function T = build_all_results_table(allResults)
    Mode = strings(0,1);
    Group = zeros(0,1);
    Label = strings(0,1);
    NumFiles = zeros(0,1);
    NumSelectedPeptides = zeros(0,1);

    for i = 1:numel(allResults)
        Mode(end+1,1) = string(allResults(i).mode); %#ok<AGROW>
        Group(end+1,1) = allResults(i).group; %#ok<AGROW>
        Label(end+1,1) = string(allResults(i).label); %#ok<AGROW>
        NumFiles(end+1,1) = numel(allResults(i).files); %#ok<AGROW>
        NumSelectedPeptides(end+1,1) = numel(extract_peptides_from_all_result(allResults(i))); %#ok<AGROW>
    end

    T = table(Mode, Group, Label, NumFiles, NumSelectedPeptides);
end


function peptides = extract_peptides_from_all_result(oneResult)
    groupResult = empty_group_result();
    groupResult.didRun = true;
    groupResult.result = oneResult.result;
    peptides = extract_peptides_from_group_result(groupResult);
end


function summary = build_suite_summary(suite)
    summary = struct();
    summary.numGlobalSelectedPeptides = numel(suite.globalSelectedPeptides);
    summary.numHonestSelectedPeptides = numel(suite.honestSelectedPeptideUnion);
    summary.numPracticalSelectedPeptides = numel(suite.practicalSelectedPeptideUnion);
    summary.numSelectedPeptidesAcrossModes = numel(suite.selectedPeptideUnionAcrossModes);
    summary.numHonestGroupsRun = sum([suite.honestGroups.didRun]);
    summary.numPracticalGroupsRun = sum([suite.practicalGroups.didRun]);
end


function result = attach_reference_metrics_if_available(result, optCfg)
    if ~isstruct(result)
        return;
    end

    if ~isfield(optCfg, 'referenceEvaluation') || ~isstruct(optCfg.referenceEvaluation)
        return;
    end

    refEval = optCfg.referenceEvaluation;
    if ~isfield(refEval, 'enabled') || ~refEval.enabled
        return;
    end

    if ~isfield(refEval, 'referenceTable') || isempty(refEval.referenceTable)
        return;
    end

    selectedPeptides = strings(0,1);
    if isfield(result, 'selectedPeptides') && istable(result.selectedPeptides) && ...
            any(strcmp(result.selectedPeptides.Properties.VariableNames, 'Peptide'))
        selectedPeptides = string(result.selectedPeptides.Peptide);
    end

    if ~isfield(result, 'combinedTable') || ~istable(result.combinedTable)
        return;
    end

    refOptions = struct();
    if isfield(refEval, 'populationWeightMode')
        refOptions.populationWeightMode = refEval.populationWeightMode;
    end
    if isfield(refEval, 'hitCountMode')
        refOptions.hitCountMode = refEval.hitCountMode;
    end

    result.referenceMetrics = Mobi_reference_population_metrics( ...
        selectedPeptides, result.combinedTable, refEval.referenceTable, refOptions);
end
