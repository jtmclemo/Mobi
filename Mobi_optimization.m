function result = Mobi_optimization(groupCSVFiles, maxPeptides, WEIGHTS, optCfg)

    cfg = Mobi_config();

    % DEFAULT INPUTS
    if nargin < 2 || isempty(maxPeptides)
        maxPeptides = cfg.optimization.maxPeptides;
    end

    if nargin < 3 || isempty(WEIGHTS)
        WEIGHTS = cfg.optimization.weights;
    end

    if nargin < 4 || isempty(optCfg)
        optCfg = struct();
    end

    % GLOBAL SETTINGS
    if isfield(optCfg, 'elRankThreshold') && ~isempty(optCfg.elRankThreshold)
        ELRANK_THRESHOLD = optCfg.elRankThreshold;
    else
        ELRANK_THRESHOLD = cfg.optimization.elRankThreshold;
    end

    if isfield(optCfg, 'familyMinOverlap') && ~isempty(optCfg.familyMinOverlap)
        FAMILY_MIN_OVERLAP = optCfg.familyMinOverlap;
    else
        FAMILY_MIN_OVERLAP = cfg.optimization.familyMinOverlap;
    end

    if isfield(optCfg, 'redundancyMinOverlap') && ~isempty(optCfg.redundancyMinOverlap)
        REDUNDANCY_MIN_OVERLAP = optCfg.redundancyMinOverlap;
    else
        REDUNDANCY_MIN_OVERLAP = cfg.optimization.redundancyMinOverlap;
    end

    optCfg.elRankThreshold = ELRANK_THRESHOLD;
    optCfg.familyMinOverlap = FAMILY_MIN_OVERLAP;
    optCfg.redundancyMinOverlap = REDUNDANCY_MIN_OVERLAP;
    Mobi_validate.optimization_input(groupCSVFiles, maxPeptides, WEIGHTS, optCfg);

    % VALIDATE INPUTS
    if isempty(groupCSVFiles)
        error('groupCSVFiles is empty.');
    end

    if ~iscell(groupCSVFiles) && ~isstring(groupCSVFiles)
        error('groupCSVFiles must be a cell array or string array of file paths.');
    end

    groupCSVFiles = cellstr(string(groupCSVFiles(:)));

    if isnan(maxPeptides) || maxPeptides <= 0 || mod(maxPeptides,1) ~= 0
        error('maxPeptides must be a positive integer.');
    end

    if isnan(ELRANK_THRESHOLD) || ELRANK_THRESHOLD <= 0
        error('EL rank threshold must be a positive numeric value.');
    end

    if isnan(FAMILY_MIN_OVERLAP) || FAMILY_MIN_OVERLAP <= 0 || mod(FAMILY_MIN_OVERLAP,1) ~= 0
        error('familyMinOverlap must be a positive integer.');
    end

    if isnan(REDUNDANCY_MIN_OVERLAP) || REDUNDANCY_MIN_OVERLAP <= 0 || mod(REDUNDANCY_MIN_OVERLAP,1) ~= 0
        error('redundancyMinOverlap must be a positive integer.');
    end

    requiredWeightFields = [ ...
        "wPatientCoverage", ...
        "wHLACoverage", ...
        "wBinding", ...
        "wPrevalence", ...
        "wFamilyNovelty", ...
        "wRedundancy"];

    for k = 1:numel(requiredWeightFields)
        if ~isfield(WEIGHTS, requiredWeightFields(k))
            error('WEIGHTS is missing required field: %s', requiredWeightFields(k));
        end
    end

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

    % READ AND COMBINE ALL GROUP CSV FILES
    combinedTable = table();

    for k = 1:length(groupCSVFiles)

        thisFile = groupCSVFiles{k};

        T = readtable(thisFile, ...
            'TextType', 'string', ...
            'VariableNamingRule', 'preserve');

        requiredCols = {'Peptide','MHC','%Rank_EL'};
        missingCols = requiredCols(~ismember(requiredCols, T.Properties.VariableNames));

        if ~isempty(missingCols)
            error('CSV %s is missing required columns: %s', ...
                thisFile, strjoin(missingCols, ', '));
        end

        % Force %Rank_EL numeric if needed
        if ~isnumeric(T.("%Rank_EL"))
            T.("%Rank_EL") = str2double(T.("%Rank_EL"));
        end

        % Keep only columns used by optimizer
        T = T(:, {'Peptide','MHC','%Rank_EL'});

        % Normalize HLA strings
        T.MHC = Mobi_utils.normalize_hla_list(T.MHC);

        % Store source file information
        T.SourceFile = repmat(string(thisFile), height(T), 1);

        combinedTable = [combinedTable; T]; %#ok<AGROW>
    end

    % FILTER TO CANDIDATE BINDERS
    combinedTable = combinedTable(combinedTable.("%Rank_EL") < ELRANK_THRESHOLD, :);

    if isempty(combinedTable)
        warning('No peptide-HLA rows passed the EL rank threshold in this group.');

        result = struct();
        result.combinedTable       = combinedTable;
        result.peptideFeatureTable = table();
        result.selectedPeptides    = table();
        result.selectionHistory    = table();
        result.weights             = WEIGHTS;
    result.maxPeptides         = maxPeptides;
    result.elrankThreshold     = ELRANK_THRESHOLD;
    result.familyMinOverlap    = FAMILY_MIN_OVERLAP;
    result.redundancyMinOverlap = REDUNDANCY_MIN_OVERLAP;
        result.coveredPatients     = strings(0,1);
        result.coveredHLAs         = strings(0,1);
        result.peptideFamily       = [];
        result.familyList          = [];
        Mobi_validate.optimization_output(result);
        return;
    end

    % BASIC SETS
    peptideList = unique(string(combinedTable.Peptide), 'stable');
    patientList = unique(string(combinedTable.SourceFile), 'stable');
    hlaList     = unique(string(combinedTable.MHC), 'stable');

    nPatients = length(patientList);
    nHLAs     = length(hlaList);
    nPeptides = length(peptideList);

    % PREALLOCATE PEPTIDE-LEVEL FEATURES
    StaticPatientCoverage = zeros(nPeptides,1);
    StaticHLACoverage     = zeros(nPeptides,1);
    BindingScore          = zeros(nPeptides,1);
    PrevalenceScore       = zeros(nPeptides,1);

    peptidePatients = cell(nPeptides,1);
    peptideHLAs     = cell(nPeptides,1);

    % BUILD PEPTIDE FAMILIES
    peptideFamily = build_peptide_families(peptideList, FAMILY_MIN_OVERLAP);
    familyList = unique(peptideFamily, 'stable');

    % PEPTIDE-LEVEL FEATURE CONSTRUCTION
    for i = 1:nPeptides

        pep = peptideList(i);

        rows = combinedTable(string(combinedTable.Peptide) == pep, :);

        filesWithPep = unique(string(rows.SourceFile), 'stable');
        hlasWithPep  = unique(string(rows.MHC), 'stable');

        peptidePatients{i} = filesWithPep;
        peptideHLAs{i}     = hlasWithPep;

        % Fraction of patients/files covered by this peptide
        StaticPatientCoverage(i) = length(filesWithPep) / max(nPatients, 1);

        % Fraction of HLAs covered by this peptide
        StaticHLACoverage(i) = length(hlasWithPep) / max(nHLAs, 1);

        % Binding score:
        % best-per-HLA transformed binding quality
        uniquePepHLAs = unique(string(rows.MHC), 'stable');
        bestPerHLA = nan(length(uniquePepHLAs),1);

        for h = 1:length(uniquePepHLAs)
            hlaRows = rows(string(rows.MHC) == uniquePepHLAs(h), :);
            ranks = hlaRows.("%Rank_EL");
            bindVals = max(0, 1 - ranks / ELRANK_THRESHOLD);
            bestPerHLA(h) = max(bindVals, [], 'omitnan');
        end

        BindingScore(i) = mean(bestPerHLA, 'omitnan');

        % Prevalence score:
        % how many distinct files contain this peptide
        PrevalenceScore(i) = length(filesWithPep) / max(nPatients, 1);
    end

    % Normalize static features to [0,1]
    StaticPatientCoverage = localNormalize01(StaticPatientCoverage);
    StaticHLACoverage     = localNormalize01(StaticHLACoverage);
    BindingScore          = localNormalize01(BindingScore);
    PrevalenceScore       = localNormalize01(PrevalenceScore);

    % ASSEMBLE PEPTIDE FEATURE TABLE

    peptideFeatureTable = table();
    peptideFeatureTable.Peptide               = peptideList;
    peptideFeatureTable.FamilyID              = peptideFamily;
    peptideFeatureTable.StaticPatientCoverage = StaticPatientCoverage;
    peptideFeatureTable.StaticHLACoverage     = StaticHLACoverage;
    peptideFeatureTable.BindingScore          = BindingScore;
    peptideFeatureTable.PrevalenceScore       = PrevalenceScore;

    % GREEDY SELECTION INITIALIZATION
    selectedMask     = false(nPeptides,1);
    selectedIdx      = [];
    coveredPatients  = strings(0,1);
    coveredHLAs      = strings(0,1);
    selectedFamilies = [];

    selectionHistory = table();

    % GREEDY PEPTIDE SELECTION
    for step = 1:min(maxPeptides, nPeptides)

        bestIdx = NaN;
        bestScore = -inf;
        bestBreakdown = [];

        for i = 1:nPeptides

            if selectedMask(i)
                continue;
            end

            % Incremental coverage terms
            newPatients = setdiff(peptidePatients{i}, coveredPatients);
            newHLAs     = setdiff(peptideHLAs{i}, coveredHLAs);

            IncrementalPatientCoverage = length(newPatients) / max(nPatients, 1);
            IncrementalHLACoverage     = length(newHLAs) / max(nHLAs, 1);

            % Family novelty
            thisFamily = peptideFamily(i);

            if ismember(thisFamily, selectedFamilies)
                FamilyNoveltyScore = 0;
            else
                FamilyNoveltyScore = 1;
            end

            % Redundancy penalty
            if isempty(selectedIdx)
                RedundancyPenalty = 0;
            else
                RedundancyPenalty = 0;
                for s = 1:length(selectedIdx)
                    selPep = peptideList(selectedIdx(s));
                    if Mobi_utils.has_min_contiguous_overlap(peptideList(i), selPep, REDUNDANCY_MIN_OVERLAP)
                        RedundancyPenalty = 1;
                        break;
                    end
                end
            end

            % Weighted score terms
            PatientTerm     = WEIGHTS.wPatientCoverage * IncrementalPatientCoverage;
            HLATerm         = WEIGHTS.wHLACoverage     * IncrementalHLACoverage;
            BindingTerm     = WEIGHTS.wBinding         * peptideFeatureTable.BindingScore(i);
            PrevalenceTerm  = WEIGHTS.wPrevalence      * peptideFeatureTable.PrevalenceScore(i);
            FamilyNovelTerm = WEIGHTS.wFamilyNovelty   * FamilyNoveltyScore;
            RedundancyTerm  = WEIGHTS.wRedundancy      * RedundancyPenalty;

            TotalScore = ...
                PatientTerm + ...
                HLATerm + ...
                BindingTerm + ...
                PrevalenceTerm + ...
                FamilyNovelTerm - ...
                RedundancyTerm;

            if TotalScore > bestScore
                bestScore = TotalScore;
                bestIdx = i;

                bestBreakdown = [ ...
                    IncrementalPatientCoverage, ...
                    IncrementalHLACoverage, ...
                    FamilyNoveltyScore, ...
                    RedundancyPenalty, ...
                    PatientTerm, ...
                    HLATerm, ...
                    BindingTerm, ...
                    PrevalenceTerm, ...
                    FamilyNovelTerm, ...
                    RedundancyTerm, ...
                    TotalScore];
            end
        end

        if isnan(bestIdx)
            break;
        end

        % Commit chosen peptide
        selectedMask(bestIdx) = true;
        selectedIdx(end+1,1) = bestIdx; %#ok<AGROW>

        coveredPatients = union(coveredPatients, peptidePatients{bestIdx}, 'stable');
        coveredHLAs     = union(coveredHLAs, peptideHLAs{bestIdx}, 'stable');

        if ~ismember(peptideFamily(bestIdx), selectedFamilies)
            selectedFamilies(end+1,1) = peptideFamily(bestIdx); %#ok<AGROW>
        end

        newRow = table( ...
            step, ...
            string(peptideFeatureTable.Peptide(bestIdx)), ...
            peptideFamily(bestIdx), ...
            bestBreakdown(1), ...
            bestBreakdown(2), ...
            bestBreakdown(3), ...
            bestBreakdown(4), ...
            bestBreakdown(5), ...
            bestBreakdown(6), ...
            bestBreakdown(7), ...
            bestBreakdown(8), ...
            bestBreakdown(9), ...
            bestBreakdown(10), ...
            bestBreakdown(11), ...
            length(coveredPatients) / max(nPatients, 1), ...
            length(coveredHLAs) / max(nHLAs, 1), ...
            length(selectedFamilies), ...
            'VariableNames', { ...
                'Step', ...
                'Peptide', ...
                'FamilyID', ...
                'IncrementalPatientCoverage', ...
                'IncrementalHLACoverage', ...
                'FamilyNoveltyScore', ...
                'RedundancyPenalty', ...
                'PatientTerm', ...
                'HLATerm', ...
                'BindingTerm', ...
                'PrevalenceTerm', ...
                'FamilyNoveltyTerm', ...
                'RedundancyTerm', ...
                'TotalScore', ...
                'CumulativePatientCoverage', ...
                'CumulativeHLACoverage', ...
                'NumFamiliesRepresented'});

        selectionHistory = [selectionHistory; newRow]; %#ok<AGROW>
    end

    % OUTPUTS

    selectedPeptides = peptideFeatureTable(selectedIdx,:);

    result = struct();
    result.combinedTable       = combinedTable;
    result.peptideFeatureTable = peptideFeatureTable;
    result.selectedPeptides    = selectedPeptides;
    result.selectionHistory    = selectionHistory;
    result.weights             = WEIGHTS;
    result.maxPeptides         = maxPeptides;
    result.elrankThreshold     = ELRANK_THRESHOLD;
    result.familyMinOverlap    = FAMILY_MIN_OVERLAP;
    result.redundancyMinOverlap = REDUNDANCY_MIN_OVERLAP;
    result.coveredPatients     = coveredPatients;
    result.coveredHLAs         = coveredHLAs;
    result.peptideFamily       = peptideFamily;
    result.familyList          = familyList;
    Mobi_validate.optimization_output(result);
end


% HELPER: normalize to [0,1]
function x = localNormalize01(x)

    x = double(x);

    if isempty(x) || max(x) == min(x)
        x = zeros(size(x));
    else
        x = (x - min(x)) / (max(x) - min(x));
    end
end


% HELPER: peptide family assignment
function familyID = build_peptide_families(peptideList, minOverlap)

    peptideList = string(peptideList);
    m = numel(peptideList);

    A = false(m,m);
    A(1:m+1:end) = true;

    for i = 1:m-1
        for j = i+1:m
            if Mobi_utils.has_min_contiguous_overlap(peptideList(i), peptideList(j), minOverlap)
                A(i,j) = true;
                A(j,i) = true;
            end
        end
    end

    familyID = Mobi_utils.connected_components_from_adjacency(A);
end
