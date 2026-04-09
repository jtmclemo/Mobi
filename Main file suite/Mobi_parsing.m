function [X, feature_names, ids, summaryTbl] = Mobi_parsing(filename, id_col, typedHLAs)

    cfg = Mobi_config();

    % DEFAULT INPUTS
    if nargin < 2 || isempty(id_col)
        id_col = cfg.parsing.idColumn;
    end

    if nargin < 3 || isempty(typedHLAs)
        typedHLAs = strings(0,1);
    end

    id_col = string(id_col);
    typedHLAs = string(typedHLAs);

    % CONFIG-DRIVEN SETTINGS
    min_overlap_family = cfg.parsing.minOverlapFamily;

    peptide_col = cfg.parsing.peptideColumn;
    hla_col     = cfg.parsing.hlaColumn;
    rank_col    = cfg.parsing.rankColumn;
    binder_col  = cfg.parsing.binderColumn; %#ok<NASGU>
    elscore_col = cfg.parsing.scoreColumn;  %#ok<NASGU>

    % READ INPUT FILE
    T = readtable(filename, ...
        'TextType', 'string', ...
        'VariableNamingRule', 'preserve');

    vars = string(T.Properties.VariableNames);

    required_cols = [id_col, peptide_col, hla_col, rank_col];
    for k = 1:numel(required_cols)
        if ~any(vars == required_cols(k))
            error('Missing required column: %s', required_cols(k));
        end
    end

    % Force %Rank_EL numeric
    if ~isnumeric(T.(rank_col))
        T.(rank_col) = str2double(T.(rank_col));
    end

    % Normalize HLA strings
    T.(hla_col) = Mobi_utils.normalize_hla_list(T.(hla_col));

    % OPTIONAL HLA FILTERING
    if ~isempty(typedHLAs)
        typedHLAs = Mobi_utils.normalize_hla_list(typedHLAs);
        keepIdx = ismember(T.(hla_col), typedHLAs);
        T = T(keepIdx, :);
    end

    if isempty(T)
        error('No NetMHCpan rows remained after HLA filtering.');
    end

    % UNIQUE SAMPLE IDENTITIES
    ids = unique(T.(id_col), 'stable');
    n = numel(ids);

    % PREALLOCATE SUMMARY FEATURES
    num_rows                = zeros(n,1);
    num_unique_peptides     = zeros(n,1);
    num_unique_HLAs_hit     = zeros(n,1);

    best_rank               = nan(n,1);
    mean_rank               = nan(n,1);
    median_rank             = nan(n,1);
    mean_top3_rank          = nan(n,1);
    mean_top5_rank          = nan(n,1);

    peptide_promiscuity_max = zeros(n,1);
    mean_peptide_length     = nan(n,1);

    rank_std                         = nan(n,1);
    peptide_promiscuity_mean         = nan(n,1);
    mean_best_rank_per_peptide       = nan(n,1);
    hla_entropy                      = nan(n,1);
    rows_per_unique_peptide          = nan(n,1);

    num_peptide_families              = zeros(n,1);
    largest_peptide_family_size       = zeros(n,1);
    frac_peptides_in_family_size_gt_1 = nan(n,1);
    mean_best_rank_per_family         = nan(n,1);

    % MAIN LOOP OVER IDENTITIES
    for i = 1:n

        thisID = ids(i);
        idx = T.(id_col) == thisID;
        S = T(idx,:);

        % BASIC COUNTS
        num_rows(i)            = height(S);
        num_unique_peptides(i) = numel(unique(S.(peptide_col)));
        num_unique_HLAs_hit(i) = numel(unique(S.(hla_col)));

        rows_per_unique_peptide(i) = ...
            num_rows(i) / max(num_unique_peptides(i), 1);

        % RANK-LEVEL SUMMARIES
        ranks = S.(rank_col);
        ranks = ranks(~isnan(ranks));

        if ~isempty(ranks)

            best_rank(i)   = min(ranks);
            mean_rank(i)   = mean(ranks);
            median_rank(i) = median(ranks);

            ranks_sorted = sort(ranks, 'ascend');
            mean_top3_rank(i) = mean(ranks_sorted(1:min(3,end)));
            mean_top5_rank(i) = mean(ranks_sorted(1:min(5,end)));

            rank_std(i) = std(ranks);
        end

        % PEPTIDE-LEVEL STRUCTURE
        peptides = unique(S.(peptide_col), 'stable');

        perPepCounts    = zeros(numel(peptides),1);
        pepLens         = zeros(numel(peptides),1);
        bestRankPerPep  = nan(numel(peptides),1);

        for p = 1:numel(peptides)

            pidx = S.(peptide_col) == peptides(p);

            perPepCounts(p) = numel(unique(S.(hla_col)(pidx)));
            pepLens(p)      = strlength(peptides(p));

            thisPepRanks = S.(rank_col)(pidx);
            thisPepRanks = thisPepRanks(~isnan(thisPepRanks));

            if ~isempty(thisPepRanks)
                bestRankPerPep(p) = min(thisPepRanks);
            end
        end

        if ~isempty(perPepCounts)
            peptide_promiscuity_max(i)  = max(perPepCounts);
            peptide_promiscuity_mean(i) = mean(perPepCounts);
            mean_peptide_length(i)      = mean(pepLens);
        end

        if any(~isnan(bestRankPerPep))
            mean_best_rank_per_peptide(i) = ...
                mean(bestRankPerPep(~isnan(bestRankPerPep)));
        end

        % HLA ENTROPY / EVENNESS
        uniqueHLAs = unique(S.(hla_col), 'stable');
        hlaCounts  = zeros(numel(uniqueHLAs),1);

        for h = 1:numel(uniqueHLAs)
            hlaCounts(h) = sum(S.(hla_col) == uniqueHLAs(h));
        end

        if sum(hlaCounts) > 0
            p_hla = hlaCounts / sum(hlaCounts);
            hla_entropy(i) = -sum(p_hla .* log(p_hla + eps));
        end

        % PEPTIDE FAMILY GRAPH
        m = numel(peptides);

        if m > 0

            comp = build_peptide_family_components(peptides, min_overlap_family);
            families = unique(comp, 'stable');
            num_peptide_families(i) = numel(families);

            famSizes = zeros(numel(families),1);
            familyBestRanks = nan(numel(families),1);

            for f = 1:numel(families)

                famIdx = (comp == families(f));
                famSizes(f) = sum(famIdx);

                theseBest = bestRankPerPep(famIdx);
                theseBest = theseBest(~isnan(theseBest));

                if ~isempty(theseBest)
                    familyBestRanks(f) = min(theseBest);
                end
            end

            if ~isempty(famSizes)
                largest_peptide_family_size(i) = max(famSizes);

                frac_peptides_in_family_size_gt_1(i) = ...
                    sum(famSizes(famSizes > 1)) / max(m,1);
            end

            if any(~isnan(familyBestRanks))
                mean_best_rank_per_family(i) = ...
                    mean(familyBestRanks(~isnan(familyBestRanks)));
            end
        end
    end

    % BUILD SUMMARY TABLE
    summaryTbl = table(ids, ...
        num_rows, ...
        num_unique_peptides, ...
        num_unique_HLAs_hit, ...
        best_rank, ...
        mean_rank, ...
        median_rank, ...
        mean_top3_rank, ...
        mean_top5_rank, ...
        peptide_promiscuity_max, ...
        mean_peptide_length, ...
        rank_std, ...
        peptide_promiscuity_mean, ...
        mean_best_rank_per_peptide, ...
        hla_entropy, ...
        rows_per_unique_peptide, ...
        num_peptide_families, ...
        largest_peptide_family_size, ...
        frac_peptides_in_family_size_gt_1, ...
        mean_best_rank_per_family, ...
        'VariableNames', {char(id_col), ...
        'num_rows', ...
        'num_unique_peptides', ...
        'num_unique_HLAs_hit', ...
        'best_rank', ...
        'mean_rank', ...
        'median_rank', ...
        'mean_top3_rank', ...
        'mean_top5_rank', ...
        'peptide_promiscuity_max', ...
        'mean_peptide_length', ...
        'rank_std', ...
        'peptide_promiscuity_mean', ...
        'mean_best_rank_per_peptide', ...
        'hla_entropy', ...
        'rows_per_unique_peptide', ...
        'num_peptide_families', ...
        'largest_peptide_family_size', ...
        'frac_peptides_in_family_size_gt_1', ...
        'mean_best_rank_per_family'});

    % BUILD FEATURE MATRIX
    X = [ ...
        summaryTbl.num_unique_peptides, ...
        summaryTbl.num_unique_HLAs_hit, ...
        summaryTbl.mean_rank, ...
        summaryTbl.best_rank, ...
        summaryTbl.rank_std, ...
        summaryTbl.peptide_promiscuity_mean, ...
        summaryTbl.mean_best_rank_per_peptide, ...
        summaryTbl.hla_entropy, ...
        summaryTbl.rows_per_unique_peptide, ...
        summaryTbl.num_peptide_families, ...
        summaryTbl.largest_peptide_family_size, ...
        summaryTbl.frac_peptides_in_family_size_gt_1, ...
        summaryTbl.mean_best_rank_per_family ...
    ];

    feature_names = [ ...
        "num_unique_peptides", ...
        "num_unique_HLAs_hit", ...
        "mean_rank", ...
        "best_rank", ...
        "rank_std", ...
        "peptide_promiscuity_mean", ...
        "mean_best_rank_per_peptide", ...
        "hla_entropy", ...
        "rows_per_unique_peptide", ...
        "num_peptide_families", ...
        "largest_peptide_family_size", ...
        "frac_peptides_in_family_size_gt_1", ...
        "mean_best_rank_per_family" ...
    ];
    Mobi_validate.parser_output(X, feature_names, ids, summaryTbl, filename);
end


function comp = build_peptide_family_components(peptides, minOverlap)

    m = numel(peptides);

    if m == 0
        comp = zeros(0,1);
        return;
    end

    parents = 1:m;
    windowToIndices = containers.Map('KeyType', 'char', 'ValueType', 'any');

    for idx = 1:m
        peptide = char(string(peptides(idx)));

        if length(peptide) < minOverlap
            continue;
        end

        nWindows = length(peptide) - minOverlap + 1;
        windows = strings(nWindows, 1);

        for startIdx = 1:nWindows
            windows(startIdx) = string(peptide(startIdx:startIdx + minOverlap - 1));
        end

        windows = unique(windows, 'stable');

        for w = 1:numel(windows)
            key = char(windows(w));
            if isKey(windowToIndices, key)
                windowToIndices(key) = [windowToIndices(key), idx];
            else
                windowToIndices(key) = idx;
            end
        end
    end

    overlapKeys = windowToIndices.keys;
    for keyIdx = 1:numel(overlapKeys)
        indices = unique(windowToIndices(overlapKeys{keyIdx}), 'stable');
        if numel(indices) < 2
            continue;
        end

        anchor = indices(1);
        for idxPos = 2:numel(indices)
            parents = union_component_roots(parents, anchor, indices(idxPos));
        end
    end

    comp = zeros(m,1);
    rootLabels = containers.Map('KeyType', 'double', 'ValueType', 'double');
    nextLabel = 0;

    for idx = 1:m
        root = find_component_root(parents, idx);
        parents(idx) = root;

        if ~isKey(rootLabels, root)
            nextLabel = nextLabel + 1;
            rootLabels(root) = nextLabel;
        end

        comp(idx) = rootLabels(root);
    end
end


function parents = union_component_roots(parents, leftIdx, rightIdx)

    leftRoot = find_component_root(parents, leftIdx);
    rightRoot = find_component_root(parents, rightIdx);

    if leftRoot ~= rightRoot
        parents(rightRoot) = leftRoot;
    end
end


function root = find_component_root(parents, idx)

    root = idx;
    while parents(root) ~= root
        root = parents(root);
    end

    while parents(idx) ~= idx
        nextIdx = parents(idx);
        parents(idx) = root;
        idx = nextIdx;
    end
end
