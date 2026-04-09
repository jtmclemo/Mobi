function cfg = Mobi_config()

    % General / version
    cfg.general.programName = "Mobi";
    cfg.general.version     = "v1.0.0";
    cfg.general.dateLabel   = "2026-04-07";

    % Parsing defaults
    cfg.parsing.idColumn = "Identity";

    % Minimum contiguous overlap length used to declare peptide-family
    % similarity in parsing summaries.
    cfg.parsing.minOverlapFamily = 7;

    % Expected NetMHCpan-style columns
    cfg.parsing.peptideColumn = "Peptide";
    cfg.parsing.hlaColumn     = "MHC";
    cfg.parsing.rankColumn    = "%Rank_EL";
    cfg.parsing.binderColumn  = "BindLevel";
    cfg.parsing.scoreColumn   = "Score_EL";

    % TDA defaults
    cfg.tda.makePlots = true;
    cfg.tda.distanceMetric = "euclidean";

    % Default TDA feature weights for 13-feature input
    cfg.tda.defaultWeights13 = [ ...
        0.80, ... % 1  num_unique_peptides
        0.90, ... % 2  num_unique_HLAs_hit
        1.20, ... % 3  mean_rank
        1.10, ... % 4  best_rank
        0.90, ... % 5  rank_std
        0.90, ... % 6  peptide_promiscuity_mean
        1.10, ... % 7  mean_best_rank_per_peptide
        0.85, ... % 8  hla_entropy
        0.75, ... % 9  rows_per_unique_peptide
        1.00, ... % 10 num_peptide_families
        1.00, ... % 11 largest_peptide_family_size
        0.95, ... % 12 frac_peptides_in_family_size_gt_1
        1.15  ... % 13 mean_best_rank_per_family
    ];

    % Thresholding rule
    cfg.tda.thresholdRule = "percentile";
    cfg.tda.defaultPercentile = 75;

    % Optimization defaults
    cfg.optimization.maxPeptides = 20;

    % Permissive NetMHCpan inclusion threshold
    cfg.optimization.elRankThreshold = 32.7;

    % Peptide family / redundancy logic
    cfg.optimization.familyMinOverlap     = 7;
    cfg.optimization.redundancyMinOverlap = 7;

    % Default optimization weights
    cfg.optimization.weights.wPatientCoverage = 0.30;
    cfg.optimization.weights.wHLACoverage     = 0.20;
    cfg.optimization.weights.wBinding         = 0.20;
    cfg.optimization.weights.wPrevalence      = 0.10;
    cfg.optimization.weights.wFamilyNovelty   = 0.15;
    cfg.optimization.weights.wRedundancy      = 0.15;

    % Whether frontend should normalize user-entered optimization weights
    % to sum to 1
    cfg.optimization.normalizeWeights = true;

    % Reference panel defaults
    cfg.reference.defaultTopN = 10;
    cfg.reference.defaultLoci = ["A","B","C"];

    cfg.reference.regionOptions = [ ...
        "Australia";
        "Central Asia";
        "Europe";
        "North Africa";
        "North America";
        "North-East Asia";
        "Oceania";
        "South and Central America";
        "South Asia";
        "South-East Asia";
        "Sub-Saharan Africa";
        "Western Asia";
        "All Regions"];

    % Frontend / output defaults
    cfg.frontend.useDiaryLog = true;
    cfg.frontend.logPrefix   = "Mobi_DetailedView_";
    cfg.frontend.logExt      = ".txt";
    cfg.frontend.timestampFormat = "yyyymmdd_HHMMSS";
    cfg.frontend.outputPrefix = "Mobi_outputs_";
    cfg.frontend.enableStructuredCSVExport = false;

    % Exclude auto-generated CSVs when scanning the NetMHCpan folder
    cfg.frontend.excludePrefixes = ["bchs4397_", "group_", "Mobi_"];

    % Prompt defaults
    cfg.frontend.defaultContinueOptimization = true;


    % Future switches
    cfg.future.enableHLARobustness       = false;
    cfg.future.enableSubgroupWeighting   = false;
end
