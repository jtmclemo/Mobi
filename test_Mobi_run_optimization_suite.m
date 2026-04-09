function tests = test_Mobi_run_optimization_suite
    tests = functiontests(localfunctions);
end


function setupOnce(testCase)
    addpath(fileparts(fileparts(mfilename('fullpath'))));
end


function testOptimizationSuiteConsistency(testCase)
    [files, cleanupObj] = make_suite_csvs(); %#ok<ASGLU>

    tdaResult = struct();
    tdaResult.rawClusters = [1; 2];
    tdaResult.groupLabels = [1; 2];
    tdaResult.clusters = [1; 2];
    tdaResult.optimizationGroups = [1; 1];

    weights = default_test_weights();
    optCfg = struct();
    optCfg.maxPeptides = 1;
    optCfg.elRankThreshold = 32.7;
    optCfg.familyMinOverlap = 7;
    optCfg.redundancyMinOverlap = 7;
    optCfg.weights = weights;
    optCfg.minGroupSizeForPracticalOptimization = 1;

    baseInput = struct();
    baseInput.sampleFiles = string(files(:));
    baseInput.sampleNames = ["S1"; "S2"];

    suite = Mobi_run_optimization_suite(baseInput, optCfg, tdaResult, "all");

    verifyTrue(testCase, suite.didRun);
    verifyEqual(testCase, numel(suite.honestGroups), 2);
    verifyEqual(testCase, numel(suite.practicalGroups), 1);
    verifyEqual(testCase, numel(suite.globalSelectedPeptides), 1);
    verifyFalse(testCase, isempty(suite.honestSelectedPeptidesByGroup));
    verifyFalse(testCase, isempty(suite.practicalSelectedPeptidesByGroup));

    expectedUnion = unique([ ...
        suite.globalSelectedPeptides(:); ...
        suite.honestSelectedPeptideUnion(:); ...
        suite.practicalSelectedPeptideUnion(:)], 'stable');

    verifyEqual(testCase, suite.selectedPeptideUnionAcrossModes, expectedUnion);
end


function [files, cleanupObj] = make_suite_csvs()
    file1 = fullfile(tempdir, ['mobi_suite_test_1_' char(java.util.UUID.randomUUID) '.csv']);
    file2 = fullfile(tempdir, ['mobi_suite_test_2_' char(java.util.UUID.randomUUID) '.csv']);
    write_opt_csv(file1, ["PLOW,HLA-A*01:01,1"; "PHIGH,HLA-A*01:01,10"]);
    write_opt_csv(file2, ["POTHER,HLA-B*07:02,5"; "PLOW,HLA-A*01:01,2"]);
    files = {file1; file2};
    cleanupObj = onCleanup(@() cleanup_files(file1, file2));
end


function weights = default_test_weights()
    weights = struct();
    weights.wPatientCoverage = 0;
    weights.wHLACoverage = 0;
    weights.wBinding = 1;
    weights.wPrevalence = 0;
    weights.wFamilyNovelty = 0;
    weights.wRedundancy = 0;
end


function write_opt_csv(path, rows)
    fid = fopen(path, 'w');
    fprintf(fid, 'Peptide,MHC,%%Rank_EL\n');
    for i = 1:numel(rows)
        fprintf(fid, '%s\n', rows(i));
    end
    fclose(fid);
end


function cleanup_files(varargin)
    for i = 1:nargin
        if exist(varargin{i}, 'file')
            delete(varargin{i});
        end
    end
end
