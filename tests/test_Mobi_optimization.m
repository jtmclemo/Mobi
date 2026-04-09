function tests = test_Mobi_optimization
    tests = functiontests(localfunctions);
end


function setupOnce(testCase)
    addpath(fileparts(fileparts(mfilename('fullpath'))));
end


function testOptimizationDeterministicSelection(testCase)
    file1 = fullfile(tempdir, ['mobi_opt_test_1_' char(java.util.UUID.randomUUID) '.csv']);
    file2 = fullfile(tempdir, ['mobi_opt_test_2_' char(java.util.UUID.randomUUID) '.csv']);
    cleanupObj = onCleanup(@() cleanup_files(file1, file2)); %#ok<NASGU>

    write_opt_csv(file1, [
        "PLOW,HLA-A*01:01,1"
        "PHIGH,HLA-A*01:01,10"]);
    write_opt_csv(file2, [
        "PLOW,HLA-A*01:01,2"
        "POTHER,HLA-B*07:02,5"]);

    weights = struct();
    weights.wPatientCoverage = 0;
    weights.wHLACoverage = 0;
    weights.wBinding = 1;
    weights.wPrevalence = 0;
    weights.wFamilyNovelty = 0;
    weights.wRedundancy = 0;

    optCfg = struct();
    optCfg.elRankThreshold = 32.7;
    optCfg.familyMinOverlap = 7;
    optCfg.redundancyMinOverlap = 7;

    result = Mobi_optimization(string([file1; file2]), 1, weights, optCfg);

    verifyEqual(testCase, height(result.selectedPeptides), 1);
    verifyEqual(testCase, string(result.selectedPeptides.Peptide(1)), "PLOW");
    verifyLessThanOrEqual(testCase, height(result.selectedPeptides), result.maxPeptides);
    verifyEqual(testCase, result.selectionHistory.Step, 1);
    verifyGreaterThanOrEqual(testCase, result.selectionHistory.CumulativePatientCoverage, 0);
    verifyLessThanOrEqual(testCase, result.selectionHistory.CumulativePatientCoverage, 1);
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
