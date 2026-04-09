function tests = test_Mobi_robustness
    tests = functiontests(localfunctions);
end


function setupOnce(testCase)
    addpath(fileparts(fileparts(mfilename('fullpath'))));
end


function testRobustnessSweepSizes(testCase)
    cfg = Mobi_config();
    baseInput = struct();
    baseInput.X = [0; 1];
    baseInput.sampleNames = ["A"; "B"];
    baseInput.sampleFiles = ["A.csv"; "B.csv"];

    runFns = struct();
    runFns.runTDA = @(baseInput, tdaCfg) fake_tda(baseInput, tdaCfg);
    runFns.runOptimization = @(baseInput, optCfg, tdaResult) fake_optimization(baseInput, optCfg, tdaResult);

    fast = Mobi_robustness(runFns, baseInput, cfg, 1);
    standard = Mobi_robustness(runFns, baseInput, cfg, 2);
    full = Mobi_robustness(runFns, baseInput, cfg, 3);

    verifyEqual(testCase, height(fast.summaryTable), 3);
    verifyEqual(testCase, fast.overall.totalRuns, 3);
    verifyEqual(testCase, string(fast.mode.label), "fast");
    verifyEqual(testCase, string(fast.signaturePeptideMode), "selectedPeptideUnionAcrossModes");

    verifyEqual(testCase, height(standard.summaryTable), 10);
    verifyEqual(testCase, standard.overall.totalRuns, 10);
    verifyEqual(testCase, string(standard.mode.label), "standard");

    verifyEqual(testCase, height(full.summaryTable), 84);
    verifyEqual(testCase, full.overall.totalRuns, 84);
    verifyEqual(testCase, string(full.mode.label), "full");
end


function result = fake_tda(baseInput, tdaCfg)
    n = numel(baseInput.sampleNames);
    result = struct();
    result.clusters = ones(n,1);
    result.rawClusters = ones(n,1);
    result.groupLabels = ones(n,1);
    result.optimizationGroups = ones(n,1);
    result.optimizationGroupSizes = n;
    result.thresholdRule = string(tdaCfg.thresholdRule);
    result.chosenPercentile = tdaCfg.percentile;
end


function result = fake_optimization(~, optCfg, ~)
    result = struct();
    result.selectedPeptideUnionAcrossModes = "PEP" + string(optCfg.maxPeptides);
    result.practicalSelectedPeptideUnion = result.selectedPeptideUnionAcrossModes;
    result.globalSelectedPeptides = result.selectedPeptideUnionAcrossModes;
    result.honestSelectedPeptideUnion = strings(0,1);
    result.signaturePeptideMode = "selectedPeptideUnionAcrossModes";
end
