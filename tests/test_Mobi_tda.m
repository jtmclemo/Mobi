function tests = test_Mobi_tda
    tests = functiontests(localfunctions);
end


function setupOnce(testCase)
    addpath(fileparts(fileparts(mfilename('fullpath'))));
end


function testTDADeterministicGeometry(testCase)
    X = [0; 1; 10];
    ids = ["A"; "B"; "C"];

    tdaCfg = struct();
    tdaCfg.makePlots = false;
    tdaCfg.distanceMetric = "euclidean";
    tdaCfg.thresholdRule = "percentile";
    tdaCfg.percentile = 50;
    tdaCfg.mergeSingletons = false;
    tdaCfg.maxOptimizationGroupSize = Inf;
    tdaCfg.minOptimizationGroupSize = 1;

    result = Mobi_tda(X, ids, tdaCfg);

    expectedXz = zscore(X);
    expectedD = pdist2(expectedXz, expectedXz, 'euclidean');
    expectedDeaths = sort([expectedD(1,2); expectedD(2,3)]);
    expectedThreshold = prctile(expectedDeaths, 50);

    verifyEqual(testCase, result.D, expectedD, 'AbsTol', 1e-12);
    verifyEqual(testCase, result.finiteDeaths, expectedDeaths, 'AbsTol', 1e-12);
    verifyEqual(testCase, result.threshold, expectedThreshold, 'AbsTol', 1e-12);
    verifyEqual(testCase, numel(result.rawClusters), 3);
    verifyEqual(testCase, numel(result.optimizationGroups), 3);
    verifyTrue(testCase, all(result.rawClusters >= 1));
    verifyTrue(testCase, isfield(result, 'practicalGroupingAudit'));
    verifyTrue(testCase, isfield(result.practicalGroupingAudit, 'events'));
    verifyTrue(testCase, isfield(result.practicalGroupingAudit, 'rawToPracticalMapping'));
    verifyTrue(testCase, isfield(result.practicalGroupingAudit, 'summary'));
    verifyGreaterThanOrEqual(testCase, numel(result.practicalGroupingAudit.events), 2);
    verifyEqual(testCase, result.practicalGroupingAudit.events(1).Stage, "initial");
    verifyEqual(testCase, result.practicalGroupingAudit.events(end).Stage, "final");
    verifyEqual(testCase, height(result.practicalGroupingAudit.rawToPracticalMapping), numel(ids));
end
