function tests = test_Mobi_parsing
    tests = functiontests(localfunctions);
end


function setupOnce(testCase)
    addpath(fileparts(fileparts(mfilename('fullpath'))));
end


function testParserFeatureRegression(testCase)
    csvPath = fullfile(tempdir, ['mobi_parser_test_' char(java.util.UUID.randomUUID) '.csv']);
    fid = fopen(csvPath, 'w');
    cleanupObj = onCleanup(@() delete_if_exists(csvPath)); %#ok<NASGU>

    fprintf(fid, 'Identity,Peptide,MHC,%%Rank_EL\n');
    fprintf(fid, 'S1,ABCDEFG,HLA-A*01:01,1\n');
    fprintf(fid, 'S1,XXABCDEFGY,HLA-A*02:01,3\n');
    fprintf(fid, 'S1,CCCCCCC,HLA-A*01:01,5\n');
    fclose(fid);

    [X, featureNames, ids, summaryTbl] = Mobi_parsing(csvPath, "Identity", strings(0,1));

    verifyEqual(testCase, size(X), [1 13]);
    verifyEqual(testCase, numel(featureNames), 13);
    verifyEqual(testCase, string(ids), "S1");
    verifyEqual(testCase, height(summaryTbl), 1);
    verifyEqual(testCase, summaryTbl.num_unique_peptides, 3);
    verifyEqual(testCase, summaryTbl.num_unique_HLAs_hit, 2);
    verifyEqual(testCase, summaryTbl.best_rank, 1);
    verifyEqual(testCase, summaryTbl.mean_rank, 3);
    verifyEqual(testCase, summaryTbl.rows_per_unique_peptide, 1);
    verifyEqual(testCase, summaryTbl.num_peptide_families, 2);
    verifyEqual(testCase, summaryTbl.largest_peptide_family_size, 2);
end


function delete_if_exists(path)
    if exist(path, 'file')
        delete(path);
    end
end
