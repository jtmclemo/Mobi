classdef Mobi_validate

    methods (Static)

        function parser_output(X, featureNames, ids, summaryTbl, sourceLabel)
            if nargin < 5 || isempty(sourceLabel)
                sourceLabel = "parser output";
            end
            sourceLabel = string(sourceLabel);

            if ~isnumeric(X)
                error('%s: X must be numeric.', sourceLabel);
            end

            ids = string(ids(:));
            featureNames = string(featureNames(:));

            if isempty(ids)
                error('%s: ids must be nonempty.', sourceLabel);
            end

            if size(X,1) ~= numel(ids)
                error('%s: size(X,1) must match numel(ids).', sourceLabel);
            end

            if size(X,2) ~= numel(featureNames)
                error('%s: size(X,2) must match numel(feature_names).', sourceLabel);
            end

            if ~istable(summaryTbl)
                error('%s: summaryTbl must be a table.', sourceLabel);
            end

            if height(summaryTbl) ~= numel(ids)
                error('%s: summaryTbl height must match numel(ids).', sourceLabel);
            end

            if any(all(~isfinite(X), 1))
                error('%s: no feature column may be entirely non-finite.', sourceLabel);
            end

            hlaLikeVars = intersect(["MHC","HLA","HLAsUsed"], string(summaryTbl.Properties.VariableNames));
            for i = 1:numel(hlaLikeVars)
                vals = string(summaryTbl.(hlaLikeVars(i)));
                if any(ismissing(vals))
                    error('%s: HLA-like summary field %s contains missing values.', ...
                        sourceLabel, hlaLikeVars(i));
                end
            end
        end


        function tda_input(X, ids, tdaCfg)
            if nargin < 3 || isempty(tdaCfg)
                tdaCfg = struct();
            end

            if ~isnumeric(X)
                error('TDA input: X must be numeric.');
            end

            if size(X,1) < 2
                error('TDA input: X must contain at least 2 rows.');
            end

            if numel(string(ids(:))) ~= size(X,1)
                error('TDA input: numel(ids) must match size(X,1).');
            end

            cfg = Mobi_config();
            thresholdRule = string(cfg.tda.thresholdRule);
            if isfield(tdaCfg, 'thresholdRule') && ~isempty(tdaCfg.thresholdRule)
                thresholdRule = string(tdaCfg.thresholdRule);
            end

            if ~ismember(lower(thresholdRule), ["percentile", "largest_gap"])
                error('TDA input: unsupported threshold rule %s.', thresholdRule);
            end

            if lower(thresholdRule) == "percentile"
                if isfield(tdaCfg, 'percentile') && ~isempty(tdaCfg.percentile)
                    percentile = tdaCfg.percentile;
                else
                    percentile = cfg.tda.defaultPercentile;
                end

                if ~isnumeric(percentile) || ~isscalar(percentile) || ...
                        isnan(percentile) || percentile < 0 || percentile > 100
                    error('TDA input: percentile must be between 0 and 100.');
                end
            end

            if isfield(tdaCfg, 'weights13') && ~isempty(tdaCfg.weights13)
                weights = tdaCfg.weights13;
                if numel(weights) ~= size(X,2) && ~(size(X,2) == 13 && numel(cfg.tda.defaultWeights13) == 13)
                    error('TDA input: weights13 length must match feature count or allow the 13-feature fallback.');
                end
            end
        end


        function tda_output(result, X, ids)
            n = size(X,1);
            requiredFields = ["D","rawClusters","optimizationGroups","optimizationGroupSizes", ...
                "practicalGroupingAudit","threshold","thresholdRule"];
            Mobi_validate.require_fields(result, requiredFields, "TDA output");

            if ~ismatrix(result.D) || size(result.D,1) ~= n || size(result.D,2) ~= n
                error('TDA output: result.D must be square with one row/column per sample.');
            end

            if numel(string(ids(:))) ~= n
                error('TDA output: ids must match X rows.');
            end

            if numel(result.rawClusters) ~= n
                error('TDA output: rawClusters length must match X rows.');
            end

            if numel(result.optimizationGroups) ~= n
                error('TDA output: optimizationGroups length must match X rows.');
            end

            Mobi_validate.require_positive_integer_labels(result.rawClusters, "TDA output rawClusters");
            Mobi_validate.require_positive_integer_labels(result.optimizationGroups, "TDA output optimizationGroups");
            Mobi_validate.require_consecutive_labels(result.optimizationGroups, "TDA output optimizationGroups");

            if numel(result.optimizationGroupSizes) ~= max(result.optimizationGroups)
                error('TDA output: optimizationGroupSizes length must match max(optimizationGroups).');
            end

            if ~(isnumeric(result.threshold) && isscalar(result.threshold) && isfinite(result.threshold))
                error('TDA output: threshold must be a finite numeric scalar.');
            end

            if lower(string(result.thresholdRule)) == "percentile"
                if ~isfield(result, 'chosenPercentile') || ...
                        ~isnumeric(result.chosenPercentile) || ...
                        ~isscalar(result.chosenPercentile) || ...
                        isnan(result.chosenPercentile) || ...
                        result.chosenPercentile < 0 || result.chosenPercentile > 100
                    error('TDA output: chosenPercentile must be valid for percentile thresholding.');
                end
            end

            audit = result.practicalGroupingAudit;
            Mobi_validate.require_fields(audit, ["events","rawToPracticalMapping","summary"], ...
                "TDA practicalGroupingAudit");

            if numel(audit.events) < 2
                error('TDA output: practicalGroupingAudit must include at least initial and final events.');
            end

            if ~istable(audit.rawToPracticalMapping) || height(audit.rawToPracticalMapping) ~= n
                error('TDA output: rawToPracticalMapping must be a table with one row per input sample.');
            end
        end


        function optimization_input(groupCSVFiles, maxPeptides, weights, optCfg)
            if nargin < 4 || isempty(optCfg)
                optCfg = struct();
            end

            if isempty(groupCSVFiles)
                error('Optimization input: groupCSVFiles must be nonempty.');
            end

            if ~iscell(groupCSVFiles) && ~isstring(groupCSVFiles)
                error('Optimization input: groupCSVFiles must be a cell array or string array.');
            end

            if ~isnumeric(maxPeptides) || ~isscalar(maxPeptides) || ...
                    isnan(maxPeptides) || maxPeptides <= 0 || mod(maxPeptides,1) ~= 0
                error('Optimization input: maxPeptides must be a positive integer.');
            end

            requiredWeightFields = [ ...
                "wPatientCoverage", ...
                "wHLACoverage", ...
                "wBinding", ...
                "wPrevalence", ...
                "wFamilyNovelty", ...
                "wRedundancy"];
            Mobi_validate.require_fields(weights, requiredWeightFields, "Optimization weights");

            weightValues = zeros(1, numel(requiredWeightFields));
            for i = 1:numel(requiredWeightFields)
                weightValues(i) = weights.(requiredWeightFields(i));
            end

            if any(~isfinite(weightValues)) || any(weightValues < 0)
                error('Optimization input: all weights must be finite and nonnegative.');
            end

            cfg = Mobi_config();
            elRankThreshold = Mobi_validate.get_field_or_default(optCfg, 'elRankThreshold', cfg.optimization.elRankThreshold);
            familyMinOverlap = Mobi_validate.get_field_or_default(optCfg, 'familyMinOverlap', cfg.optimization.familyMinOverlap);
            redundancyMinOverlap = Mobi_validate.get_field_or_default(optCfg, 'redundancyMinOverlap', cfg.optimization.redundancyMinOverlap);

            if ~isnumeric(elRankThreshold) || ~isscalar(elRankThreshold) || ...
                    isnan(elRankThreshold) || elRankThreshold <= 0
                error('Optimization input: elRankThreshold must be a positive numeric scalar.');
            end

            if ~Mobi_validate.is_positive_integer_scalar(familyMinOverlap)
                error('Optimization input: familyMinOverlap must be a positive integer.');
            end

            if ~Mobi_validate.is_positive_integer_scalar(redundancyMinOverlap)
                error('Optimization input: redundancyMinOverlap must be a positive integer.');
            end
        end


        function optimization_output(result)
            requiredFields = [ ...
                "peptideFeatureTable", ...
                "selectedPeptides", ...
                "selectionHistory", ...
                "maxPeptides", ...
                "coveredPatients", ...
                "coveredHLAs", ...
                "peptideFamily", ...
                "familyList"];
            Mobi_validate.require_fields(result, requiredFields, "Optimization output");

            if ~istable(result.peptideFeatureTable) || ~istable(result.selectedPeptides)
                error('Optimization output: peptideFeatureTable and selectedPeptides must be tables.');
            end

            if ~isempty(result.selectedPeptides)
                if ~any(strcmp(result.peptideFeatureTable.Properties.VariableNames, 'Peptide')) || ...
                        ~any(strcmp(result.selectedPeptides.Properties.VariableNames, 'Peptide'))
                    error('Optimization output: peptide tables must contain Peptide columns.');
                end

                selected = string(result.selectedPeptides.Peptide);
                universe = string(result.peptideFeatureTable.Peptide);
                if any(~ismember(selected, universe))
                    error('Optimization output: selectedPeptides must be a subset of peptideFeatureTable.Peptide.');
                end
            end

            if height(result.selectedPeptides) > result.maxPeptides
                error('Optimization output: selected peptide count exceeds maxPeptides.');
            end

            if ~istable(result.selectionHistory)
                error('Optimization output: selectionHistory must be a table.');
            end

            if ~isempty(result.selectionHistory)
                if ~any(strcmp(result.selectionHistory.Properties.VariableNames, 'Step'))
                    error('Optimization output: selectionHistory must contain Step.');
                end
                expectedSteps = (1:height(result.selectionHistory))';
                if ~isequal(result.selectionHistory.Step(:), expectedSteps)
                    error('Optimization output: selectionHistory.Step must increase strictly from 1.');
                end

                vars = string(result.selectionHistory.Properties.VariableNames);
                coverageVars = vars(contains(vars, "Coverage"));
                for i = 1:numel(coverageVars)
                    vals = result.selectionHistory.(coverageVars(i));
                    if any(vals < 0 | vals > 1 | ~isfinite(vals))
                        error('Optimization output: %s must stay in [0,1].', coverageVars(i));
                    end
                end
            end

            Mobi_validate.require_unique_string_list(result.coveredPatients, "Optimization output coveredPatients");
            Mobi_validate.require_unique_string_list(result.coveredHLAs, "Optimization output coveredHLAs");

            if ~isempty(result.peptideFamily)
                if numel(result.peptideFamily) ~= height(result.peptideFeatureTable)
                    error('Optimization output: peptideFamily length must match peptideFeatureTable height.');
                end

                if ~isequal(unique(result.peptideFamily, 'stable'), result.familyList(:))
                    error('Optimization output: familyList must match unique peptideFamily values.');
                end
            end
        end


        function optimization_suite_output(suite)
            requiredFields = [ ...
                "global", ...
                "honestGroups", ...
                "practicalGroups", ...
                "allResultsTable", ...
                "globalSelectedPeptides", ...
                "honestSelectedPeptidesByGroup", ...
                "practicalSelectedPeptidesByGroup", ...
                "honestSelectedPeptideUnion", ...
                "practicalSelectedPeptideUnion", ...
                "selectedPeptideUnionAcrossModes", ...
                "didRun", ...
                "optCfg", ...
                "summary"];
            Mobi_validate.require_fields(suite, requiredFields, "Optimization suite output");

            unionExpected = unique([ ...
                string(suite.globalSelectedPeptides(:)); ...
                string(suite.honestSelectedPeptideUnion(:)); ...
                string(suite.practicalSelectedPeptideUnion(:))], 'stable');
            unionExpected = unionExpected(unionExpected ~= "");

            actual = string(suite.selectedPeptideUnionAcrossModes(:));
            actual = actual(actual ~= "");

            if ~isequal(actual, unionExpected)
                error('Optimization suite output: selectedPeptideUnionAcrossModes does not match the mode-specific union.');
            end

            if ~istable(suite.allResultsTable)
                error('Optimization suite output: allResultsTable must be a table.');
            end
        end

    end


    methods (Static, Access = private)

        function require_fields(s, fields, label)
            for i = 1:numel(fields)
                if ~isfield(s, fields(i))
                    error('%s: missing required field %s.', string(label), fields(i));
                end
            end
        end


        function require_positive_integer_labels(labels, label)
            labels = labels(:);
            if isempty(labels) || any(labels < 1) || any(mod(labels,1) ~= 0) || any(~isfinite(labels))
                error('%s: labels must be positive integers.', string(label));
            end
        end


        function require_consecutive_labels(labels, label)
            u = unique(labels(:), 'stable');
            if ~isequal(u(:), (1:numel(u))')
                error('%s: labels must be consecutive after relabeling.', string(label));
            end
        end


        function tf = is_positive_integer_scalar(x)
            tf = isnumeric(x) && isscalar(x) && isfinite(x) && x > 0 && mod(x,1) == 0;
        end


        function val = get_field_or_default(s, fieldName, defaultVal)
            if isfield(s, fieldName) && ~isempty(s.(fieldName))
                val = s.(fieldName);
            else
                val = defaultVal;
            end
        end


        function require_unique_string_list(vals, label)
            vals = string(vals(:));
            if any(ismissing(vals))
                error('%s must not contain missing values.', string(label));
            end
            if numel(unique(vals, 'stable')) ~= numel(vals)
                error('%s must contain unique values.', string(label));
            end
        end

    end
end
