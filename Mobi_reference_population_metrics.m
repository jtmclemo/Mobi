function metrics = Mobi_reference_population_metrics(selectedPeptides, combinedTable, referenceTable, options)

    if nargin < 4 || isempty(options)
        options = struct();
    end

    options = local_normalize_options(options);

    metrics = struct();
    metrics.coveragePercent = 0;
    metrics.expectedHits = 0;
    metrics.conditionalMeanHits = 0;
    metrics.pc90 = 0;
    metrics.pc90Semantics = "weighted mixture across reference populations";
    metrics.populationWeightMode = string(options.populationWeightMode);
    metrics.hitCountMode = string(options.hitCountMode);
    metrics.numSelectedPeptides = 0;
    metrics.numHitHLAs = 0;
    metrics.numReferencePopulations = 0;
    metrics.populationMetrics = table();
    metrics.hitAlleleSummary = table();
    metrics.weightedMixturePMF = [];
    metrics.didEvaluate = false;

    selectedPeptides = unique(string(selectedPeptides(:)), 'stable');
    selectedPeptides = selectedPeptides(selectedPeptides ~= "");
    metrics.numSelectedPeptides = numel(selectedPeptides);

    if isempty(selectedPeptides) || ~istable(combinedTable) || isempty(combinedTable) || ...
            ~istable(referenceTable) || isempty(referenceTable)
        return;
    end

    requiredCombined = ["Peptide","MHC"];
    requiredReference = ["Population","Locus","Allele","AlleleFrequency"];

    if ~all(ismember(requiredCombined, string(combinedTable.Properties.VariableNames)))
        return;
    end

    if ~all(ismember(requiredReference, string(referenceTable.Properties.VariableNames)))
        return;
    end

    hitRows = combinedTable(ismember(string(combinedTable.Peptide), selectedPeptides), :);
    if isempty(hitRows)
        return;
    end

    hitRows = unique(hitRows(:, {'Peptide','MHC'}), 'rows', 'stable');
    hitRows.MHC = Mobi_utils.normalize_hla_list(hitRows.MHC);

    [hitAlleles, hitCounts] = local_build_hit_counts(hitRows, options.hitCountMode);
    if isempty(hitAlleles)
        return;
    end

    metrics.hitAlleleSummary = table(hitAlleles, hitCounts, ...
        'VariableNames', {'Allele','HitCount'});
    metrics.numHitHLAs = numel(hitAlleles);

    ref = local_prepare_reference_table(referenceTable, options.tolerance);
    if isempty(ref)
        return;
    end

    pops = unique(ref.Population, 'stable');
    loci = unique(ref.Locus, 'stable');

    popCoverage = zeros(numel(pops), 1);
    popExpectedHits = zeros(numel(pops), 1);
    popConditionalMeanHits = zeros(numel(pops), 1);
    popPC90 = zeros(numel(pops), 1);
    popWeights = zeros(numel(pops), 1);
    popPMFs = cell(numel(pops), 1);

    for p = 1:numel(pops)
        popName = pops(p);
        S = ref(ref.Population == popName, :);

        totalPMF = 1;
        observedLoci = 0;

        for l = 1:numel(loci)
            locusName = loci(l);
            L = S(S.Locus == locusName, :);
            if isempty(L)
                continue;
            end

            observedLoci = observedLoci + 1;
            locusPMF = local_population_locus_pmf(L, hitAlleles, hitCounts, options);
            totalPMF = conv(totalPMF, locusPMF);
        end

        if observedLoci < options.minLociForStableEstimate
            warning('Mobi:ReferenceMetrics:SparsePopulation', ...
                'Population "%s" has only %d locus block(s); estimates may be unstable.', ...
                popName, observedLoci);
        end

        totalPMF = totalPMF(:)';
        totalPMF = totalPMF / sum(totalPMF);
        support = 0:(numel(totalPMF)-1);

        popCoverage(p) = 1 - totalPMF(1);
        popExpectedHits(p) = sum(support .* totalPMF);

        if popCoverage(p) > 0
            popConditionalMeanHits(p) = popExpectedHits(p) / popCoverage(p);
        else
            popConditionalMeanHits(p) = 0;
        end

        popPC90(p) = local_pc90_from_pmf(totalPMF);
        popWeights(p) = local_population_weight(S, options.populationWeightMode);
        popPMFs{p} = totalPMF;
    end

    if all(popWeights <= 0 | isnan(popWeights))
        popWeights = ones(size(popWeights));
    end

    popWeights = popWeights(:);
    popWeights = popWeights / sum(popWeights);
    mixturePMF = local_weighted_mixture_pmf(popPMFs, popWeights);

    weightedCoverage = sum(popWeights .* popCoverage);
    weightedExpectedHits = sum(popWeights .* popExpectedHits);

    metrics.coveragePercent = 100 * weightedCoverage;
    metrics.expectedHits = weightedExpectedHits;
    if weightedCoverage > 0
        metrics.conditionalMeanHits = weightedExpectedHits / weightedCoverage;
    else
        metrics.conditionalMeanHits = 0;
    end
    metrics.pc90 = local_pc90_from_pmf(mixturePMF);
    metrics.weightedMixturePMF = mixturePMF;
    metrics.numReferencePopulations = numel(pops);
    metrics.populationMetrics = table( ...
        pops(:), ...
        popWeights(:), ...
        popCoverage(:) * 100, ...
        popExpectedHits(:), ...
        popConditionalMeanHits(:), ...
        popPC90(:), ...
        'VariableNames', {'Population','Weight','CoveragePercent','ExpectedHits','ConditionalMeanHits','PC90'});
    metrics.didEvaluate = true;
end


function options = local_normalize_options(options)
    if ~isfield(options, 'populationWeightMode') || strlength(string(options.populationWeightMode)) == 0
        options.populationWeightMode = "sample_size";
    else
        options.populationWeightMode = lower(strtrim(string(options.populationWeightMode)));
    end

    if ~ismember(options.populationWeightMode, ["sample_size","equal"])
        error('populationWeightMode must be "sample_size" or "equal".');
    end

    if ~isfield(options, 'hitCountMode') || strlength(string(options.hitCountMode)) == 0
        options.hitCountMode = "multiplicity";
    else
        options.hitCountMode = lower(strtrim(string(options.hitCountMode)));
    end

    if ~ismember(options.hitCountMode, ["multiplicity","binary_per_allele"])
        error('hitCountMode must be "multiplicity" or "binary_per_allele".');
    end

    if ~isfield(options, 'tolerance') || isempty(options.tolerance)
        options.tolerance = 1e-9;
    end

    if ~isfield(options, 'minLociForStableEstimate') || isempty(options.minLociForStableEstimate)
        options.minLociForStableEstimate = 2;
    end
end


function [hitAlleles, hitCounts] = local_build_hit_counts(hitRows, hitCountMode)
    hitAlleles = unique(string(hitRows.MHC), 'stable');
    hitAlleles = hitAlleles(hitAlleles ~= "");

    hitCounts = zeros(numel(hitAlleles), 1);
    for i = 1:numel(hitAlleles)
        switch string(hitCountMode)
            case "multiplicity"
                hitCounts(i) = sum(string(hitRows.MHC) == hitAlleles(i));
            case "binary_per_allele"
                hitCounts(i) = any(string(hitRows.MHC) == hitAlleles(i));
        end
    end
end


function ref = local_prepare_reference_table(referenceTable, tolerance)
    ref = referenceTable;
    ref.Population = strtrim(string(ref.Population));
    ref.Locus = upper(strtrim(string(ref.Locus)));
    ref.Allele = Mobi_utils.normalize_hla_list(ref.Allele);

    ref.AlleleFrequency = local_to_fraction(ref.AlleleFrequency);
    if any(strcmp(string(ref.Properties.VariableNames), 'IndividualsFrequency'))
        ref.IndividualsFrequency = local_to_fraction(ref.IndividualsFrequency);
    else
        ref.IndividualsFrequency = nan(height(ref), 1);
    end

    if any(strcmp(string(ref.Properties.VariableNames), 'SampleSize'))
        ref.SampleSize = local_to_numeric(ref.SampleSize);
        badWeights = ~isnan(ref.SampleSize) & ref.SampleSize <= 0;
        if any(badWeights)
            warning('Mobi:ReferenceMetrics:NonpositiveSampleSize', ...
                'Nonpositive sample sizes were detected; equal-weight fallback may be used.');
        end
    else
        ref.SampleSize = nan(height(ref), 1);
    end

    ref = ref(ref.Population ~= "" & ref.Locus ~= "" & ref.Allele ~= "", :);
    if isempty(ref)
        return;
    end

    missingAlleleFreq = isnan(ref.AlleleFrequency) & ~isnan(ref.IndividualsFrequency);
    ref.AlleleFrequency(missingAlleleFreq) = 1 - sqrt(max(0, 1 - ref.IndividualsFrequency(missingAlleleFreq)));

    badFreq = ~isnan(ref.AlleleFrequency) & ...
        (ref.AlleleFrequency < -tolerance | ref.AlleleFrequency > 1 + tolerance);
    if any(badFreq)
        error('Mobi:ReferenceMetrics:InvalidAlleleFrequency', ...
            'Reference allele frequencies must lie in [0,1] after normalization.');
    end

    ref.AlleleFrequency = min(1, max(0, ref.AlleleFrequency));
    ref = ref(~isnan(ref.AlleleFrequency), :);
end


function locusPMF = local_population_locus_pmf(L, hitAlleles, hitCounts, options)
    [alleles, ~, ia] = unique(string(L.Allele), 'stable');
    alleleFreqs = accumarray(ia, L.AlleleFrequency, [], @max);
    alleleFreqs = double(alleleFreqs(:));

    [tfHit, locHitIdx] = ismember(alleles, hitAlleles);
    alleleHitCounts = zeros(numel(alleles), 1);
    alleleHitCounts(tfHit) = hitCounts(locHitIdx(tfHit));

    pKnown = sum(alleleFreqs, 'omitnan');
    if pKnown > 1 + options.tolerance
        warning('Mobi:ReferenceMetrics:AlleleMassExceedsOne', ...
            'Observed allele mass %.6f exceeds 1 for a population-locus block; truncating residual mass to zero.', ...
            pKnown);
    end

    if numel(alleles) < 2
        warning('Mobi:ReferenceMetrics:SparseLocusBlock', ...
            'A population-locus block has fewer than 2 observed alleles; estimates may be unstable.');
    end

    pOther = max(0, 1 - pKnown);
    locusAlleleFreqs = [alleleFreqs; pOther];
    locusHitCounts = [alleleHitCounts; 0];

    locusPMF = local_locus_hit_pmf(locusAlleleFreqs, locusHitCounts);
end


function values = local_to_numeric(values)
    if ~isnumeric(values)
        values = str2double(string(values));
    else
        values = double(values);
    end
end


function values = local_to_fraction(values)
    values = local_to_numeric(values);
    idxPercent = values > 1;
    values(idxPercent) = values(idxPercent) / 100;
end


function locusPMF = local_locus_hit_pmf(alleleFreqs, alleleHitCounts)
    alleleFreqs = double(alleleFreqs(:));
    alleleHitCounts = double(alleleHitCounts(:));

    alleleFreqs = max(alleleFreqs, 0);
    total = sum(alleleFreqs);
    if total <= 0
        locusPMF = 1;
        return;
    end

    alleleFreqs = alleleFreqs / total;

    maxHits = max(alleleHitCounts);
    locusPMF = zeros(1, 2*maxHits + 1);

    n = numel(alleleFreqs);
    for i = 1:n
        for j = i:n
            if i == j
                prob = alleleFreqs(i)^2;
            else
                prob = 2 * alleleFreqs(i) * alleleFreqs(j);
            end

            hitTotal = alleleHitCounts(i);
            if i ~= j
                hitTotal = hitTotal + alleleHitCounts(j);
            end

            locusPMF(hitTotal + 1) = locusPMF(hitTotal + 1) + prob;
        end
    end
end


function pc90 = local_pc90_from_pmf(pmf)
    coverageTail = fliplr(cumsum(fliplr(pmf)));
    hits = 0:(numel(pmf)-1);
    idx = find(coverageTail >= 0.90, 1, 'last');

    if isempty(idx)
        pc90 = 0;
    else
        pc90 = hits(idx);
    end
end


function weight = local_population_weight(S, populationWeightMode)
    switch string(populationWeightMode)
        case "equal"
            weight = 1;
        case "sample_size"
            weight = 1;
            if any(~isnan(S.SampleSize))
                weight = max(S.SampleSize, [], 'omitnan');
                if ~isfinite(weight) || weight <= 0
                    warning('Mobi:ReferenceMetrics:FallbackEqualWeight', ...
                        'Population weight fell back to equal weighting because sample size was invalid.');
                    weight = 1;
                end
            end
    end
end


function mixturePMF = local_weighted_mixture_pmf(pmfs, weights)
    maxLen = 0;
    for i = 1:numel(pmfs)
        maxLen = max(maxLen, numel(pmfs{i}));
    end

    mixturePMF = zeros(1, maxLen);
    for i = 1:numel(pmfs)
        thisPMF = pmfs{i};
        mixturePMF(1:numel(thisPMF)) = mixturePMF(1:numel(thisPMF)) + weights(i) * thisPMF;
    end

    mixturePMF = mixturePMF / sum(mixturePMF);
end
