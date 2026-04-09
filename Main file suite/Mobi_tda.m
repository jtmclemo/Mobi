function result = Mobi_tda(X, ids, tdaInput)

    cfg = Mobi_config();

    % INPUT HANDLING
    if nargin < 2 || isempty(ids)
        ids = "case_" + string((1:size(X,1))');
    end

    if nargin < 3 || isempty(tdaInput)
        tdaCfg = struct();
    elseif isstruct(tdaInput)
        tdaCfg = tdaInput;
    else
        tdaCfg = struct();
        tdaCfg.makePlots = logical(tdaInput);
    end

    if isfield(tdaCfg, 'makePlots') && ~isempty(tdaCfg.makePlots)
        makePlots = tdaCfg.makePlots;
    else
        makePlots = cfg.tda.makePlots;
    end

    if isfield(tdaCfg, 'distanceMetric') && ~isempty(tdaCfg.distanceMetric)
        distanceMetric = char(tdaCfg.distanceMetric);
    else
        distanceMetric = char(cfg.tda.distanceMetric);
    end

    if isfield(tdaCfg, 'thresholdRule') && ~isempty(tdaCfg.thresholdRule)
        thresholdRule = string(tdaCfg.thresholdRule);
    else
        thresholdRule = string(cfg.tda.thresholdRule);
    end

    if isfield(tdaCfg, 'percentile') && ~isempty(tdaCfg.percentile)
        configuredPercentile = tdaCfg.percentile;
    else
        configuredPercentile = cfg.tda.defaultPercentile;
    end

    if isfield(tdaCfg, 'maxOptimizationGroupSize') && ~isempty(tdaCfg.maxOptimizationGroupSize)
        maxOptimizationGroupSize = tdaCfg.maxOptimizationGroupSize;
    else
        maxOptimizationGroupSize = Inf;
    end

    if isfield(tdaCfg, 'minOptimizationGroupSize') && ~isempty(tdaCfg.minOptimizationGroupSize)
        minOptimizationGroupSize = tdaCfg.minOptimizationGroupSize;
    else
        minOptimizationGroupSize = 1;
    end

    if isfield(tdaCfg, 'mergeSingletons') && ~isempty(tdaCfg.mergeSingletons)
        mergeSingletons = logical(tdaCfg.mergeSingletons);
    else
        mergeSingletons = false;
    end

    ids = string(ids(:));
    Mobi_validate.tda_input(X, ids, tdaCfg);

    if ~isnumeric(X)
        error('X must be numeric.');
    end

    if size(X,1) < 2
        error('Need at least 2 rows in X for TDA.');
    end

    if minOptimizationGroupSize < 1 || floor(minOptimizationGroupSize) ~= minOptimizationGroupSize
        error('minOptimizationGroupSize must be a positive integer.');
    end

    if isfinite(maxOptimizationGroupSize)
        if maxOptimizationGroupSize < 1 || floor(maxOptimizationGroupSize) ~= maxOptimizationGroupSize
            error('maxOptimizationGroupSize must be a positive integer or Inf.');
        end
        if minOptimizationGroupSize > maxOptimizationGroupSize
            error('minOptimizationGroupSize cannot exceed maxOptimizationGroupSize.');
        end
    end

    % Replace non-finite entries conservatively before normalization
    X(~isfinite(X)) = NaN;
    X = fillmissing(X, 'constant', 0);

    % Z-SCORE NORMALIZATION

    Xz = zscore(X);
    Xz(:, any(~isfinite(Xz),1)) = 0;

    numFeatures = size(Xz,2);

    % WEIGHT SELECTION
    if isfield(tdaCfg, 'weights13') && ~isempty(tdaCfg.weights13)
        candidateWeights = tdaCfg.weights13;
    else
        candidateWeights = cfg.tda.defaultWeights13;
    end

    if numFeatures == numel(candidateWeights)
        weights = candidateWeights;
    elseif numFeatures == 13 && numel(cfg.tda.defaultWeights13) == 13
        weights = cfg.tda.defaultWeights13;
        warning(['Provided TDA weight vector length did not match X. ', ...
                 'Falling back to cfg.tda.defaultWeights13.']);
    else
        warning('X has %d columns, but no matching weight vector was available. Using equal weights.', numFeatures);
        weights = ones(1, numFeatures);
    end

    if numel(weights) ~= numFeatures
        error('Weight vector length must match the number of columns in X.');
    end

    % Apply weights columnwise
    Xweighted = Xz .* weights;

    % PAIRWISE DISTANCE MATRIX
    D = pdist2(Xweighted, Xweighted, distanceMetric);
    n = size(D,1);

    % EDGE LIST CONSTRUCTION
    numEdges = n*(n-1)/2;
    edges = zeros(numEdges, 3);
    c = 1;

    for i = 1:n-1
        for j = i+1:n
            edges(c,:) = [i, j, D(i,j)];
            c = c + 1;
        end
    end

    % Sort edges by distance so merges occur from smallest epsilon upward
    edges = sortrows(edges, 3);

    % BUILD H0 BARCODE
    parent = 1:n;
    compBirth = zeros(n,1);

    barcode = zeros(n,4);
    for i = 1:n
        barcode(i,:) = [0, Inf, i, i];
    end

    barcodeRowOfRoot = (1:n)';

    for e = 1:size(edges,1)

        i = edges(e,1);
        j = edges(e,2);
        epsVal = edges(e,3);

        ri = uf_find(parent, i);
        rj = uf_find(parent, j);

        if ri ~= rj

            bi = compBirth(ri);
            bj = compBirth(rj);

            % Earlier birth survives; ties broken by smaller root index
            if bi < bj
                survivor = ri;
                dead = rj;
            elseif bj < bi
                survivor = rj;
                dead = ri;
            else
                survivor = min(ri, rj);
                dead = max(ri, rj);
            end

            % Record death time of the component that disappears
            deadRow = barcodeRowOfRoot(dead);
            barcode(deadRow,2) = epsVal;

            % Merge dead component into survivor
            parent(dead) = survivor;
        end
    end

    % Keep only birth/death interval columns
    barcodeIntervals = barcode(:,1:2);

    % Sort finite intervals by death time, then append infinite bars
    finiteMask = isfinite(barcodeIntervals(:,2));
    finitePart = sortrows(barcodeIntervals(finiteMask,:), 2);
    infPart = barcodeIntervals(~finiteMask,:);
    barcodeIntervals = [finitePart; infPart];

    finiteDeaths = barcodeIntervals(isfinite(barcodeIntervals(:,2)), 2);

    % THRESHOLD SELECTION
    if isempty(finiteDeaths)
        threshold = 0;
        chosenPercentile = NaN;
        thresholdInfo = struct();
        thresholdInfo.Rule = thresholdRule;
        thresholdInfo.FiniteDeaths = finiteDeaths;
        thresholdInfo.Note = "No finite H0 deaths were available.";
    else
        switch lower(thresholdRule)

            case "percentile"
                chosenPercentile = configuredPercentile;

                if isnan(chosenPercentile) || chosenPercentile < 0 || chosenPercentile > 100
                    error('Percentile must be a numeric value between 0 and 100.');
                end

                threshold = prctile(finiteDeaths, chosenPercentile);

                thresholdInfo = struct();
                thresholdInfo.Rule = "percentile";
                thresholdInfo.ChosenPercentile = chosenPercentile;
                thresholdInfo.FiniteDeaths = finiteDeaths;

            case "largest_gap"
                sortedDeaths = sort(finiteDeaths(:), 'ascend');

                if isscalar(sortedDeaths)
                    threshold = sortedDeaths(1);
                    gapIndex = 1;
                else
                    gaps = diff(sortedDeaths);
                    [~, gapIndex] = max(gaps);
                    threshold = sortedDeaths(gapIndex);
                end

                chosenPercentile = NaN;

                thresholdInfo = struct();
                thresholdInfo.Rule = "largest_gap";
                thresholdInfo.FiniteDeaths = sortedDeaths;
                thresholdInfo.LargestGapIndex = gapIndex;

            otherwise
                error('Unsupported cfg.tda.thresholdRule: %s', thresholdRule);
        end
    end

    % CLUSTER EXTRACTION AT CHOSEN THRESHOLD

    [clusters, componentSizes] = clusters_at_threshold(D, threshold);

    % PRACTICAL OPTIMIZATION GROUPING
    [optimizationGroups, practicalGroupingAudit] = build_optimization_groups( ...
        D, clusters, maxOptimizationGroupSize, minOptimizationGroupSize, mergeSingletons, ids);

    optimizationGroupSizes = accumarray(optimizationGroups, 1);
    optimizationGroupMembers = group_members_from_labels(optimizationGroups, ids);
    optimizationGroupFileNumbers = group_file_numbers_from_labels(optimizationGroups);

    % PACKAGE OUTPUT
    result = struct();
    result.X_original                 = X;
    result.X_zscore                   = Xz;
    result.X_weighted                 = Xweighted;
    result.weights                    = weights;
    result.distanceMetric             = distanceMetric;
    result.D                          = D;
    result.edges                      = edges;
    result.barcode                    = barcodeIntervals;
    result.finiteDeaths               = finiteDeaths;
    result.ids                        = ids;
    result.threshold                  = threshold;
    result.thresholdRule              = thresholdRule;
    result.thresholdInfo              = thresholdInfo;
    result.chosenPercentile           = chosenPercentile;

    % Raw TDA clustering
    result.clusters                   = clusters;
    result.groupLabels                = clusters;
    result.componentSizes             = componentSizes;
    result.rawClusters                = clusters;
    result.rawComponentSizes          = componentSizes;

    % Practical optimization grouping
    result.optimizationGroups         = optimizationGroups;
    result.optimizationGroupLabels    = optimizationGroups;
    result.optimizationGroupSizes     = optimizationGroupSizes;
    result.optimizationGroupMembers   = optimizationGroupMembers;
    result.optimizationGroupFileNumbers = optimizationGroupFileNumbers;
    result.practicalGroupingAudit = practicalGroupingAudit;
    result.mergeSingletons            = mergeSingletons;
    result.maxOptimizationGroupSize   = maxOptimizationGroupSize;
    result.minOptimizationGroupSize   = minOptimizationGroupSize;
    Mobi_validate.tda_output(result, X, ids);

    % OPTIONAL PLOTS
    if makePlots

        % Distance matrix heatmap
        figure;
        imagesc(D);
        colormap([ones(256,1), linspace(0,1,256)', zeros(256,1)]);
        colorbar;
        axis square;
        title(sprintf('Pairwise distance matrix (%s, weighted z-score feature space)', distanceMetric));
        xlabel('Cases');
        ylabel('Cases');

        % H0 barcode plot
        figure;
        hold on;

        for k = 1:size(barcodeIntervals,1)

            b = barcodeIntervals(k,1);
            d = barcodeIntervals(k,2);

            if isinf(d)
                dplot = max(edges(:,3))*1.05;
            else
                dplot = d;
            end

            plot([b dplot], [k k], 'LineWidth', 2);
        end

        hold off;
        xlabel('\epsilon');
        ylabel('Connected component interval');
        title('H_0 barcode (0D persistence)');
        grid on;

        % Finite death plot
        if ~isempty(finiteDeaths)
            figure;
            plot(sort(finiteDeaths, 'ascend'), 'LineWidth', 2);
            hold on;
            yline(threshold, '--', 'LineWidth', 1.5);
            hold off;
            xlabel('Finite death index');
            ylabel('Death value');
            title(sprintf('Finite H_0 death times with chosen threshold = %.4f', threshold));
            grid on;
        end

        % PCA projection for visualization only
        [~, score, ~] = pca(Xweighted);

        figure;
        gscatter(score(:,1), score(:,2), clusters);
        text(score(:,1), score(:,2), "  " + ids);
        xlabel('PC1');
        ylabel('PC2');
        title(sprintf('PCA view colored by raw TDA clusters at \\epsilon = %.3f', threshold));
        grid on;
    end
end


% UNION-FIND HELPER
function r = uf_find(parent, x)

    r = x;
    while parent(r) ~= r
        r = parent(r);
    end
end


% CLUSTERS AT THRESHOLD
function [clusters, componentSizes] = clusters_at_threshold(D, threshold)

    n = size(D,1);

    % Build threshold graph
    A = (D <= threshold) - eye(n);
    A = A > 0;

    clusters = zeros(n,1);
    visited = false(n,1);
    cid = 0;

    for i = 1:n

        if ~visited(i)

            cid = cid + 1;

            queue = i;
            visited(i) = true;
            clusters(i) = cid;

            while ~isempty(queue)

                v = queue(1);
                queue(1) = [];

                nbrs = find(A(v,:));

                for u = nbrs
                    if ~visited(u)
                        visited(u) = true;
                        clusters(u) = cid;
                        queue(end+1) = u; %#ok<AGROW>
                    end
                end
            end
        end
    end

    componentSizes = accumarray(clusters, 1);
end


% PRACTICAL OPTIMIZATION GROUP BUILDER
function [optimizationGroups, audit] = build_optimization_groups( ...
    D, rawClusters, maxGroupSize, minGroupSize, mergeSingletons, ids)

    optimizationGroups = relabel_groups_consecutively(rawClusters(:));
    initialGroups = optimizationGroups;
    auditEvents = empty_audit_events();
    stepNumber = 0;
    [auditEvents, stepNumber] = add_audit_event( ...
        auditEvents, stepNumber, "initial", "noop", rawClusters, ...
        optimizationGroups, optimizationGroups, (1:numel(optimizationGroups))', ids, ...
        "raw TDA clusters used as initial practical groups", NaN, ...
        minGroupSize, maxGroupSize);

    n = numel(optimizationGroups);

    if n == 0
        audit = build_practical_grouping_audit(rawClusters, initialGroups, ...
            optimizationGroups, ids, minGroupSize, maxGroupSize, auditEvents);
        return;
    end

    % If the minimum group size is impossible to satisfy globally, fall back
    % to a single group. This avoids stranded leftovers.
    if minGroupSize > n
        minGroupSize = n;
    end

    % optionally merge singleton groups immediately
    if mergeSingletons
        [optimizationGroups, auditEvents, stepNumber] = repair_small_groups( ...
            D, rawClusters, optimizationGroups, maxGroupSize, 2, ids, ...
            auditEvents, stepNumber, "merge_small_group", ...
            "singleton merged to nearest compatible group");
    end

    % split oversized groups, repairing after each split
    if isfinite(maxGroupSize)
        safetyCounter = 0;
        maxIterations = max(50, 10 * n);

        while true
            safetyCounter = safetyCounter + 1;
            if safetyCounter > maxIterations
                warning('Reached maximum iterations while building optimization groups.');
                break;
            end

            groupSizes = accumarray(optimizationGroups, 1);
            oversizedGroups = find(groupSizes > maxGroupSize);

            if isempty(oversizedGroups)
                break;
            end

            % Split the largest oversized group first
            [~, idxLargest] = max(groupSizes(oversizedGroups));
            g = oversizedGroups(idxLargest);

            thisIdx = find(optimizationGroups == g);
            Dsub = D(thisIdx, thisIdx);

            [splitLabels, splitScore] = split_group_by_distance_balanced(Dsub, minGroupSize, maxGroupSize);

            % If the split failed to produce two actual parts, stop trying on this group
            if numel(unique(splitLabels)) < 2
                break;
            end

            group1Idx = thisIdx(splitLabels == 1);
            group2Idx = thisIdx(splitLabels == 2);

            if isempty(group1Idx) || isempty(group2Idx)
                break;
            end

            labelsBefore = optimizationGroups;
            newLabel = max(optimizationGroups) + 1;
            optimizationGroups(group2Idx) = newLabel;
            optimizationGroups = relabel_groups_consecutively(optimizationGroups);
            [auditEvents, stepNumber] = add_audit_event( ...
                auditEvents, stepNumber, "split_large_group", "split", rawClusters, ...
                labelsBefore, optimizationGroups, thisIdx, ids, ...
                "group exceeded maxOptimizationGroupSize", splitScore, ...
                minGroupSize, maxGroupSize);

            % Immediately repair any new undersized groups created by splitting
            [optimizationGroups, auditEvents, stepNumber] = repair_small_groups( ...
                D, rawClusters, optimizationGroups, maxGroupSize, minGroupSize, ids, ...
                auditEvents, stepNumber, "repair_small_group", ...
                "repair after split created undersized child");
        end
    end

    % final repair to enforce the true minimum size
    [optimizationGroups, auditEvents, stepNumber] = repair_small_groups( ...
        D, rawClusters, optimizationGroups, maxGroupSize, minGroupSize, ids, ...
        auditEvents, stepNumber, "repair_small_group", ...
        "final minimum-size repair pass");

    % Final cleanup
    labelsBefore = optimizationGroups;
    optimizationGroups = relabel_groups_consecutively(optimizationGroups);
    if ~isequal(labelsBefore, optimizationGroups)
        [auditEvents, stepNumber] = add_audit_event( ...
            auditEvents, stepNumber, "final", "relabel", rawClusters, ...
            labelsBefore, optimizationGroups, (1:n)', ids, ...
            "final relabeling to consecutive practical group labels", NaN, ...
            minGroupSize, maxGroupSize);
    end

    [auditEvents, ~] = add_audit_event( ...
        auditEvents, stepNumber, "final", "noop", rawClusters, ...
        optimizationGroups, optimizationGroups, (1:n)', ids, ...
        "final practical grouping summary", NaN, ...
        minGroupSize, maxGroupSize);
    audit = build_practical_grouping_audit(rawClusters, initialGroups, ...
        optimizationGroups, ids, minGroupSize, maxGroupSize, auditEvents);
end


function [splitLabels, bestScore] = split_group_by_distance_balanced(Dsub, minGroupSize, maxGroupSize)

    n = size(Dsub,1);

    if n <= 1
        splitLabels = ones(n,1);
        bestScore = NaN;
        return;
    elseif n == 2
        splitLabels = [1; 2];
        bestScore = 0;
        return;
    end

    Y = squareform(Dsub, 'tovector');
    Z = linkage(Y, 'average');

    bestLabels = [];
    bestScore = Inf;

    % Try all possible 2-way cuts induced by the dendrogram heights
    for k = 1:(n-1)
        candidate = cluster(Z, 'maxclust', k+1);

        % Only keep candidates with exactly 2 groups
        u = unique(candidate);
        if numel(u) ~= 2
            continue;
        end

        labels2 = zeros(size(candidate));
        labels2(candidate == u(1)) = 1;
        labels2(candidate == u(2)) = 2;

        sizes = [sum(labels2 == 1), sum(labels2 == 2)];

        if all(sizes >= 1)
            penalty = 0;

            if any(sizes < minGroupSize)
                penalty = penalty + 1e6 * sum(max(0, minGroupSize - sizes));
            end

            if any(sizes > maxGroupSize)
                penalty = penalty + 1e6 * sum(max(0, sizes - maxGroupSize));
            end

            % Prefer more balanced splits when penalties tie
            balanceScore = abs(sizes(1) - sizes(2));
            score = penalty + balanceScore;

            if score < bestScore
                bestScore = score;
                bestLabels = labels2;
            end
        end
    end

    % Fallback if no better split found
    if isempty(bestLabels)
        bestLabels = cluster(Z, 'maxclust', 2);

        u = unique(bestLabels);
        labels2 = zeros(size(bestLabels));
        labels2(bestLabels == u(1)) = 1;
        if numel(u) >= 2
            labels2(bestLabels == u(2)) = 2;
        else
            labels2(1:floor(n/2)) = 1;
            labels2(floor(n/2)+1:end) = 2;
        end
        bestLabels = labels2;
    end

    splitLabels = bestLabels(:);
end


function [labels, auditEvents, stepNumber] = repair_small_groups( ...
    D, rawClusters, labels, maxGroupSize, minGroupSize, ids, ...
    auditEvents, stepNumber, stage, reason)

    labels = relabel_groups_consecutively(labels(:));

    if minGroupSize <= 1
        return;
    end

    safetyCounter = 0;
    n = numel(labels);
    maxIterations = max(50, 10 * n);

    while true
        safetyCounter = safetyCounter + 1;
        if safetyCounter > maxIterations
            warning('Reached maximum iterations while repairing small groups.');
            break;
        end

        groupSizes = accumarray(labels, 1);
        smallGroups = find(groupSizes < minGroupSize);

        if isempty(smallGroups)
            break;
        end

        % Repair the smallest group first
        [~, idxSmallest] = min(groupSizes(smallGroups));
        gSmall = smallGroups(idxSmallest);
        smallIdx = find(labels == gSmall);
        smallSize = numel(smallIdx);

        candidateGroups = setdiff(unique(labels), gSmall);

        if isempty(candidateGroups)
            break;
        end

        bestGroup = NaN;
        bestScore = Inf;

        for g = candidateGroups(:)'
            memberIdx = find(labels == g);
            candidateSize = numel(memberIdx);

            if isfinite(maxGroupSize) && (candidateSize + smallSize > maxGroupSize)
                continue;
            end

            meanDist = mean(D(smallIdx, memberIdx), 'all', 'omitnan');

            % Prefer smaller distance; slight preference for landing near min size
            resultingSize = candidateSize + smallSize;
            sizePenalty = abs(resultingSize - max(minGroupSize, min(resultingSize, maxGroupSize)));
            score = meanDist + 1e-6 * sizePenalty;

            if score < bestScore
                bestScore = score;
                bestGroup = g;
            end
        end

        % If no candidate respects max size, merge into nearest anyway.
        % This avoids impossible stranded groups when constraints conflict.
        if isnan(bestGroup)
            for g = candidateGroups(:)'
                memberIdx = find(labels == g);
                meanDist = mean(D(smallIdx, memberIdx), 'all', 'omitnan');

                if meanDist < bestScore
                    bestScore = meanDist;
                    bestGroup = g;
                end
            end
        end

        if isnan(bestGroup)
            break;
        end

        labelsBefore = labels;
        labels(smallIdx) = bestGroup;
        labels = relabel_groups_consecutively(labels);
        [auditEvents, stepNumber] = add_audit_event( ...
            auditEvents, stepNumber, stage, "merge", rawClusters, ...
            labelsBefore, labels, smallIdx, ids, reason, bestScore, ...
            minGroupSize, maxGroupSize);
    end
end


function auditEvents = empty_audit_events()

    auditEvents = repmat(struct( ...
        'StepNumber', 0, ...
        'Stage', "", ...
        'ActionType', "", ...
        'SourceRawGroups', [], ...
        'SourcePracticalGroupsBefore', [], ...
        'TargetPracticalGroupsAfter', [], ...
        'AffectedFileIndices', [], ...
        'AffectedFileNames', strings(0,1), ...
        'PreviousGroupSizes', [], ...
        'NewGroupSizes', [], ...
        'Reason', "", ...
        'DistanceEvidence', NaN, ...
        'MinGroupSize', NaN, ...
        'MaxGroupSize', NaN), 0, 1);
end


function [auditEvents, stepNumber] = add_audit_event( ...
    auditEvents, stepNumber, stage, actionType, rawClusters, labelsBefore, ...
    labelsAfter, affectedIdx, ids, reason, distanceEvidence, minGroupSize, maxGroupSize)

    stepNumber = stepNumber + 1;
    affectedIdx = affectedIdx(:);

    event = struct();
    event.StepNumber = stepNumber;
    event.Stage = string(stage);
    event.ActionType = string(actionType);
    event.SourceRawGroups = unique(rawClusters(affectedIdx), 'stable')';
    event.SourcePracticalGroupsBefore = unique(labelsBefore(affectedIdx), 'stable')';
    event.TargetPracticalGroupsAfter = unique(labelsAfter(affectedIdx), 'stable')';
    event.AffectedFileIndices = affectedIdx';
    event.AffectedFileNames = string(ids(affectedIdx));
    event.PreviousGroupSizes = group_size_vector(labelsBefore);
    event.NewGroupSizes = group_size_vector(labelsAfter);
    event.Reason = string(reason);
    event.DistanceEvidence = distanceEvidence;
    event.MinGroupSize = minGroupSize;
    event.MaxGroupSize = maxGroupSize;

    auditEvents(end+1,1) = event; %#ok<AGROW>
end


function audit = build_practical_grouping_audit( ...
    rawClusters, initialGroups, finalGroups, ids, minGroupSize, maxGroupSize, events)

    rawClusters = rawClusters(:);
    initialGroups = initialGroups(:);
    finalGroups = finalGroups(:);
    ids = string(ids(:));

    FileNumber = (1:numel(ids))';
    FileName = ids;
    RawCluster = rawClusters;
    InitialPracticalGroup = initialGroups;
    PracticalGroup = finalGroups;
    MovedFromInitialPracticalGroup = InitialPracticalGroup ~= PracticalGroup;
    AuditReasonCode = strings(numel(ids),1);

    for i = 1:numel(ids)
        reasons = strings(0,1);
        for e = 1:numel(events)
            isFinalNoop = events(e).Stage == "final" && events(e).ActionType == "noop";
            if any(events(e).AffectedFileIndices == i) && events(e).Stage ~= "initial" && ~isFinalNoop
                reasons(end+1,1) = events(e).Stage + ":" + events(e).ActionType; %#ok<AGROW>
            end
        end
        if isempty(reasons)
            AuditReasonCode(i) = "unchanged";
        else
            AuditReasonCode(i) = strjoin(unique(reasons, 'stable'), ";");
        end
    end

    audit = struct();
    audit.events = events;
    audit.rawToPracticalMapping = table( ...
        FileNumber, FileName, RawCluster, InitialPracticalGroup, PracticalGroup, ...
        MovedFromInitialPracticalGroup, AuditReasonCode);
    audit.summary = struct();
    audit.summary.NumEvents = numel(events);
    audit.summary.NumMerges = sum([events.ActionType] == "merge");
    audit.summary.NumSplits = sum([events.ActionType] == "split");
    audit.summary.NumRepairSteps = sum([events.Stage] == "repair_small_group");
    audit.summary.NumRelabels = sum([events.ActionType] == "relabel");
    audit.summary.NumFilesChanged = sum(MovedFromInitialPracticalGroup);
    audit.summary.FinalPracticalGroupSizes = group_size_vector(finalGroups);
    audit.summary.MinGroupSize = minGroupSize;
    audit.summary.MaxGroupSize = maxGroupSize;
end


function sizes = group_size_vector(labels)

    labels = labels(:);
    if isempty(labels)
        sizes = [];
    else
        sizes = accumarray(labels, 1)';
    end
end


function memberCells = group_members_from_labels(labels, ids)

    labels = labels(:);
    ids = string(ids(:));
    uniqueGroups = unique(labels, 'stable');

    memberCells = cell(numel(uniqueGroups), 1);

    for i = 1:numel(uniqueGroups)
        memberCells{i} = ids(labels == uniqueGroups(i));
    end
end


function fileNumberCells = group_file_numbers_from_labels(labels)

    labels = labels(:);
    uniqueGroups = unique(labels, 'stable');

    fileNumberCells = cell(numel(uniqueGroups), 1);

    for i = 1:numel(uniqueGroups)
        fileNumberCells{i} = find(labels == uniqueGroups(i))';
    end
end


function labels = relabel_groups_consecutively(labels)

    u = unique(labels(:), 'stable');
    newLabels = zeros(size(labels));

    for i = 1:numel(u)
        newLabels(labels == u(i)) = i;
    end

    labels = newLabels;
end
