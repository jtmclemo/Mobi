classdef Mobi_utils

    methods (Static)

        function h = normalize_hla_list(h)
            h = string(h);
            h = strtrim(h);
            h(ismissing(h)) = "";

            h = replace(h, ",", "");
            h = replace(h, " ", "");
            h = upper(h);

            needsPrefix = startsWith(h, ["A*" "B*" "C*"]);
            h(needsPrefix) = "HLA-" + h(needsPrefix);

            h = h(h ~= "");
        end


        function tf = has_min_contiguous_overlap(p1, p2, minLen)
            p1 = char(string(p1));
            p2 = char(string(p2));

            tf = false;

            if length(p1) < minLen || length(p2) < minLen
                return;
            end

            for startIdx = 1:(length(p1) - minLen + 1)
                sub = p1(startIdx:startIdx + minLen - 1);
                if contains(string(p2), string(sub))
                    tf = true;
                    return;
                end
            end
        end


        function comp = connected_components_from_adjacency(A)
            n = size(A,1);
            comp = zeros(n,1);
            currentComp = 0;

            for node = 1:n
                if comp(node) ~= 0
                    continue;
                end

                currentComp = currentComp + 1;
                stack = node;

                while ~isempty(stack)
                    v = stack(end);
                    stack(end) = [];

                    if comp(v) ~= 0
                        continue;
                    end

                    comp(v) = currentComp;

                    nbrs = find(A(v,:));
                    nbrs = nbrs(comp(nbrs) == 0);
                    stack = [stack, nbrs]; %#ok<AGROW>
                end
            end
        end

    end
end
