function [pKW,pRS] = KWandRSTests(data,conditionNames,alpha)
    % data - matric nTrials x nConditions
    % conditionNames - string array of names of each condition
    % alpha - (optional) p val threshold to test (default: 0.05)
    % dataName - (optional) specifier of data to ease ready of printed results

    if nargin < 3
        alpha = 0.05;
    end

    % Kruskal Wallis test w/ Multiple Comparisons
    [pKW,tbl,stat] = kruskalwallis(data,[],'off');
    fprintf('\nKruskal-Wallis Test:\n p = %0.3g\n', pKW)
    
    % Multiple Comparisons
    fprintf('\nmultiple comparisons, bonferroni corrected:\n')
    mc = multcompare(stat,'Alpha',alpha,'CType','bonferroni',"Display","off");
    multcomp = cat(2,conditionNames(mc(:,1))',conditionNames(mc(:,2))',num2cell(round(mc(:,6),3,'significant')));
    MCtable = table(multcomp(:,1),multcomp(:,2),multcomp(:,3),'VariableNames',{'TrialType1','TrialType2','pval'});
    disp(MCtable)

    % Pairwise Wilcoxon Rank Sum tests
    nCond = length(conditionNames); %number of conditions /groups
    combinations = nchoosek(1:nCond,2);

    pRS = ones(size(combinations,1),1);
    h = zeros(size(combinations,1),1);

    fprintf('Pairwise Wilcoxon Rank Sum tests:\n')
    for i = 1:size(combinations,1)
        c1 = combinations(i,1);
        c2 = combinations(i,2);
        [pRS(i),h(i)] = ranksum(data(:,c1),data(:,c2),'alpha',alpha);     
%         disp([conditionNames{c1} ' vs ' conditionNames{c2} ' : p = ' num2str(pRS(i)) ]);
        fprintf(' - %s vs %s: p = %0.3g\n',conditionNames{c1}, conditionNames{c2}, pRS(i))
    end

end