function [chanceCI, ciAlpha] = GetChanceCI(labels,ciAlpha, N)
    if nargin < 3
        N = 10e3;
    end
    if nargin < 2
        ciAlpha = 5;
    end
    if nargin == 0
        chanceCI = ciAlpha;
        return
    end

    trueLabels = labels;
    shuffledLabels = trueLabels;
    nLblTrials = length(shuffledLabels);
    prctScale = 100/nLblTrials;
    chanceCorrect = nan(N,1);

    for ii = 1:N
        shuffledLabels = trueLabels(randi(nLblTrials,nLblTrials,1));
        chanceCorrect(ii) = sum(trueLabels == shuffledLabels)*prctScale;
    end


    prctileVals = [ciAlpha/2 100-ciAlpha/2];
    chanceCI = prctile(chanceCorrect, prctileVals);
end