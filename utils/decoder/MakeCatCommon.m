function newCatWithOrgCatOrder = MakeCatCommon(orgCat,newCat)
% Makes categories of newCat be in the same order as orgCat
% 
% Order here is defined if you specify double(newCat) any common categories
% with orgCat will match double(orgCat)'s value.
% 
% Example
% c1 = categorical({'d','b','c','z'});
% c2 = categorical({'f','c','g','a'});
% double(c1)
% double(c2)
% c3 = MakeCatCommon(c1, c2);
% double(c3)
if iscategorical(orgCat) && iscategorical(newCat)
    orgCatStr = categories(orgCat);
    newCatStr = categories(newCat);
    commonCat = union(orgCatStr, newCatStr, 'stable');
    newCatWithOrgCatOrder = setcats(newCat, commonCat ); % Make test label match trained data
else
    % Not categorical, just pass newCat to the output
    newCatWithOrgCatOrder = newCat;
end
end