function [c,ia] = LIGHT_setdiff(a,b)
% [c,ia] = LIGHT_setdiff(a,b) Light version of sediff extracted
% from the original MatLab SETDIFF code.
%

% Antti Pulkkinen, April 2017.

% % Convert to columns.
% a = a(:);
% b = b(:);

% Call ISMEMBER to determine list of non-matching elements of A.
logUA = ~(ismember(a,b));
c = a(logUA);

% % Call UNIQUE to remove duplicates from list of non-matches. Do not sort
% % to reduce operations.
% [c,ndx] = unique(c,'stable');
% 
% % Find indices by using logUA and NDX.
% indlogUA = find(logUA);
% ia = indlogUA(ndx);

% Find indices by using logUA. NO REMOVAL OF POSSIBLY REPEATED INDICES.
ia = find(logUA);


