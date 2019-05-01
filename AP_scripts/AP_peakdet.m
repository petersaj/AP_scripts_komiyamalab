function [maxtab, mintab]=ap_peakdet(v, delta, x)
%PEAKDET Detect peaks in a vector
%        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
%        maxima and minima ("peaks") in the vector V.
%        MAXTAB and MINTAB consists of two columns. Column 1
%        contains indices in V, and column 2 the found values.
%      
%        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
%        in MAXTAB and MINTAB are replaced with the corresponding
%        X-values.
%
%        A point is considered a maximum peak if it has the maximal
%        value, and was preceded (to the left) by a value lower by
%        DELTA. DELTA can ve a vector, to only find peaks which satistfy
%        given order of criterion in Schmitt-trigger style (AP)
%        Negative values mean it has to be UNDER that amount (12/2/11 AP)


% Eli Billauer, 3.4.05 (Explicitly not copyrighted).
% This function is released to the public domain; Any use is allowed.


maxtab = [];
mintab = [];

if nargin < 3
  x = (1:length(v))';
else 
  x = x(:);
  if length(v)~= length(x)
    error('Input vectors v and x must have same length');
  end
end

if delta <= 0
  error('Input argument DELTA must be positive');
end

trigger_length = length(delta);
delta_pos = find(delta > 0);
delta_neg = find(delta < 0);
trigger_length_forward = length(delta_pos);
trigger_length_backward = length(delta_neg);

mn = Inf; mx = -Inf;
mnpos = NaN; mxpos = NaN;

lookformax = 1;

for i=1+trigger_length_backward:length(v)-trigger_length_forward
    this = v(i);
    if this > mx, mx = this; mxpos = x(i); end
    if this < mn, mn = this; mnpos = x(i); end
    
    if lookformax
        if this < mx-delta(delta_pos(1))
            maxtab = [maxtab ; mxpos mx];
            mn = this; mnpos = x(i);
            lookformax = 0;
        end
    else
        if ~isempty(delta_neg)
            if v(i:i+trigger_length_forward-1) >= mn+delta(delta_pos) &&  v(i-trigger_length_backward:i-1) < mn+delta(delta_neg)
                mintab = [mintab ; mnpos mn];
                mx = this; mxpos = x(i);
                lookformax = 1;
            end
        elseif isempty(delta_neg)
            if v(i:i+trigger_length_forward-1) >= mn+delta(delta_pos)
                mintab = [mintab ; mnpos mn];
                mx = this; mxpos = x(i);
                lookformax = 1;
            end
        end
    end
end

