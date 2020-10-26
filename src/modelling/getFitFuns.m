function fun = getFitFuns(expmnt, md, metric, noOff)

if expmnt == "affinities" && md=="simpleweak" && metric=="fraction"
    fun = @(x, p)subfun_simplebinding_weak_fraction_std2(x, p);
elseif expmnt == "affinities" && md=="simpleweak" && metric=="fluo"
    if noOff
        fun = @(x, p)subfun_simplebinding_weak_fluo_std2_nooff(x, p);
    else
        fun = @(x, p)subfun_simplebinding_weak_fluo_std2(x, p);
    end
elseif expmnt == "phases" && md=="simpleweak" && metric=="fraction"
    fun = @(x, p)subfun_simplebinding_weak_fraction_std2_phases(x, p);
elseif expmnt == "phases" && md=="simpleweak" && metric=="fluo"
    if noOff
        fun = @(x, p)subfun_simplebinding_weak_fluo_std2_phases_nooff(x, p);
    else
        fun = @(x, p)subfun_simplebinding_weak_fluo_std2_phases(x, p);
    end
elseif expmnt == "phaff" && md=="simpleweak" && metric=="fraction"
    fun = @(x, p)subfun_simplebinding_weak_fraction_std2_phaff(x, p);
elseif expmnt == "phaff" && md=="simpleweak" && metric=="fluo"
  fun = @(x, p)subfun_simplebinding_weak_fluo_std2_phaff_nooff(x, p);
end


end


%% simple weak fraction



function yfit = subfun_simplebinding_weak_fraction_std2(x, params)
%simplebinding in the weak promoter limit.

X = [];
n = 1;
X(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X(k, 1) = x(k);
    X(k, 2) = n;
end

dsid = X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec
omegaDP = params(1);
KD = params(2:max(dsid)+1)';
yfit = (((x./KD(dsid)).*omegaDP)./(1+ (x./KD(dsid))+ ((x./KD(dsid)).*omegaDP)));

end

%% simple weak fluo



function yfit = subfun_simplebinding_weak_fluo_std2(x, params)
%simplebinding in the weak promoter limit.

X = [];
n = 1;
X(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X(k, 1) = x(k);
    X(k, 2) = n;
end

dsid = X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec
omegaDP = params(1);
KD = params(2:max(dsid)+1)';
amp = params(max(dsid)+2);
offset = params(max(dsid)+3);
yfit = amp.*(((x./KD(dsid)).*omegaDP)./(1+ (x./KD(dsid))+ ((x./KD(dsid)).*omegaDP))) + offset;

end


function yfit = subfun_simplebinding_weak_fluo_std2_nooff(x, params)
%simplebinding in the weak promoter limit.

X = [];
n = 1;
X(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X(k, 1) = x(k);
    X(k, 2) = n;
end

dsid = X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec
omegaDP = params(1);
KD = params(2:max(dsid)+1)';
amp = params(max(dsid)+2);
yfit = amp.*(((x./KD(dsid)).*omegaDP)./(1+ (x./KD(dsid))+ ((x./KD(dsid)).*omegaDP)));

end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%phases


function yfit = subfun_simplebinding_weak_fraction_std2_phases(x, params)
%simplebinding in the weak promoter limit.

X = [];
n = 1;
X(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X(k, 1) = x(k);
    X(k, 2) = n;
end

dsid = X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec
KD = params(1);
omegaDP = params(2:max(dsid)+1)';
yfit = (((x./KD).*omegaDP(dsid))./(1+ (x./KD)+((x./KD).*omegaDP(dsid))));

end

%% simple weak fluo



function yfit = subfun_simplebinding_weak_fluo_std2_phases(x, params)
%simplebinding in the weak promoter limit.

X = [];
n = 1;
X(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X(k, 1) = x(k);
    X(k, 2) = n;
end

dsid = X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec
KD = params(1);
omegaDP = params(2:max(dsid)+1)';
amp = params(max(dsid)+2);
offset = params(max(dsid)+3);
yfit = amp.*(((x./KD).*omegaDP(dsid))./(1+ (x./KD)+ ((x./KD).*omegaDP(dsid)))) + offset;

end


function yfit = subfun_simplebinding_weak_fluo_std2_phases_nooff(x, params)
%simplebinding in the weak promoter limit.

X = [];
n = 1;
X(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X(k, 1) = x(k);
    X(k, 2) = n;
end

dsid = X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec
KD = params(1);
omegaDP = params(2:max(dsid)+1)';
amp = params(max(dsid)+2);
yfit = amp.*(((x./KD).*omegaDP(dsid))./(1+ (x./KD)+ ((x./KD).*omegaDP(dsid))));

end




function yfit = subfun_simplebinding_weak_fluo_std2_phaff_nooff(x, params)
%simplebinding in the weak promoter limit.
%note that this could probably be cleaner if dsid was a 2d matrix and
%params were referenced in that matrix like data set 1 has kd1, omega1,
%etc.

X = [];
n = 1;
X(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X(k, 1) = x(k);
    X(k, 2) = n;
end

nph = 3;
naff = 7;
dsid = X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec
omegaDP = [params(1)*ones(1, naff), params(2:nph)]';
KD = [params(1:naff), params(1)*ones(1, nph-1)]';
amp = params(nph+naff+1);
yfit = amp.*(((x./KD(dsid)).*omegaDP(dsid))./(1+ (x./KD(dsid))+ ((x./KD(dsid)).*omegaDP(dsid))));

end



function yfit = subfun_simplebinding_weak_fraction_std2_phaff(x, params)
%simplebinding in the weak promoter limit.

nph = 3;
naff = 7;

X = [];
n = 1;
X(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X(k, 1) = x(k);
    X(k, 2) = n;
end

dsid = X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec
omegaDP = [params(1)*ones(1, naff), params(2:nph)]';
KD = [params(1:naff), params(1)*ones(1, nph-1)]';
yfit = (((x./KD(dsid)).*omegaDP(dsid))./(1+ (x./KD(dsid))+ ((x./KD(dsid)).*omegaDP(dsid))));

end


function yfit = subfun_simplebinding_weak_fraction_std2_phases_fixKD(x, params)
%simplebinding in the weak promoter limit.

X = [];
n = 1;
X(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X(k, 1) = x(k);
    X(k, 2) = n;
end

dsid = X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec
KD = params(1);
omegaDP = params(2:max(dsid)+1)';
yfit = (((x./KD).*omegaDP(dsid))./(1+ (x./KD)+((x./KD).*omegaDP(dsid))));

end

function yfit = subfun_simplebinding_weak_fluo_std2_phases_nooff_fixKD(x, params)
%simplebinding in the weak promoter limit.

X = [];
n = 1;
X(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X(k, 1) = x(k);
    X(k, 2) = n;
end

dsid = X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec
KD = params(1);
omegaDP = params(2:max(dsid)+1)';
amp = params(max(dsid)+2);
yfit = amp.*(((x./KD).*omegaDP(dsid))./(1+ (x./KD)+ ((x./KD).*omegaDP(dsid))));

end