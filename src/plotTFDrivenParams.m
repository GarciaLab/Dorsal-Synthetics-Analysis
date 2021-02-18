% load('C:\Users\Armando\Dropbox\DorsalSyntheticsDropbox\tfentryexit_paramsearch.mat')
function goodMatrixIndices = plotTFDrivenParams(factive, dt, mfpts, varargin)

nPoints = []; %if this is an integer, plotting will subsample using this many points. Otherwise, no subsampling
shouldRound = false; %if true, subsample by discretizing data points and rounding nearby values
fig = [];
goodMatrixIndices = [];
dim = 3;

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

mfpts0 = mfpts;
dt0 = dt;
factive0 = factive;

factive(isnan(mfpts)) = [];
dt(isnan(mfpts)) = [];
mfpts(isnan(mfpts)) = [];

factive(isnan(dt)) = [];
mfpts(isnan(dt)) = [];
dt(isnan(dt)) = [];

x = factive(:);
y = mfpts(:);
z = dt(:);

%% make  convex hull out of the allowable region
% fmin = .1;
fmin = 0;
fmax = 1;
tmin = 3.5;
tmax = 7;
deltatmin = 0;
deltatmax = 3;

box = makeBox([fmin, tmin, deltatmin], [fmax, tmax, deltatmax]);
x_hull = box(:, 1);
y_hull = box(:, 2);
z_hull = box(:, 3);


hull = alphaShape(x_hull, y_hull, z_hull,Inf,'HoleThreshold',1E30 );
%%

%% Subsample and/or round
if ~isempty(nPoints)
    subSample = randsample(length(x),nPoints);
    x = x(subSample);
    y = y(subSample);
    z = z(subSample);
end

if shouldRound
    r_rounded = unique([round(x,2) round(y,1) round(z,-1)], 'rows'); %factive, mfpts, dt
    x = r_rounded(:, 1);
    y = r_rounded(:, 2);
    z = r_rounded(:, 3);
end

in = inShape(hull,x, y, z);

%%
if ~isempty(fig)
    gcf;
end


if nargout == 0
    
    if dim == 3
        plot(hull)
        hold on
        
        scatter3(x(in),y(in),z(in),'r.')
        scatter3(x(~in),y(~in), z(~in),'b.')
        %
        % xlabel('factive')
        % ylabel('mean turn on (min)')
        % zlabel('delta t')
        % legend('viable region', 'viable parameters', 'unphysical parameters');
        
        ax = gca;
        % ax.Children(3).EdgeColor = 'none';
        ax.Children(3).FaceAlpha = .5;
        
        % title('Parameter space for TF Driven model')
        
        axis square;
        
        xlim([0, 1]);
        ylim([0, 10]);
        fig = gcf;
        fig.Renderer='Painters';
        
        view(0, -90)
        
    elseif dim == 2
        
        colormap(brewermap(20,'Blues'));
        
        [dat_y, dat_x] = plotGreenBoxWithData;
        
        nBins = [15, 5];
        h = binscatter(dat_x,dat_y, nBins);
%         h.ShowEmptyBins = 'on';
        colormap(brewermap(20,'Blues'));
        xlim([0, 1])
        

        
        hold on
        
        scatter(x(in),y(in),'r.')
        
        hold on
        
        scatter(x(~in),y(~in),'b.')
        
%         set (gca,'Ydir','reverse')
        
   
        
    end
    
    %%% Let's return the good parameters
    if ~shouldRound && isempty(nPoints)
        
        factive0(isnan(mfpts0)) = -1;
        dt0(isnan(mfpts0)) = 100;
        mfpts0(isnan(mfpts0)) = 100;
        
        factive0(isnan(dt0)) = -1;
        mfpts0(isnan(dt0)) = 100;
        dt0(isnan(dt0)) = 100;
        
        in2 = inShape(hull,factive0(:), mfpts0(:), dt0(:));
        goodLinearIndices = find(in2);
        [in1, in2, in3, in4, in5] = ind2sub(size(mfpts0), goodLinearIndices);
        goodMatrixIndices = [in1 in2 in3 in4 in5];
        
    end
    
end