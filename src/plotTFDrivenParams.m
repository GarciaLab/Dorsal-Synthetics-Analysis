% load('C:\Users\Armando\Dropbox\DorsalSyntheticsDropbox\tfentryexit_paramsearch.mat')
function plotTFDrivenParams(factive, dt, mfpts, varargin)

nPoints = []; %if this is an integer, plotting will subsample using this many points. Otherwise, no subsampling
shouldRound = false; %if true, subsample by discretizing data points and rounding nearby values

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

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
fmin = .1;
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
    r_rounded = unique([round(x,1) round(y,1) round(z,1)], 'rows');
    x = r_rounded(:, 1);
    y = r_rounded(:, 2);
    z = r_rounded(:, 3);
end

%%
figure;
plot(hull)
hold on

in = inShape(hull,x, y, z);
scatter3(x(in),y(in),z(in),'r.')
scatter3(x(~in),y(~in), z(~in),'b.')

xlabel('factive')
ylabel('mean turn on (min)')
zlabel('delta t')
legend('viable region', 'viable parameters', 'unphysical parameters');

ax = gca;
% ax.Children(3).EdgeColor = 'none';
ax.Children(3).FaceAlpha = .5;

title('Parameter space for TF Driven model')

axis square;

fig = gcf;
fig.Renderer='Painters';

view(0, -90)

end