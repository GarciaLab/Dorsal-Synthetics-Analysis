% load('C:\Users\Armando\Dropbox\DorsalSyntheticsDropbox\tfentryexit_paramsearch.mat')
function plotTFDrivenParams(factive, dt, mfpts)

    factive(isnan(mfpts)) = [];
    dt(isnan(mfpts)) = [];
    mfpts(isnan(mfpts)) = [];

    factive(isnan(dt)) = [];
    mfpts(isnan(dt)) = [];
    dt(isnan(dt)) = [];

    fmin = .1;
    fmax = 1;
    tmin = 3.5;
    tmax = 7;
    deltatmin = 0;
    deltatmax = 3;

    box = makeBox([fmin, tmin, deltatmin], [fmax, tmax, deltatmax]);
    qx = box(:, 1);
    qy = box(:, 2);
    qz = box(:, 3);

    qx2 = factive(:);
    qy2 = mfpts(:);
    qz2 = dt(:);
    shp = alphaShape(qx, qy, qz,Inf,'HoleThreshold',1E30 );

    subsample = randsample(length(qx2),10);

    qx2down = qx2(subsample);
    qy2down = qy2(subsample);
    qz2down = qz2(subsample);

    q = [round(qx2,1) round(qy2,1) round(qz2,1)];
    qu = unique(q,'rows');
    qx2u = qu(:, 1);
    qy2u = qu(:, 2);
    qz2u = qu(:, 3);
    figure;
    plot(shp)
    hold on

    qx2u = qx2;
    qy2u = qy2;
    qz2u = qz2;

    % in = inShape(shp,qx2u, qy2u, qz2u);
    % scatter3(qx2u(in),qy2u(in),qz2u(in),'r.')
    % scatter3(qx2u(~in),qy2u(~in), qz2u(~in),'b.')

    in = inShape(shp,qx2down, qy2down, qz2down);
    scatter3(qx2down(in),qy2down(in),qz2down(in),'r.')
    scatter3(qx2down(~in),qy2down(~in), qz2down(~in),'b.')

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

end