close all;
clear;
% clc

addpath('./RiemannianTools');

%%
[Data, Covs, vF, vC] = MakeData();
ax1 = PlotCovs(Covs, vF, vC);

%%
CovsPT  = PT(Covs);
ax2     = PlotCovs(CovsPT, vF, vC);
linkprop([ax1, ax2], {'CameraPosition','CameraUpVector'});

%======End=Sciprt=======================================================================================%

%%
function [Data, Covs, vF, vC] = MakeData()
    T = 500;
    N = 200;
    t = linspace(0, 1, T + 1); t(end) = [];
    
    f1   = 10; 
    A1 = [-0.5845, 2.0633;
          -0.0741, 0.6288];
    
    U = [-1, 0;
         0,  1];

    A2 = 1.5 * U * A1;

    
    vF        = [];
    A         = {A1, A2};
    Covs{N,2} = [];
    Data{N,2} = [];
    vC        = [1 * ones(N, 1);
                 2 * ones(N, 1)];
    for ss = 1 : 2
        for ii = 1 : N
            p0        = pi/2 * (rand(1) - 1);
            mX        = [sin(2 * pi * f1 * t);
                         cos(2 * pi * f1 * t + p0)];
            mX        = mX + randn(size(mX)) / 10;
            vF(end+1) = p0;
            mX        = A{ss} * mX;
            
            Covs{ii,ss} = cov(mX');
            Data{ii,ss} = mX;
        end
    end
end

%%
function ax = PlotCovs(Covs, vF, vC)
        
    mCov      = cat(3, Covs{:});
    mCov      = reshape(mCov, 4, []);
    
    K = max(mCov(:));
    vK = linspace(0, K, 100);
    figure; PlotPositiveMatrix(K); hold on; grid on; 
    scatter3(mCov(1,:), mCov(2,:), mCov(4,:),   100, vF,  'Fill');
    scatter3(vK,        0*vK,      vK, 10, 'r',           'Fill');
    plot3   (1,         0,         1,  '.b', 'MarkerSize', 20);
    ax(1) = gca;
    
    hold on; set(gca, 'FontSize', 16);
    h = colorbar;
    title(h, '$\phi$', 'Interpreter', 'Latex');

    %--
    figure; PlotPositiveMatrix(K); hold on; grid on; 
    scatter3(mCov(1,:), mCov(2,:), mCov(4,:),   100, vC,  'Fill');
    scatter3(vK,        0*vK,      vK, 10, 'r',           'Fill');
    plot3   (1,         0,         1,  '.b', 'MarkerSize', 20);
    ax(2) = gca;
    
    hold on; set(gca, 'FontSize', 16);
    h1 = plot3(nan, nan, nan, 'y.', 'MarkerSize', 20);
    h2 = plot3(nan, nan, nan, 'b.', 'MarkerSize', 20);
    legend([h1, h2], 'M_1', 'M_2');
end

%%
function CovsPT = PT(Covs)

M{1} = RiemannianMean(cat(3, Covs{:,1}));
M{2} = RiemannianMean(cat(3, Covs{:,2}));
D    = RiemannianMean(cat(3, M{1}, M{2}));
% D    = RiemannianMean(cat(3, Covs{:}));

CovsPT = Covs;
for ss = 1 : 2
    E      = (D * M{ss}^(-1))^(1/2);
    for ii = 1 : size(Covs, 1)
        CovsPT{ii,ss} = E * Covs{ii,ss} * E';
    end
end

end

%%
function CovsMT = MeanTransport(Covs)

M{1} = RiemannianMean(cat(3, Covs{:,1}));
M{2} = RiemannianMean(cat(3, Covs{:,2}));
D    = RiemannianMean(cat(3, M{1}, M{2}));

CovsMT = Covs;
for ss = 1 : 2
    for ii = 1 : size(Covs, 1)
        CovsMT{ii,ss} = ExpMap(D, LogMap(M{ss}, Covs{ii,ss}));
    end
end

end
