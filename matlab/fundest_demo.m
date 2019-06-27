matchesfile='../test/matches.txt';

% read in matches
[x1, y1, x2, y2]=textread(matchesfile, '%f%f%f%f', 'commentstyle', 'shell');
pts0(:, :)=[x1'; y1']';
pts1(:, :)=[x2'; y2']';
%pts0, pts1

%[F, idxOutliers]=fundest(pts0, pts1, 0.7, 1, 'epip_dist'); % minimize symmetric distances to epipolar lines
%[F, idxOutliers]=fundest(pts0, pts1, 0.7, 1, 'sampson_err', 2); % verbose
[F, idxOutliers]=fundest(pts0, pts1, 0.7, 1, 'sampson_err');
nbOutliers=max(size(idxOutliers));

format long g
F, nbOutliers

Fini=[
-1.e-07     -2.e-06       0.002
 5.e-06      3.e-07      -0.024
-0.003       0.022        0.975
];
[Fi, idxOutliers]=fundest(pts0, pts1, -0.8, Fini, 'sampson_err'); % 80% outliers; minus sign necessary!
Fi
