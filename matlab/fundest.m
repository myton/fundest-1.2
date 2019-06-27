%FUNDEST Estimate the Fundamental matrix of two images from corresponding points.
%
% [F, idxOutliers]=fundest(pts0, pts1, inlPcent, normalize, howto, verbose);
%        pts0, pts1 are the input matching point pairs
%        inlPcent is the expected fraction of inliers in the input points
%        normalize specifies whether Hartley's normalization should be applied
%             to input points prior to DLT estimation (default)
%        howto controls how estimation is performed; can be one of:
%             '8pt', 'algmin', 'epip_dist', 'sampson_err'or empty (implying '8pt', default)
%             A synonym for '8pt' is 'norefine'
%        verbose is optional
%
% If an initial, approximate estimate Fini of the fund. matrix is available,
% fundest can be also invoked as
%
% [F, idxOutliers]=fundest(pts0, pts1, inlThresh, Fini, howto, verbose);
%        inlThresh is a threshold used to determine inliers; see fundest_wie() for details.
%
% Returns the 3x3 fundamental matrix F and the indices of outlying point pairs from pts0, pts1.
%
% Notes::
%   MEX file.
%

if ~exist('fundest', 'file')
  error('you need to build the MEX version of fundest, see mex');
end

