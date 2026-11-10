function z = rsmooth(y)
% RSmooth - Robust smoothing of gridded data in one and higher dimensions.
%
%   This function implements the algorithm described by:
%   Garcia D. Robust smoothing of gridded data in one and higher dimensions with missing values.
%   Comput Stat Data Anal. 2010 Apr 1;54(4):1167-1178. doi: 10.1016/j.csda.2009.09.020.
%   PMID: 24795488; PMCID: PMC4008475.
%
%   y - Input data array to be smoothed.
%   z - Smoothed output array.
%
%   Note: Please ensure the proper acknowledgment of the original author
%         and source when using this code.
% -----------------------------------------------------------------------------

if isempty(y)
    z = y;
    return
end
if ~ismatrix(y) || any(~isfinite(y(:)))
    error('Input must be a finite vector or matrix; missing data are unsupported.')
end
if numel(y) < 3 || max(y(:))-min(y(:)) <= 100*eps*max(1,max(abs(y(:))))
    z = y;
    return
end
[n1,n2] = size(y); n = n1*n2;
N = sum([n1,n2]~=1);
Lambda = bsxfun(@plus,repmat(-2+2*cos((0:n2-1)*pi/n2),n1,1),...
    -2+2*cos((0:n1-1).'*pi/n1));
W = ones(n1,n2);
zz = y;
s = [];

for k = 1:6
    tol = Inf;
    iter = 0;
    while tol>1e-5
        iter = iter+1;
        if iter > 200
            error('Smoothing iteration failed to converge.')
        end
        DCTy = dct2(W.*(y-zz)+zz);
        p = fminbnd(@GCVscore,-15,38);
        GCVscore(p);
        tol = norm(zz(:)-z(:))/max(norm(z(:)),eps);
        zz = z;
    end
    tmp = sqrt(1+16*s);
    h = (sqrt(1+tmp)/sqrt(2)/tmp)^N;
    W = bisquare(y-z,h);
end

    function GCVs = GCVscore(p)
        s = 10^p;
        Gamma = 1./(1+s*Lambda.^2);
        z = idct2(Gamma.*DCTy);
        RSS = norm(sqrt(W(:)).*(y(:)-z(:)))^2;
        TrH = sum(Gamma(:));
        GCVs = RSS/n/max(1-TrH/n,eps)^2;
    end
end

function W = bisquare(r,h)
    MAD = median(abs(r(:)-median(r(:))));
    if MAD <= eps*max(1,max(abs(r(:))))
        W = ones(size(r));
        return
    end
    u = abs(r/(1.4826*MAD)/sqrt(max(1-h,eps)));
    W = (1-(u/4.685).^2).^2.*((u/4.685)<1);
end