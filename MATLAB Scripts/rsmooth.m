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

[n1,n2] = size(y); n = n1*n2;
N = sum([n1,n2]~=1);
Lambda = bsxfun(@plus,repmat(-2+2*cos((0:n2-1)*pi/n2),n1,1),...
    -2+2*cos((0:n1-1).'*pi/n1));
W = ones(n1,n2);
zz = y;
s = [];

for k = 1:6
    tol = Inf;
    while tol>1e-5
        DCTy = dct2(W.*(y-zz)+zz);
        fminbnd(@GCVscore,-15,38);
        tol = norm(zz(:)-z(:))/norm(z(:));
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
        GCVs = RSS/n/(1-TrH/n)^2;
    end
end

function W = bisquare(r,h)
    MAD = median(abs(r(:)-median(r(:))));
    u = abs(r/(1.4826*MAD)/sqrt(1-h));
    W = (1-(u/4.685).^2).^2.*((u/4.685)<1);
end
