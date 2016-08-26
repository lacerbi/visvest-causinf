function [postright_c2,likec2,err,fevals] = VestBMS_c2corrpostandlikec2adapt(priorpdf2d,like_vis,like_vest,TolErr)
%VESTBMS_C2CORRPOSTANDLIKEC2QTRAPZ Multiple computations for C=2 (correlated) 
%
% ================ INPUT VARIABLES ====================
% PRIORPDF2D: 2d prior over s_vis, s_vest. [S-by-S] (double)
% LIKE_VIS: visual likelihood. [S-by-K] (double)
% LIKE_VEST: vestibular likelihood. [S-by-1-by-K] (double)
% TOLERR: relative error tolerance. [scalar] (double)
%
% ================ OUTPUT VARIABLES ==================
% POSTRIGHT_C2: p(right|x_vis,x_vest,C=2). [K-by-K] (double)
% LIKEC2: p(x_vis,x_vest|C=2). [K-by-K] (double)
% ERR: Simple estimate. [K-by-K] (double)
% FEVALS: Function evaluations. [K-by-K] (int)

K = size(like_vis,2);

tmp_vis(:,1,:) = like_vis;
tmp_vest(1,:,1,:) = squeeze(like_vest);
postpdf_c2(:,:,:,:) = bsxfun(@times,bsxfun(@times,priorpdf2d,tmp_vis),tmp_vest);

NGRID = 5;
M = NGRID^3;
col = size(postpdf_c2,1);
nx = size(postpdf_c2,1) / NGRID^3;
ny = size(postpdf_c2,2) / NGRID^3;

postright_c2 = zeros(K,K);
likec2 = zeros(K,K);
err = zeros(K,K);
fevals = zeros(K,K);

% Outer loop
for k = 1:K
    for m = 1:K
        Z = squeeze(postpdf_c2(:,:,k,m));
        
        [y1,err1,fevals1] = int2dadapt(Z(:,1:end/2),col,nx,ny/2,M,TolErr);
        [y2,err2,fevals2] = int2dadapt(Z(:,end/2+1:end),col,nx,ny/2,M,TolErr);
        
        err(k,m) = err1 + err2;
        fevals(k,m) = fevals1 + fevals2;
        
        likec2(k,m) = y1 + y2 + realmin;
        postright_c2(k,m) = (y2 + realmin)./(y1 + y2 + realmin);
        
    end
end


end


function [y,err,fevals] = int2dadapt(Z,col,nx,ny,M,TolErr)

idx = 1 + ((M+1)/2-1)*(col+1) + M*bsxfun(@plus, (0:ny-1)*col, (0:nx-1)'); idx = idx(:);
err = sum(Z(idx))*M^2;
tol = TolErr*err;   % Maximum error

[y,fevals] = evalint(Z,1,M,nx,ny,col,tol);

end


function [z,fevals] = evalint(Z,idx0,M,nx,ny,col,tol)

NGRID = 5;

idx = idx0 + ((M+1)/2-1)*(col+1) + M*bsxfun(@plus, (0:ny-1)*col, (0:nx-1)');
idx = idx(:);

y = Z(idx)*M^2;
fevals = ones(1,numel(y));

% Estimate of the error
if M > 1
    dy = max(y) - min(y);
    if dy > tol
        Snew = M/NGRID;
        for j = 1:ny
            for i = 1:nx
                idx0sub = idx0 + (i-1)*M + (j-1)*col*M;
                [y(i+(j-1)*nx),fevals(i+(j-1)*nx)] = evalint(Z,idx0sub,Snew,NGRID,NGRID,col,tol);
            end
        end
    end
    
end

z = sum(y);
fevals = sum(fevals);

end