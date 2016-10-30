function [gmm,isactive] = merging(gmm,isactive,cst)

% MERGING merges components of a Gaussian mixture that are similar w.r.t. a
%         certain distance (either Hellinger or Mahalanobis).
%
% IN        gmm         Gaussian mixture
%           cst         collection of constants
%
% OUT       gmm         merged mixture
%
% AUTHOR    Isabel Schlangen, (c) 2016

active = find(isactive);     % collect all components
for ii=1:length(active)
    if gmm(active(ii)).i %i.e. if the component is still active
        indcollect = active(ii);
        for jj=ii+1:length(active)
            if gmm(active(jj)).i
                diff = gmm(active(ii)).m - gmm(active(jj)).m;
                switch cst.dist
                    case 1 % Mahalanobis (that does not work with the fat Gaussian!)
                        dist = diff'/gmm(active(ii)).C*diff;
                    case 2 % Hellinger
                        Ctmp = 0.5*(gmm(active(ii)).C + gmm(active(jj)).C);
                        dist = 1-sqrt(sqrt(det(gmm(active(ii)).C*gmm(active(jj)).C))...
                            /det(Ctmp)) * exp(-0.125 * diff'/Ctmp*diff);
                    case 3 % Bhattacharyya
                        Ctmp = 0.5*(gmm(active(ii)).C + gmm(active(jj)).C);
                        dist = diff'/Ctmp*diff*0.25 + 0.5*log(det(Ctmp)...
                            / sqrt(det(gmm(active(ii)).C)+det(gmm(active(ii)))));
                    otherwise
                        disp('invalid distance measure.')
                        keyboard;
                end
                if dist<cst.tmerge
                    indcollect = [indcollect, active(jj)];
                    gmm(active(jj)).i = 0; %deactivate this component
                    isactive(active(jj)) = 0;
                end
            end
        end
        if length(indcollect) > 1
            gmmnew.w = sum([gmm(indcollect).w]);
            gmmnew.m = sum(repmat([gmm(indcollect).w],cst.Nx,1).*[gmm(indcollect).m],2)/gmmnew.w;
            gmmnew.C = zeros(cst.Nx,cst.Nx);
            for kk=1:length(indcollect)
                gmmnew.C = gmmnew.C+gmm(indcollect(kk)).w*(gmm(indcollect(kk)).C...
                    +(gmmnew.m-gmm(indcollect(kk)).m)*(gmmnew.m-gmm(indcollect(kk)).m)');
            end
            gmmnew.i = 1;
            gmm(active(ii)) = gmmnew;
            isactive(active(ii)) = 1;
        end
    end
end