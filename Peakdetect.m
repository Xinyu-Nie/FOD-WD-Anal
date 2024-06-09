function [DirVol,FracVol,peaks,peakres] = Peakdetect(y,yref,THD,NumPeaks,coord,VertNbrCellArr)
% Contribution
%  Author: Xinyu Nie
%  Created: 2024/6/1
%  Copyright:The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu
%  FOD peak lobes detection
DirVol = zeros(3,NumPeaks);
FracVol = zeros(1,NumPeaks);
peaks=zeros(1,NumPeaks);
peakres=peaks;
                flag = zeros(length(y),1);
                ind = find(yref>THD);
                N = length(ind);
                for m = 1:N
                    if y(ind(m))>max(y(VertNbrCellArr{ind(m)})) %peak detected
                        flag(ind(m)) = yref(ind(m));
                    end
                end
                [a,b] = sort(flag);
                %evenly sample the arrays since FODs are symmetric
                for m=1:NumPeaks
                    M=2*m;
                    if M<=N
                    FracVol(m) = yref(b(end-2*m+2));
                    peak_index = b(end-2*m+2);
                    DirVol(:,m) = coord(peak_index,:)';
                    peaks(m)=peak_index;
                    peakres(m)=b(end-2*m+1);
                    end
                end
        

end