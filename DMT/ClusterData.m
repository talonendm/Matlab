function [idx, cc, sumd] = ClusterData(klusteri)
% *************************************************************************
% 20.12.2006 - talonen (& Parviainen)

P2 = double(klusteri.P);
X  = P2 .* repmat(klusteri.L', size(P2,1), 1);


if klusteri.maara == 0
    
%alustukset
% iteroi .. t‰m‰n j‰lkeen k‰ytt‰j‰n helpompi valita klusterien m‰‰r‰

koko = 7;
if length(P2)<koko 
    koko=length(P2)-1
end
    
AIDX = cell(koko,1);
ACC  = cell(koko,1);
ASC = cell(koko,1);
for k=2:koko+1
    k
    [idx, cc, sumd] = kmeans(X, k, 'display','iter');
    AIDX{k-1} = idx;
    ACC{k-1} = cc;
    ASC{k-1} = sumd;
    pause;
end;
%[[1:length(klusteri.P)]' AIDX{3}]
%[[1:59]' AIDX{3}]
%find(AIDX{3}==1)
%find(AIDX{3}==2)
%find(AIDX{3}==3)
%find(AIDX{3}==4)


else
    [idx, cc, sumd] = kmeans(X, klusteri.maara, 'display','iter');
    %idx
    %idx = muuttujat(idx)
end


% *************************************************************************