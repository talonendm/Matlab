function Xsmooth = SmoothData(X,l,s)
% *************************************************************************
% 10.10.2006 - talonen
% Data preparation - Kohinan suodatus
% X = k‰sitelt‰v‰ data
% l = tasoituksen pituus
% s = tyyli [1=liukuva ka]
% *************************************************************************
if s==1
    % liukuvakeskiarvo
    smoothness = l;
    [n,m] = size(X);
    Xsmooth2 = zeros(n-smoothness+1,m);
    for i=1:smoothness
        Xsmooth2 = Xsmooth2 + X(i:n-smoothness+i,:);
        %Xsmooth = Xsmooth + X(i:n-smoothness+i,:);
    end
    Xsmooth = [zeros(smoothness-1,m); Xsmooth2/smoothness];
    %Xsmooth =  Xsmooth2/smoothness;
    for i=1:smoothness-1
        for j=1:i
            Xsmooth(i,:)=Xsmooth(i,:)+X(j,:);
        end
        Xsmooth(i,:)=Xsmooth(i,:)/i;
    end
end
if s==2
    % ennustaa tulevaan.. joten ei kovinkaan loogista k‰ytt‰‰ t‰t‰?!
    b = ones(1,l)/l;
    Xsmooth = filtfilt(b,1,X);
    % voisi olla painotettu uusimpaan dataan.
end
% *************************************************************************