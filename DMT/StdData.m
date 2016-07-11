function Xstd = StdData(X,l,s)
% *************************************************************************
% 7.12.2006 - talonen
% Data preparation - Liukuvan virherajan määritys
% X = käsiteltävä data
% l = aikaikkunan pituus
% s = ennustavuus (=1 normaali, 0=ei siirtoa, aina rajojen sisällä)
% *************************************************************************

    % liukuvakeskiarvo
    [n,m] = size(X);
    Xstd2 = zeros(n-l+1,m);
    for i=2+s:l+s
        Xstd2(i,:) = std(X(1:i-s,:));
        %Xsmooth2 = Xsmooth2 + X(i:n-l+i,:);
        %Xsmooth = Xsmooth + X(i:n-smoothness+i,:);
    end
    Xstd = [zeros(s,m); Xstd2; zeros(l-1-s,m)];
    %Xsmooth =  Xsmooth2/smoothness;
    
 
        for i=l+s+1:n
        
            Xstd(i,:)=std(X(i-l-s:i-s,:));
        
        end

    
% 14.9. HUOM! esim. StdData(a,1,0): antaa a vektorista std tämän ja
% edellisen askeleen MSDV:nn

%if s>0
 %   Xstd = [zeros(s,m); Xstd(1:n-s,m)];
%end
    
end
% *************************************************************************