function b=databinaariseksi(hs,m,p)
% *************************************************************************
% DATA BINAARISEKSI
% *************************************************************************
n = length(hs.Ika);

clear apu;
k=0;
% *************************************************************************
% 31 questions
kymyksyksiaotetaan = length(p.mitkakyssaritmukana); %31

% *************************************************************************
for iA=1:kymyksyksiaotetaan;%31
    
    i = p.mitkakyssaritmukana(iA);
    
    % ************************** multiple choice
    if i==21 || i==31
        clear fb
        clear maa
        apu = (hs.k.Vastaus(:,i));
        k4=0;
        for i4=1:length(apu)
            pilkut = findstr(apu{i4},',');
            if length(pilkut) == 0
                k4=k4+1;
                maa{k4} = apu{i4};
            else
                if length(pilkut)==1
                    k4=k4+1;
                    maa{k4} = apu{i4}(1:pilkut-1);
                    k4=k4+1;
                    maa{k4} = apu{i4}(pilkut+1:end);
                else
                    k4=k4+1;
                    maa{k4} = apu{i4}(1:pilkut(1)-1);
                    for j=1:length(pilkut)-1
                        k4=k4+1;
                        maa{k4} = apu{i4}(pilkut(j)+1:pilkut(j+1)-1);
                    end
                    k4=k4+1;
                    maa{k4} = apu{i4}(pilkut(end)+1:end);
                end
            end
    
    
    
        end

        fb.ehdokkaat = unique(maa);

        disp(' ');
        disp(num2str(i))
        disp(hs.k.Kysymys{i});
        
        for j=1:length(fb.ehdokkaat)
            disp(fb.ehdokkaat{j});
            if strcmp(fb.ehdokkaat{j},'-')
                disp('pois');
            else
                k=k+1;
                b.vastaus{k} = fb.ehdokkaat{j};
                b.kysymysnumero(k) = i;
            end
        end
        
    % ************************** single choice    
    else
        disp(' ');
        disp(num2str(i))
        disp(hs.k.Kysymys{i});
        apu = unique(hs.k.Vastaus(:,i));
        for j=1:length(apu)
            disp(apu{j});
            if strcmp(apu{j},'-')
                disp('pois');
            else
                k=k+1;
                b.vastaus{k} = apu{j};
                b.kysymysnumero(k) = i;
            end
        end

    end
end

k=0;
k2=0;
k3=0;
b.binaaritaulu = zeros(n,length(b.vastaus));
b.binaarimassa = zeros(n,length(b.vastaus));
b.vastaustaulu = zeros(n,length(b.vastaus));
for iA=1:kymyksyksiaotetaan
    
    
    i = p.mitkakyssaritmukana(iA);
    
    
    if i==21 || i==31
        disp(' ');
        disp(num2str(i));
        disp(hs.k.Kysymys{i});
        
        apu2 = b.vastaus(b.kysymysnumero==i);
        kysymysvaihtoehtoja = length(apu2);
       
        
        %clear fb
        
        apu = (hs.k.Vastaus(:,i));
        
        for i4=1:length(apu)
            k4=0;
            clear maa
            pilkut = findstr(apu{i4},',');
            if length(pilkut) == 0
                k4=k4+1;
                maa{k4} = apu{i4};
            else
                if length(pilkut)==1
                    k4=k4+1;
                    maa{k4} = apu{i4}(1:pilkut-1);
                    k4=k4+1;
                    maa{k4} = apu{i4}(pilkut+1:end);
                else
                    k4=k4+1;
                    maa{k4} = apu{i4}(1:pilkut(1)-1);
                    for j=1:length(pilkut)-1
                        k4=k4+1;
                        maa{k4} = apu{i4}(pilkut(j)+1:pilkut(j+1)-1);
                    end
                    k4=k4+1;
                    maa{k4} = apu{i4}(pilkut(end)+1:end);
                end
            end
            for jjj = 1:length(maa)
                missa = find(strcmp(apu2,maa{jjj}));
    
            
                if length(missa)==0
                    %disp('pois');
                else
                    %k
                    kkkk=missa+k;
                    b.binaarimassa(i4,kkkk) = hs.k.Merkitys(i4,i); % kyssariin siis ainasama
                    b.binaaritaulu(i4,kkkk) = 1;%./length(b.v{i}.v);%./length(b.v{i}.v);%/length(b.v{i}.v);
                    b.vastaustaulu(i4,kkkk) = b.binaaritaulu(i4,kkkk);
                    
                    % TÄNNE 110511 a b painotus
                    %b.binaaritaulu = b.binaaritaulu.*(1+(b.binaarimassa==1)*p.kysymyksenpainotusB - (b.binaarimassa==-1)*p.kysymyksenpainotusA);
                    %110511 - ennen kysymysmäärällä painottamista
                    b.binaaritaulu(i4,kkkk) = b.binaaritaulu(i4,kkkk) + b.binaaritaulu(i4,kkkk)*((b.binaarimassa(i4,kkkk)==1)*p.kysymyksenpainotusB - (b.binaarimassa(i4,kkkk)==-1)*p.kysymyksenpainotusA);
                    
                    % tänne kans!
                    % ***************** kyllä ei paino suurempi kuin monta eri?
                    if p.vastausmaarallapainottaminen == 1
                        % '-' otettu pois: kyllä ei = 0.5... 
                        % neljävaihtoehtoa = 0.25, jne..
                        if p.vastausmaaraneliojuuressa == 0
                            b.binaaritaulu(i4,kkkk) = b.binaaritaulu(i4,kkkk)/((kysymysvaihtoehtoja));%length(b.v{i}.v);
                        else
                            b.binaaritaulu(i4,kkkk) = b.binaaritaulu(i4,kkkk)/sqrt(kysymysvaihtoehtoja); %-1 poistettu silla on tossa otettu jo
                        end
                    end
                    
                    
                    
                end
            end
            
        end
        
        % ***************** pakotettu kompo, tarkastellaan ekaa
                    if p.kompo2maaritys == 1
                        if sum(p.X1==i)
                            for kkkk2 = 1:length(b.vastaus(b.kysymysnumero==i))
                                k2= k2+1;
                                b.Xalue(k2) = k+kkkk2;
                            end
                        end
                        if sum(p.X2==i)
                            for kkkk2 = 1:length(b.vastaus(b.kysymysnumero==i))
                                k3= k3+1;
                                b.Xalue2(k3) = k+kkkk2;
                            end
                        end
                    end
        
        k=k+kysymysvaihtoehtoja;
        
        
    else
        % ************************** single choice 
        disp(' ');
        disp(num2str(i));
        disp(hs.k.Kysymys{i});
        apu = unique(hs.k.Vastaus(:,i));
        kysymysvaihtoehtoja = length(apu)-1;
        for j=1:length(apu)
            %disp(apu{j});
        
            if strcmp(apu{j},'-')
                %disp('pois');
            else
                disp(apu{j});
                k=k+1;
                b.binaarimassa(:,k) = hs.k.Merkitys(:,i); % kyssariin siis ainasama
                b.binaaritaulu(:,k) = strcmp(hs.k.Vastaus(:,i),b.vastaus{k});%./length(b.v{i}.v);%./length(b.v{i}.v);%/length(b.v{i}.v);
                b.vastaustaulu(:,k) = b.binaaritaulu(:,k);
                
                 % TÄNNE 110511 a b painotus
                 b.binaaritaulu(:,k) = b.binaaritaulu(:,k).*(1+(b.binaarimassa(:,k)==1)*p.kysymyksenpainotusB - (b.binaarimassa(:,k)==-1)*p.kysymyksenpainotusA);
                 %110511 - ennen kysymysmäärällä painottamista
                    
                
                % ***************** kyllä ei paino suurempi kuin monta eri?
                if p.vastausmaarallapainottaminen == 1
                    % '-' otettu pois: kyllä ei = 0.5... 
                    % neljävaihtoehtoa = 0.25, jne..
                    if p.vastausmaaraneliojuuressa == 0
                        b.binaaritaulu(:,k) = b.binaaritaulu(:,k)/(length(apu)-1);%length(b.v{i}.v);
                    else
                        b.binaaritaulu(:,k) = b.binaaritaulu(:,k)/sqrt(length(apu)-1);  % -1 koska apussa myös vastaus '-'
                    end
                end
                % ***************** pakotettu kompo, tarkastellaan ekaa
                if p.kompo2maaritys == 1
                    if sum(p.X1==i)
                        k2= k2+1;
                        b.Xalue(k2) = k;
                    end
                    if sum(p.X2==i)
                        k3= k3+1;
                        b.Xalue2(k3) = k;
                    end
                end
            end
        end
    end    
end
% plot(sum(b.binaaritaulu(:,:)'))  %tarkistus, kaikki vastannut ainakin-
%2  16 20 22 24 25
%%
% TARKISTUSKUVAT % ihan ok.. 
tarkistus = 0;
if tarkistus == 1
figure(10);
bar(sum(b.binaaritaulu(:,:)))
for i=1:length(b.vastaus)
    text(i,sum(b.binaaritaulu(:,i)),b.vastaus{i});
end

figure(16);
kyss = 2;
kyssareita = sum(b.kysymysnumero==kyss);
barh(sum(b.binaaritaulu(:,find(b.kysymysnumero==kyss))))
colormap(cool)
%axis tight


k=0;
for i=find(b.kysymysnumero==kyss)
    k=k+1;
    apusum = sum(b.binaaritaulu(:,i));
    text(apusum(1)+3,k-0.1,b.vastaus{i},'verticalalignment','middle');
end
axis([0.5 2200 0 kyssareita+1]);
title(hs.k.Kysymys{kyss});
end

% NANit massasta pois IDA test
for i = 1:2306
    for j=1:size(b.binaarimassa,2)
        % 110509
        if isnan(b.binaarimassa(i,j))
            b.binaarimassa(i,j) = 0;
        end
    end
end
            



end