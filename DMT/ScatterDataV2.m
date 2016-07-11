function correlations = ScatterDataV2(tarkastelu,parivertailu)
% *************************************************************************
% 10.10.2006 - talonen
% Data preparation - differentaalien tarkastelu
% X = käsiteltävä data
% dX = käsiteltävä data, differenssi
% t = datan aikavektori
% mn = muuttujien nimet
% m = piirrettävien/tarkasteltavien muuttujien lkm
% s = tyyli [1=perus, 2=differentiaalierilainen]
% *************************************************************************
% SIISTITTÄVÄ!
% *************************************************************************
t = tarkastelu.t;
mn = tarkastelu.nimet;
if length(parivertailu.variables)==0
    tarkasteltavat = tarkastelu.variables;
else
    tarkasteltavat = parivertailu.variables;
end

m = parivertailu.piirtotiheys;
m2 = parivertailu.varinPiirtotiheys;
s = parivertailu.tyyli;
%X = tarkastelu.X;
%dX = tarkastelu.dX;
yk = tarkastelu.yksikko;

figure(9);
% 25.9.2006
% Toimii - scatterit riippuvuuksien visuaaliseen tarkasteluun
% lävistäjään ajansuhteessa
n = length(tarkastelu.eaX);
muuttujia = length(tarkasteltavat);

if muuttujia>7 
    muuttujia = 7
    disp('error: Too many variables for ScatterPlot');
end

% 31.10.2006
X = tarkastelu.efX;   % tarkastellaanko A/B?!
dX = tarkastelu.defX;   % tarkastellaanko A/B?!

%if tarkastelu.esikasittely == 1
%    cd(tyokaluhakemisto);
%    [A barA] = regrCenter(A);
%    [A backmap] = regrScale(A);
    
%    [B barB] = regrCenter(B);
%    [B backmapB] = regrScale(B);
%end
%X = A;
%dX = B;
%cd(tyoskentelyhakemisto);
%m = 100
%m2 = 50   %väriapupisteet, kuinka usein (vaalee eka, tumma lopus)  28.9.2006
%tarkasteltavat = tarkasteltavat'
%muuttujia = 6
% mitä muuttujia haluataan tarkastella

% ******************************************************
% lisäys 28.9.2006 - saadaan valittua halutut variables
% myöhemmin voisi pyrkiä luokittelemaan prosessista datasettejä, jotka
% liittyvät toisiinsa tai on yleensä mielekästä tutkailla riippuvuuksia
%tarkasteltavat = zeros(muuttujia,1);
  
% ******************************************************
%plot(X(1:m:n,1),X(1:m:n,2),'.')
%

%k=0;


subplot(muuttujia,muuttujia,1);
% *************************************************************************
% Luodaan alakolmio (parivertailu) sekä lävistäjälle kuvaajat ajansuhteen
% *************************************************************************
%plot(t(1:m:n),X(1:m:n,1),'color',[0 0 0]);
% ajanfunktiona tärkeää piirtää kaikki pisteet:
plot(t,X(:,1),'color',[0 0 0]);
axis([t(1),t(end),min(X(:,1)),max(X(:,1))]);
otsikko = strcat(num2str(tarkasteltavat(1)),'-',mn(1),'[',yk(1),']');
title(otsikko);
%ylabel(yk(1));
% *************************************************************************
for i=2:muuttujia
    
    subplot(muuttujia,muuttujia,i+(i-1)*muuttujia);
   
    %plot(t(1:m:n),X(1:m:n,i),'color',[0 0 0]);
    plot(t,X(:,i),'color',[0 0 0]);
    axis([t(1),t(end),min(X(:,i)),max(X(:,i))]);
    otsikko = strcat(num2str(tarkasteltavat(i)),'-',mn((i)),'[',yk(i),']');
    title(otsikko);
    %ylabel(yk(i));
    
    for j=1:i-1
        k=(i-1)*muuttujia+j;
        subplot(muuttujia,muuttujia,k);
        plot(X(1:m:n,(i)),X(1:m:n,(j)),'.','MarkerSize',1,'color',[0 0 0]);
        
        %******************************************************************
        % Väripisteiden lisäys (sininen vanha, punainen uusi data)
        hold on;
        for apuj = 1:m2:n
            apuvari = apuj/(n*1.2);
            plot(X(apuj,(i)),X(apuj,(j)),'.','MarkerSize',1,'color',[apuvari 1-apuvari 1-apuvari])
        end
        hold off;
        %******************************************************************
        % lin.korrelaatio
        r = corr(X(:,(i)),X(:,(j)));
        otsikko = strcat('r=',num2str(r));
        title(otsikko);
    end
end
%%
% *************************************************************************
% yläkolmioon differentiaalit
% *************************************************************************
% scatteriin differentiaaliparit: muuttuu huomattavasti riippuen
% slidingista ja differentiaalin muodostamistyylist


if s==1
for i=1:muuttujia
    for j=1+i:muuttujia
        k=(i-1)*muuttujia+j;
        subplot(muuttujia,muuttujia,k);
        plot(dX(1:m:n,(i)),dX(1:m:n,(j)),'.','MarkerSize',1,'color',[0 0 0])
        %******************************************************************
        % Väripisteiden lisäys (sininen vanha, punainen uusi data)
        hold on;
        for apuj = 1:m2:n
            apuvari = apuj/(n*1.2);
            plot(dX(apuj,(i)),dX(apuj,(j)),'.','MarkerSize',1,'color',[apuvari 1-apuvari 1-apuvari])
        end
        hold off;
        %******************************************************************
        % lin.korrelaatio
        r = corr(dX(:,(i)),dX(:,(j)));
        %otsikko = strcat(mn(tarkasteltavat(i)),'-',mn(tarkasteltavat(j)),' r=',num2str(r));
        otsikko = strcat('r=',num2str(r));
        title(otsikko)
    end
end

end

%%
% scatteriin differentiaaliparit: muuttuu huomattavasti riippuen
% slidingista ja differentiaalin muodostamistyylist
% lisäksi muunneltu korrelaation laskutapa differentaaliin
% *************************************************************************
% Korrelaatio diff. pareissa vain diff_rajan ylittävät otetaan huomioon
% *************************************************************************

% PITÄISI TOIMIA... vaatiin vaan hieman säätöä!
if s==2
    
    
% 10.10.2006 tässä mahdollisesti nyt ongelma, sillä ei tuoda kuin
% tarkasteltavat vektorit
diff_rajat = std(dX)*2;


for i=1:muuttujia
    r1_o = abs(dX(:,(i)))>diff_rajat((i));
       
    for j=1+i:muuttujia
        k=(i-1)*muuttujia+j;
        subplot(muuttujia,muuttujia,k);
                     
        plot(dX(1:m:n,(i)),dX(1:m:n,(j)),'.','MarkerSize',1,'color',[0 0 0])
        % scatterin väritys ajansuhteen!? esimerkiksi jako 3 osaan eri
        % värisävyillä. 
        %******************************************************************
        % Väripisteiden lisäys (sininen vanha, punainen uusi data)
        hold on;
        
        %line([diff_rajat 1],[1 1]);
        
        for apuj = 1:m2:n
            apuvari = apuj/(n*1.2);
            plot(dX(apuj,(i)),dX(apuj,(j)),'.','MarkerSize',1,'color',[apuvari 1-apuvari 1-apuvari])
        end
        hold off;
        %******************************************************************
        % lin.korrelaatio
        r = corr(dX(:,(i)),dX(:,(j)));
        
        r2_o = abs(dX(:,(j)))>diff_rajat((j));
        
        r3_o=r1_o+r2_o;
        r_data = zeros(2,n)';
        l=0;    
        hold on;
        for mm=1:n
            if (r3_o(mm)==2)
                l=l+1;
                r_data(l,1) = dX(mm,(i));
                r_data(l,2) = dX(mm,(j));
                apuvari = mm/(n*1.2);
                plot(dX(mm,(i)),dX(mm,(j)),'.','MarkerSize',10,'color',[apuvari 1-apuvari 1-apuvari])
            end
        end
        hold off;
        %plot(r_data(:,1),r_data(:,2),'.')
        % aika hyvä scatteri
        
        % TÄHÄN pitäisi saada: valitsemaan vain ne pisteparit mukaan corr
        % analyysiin jotka ylittää ton rajan
        
        % lopuksi otetaan matriisista talteen kaikki yli rajan ylittäneet
        % pisteparit ja jätetään nollat veks
        if (l>0)
            r_apu1=r_data(1:l,1);
            r_apu2=r_data(1:l,2);
            r_oma = corr(r_apu1,r_apu2);
        else
            r_oma=0;
        end
        
        %otsikko = strcat(mn(tarkasteltavat(i)),'-',mn(tarkasteltavat(j)),' r=',num2str(r));
        otsikko = strcat('r=',num2str(r),' r_o=',num2str(r_oma));
        title(otsikko)
    end
end

end

if parivertailu.kuvantallennus == 1
    % hakemisto
    cd(p.kuvahakemisto);
    %kuvatiedosto = strcat('t',num2str(tarkastelu.kuvanumero),'m',num2str(muuttujia),'e',num2str(tarkastelu.esikasittely),'_',DATESTR(now,'yymmdd'),'.png');
    apu = 'm';
    for i=tarkastelu.variables
        apu = strcat(apu,num2str(i),'_');
    end
    kuvatiedosto = strcat('t',num2str(0),apu,'e',num2str(1),'_',DATESTR(now,'yymmdd'),'.png');
    print('-dpng',kuvatiedosto);
    cd(p.tyoskentelyhakemisto);
end
%correlations = r;

