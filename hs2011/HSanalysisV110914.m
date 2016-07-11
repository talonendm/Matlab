% *************************************************************************
% HSopen vaalikone ja tulos analysointi
% talonen
% 29.4.2011 (siistitty koodi)
% *************************************************************************
clear all;
close all;
% path(path,'/home/talonen/DMT2demo/extratools')
% path(path,'L:\regrtoolbox\regrtoolbox')
% *************************************************************************
% Haetaan ohjelmiston sijainti
% *************************************************************************
p.cd = cd;
p.mat = strcat(p.cd(1:end-6),'mat/');
p.fig = strcat(p.cd(1:end-6),'fig/');
p.avi = strcat(p.cd(1:end-6),'avi/');
p.f = strcat(p.cd(1:end),'/f/');                   % functioita
cd(p.mat);
load hs.mat;        % Hsopen 2011
load logot.mat;     % omat logot 29 kpl
load A.mat;         % yle data
load ehdokasraha935v110530.mat;


cd(p.cd);

% *************************************************************************
% korjataan tietokantoi
% *************************************************************************
cd(p.f);[hs A] = KorjaaTietokantojaW7(hs,A); cd ..
% *************************************************************************
%%
% *************************************************************************
% parametriasetuksissa
p.idapaperfigure = 1;
% *************************************************************************
p = loadHSparameters(p);
% *************************************************************************

p.idaposterprint = 1;


% *************************************************************************
% *************************************************************************
% purkka - muuta A.aani maara rahaksi
% 110530
% *************************************************************************
p.mitarahaatarkastellaan = 2;  %12 sanomalehti, 6 yritysrahoitus
for i=1:length(ehdokasraha.r2A)
    if isnan(ehdokasraha.euroa(i,p.mitarahaatarkastellaan))
        ehdokasraha.euroa(i,p.mitarahaatarkastellaan) = 0;
    end
    mediaaniraha = median(ehdokasraha.euroa(:,p.mitarahaatarkastellaan));
end

A.rahaa = ones(1,2315)*mediaaniraha;  % mediaaniraha on oletus
for i=1:length(ehdokasraha.r2A)
    
    if ehdokasraha.r2A(i)>0
        A.rahaa(ehdokasraha.r2A(i)) = ehdokasraha.euroa(i,p.mitarahaatarkastellaan);
        % 825 tietoa matchaa..
    end   
end
% *************************************************************************
% *************************************************************************





% *************************************************************************



%%
% *************************************************************************
% F - display answers
% *************************************************************************
if p.vertailekompoissakahta == 1
    cd(p.f);vertailevastauksia(hs,p,'Soini Timo','Katainen Jyrki'); cd ..
end
% *************************************************************************


% **********
% *************************************************************************
%%
% *************************************************************************
% Lasketaan HS-datasta apumuuttujia
p.ikahistX = 20:5:75;%p.ikahistX = [20:5:75];
m.vaalipiirit = unique(hs.Vaalipiiri);
m.ammatit = unique(hs.Ammatti);
%m.sukupuolet = unique(hs.Sukupuoli);  % 0=nainen, 1=mies , onko?
m.outlier.ika = find(hs.Ika==max(hs.Ika)); % ehdokas 110 vuotias..
hs.Ika(m.outlier.ika) = nan;
n = length(hs.Ika); %2306
% ********************************************
m.n.vaalipiirit = length(m.vaalipiirit);
for i=1:m.n.vaalipiirit
    m.ehdokkaat{i} = find(strcmp(hs.Vaalipiiri,m.vaalipiirit{i})==1);
end
% *******************************************
m.puolueet = unique(hs.Puolue);
m.puolueet{length(m.puolueet)+1} = 'MUUT';

m.puolueetYLE = m.puolueet;
m.puolueetYLE{2} = 'IPU';
m.puolueetYLE{3} = 'KÖY';
m.puolueetYLE{12} = 'M11';
m.puolueetYLE{20} = 'SSP';


m.n.puolueet = length(m.puolueet);
for i=1:m.n.puolueet
    m.puolueN{i} = find(strcmp(hs.Puolue,m.puolueet{i})==1);
    %m.puolueN{i} = apu{1};
    m.puolueEhdokkaita(i) = length(m.puolueN{i});
end
% *************************************************************************
% party color setting
cd(p.f); p = pvarit(m,p); cd ..   
% *************************************************************************
%%
% *************************************************************************
% ! FIGS ! 
% *************************************************************************
% Ikajakauma vaalipiireittäin
% FIG1
% cd(p.f);ikahistogrammi(hs,m,p,1);cd ..
% *************************************************************************
% FIGS
% Tarkastellaan arvomaailmaa (yla-anttila)
% cd(p.f);arvomaailmasplini(hs,m,p,2);cd ..
% *************************************************************************

% *************************************************************************
% kysymykset ja vastaukset
% cd(p.f);b = kysymyksetjavastaukset(hs,m,p);cd ..
% *************************************************************************
%%
% *************************************************************************
% F
% *************************************************************************
% data binaariseksi
cd(p.f);b = databinaariseksi(hs,m,p);cd ..
% *************************************************************************
% 110909 - tehdään uusimuuttuja. HS apps4finlandia varten

%
for iii=1:length(b.vastaus)
    disp(strcat(num2str(iii),':',b.vastaus{iii}));
end
%


p.konsessus1 = zeros(1,173);
p.konsessus2 = zeros(1,173);


% http://tstm.info/files/konsensuspisteet.pdf

p.konsessus1(1) = -100;
p.konsessus1(2) = -50;
p.konsessus1(4) = 50;
p.konsessus1(5) = 100;

p.konsessus2(6) = -100;
p.konsessus2(7) = 100;

p.konsessus1(8) = -50; % huom. ei ja kyllä omassa listassa toisten päin
p.konsessus1(9) = 50;

p.konsessus2(8) = 50;
p.konsessus2(9) = -50;

p.konsessus1(10) = -80;
p.konsessus2(10) = 0;

p.konsessus1(11) = -30;
p.konsessus2(11) = -25;

p.konsessus1(12) = -60;
p.konsessus2(12) = -25;

%13

p.konsessus1(14) = 50;
p.konsessus1(15) = 25;
p.konsessus1(16) = -25;

p.konsessus1(17) = -40;
p.konsessus2(17) = 50;

p.konsessus1(19) = 20;
p.konsessus2(19) = -25;

p.konsessus1(18) = 60;
p.konsessus2(18) = -50;

p.konsessus1(24) = 40;
p.konsessus2(24) = -25;

p.konsessus1(25) = 40;
p.konsessus2(25) = -25;

p.konsessus1(22) = 40;
p.konsessus2(22) = -25;

p.konsessus1(23) = 100;
p.konsessus2(23) = -100;

p.konsessus1(30) = -30;
p.konsessus1(26) = -30;
p.konsessus1(27) = -30;
p.konsessus1(28) = 30;


p.konsessus1(36) = 25;
p.konsessus2(36) = 25;

p.konsessus1(32) = 50;
p.konsessus2(32) = 40;

p.konsessus1(35) = 70;
p.konsessus2(35) = 25;

p.konsessus1(40) = 25;
p.konsessus2(40) = 25;

p.konsessus1(34) = -50;
p.konsessus2(34) = -25;

p.konsessus1(37) = 50;
p.konsessus2(37) = 0;

p.konsessus1(33) = 70;
p.konsessus2(33) = 50;

p.konsessus1(38) = 25;
p.konsessus2(38) = 25;

p.konsessus1(39) = -100;
p.konsessus2(39) = 0;


p.konsessus1(45) = 25;
p.konsessus2(45) = -50;

p.konsessus1(46) = 50;
p.konsessus2(46) = -25;

p.konsessus1(42) = 100;
p.konsessus2(42) = 0;

p.konsessus1(41) = 25;
p.konsessus2(41) = 25;

p.konsessus1(43) = -25;
p.konsessus2(43) = 40;

p.konsessus1(44) = -100;
p.konsessus2(44) = 50;


p.konsessus1(49) = 100;
p.konsessus2(49) = 50;

p.konsessus1(48) = -50;
p.konsessus2(48) = -25;

p.konsessus1(47) = -100;
p.konsessus2(47) = -50;

p.konsessus1(50) = -50;
p.konsessus1(57) = -50;
p.konsessus1(58) = -70;
p.konsessus1(51) = -10;
p.konsessus1(54) = -25;
p.konsessus1(55) = -50;

p.konsessus1(56) = -70;
p.konsessus2(56) = -25;

p.konsessus1(52) = 50;
p.konsessus1(53) = 100;


p.konsessus1(63) = -50;
p.konsessus2(63) = 30;

p.konsessus1(62) = -25;
p.konsessus2(62) = 30;

p.konsessus1(59) = 50;
p.konsessus2(59) = 20;

p.konsessus1(65) = -25;
p.konsessus2(65) = 20;

p.konsessus1(61) = -50;
p.konsessus2(61) = 0;

p.konsessus1(64) = 50;
p.konsessus2(64) = 0;

p.konsessus1(60) = 0;
p.konsessus2(60) = -50;



p.konsessus1(67) = -100;
p.konsessus2(67) = 30;

p.konsessus1(66) = -50;
p.konsessus2(66) = 15;

p.konsessus1(69) = 50;
p.konsessus2(69) = 0;



p.konsessus1(73) = 0;
p.konsessus2(73) = -100;

p.konsessus1(76) = 25;
p.konsessus2(76) = 25;

p.konsessus1(75) = 0;
p.konsessus2(75) = 50;

p.konsessus1(72) = 50;
p.konsessus2(72) = 100;


p.konsessus2(79) = -50;
p.konsessus2(81) = -25;
p.konsessus2(80) = -25;
p.konsessus2(77) = 25;
p.konsessus2(78) = 50;



p.konsessus1(84) = 0;
p.konsessus2(84) = -25;

p.konsessus1(85) = 0;
p.konsessus2(85) = -25;

p.konsessus1(83) = 50;
p.konsessus2(83) = 0;

p.konsessus1(82) = 25;
p.konsessus2(82) = 50;



p.konsessus1(88) = -100;
p.konsessus2(88) = -100;

p.konsessus1(87) = -50;
p.konsessus2(87) = -75;

p.konsessus1(86) = 50;
p.konsessus2(86) = 50;


p.konsessus1(93) = -100;
p.konsessus1(92) = -50;
p.konsessus1(91) = 50;
p.konsessus1(89) = 100;


p.konsessus2(121) = -100;
p.konsessus2(122) = -50;
p.konsessus2(120) = 50;



p.konsessus2(124) = 100;
p.konsessus2(123) = -100;




p.konsessus2(126) = 25;
p.konsessus2(127) = 25;
p.konsessus2(125) = -25;



p.konsessus1(130) = -25;
p.konsessus2(130) = 50;

p.konsessus1(129) = 25;
p.konsessus2(129) = -50;



p.konsessus2(134) = 25;
p.konsessus2(136) = -25;
p.konsessus2(135) = -50;
p.konsessus2(132) = -75;



p.konsessus1(138) = 100;
p.konsessus2(138) = 30;

p.konsessus1(140) = 75;
p.konsessus2(140) = 20;

p.konsessus1(141) = 50;
p.konsessus2(141) = 10;

p.konsessus1(139) = -25;
p.konsessus2(139) = 0;

p.konsessus1(142) = -80;
p.konsessus2(142) = -30;



p.konsessus1(152) = 25;
p.konsessus2(152) = 50;

p.konsessus1(151) = 25;
p.konsessus2(151) = -50;

p.konsessus1(149) = -50;
p.konsessus2(149) = 50;



p.konsessus1(154) = -100;
p.konsessus1(155) = 50;
p.konsessus1(156) = 100;


b.konsessustaulu1 = b.vastaustaulu.*(ones(length(b.vastaustaulu),1)*p.konsessus1);
b.konsessustaulu2 = b.vastaustaulu.*(ones(length(b.vastaustaulu),1)*p.konsessus2);


b.konsessustaulu1m = b.konsessustaulu1.*(1+b.binaarimassa/2);
b.konsessustaulu2m = b.konsessustaulu2.*(1+b.binaarimassa/2);

b.konsessusX = sum(b.konsessustaulu1')';
b.konsessusY = sum(b.konsessustaulu2')';

b.konsessusXm = sum(b.konsessustaulu1m')';  % painotettu 0.5,1.1.5 :llä
b.konsessusYm = sum(b.konsessustaulu2m')';


%%
% *************************************************************************
%baktaulu = sign(b.binaaritaulu);
vastannutkyssareihin = find((sum(b.vastaustaulu,2)>p.moneenkovastannuettamukaan))'  % virheellisesti painotusten jälkeen rajaus
% *************************************************************************
% list: number of answers
p.listaavajaastivastaanneet = 0;
if p.listaavajaastivastaanneet == 1
    apusumma = 0;
    for i=0:47
        vajaastivastanneet = find((sum(b.vastaustaulu,2)==i));
        disp(strcat(num2str(i),':',num2str(length(vajaastivastanneet))));
        for j=1:length(vajaastivastanneet)
            %disp(hs.Nimi{vajaastivastanneet(j)});
        end
        if i>20 && i<33
            apusumma = apusumma + length(vajaastivastanneet);
        end
    end
end
% *************************************************************************


% *************************************************************************
% 110504 PAINOTUS Kokoomuksen aanien perusteella vastaukselle
% idea on saada suuntaa vastaukselle kysymyksissä
% kokoomus luokittelijana siis - voi olla mikä puolue vaan periaatteessa
% *************************************************************************
if p.puoluekohtainenpainotus>0
    for i=1:length(b.vastaus)
        % hs.k.Vastaus(m.puolueN{p.puoluekohtainenpainotus},:)
        summa = 0;
        for j=1:m.puolueEhdokkaita(p.puoluekohtainenpainotus)
            if b.binaaritaulu(m.puolueN{p.puoluekohtainenpainotus}(j),i)>0
                summa = summa + b.binaaritaulu(m.puolueN{p.puoluekohtainenpainotus}(j),i);
            end
        end
        disp(summa)
        b.binaaritaulu(:,i) = b.binaaritaulu(:,i).*sqrt(1+summa);
    end
end
% *************************************************************************

% *************************************************************************
% vastaamatta jättäneiden poistaminen
% *************************************************************************



b.binaaritaulu = b.binaaritaulu(vastannutkyssareihin,:); % ks. ylempana



b.konsessusXv = b.konsessusX(vastannutkyssareihin,:);
b.konsessusYv = b.konsessusY(vastannutkyssareihin,:);

b.konsessusXvm = b.konsessusXm(vastannutkyssareihin,:);
b.konsessusYvm = b.konsessusYm(vastannutkyssareihin,:);


% TESTI!
%b.binaaritaulu = b.binaaritaulu(1:200,:);


% 110803 - aikasarjaa - pikakoodii = paskaa
p.valtiopaivat = 0;
if p.valtiopaivat == 1
    b.binaaritaulu = [b.binaaritaulu rand(length(b.binaaritaulu),3)];
    b.binaarimassa = [b.binaarimassa ones(length(b.binaarimassa),3)];
    
    b.vastaus = {b.vastaus{:}, 'fdf1', 'fdf2', 'fdf3'};
    b.kysymysnumero(174) = 31;
    b.kysymysnumero(175) = 31;
    b.kysymysnumero(176) = 31;
    
end


%p.regrKompoja = 40;
%b.binaaritaulu = b.binaaritaulu(:,1:p.regrKompoja);
% *************************************************************************


% *************************************************************************
% muu painotus KATAINEN
% kielipainotus (mieti tosiaan tuo kyllä ei kohtisuoraakselit (ei täysin
% painojen takia.. mitä jos tekis summa muuttujan?
% *************************************************************************
p.kieliyhdeksimuuttujaksitesti = 0;
% *************************************************************************
if p.kieliyhdeksimuuttujaksitesti == 1
    paikka = find(b.kysymysnumero==23)
    b.binaaritaulu(:,paikka(1)) = (b.binaaritaulu(:,paikka(1)) - b.binaaritaulu(:,paikka(2)))*3; 
    b.binaaritaulu(:,paikka(2)) = b.binaaritaulu(:,paikka(1))*0;
    % tyhjä sarake pitäis poistaa
end
% *************************************************************************
% Muut painot
% *************************************************************************
if p.kompo2maaritys ==1
    for i=1:length(p.katainen)
    %for i=p.katainen
        % ei toimi näin
        % aputaulu = TODELLA mystistä.. workspaveen ei päivity, mutta
        % toimii muuten 110504
        b.binaaritaulu(:,b.kysymysnumero==p.katainen(i)) = b.binaaritaulu(:,b.kysymysnumero==p.katainen(i)).*p.katainenpaino(i);%.*p.katainenpaino(i);
        %b.binaaritaulu(:,find(b.kysymysnumero==p.katainen(i))) = b.binaaritaulu(:,find(b.kysymysnumero==p.katainen(i))).*p.katainenpaino(i);%.*p.katainenpaino(i);
 
    end
end
% *************************************************************************

%%

b.b2hs = vastannutkyssareihin;

jiiS.vasenoikea =  ones(A.n,1)*100;  % tallenna tänh %2315
jiiS.libkon =  ones(A.n,1)*100;  % tallenna tänh %2315

% *************************************************************************
if p.kompo2maaritys == 1
    % tehdään omat akselistot
    if p.regrPCA == 0
        [COEFF,SCO,latent,tsquare] = princomp(b.binaaritaulu(:,b.Xalue));
        [COEFF2,SCO2,latent2,tsquare2] = princomp(b.binaaritaulu(:,b.Xalue2));
        SCORE(:,1) = SCO(:,1);
        SCORE(:,2) = SCO2(:,1);   
    else
        
        % center the data
        b.binaaritaulu = regrCenter(b.binaaritaulu);
    
        %[theta,lambda,Q,T2] = regrPCA(X,N)
        p.regrKompoja = length(b.Xalue);
        [COEFF,latent,Q,tsquare] = regrPCA(b.binaaritaulu(:,b.Xalue),p.regrKompoja);
        p.regrKompoja = length(b.Xalue2);
        [COEFF2,latent2,Q2,tsquare2] = regrPCA(b.binaaritaulu(:,b.Xalue2),p.regrKompoja);
        %[COEFF,latent,Q,tsquare] = regrPCA(b.binaaritaulu);
        % scores
        %Z = X*theta_pca;
        SCORE = zeros(length(b.binaaritaulu),2);
        kSCORE = b.binaaritaulu(:,b.Xalue)*COEFF;
        SCORE(:,1) = kSCORE(:,1);
        kSCORE2 = b.binaaritaulu(:,b.Xalue2)*COEFF2;
        SCORE(:,2) = kSCORE2(:,1);
    end
    
else
    % 1 ja 2 akseli mukaan
    if p.regrPCA == 0
        [COEFF,SCORE,latent,tsquare] = princomp(b.binaaritaulu);
    else
        
        % center the data
        
        mallit.ka0 = mean(b.binaaritaulu);
        
        b.binaaritaulu = regrCenter(b.binaaritaulu);

        
        mallit.ka1 = mean(b.binaaritaulu);
        
        %[theta,lambda,Q,T2] = regrPCA(X,N)
        p.regrKompoja = length(b.vastaus);
        [COEFF,latent,Q,tsquare] = regrPCA(b.binaaritaulu,p.regrKompoja);
        %[COEFF,latent,Q,tsquare] = regrPCA(b.binaaritaulu);
        % scores
        %Z = X*theta_pca;
        SCORE = b.binaaritaulu*COEFF;
    end
        
    % SCORE 1 -> 3
    % SCORE(:,1) = SCORE(:,3);
end


% *************************************************************************
NiemiTapani = 0;
disp('Selvitetään dublikaatit');
% *************************************************************************
for j=1:length(hs.Nimi)
    % KESKEN
    H.puolue(j) = find(strcmp(m.puolueet,hs.Puolue{j}));
    
    apu = find(strcmp(A.nimi,hs.Nimi{j}));
    
    
    % ONGELMA - Matti Kyllönen dublikaatti, mutta ei pelkästään
    % binääritaulussa. Tämä selvittää vain sen.
    
    if length(apu)>1
        % pyritään puolueen avulla
        disp(strcat(num2str(apu(1)),',',num2str(apu(2))))
        apu = apu(find(strcmp(A.puolue(apu),m.puolueet{H.puolue(j)})));
        if length(apu)==1
            disp(apu)
            disp('ok..');
        else
            disp('sama nimi ja sama puolue');
            disp(A.nimi(apu(1)));
            % pitäis selvittää, missä vaalipiirissä, muttei jaksa ny
            disp('arvio kumpi on kampi');
            NiemiTapani = NiemiTapani +1;
            if NiemiTapani==1
                apu = apu(1);
            else
                apu = apu(2);
            end
            disp(apu)
            
        end
    end
    
    
    if length(apu)==1
        H.ylenimi(j) = apu;
        H.yleStr{j} = hs.Nimi{j};
    else
        if length(apu)>2
            disp('kolme samaa nimea');
            disp((apu(1)));
            disp((apu(2)));
            disp((apu(3)));
           
        end
        if length(apu) == 2
            disp('SAMA NIMI - virhe koodissa'); 
            disp((apu(1)));
            disp((apu(2)));
            if find(strcmp(hs.Nimi{apu(1)},'PS'))
                H.ylenimi(j) = apu(1);
            else
                H.ylenimi(j) = apu(2); % ainakin Tossavainen
            end
        end
        H.ylenimi(j) = 0;
        H.yleStr{j} = hs.Nimi{j};
        
        disp(strcat('Eri nimi YLE vs. HS',hs.Nimi{j}));
        for jjj=1:length(A.nimi)
            if strcmp(A.nimi{jjj}(1:3),hs.Nimi{j}(1:3))
                %disp(strcat('yle:',num2str(jjj),A.nimi{jjj}))
                disp(strcat('yle:',num2str(jjj),A.nimi{jjj}))
                
            end
        end
        
        disp('Alhojärvi Lauri puuttuu');
        if strcmp(hs.Nimi{j},'Alhojärvi Lauri')
            H.yleStr{j} = 'Alhojärvi Lauri puuttuu';
            H.ylenimi(j) = nan;
        end
        
    end
    
end


b.ylenimi = H.ylenimi(b.b2hs);
b.puolue = H.puolue(b.b2hs);
%b.yleStr




%%
% *************************************************************************
% PCA:n latenttien ja coeffin tarkastelu
% *************************************************************************
p.pcatarkkailu = 1;
if p.pcatarkkailu == 1
    figure(35);
    m.pca.selitysaste = latent./sum(latent)*100;
    subplot(2,1,1);
    bar(m.pca.selitysaste(1:7));
    title('latent');
    subplot(2,1,2);
    bar(cumsum(m.pca.selitysaste(1:7)));
    title('cusum latent');
    
    figure(36);
    subplot(2,1,1);
    bar(COEFF(:,1));
    colormap(cool);
    if p.kompo2maaritys ==0
    for i=1:length(COEFF(:,1))
        text(i,COEFF(i,1)+sign(COEFF(i,1))*0.02, b.vastaus{i},'horizontalalignment','left','verticalalignment','middle','fontsize',8)
        if i>1
            if b.kysymysnumero(i-1)~=b.kysymysnumero(i)
                text(i,COEFF(i,1)+sign(COEFF(i,1))*0.1, hs.k.Kysymys{b.kysymysnumero(i)},'color',[1 0 0],'horizontalalignment','left','verticalalignment','middle','fontsize',8)
            end
        else
            text(i,COEFF(i,1)+sign(COEFF(i,1))*0.1, hs.k.Kysymys{b.kysymysnumero(i)},'color',[1 0 0],'horizontalalignment','left','verticalalignment','middle','fontsize',8)
        end
    end
    else
        subplot(2,1,1);
        bar(COEFF(:,1));
        title(p.labelX);
        
        koffi = 1;
        for i=1:length(COEFF(:,koffi))
            text(i,COEFF(i,koffi)+sign(COEFF(i,koffi))*0.02, b.vastaus{b.Xalue(i)},'horizontalalignment','left','verticalalignment','middle','fontsize',8)
            if i>1
                if b.kysymysnumero(b.Xalue(i-1))~=b.kysymysnumero(b.Xalue(i))
                    text(i,COEFF(i,koffi)+sign(COEFF(i,koffi))*0.1, hs.k.Kysymys{b.kysymysnumero(b.Xalue(i))},'color',[1 0 0],'horizontalalignment','left','verticalalignment','middle','fontsize',8)
                end
            else
                text(i,COEFF(i,koffi)+sign(COEFF(i,koffi))*0.1, hs.k.Kysymys{b.kysymysnumero(b.Xalue(i))},'color',[1 0 0],'horizontalalignment','left','verticalalignment','middle','fontsize',8)
            end
        end
        
    end
    
    
   
    if p.kompo2maaritys ==0
        
        subplot(2,1,2);
        bar(COEFF(:,2));
        
        
    koffi = 2;
    for i=1:length(COEFF(:,koffi))
        text(i,COEFF(i,koffi)+sign(COEFF(i,koffi))*0.02, b.vastaus{i},'horizontalalignment','left','verticalalignment','middle','fontsize',8)
        if i>1
            if b.kysymysnumero(i-1)~=b.kysymysnumero(i)
                text(i,COEFF(i,koffi)+sign(COEFF(i,koffi))*0.1, hs.k.Kysymys{b.kysymysnumero(i)},'color',[1 0 0],'horizontalalignment','left','verticalalignment','middle','fontsize',8)
            end
        else
            text(i,COEFF(i,koffi)+sign(COEFF(i,koffi))*0.1, hs.k.Kysymys{b.kysymysnumero(i)},'color',[1 0 0],'horizontalalignment','left','verticalalignment','middle','fontsize',8)
        end
    end
    
    else
       
        subplot(2,1,2);
        bar(COEFF2(:,1));
        title(p.labelY);
        koffi = 1;
        for i=1:length(COEFF2(:,koffi))
            text(i,COEFF2(i,koffi)+sign(COEFF2(i,koffi))*0.02, b.vastaus{b.Xalue2(i)},'horizontalalignment','left','verticalalignment','middle','fontsize',8)
        if i>1
            if b.kysymysnumero(b.Xalue2(i-1))~=b.kysymysnumero(b.Xalue2(i))
                text(i,COEFF2(i,koffi)+sign(COEFF2(i,koffi))*0.1, hs.k.Kysymys{b.kysymysnumero(b.Xalue2(i))},'color',[1 0 0],'horizontalalignment','left','verticalalignment','middle','fontsize',8)
            end
        else
            text(i,COEFF2(i,koffi)+sign(COEFF2(i,koffi))*0.1, hs.k.Kysymys{b.kysymysnumero(b.Xalue2(i))},'color',[1 0 0],'horizontalalignment','left','verticalalignment','middle','fontsize',8)
        end
    end
        
        
    end
    
    
    
end
% *************************************************************************    
%%
% ************************************************************************* 
% SEURAAVASSA tutkimuksessa. käännä esim. kok. mahd. oikealle.

if p.rotaatio == 1
    % ************************************************************************* 
    % ROTAATIOtesti
    
    
    % VASEMMISTO - OIK
    keta1 = find(b.puolue==find(strcmp(m.puolueet,p.rotaatiopuolue1)))';
    kx1 = median(SCORE(keta1,1));
    ky1 = median(SCORE(keta1,2));
    keta2 = find(b.puolue==find(strcmp(m.puolueet,p.rotaatiopuolue2)))';
    kx2 = median(SCORE(keta2,1));
    ky2 = median(SCORE(keta2,2));
    
    
    % HUOM! 3D kuvaajissa asteikko eritavalla
    if p.idapaperfigure > 0
        ast = (2*pi/360)*-90; % HUOM! 3D kuvaajissa asteikko eritavalla
    else
        ast = (2*pi/360)*180; % 2D kuviin
        ast = (2*pi/360)*0; % HUOM! 3D kuvaajissa asteikko eritavalla
    end
    ast = ast-atan((ky1-ky2)/(kx1-kx2));
    sc1 = SCORE(:,1);
    sc2 = SCORE(:,2);

    
    
    SCORE(:,1) = cos(ast)*sc1 - sin(ast)*sc2;
    SCORE(:,2) = sin(ast)*sc1 + cos(ast)*sc2;
end
% ************************************************************************* 


% *************************************************************************    
% PCA original
%figure(3);


figure(334);

subplot(2,2,1);
hold on;
for i=1:length(b.binaaritaulu)
    plot(b.konsessusXv(i),SCORE(i,1),'.','color',p.vari{(strcmp(m.puolueet,hs.Puolue{vastannutkyssareihin(i)}))}); 
end
ylabel('pc1');
%xlabel('vasemmisto-oikeisto');
xlabel('left-right');

axis tight;

subplot(2,2,2);
hold on;
for i=1:length(b.binaaritaulu)
    plot(b.konsessusYv(i),SCORE(i,2),'.','color',p.vari{(strcmp(m.puolueet,hs.Puolue{vastannutkyssareihin(i)}))}); 
end
ylabel('pc2');
xlabel('restrictive-liberal');
axis tight;
subplot(2,2,3);
hold on;
for i=1:length(b.binaaritaulu)
  %plot(SCORE(i,1),SCORE(i,2),'.','color',p.vari{find(strcmp(m.puolueet,hs.Puolue{vastannutkyssareihin(i)}))}); 
  plot(SCORE(i,1),SCORE(i,2),'.','color',p.vari{(strcmp(m.puolueet,hs.Puolue{vastannutkyssareihin(i)}))}); 
end
xlabel('pc1');
ylabel('pc2');

axis tight;
subplot(2,2,4);
hold on;
for i=1:length(b.binaaritaulu)
    plot(b.konsessusXv(i),b.konsessusYv(i),'.','color',p.vari{(strcmp(m.puolueet,hs.Puolue{vastannutkyssareihin(i)}))}); 
end
xlabel('restrictive-liberal');
xlabel('left-right');


%title('Martti Leppänen');
axis tight;



% **********************************************
% Print posterpic - ei joka kerta- kun onnistunut kuva, niin 0
p.idaposterprint =0;
% **********************************************
if p.idaposterprint == 1
    
     p.tallennusX = 10;
     p.tallennusY = 10;
    
     set(gcf,'PaperUnits','inches','PaperPosition',[0 0 p.tallennusX p.tallennusY]);
            
     % HUOM! ZEKKAA TIEDOSTO
     tallennustiedosto = strcat('posterpics/','k334','.png');
     print('-dpng','-r300',tallennustiedosto);% pngWille/testi111015.png -r150
    
end
% **********************************************






figure(335);
hold on;
for i=1:length(b.binaaritaulu)
    plot(b.konsessusXvm(i),b.konsessusYvm(i),'.','color',p.vari{(strcmp(m.puolueet,hs.Puolue{vastannutkyssareihin(i)}))}); 
end
ylabel('restriktiivi-liberaali');
xlabel('vasemmisto-oikeisto');
title('Painotettu');
axis tight;





% *************************************************************************
% Samankaltaisuus:
%tarkastele = 44;
%for i=1:length(b.binaaritaulu)
%   eroavastauksissa(i) = sum(abs(b.binaaritaulu(tarkastele,:)-b.binaaritaulu(i,:)));
%end
%figure(66);plot(eroavastauksissa);
% *************************************************************************
% Skaalaus -100..100
p.skaalaus = 1;
if p.skaalaus == 1
for i=1:2
    %SCORE(:,i) = - SCORE(:,i);
    
    m.pca.min(i) = min(SCORE(:,i));
    m.pca.max(i) = max(SCORE(:,i));
    
    SCORE(:,i) = SCORE(:,i)-m.pca.min(i);
    % huom. max min korjauksen jälkeen ok
    SCORE(:,i) = SCORE(:,i)/max(SCORE(:,i))*200-100;
    
    
end
end
% *************************************************************************

%%
figure(31);
%subplot(2,2,1);
hold on;

for i=1:length(b.binaaritaulu)
    j = vastannutkyssareihin(i);
    plot(SCORE(i,1),SCORE(i,2),'.','color',p.vari{find(strcmp(m.puolueet,hs.Puolue{j}))}); 
    
end





%%




for i=1:length(m.puolueet)
    keta = find(b.puolue==i);
    if length(keta)>2
        kx = mean(SCORE(keta,1));
        ky = mean(SCORE(keta,2));
        text(kx,ky,m.puolueet{i},'color',[0 0 0],'BackgroundColor',p.vari{i},'fontsize',8,'horizontalalignment','center','verticalalignment','middle');
    end
    
        
end

title('All candidates - voting application data');
%xlabel('vasemmisto - oikeisto');
%ylabel('liberaali - konservatiivi');

%%

%for i=1:2315
%    if strcmp(A.nimi{i}(end-4:end),'Lauri')
%        A.nimi{i}
%    end
%end


%%



% 110421 seuraavaksi
% ehdokkaat hilaan 200x200
% kerrotaan äänimäärällä

%%
% MUOKKAA. tässä randomilla
p.omacolormap = [];
for i=1:length(m.puolueet)
    p.omacolormap = [p.omacolormap;p.vari{i}];
end


%%

% ******************************************************
% TALLENNETAAN MALLIT
% 21.9.2011
% ******************************************************

mallit.kons1 = p.konsessus1';
mallit.kons2 = p.konsessus2';

mallit.pc1 = COEFF(:,1);
mallit.pc2 = COEFF(:,2);

% SKAALAUSTIETO
mallit.pca.pc.min = m.pca.min;
mallit.pca.pc.max = m.pca.max;

% Aihepiirettäin
b.aihepiiri{1} = [1 2 3 4 5];
b.aihepiiri{2} = [6 7 8];
b.aihepiiri{3} = [9 10 11];
b.aihepiiri{4} = [12 13 14 15];
b.aihepiiri{5} = [16 17];
b.aihepiiri{6} = [18 19 20 21];
b.aihepiiri{7} = [22 23 24 25 26];
b.aihepiiri{8} = [27 28 29 30];
b.aihepiiri{9} = [31];

for j=1:9

    mitkakysymykset = zeros(173,1)';
    for i=1:length(b.aihepiiri{j})
        mitkakysymykset = mitkakysymykset + (b.kysymysnumero == b.aihepiiri{j}(i));
    end

    mitkakysymykset = find(mitkakysymykset)
    [mallit.aihepiiri{j},latent11,Q11,tsquare11] = regrPCA(b.binaaritaulu(:,mitkakysymykset),1);

    malliSCORE = b.binaaritaulu(:,mitkakysymykset)*mallit.aihepiiri{j};

    mallit.pca.min(j) = min(malliSCORE(:,1));
    mallit.pca.max(j) = max(malliSCORE(:,1));
    
    %SCORE(:,i) = SCORE(:,i)-m.pca.min(i);
    
    %SCORE(:,i) = SCORE(:,i)/max(SCORE(:,i))*200-100;

    

end

mallit.rivit = b.kysymysnumero;

mallit.aihe{1} = '(yleiset) Kysymykset 1-5';
mallit.aihe{2} = 'Eläkkeet 6-8';
mallit.aihe{3} = 'Talous 9-11';
mallit.aihe{4} = 'Verot 12-15';
mallit.aihe{5} = 'Puolustus 16-17';
mallit.aihe{6} = 'Ulkomaat 18-21';
mallit.aihe{7} = 'Kotimaa 22-26';
mallit.aihe{8} = 'Kunnat 27-30';
mallit.aihe{9} = 'Hallituspohja 31';


mallit.vastaus = b.vastaus;












%%
% ******************************************************
% Ajele t�st� kuvia 110424
% ******************************************************
% surfacessa eri tavalla
%SCORE(:,1) = - SCORE(:,1);   % 110422 mieti viel� onko tarpeen
%SCORE(:,2) = - SCORE(:,2);
disp('Asetetaan ehdokkaat hilaan');
% EDITOI! 110421. jokin virhe maksimi pisteen haussa!

if p.logotkarttaan == 1
    p.hilakoko=201;%261;
    
    
    if p.puoluereunat == 1
        p.varianssi = 6;%8;  %esim 3..5
    else
        %p.varianssi = 8;%12;%8;%8;  %esim 3..5  % TAA IDASSA. p.varianssi = standard deviation in pdf
        
        % pistetään puolue- presidenttiä varten!
        p.varianssi = 12;%12;%8;%8;  %esim 3..5  % TAA IDASSA. p.varianssi = standard deviation in pdf
    end
    
    
    % HSopen2:sta varten.. IDA varteen tossa yllä
    %p.hilakoko=101;
    %p.varianssi = 5;  %esim 3..5
    
else
    p.hilakoko=101;
    p.varianssi = 5;%5;  %esim 3..5
    
    %p.hilakoko=301;
    %p.varianssi = 13;  %esim 3..5
    
end

p.hilaskaala = (p.hilakoko-1)/2;

b.hila = zeros(p.hilakoko,p.hilakoko);
b.puoluehila = zeros(p.hilakoko,p.hilakoko,length(m.puolueet));
b.rahahila = zeros(p.hilakoko,p.hilakoko);


% 110424 kysymyskentt� tarkastele erillisi� kyss�reit�
% mieti painotustasoa, pit�isk� eliminoida.. sill� kyse vain mielipiteen
% suunnasta / ehdokas?! 
p.tarkasteltavakysymys = 3;%23;
p.tarkasteltaviakysymyksia = sum(b.kysymysnumero == p.tarkasteltavakysymys);
p.tarkasteltavat = find(b.kysymysnumero == p.tarkasteltavakysymys);


b.kysymyshila= zeros(p.hilakoko,p.hilakoko,p.tarkasteltaviakysymyksia);  % 2 = vastausten m��r�
b.vastausmaarahila= zeros(p.hilakoko,p.hilakoko);  % 2 = vastausten m��r�
puhislaskuri = 0;
kansislaskuri = 0;
for i=1:length(b.binaaritaulu)
    %b.hila(round(SCORE(i,1)+101),round(SCORE(i,2)+101)) = b.hila(round(SCORE(i,1)+101),round(SCORE(i,2)+101)) + rand;

    if p.kokosuomi == 0
        if strcmp(hs.Vaalipiiri{b.b2hs(i)},p.vaalipiirivalinta)
            p.analyysiinmukaan = 1;
        else
            p.analyysiinmukaan = 0;
        end
    else
        p.analyysiinmukaan = 1;
    end
    
    
    if b.ylenimi(i)>0 && p.analyysiinmukaan == 1
    %if H.ylenimi(b.b2hs(i))>0 && p.analyysiinmukaan == 1
        
        px1 = round((SCORE(i,1)+101)/201*p.hilakoko);
        px2 = round((SCORE(i,2)+101)/201*p.hilakoko);
        % Hetkinen eihän tossa yllä olevassa ole järkee.. 1106014
        % eikus on.. ollaan skaalattu -100.. 100 data.
        
        
        
        
        p1 = pdf('Normal',1:1:p.hilakoko,px1,p.varianssi);  % p.varianssi
        p2 = pdf('Normal',1:1:p.hilakoko,px2,p.varianssi);
        
        % jos p1 tai p2 yli
        p1 = p1 / sum(p1);
        p2 = p2 / sum(p2);
   
        
        b.all.X(i) = px1;
        b.all.Y(i) = px2;
        b.all.nimi{i} = A.nimi{b.ylenimi(i)};
        b.all.aania(i) = A.aania(b.ylenimi(i));
        b.all.puolue{i} = A.puolue{b.ylenimi(i)};
        
        jiiS.vasenoikea(b.ylenimi(i)) = px1;
        jiiS.libkon(b.ylenimi(i)) = px2;
        
        
        rahaa = A.rahaa(b.ylenimi(i));
        
        if p.kansalaistarkkailu == 1
            % A.nimi(b.ylenimi(487)) esim. timosoini
            aania = A.aania(b.ylenimi(i));
        else
            
            if p.kaikkiehdokkaat == 1
                aania = 1;
            else
                % pelkastaan valitut (muille = 0)..
                aania = A.valitaan(b.ylenimi(i))==2;  % VALITTU, muuten 0
            end
        end    
            if A.valitaan(b.ylenimi(i))==2
                % otetaan ylös ehdokkaan koordinaatit
                kansislaskuri = kansislaskuri + 1;
                b.kansanedustaja.X(kansislaskuri) = px1;
                b.kansanedustaja.Y(kansislaskuri) = px2;
                b.kansanedustaja.nimi{kansislaskuri} = hs.Nimi{b.b2hs(i)};
                b.kansanedustaja.approksimoitu(kansislaskuri) = 0;
                
                b.kansanedustaja.puolue{kansislaskuri} = A.puolue{b.ylenimi(i)};
                % nämä lähinnä tarkistusta varten
                b.kansanedustaja.aania(kansislaskuri) = A.aania(b.ylenimi(i));
                b.kansanedustaja.uusi(kansislaskuri) = A.uusi(b.ylenimi(i));
                
                b.kansanedustaja.puolueN(kansislaskuri) = b.puolue(i);
                
                
                b.kansanedustaja.puhis(kansislaskuri) = 0;
                
                for ppp = 1:length(p.puhislista) 
                    if strcmp(hs.Nimi{b.b2hs(i)},p.puhislista{ppp})
                        % turhaa periaatteessa
                        puhislaskuri = puhislaskuri + 1;
                        b.puheenjohtaja.X(puhislaskuri) = px1;
                        b.puheenjohtaja.Y(puhislaskuri) = px2;
                        b.puheenjohtaja.nimi{puhislaskuri} = hs.Nimi{b.b2hs(i)};
                        
                        b.kansanedustaja.puhis(kansislaskuri) = 1;
                        
                        b.puheenjohtaja.puolue(puhislaskuri) = b.puolue(i);
                        
                    end
                end
                
                
            end
            
            
            
            
        %end
        
        %viimesinnimi = A.nimi{b.ylenimi(i)};
        %aania = 1;
        
        p.optimointihila = 0;  %110422 ei toimi, en l�yd� vikaa
        
        if p.optimointihila == 0
            % ******************************************************
            % jakautuu hilaan
            % ******************************************************
            p2e = p1'*p2;
            % pit�is tarkistaa onko sum,sum p2e = 1 
            p2d = p2e.*aania;
            b.hila = b.hila + p2d;
            b.rahahila = b.rahahila + p2e.*rahaa;
            
            b.puoluehila(:,:,b.puolue(i)) = b.puoluehila(:,:,b.puolue(i)) + p2d;
            
            
            
            % ******************************************************
            
            % ******************************************************
            % ********************** tarkastellaan samalla erikysymysten
            % painoa arvoalueella
            % ******************************************************
            for iiii=1:p.tarkasteltaviakysymyksia
                %if b.binaaritaulu(i,p.tarkasteltavat(iiii))~=0
                % 110530 - korjaus
                if b.binaaritaulu(i,p.tarkasteltavat(iiii))>0 
                    b.kysymyshila(:,:,iiii) = b.kysymyshila(:,:,iiii) + p2e;%.*b.binaaritaulu(i,6);
                end
                
            end
            
            %if b.binaaritaulu(i,6)~=0
                % painotus poistettu
             %   b.kysymyshila(:,:,1) = b.kysymyshila(:,:,1) + p2e;%.*b.binaaritaulu(i,6);
            %end
            %if b.binaaritaulu(i,7)~=0
             %   b.kysymyshila(:,:,2) = b.kysymyshila(:,:,2) + p2e;%.*b.binaaritaulu(i,7);
            %end
            % onko jokin ajatus j��nyt kesken?
            %if 
            b.vastausmaarahila = b.vastausmaarahila + p2e;
            % b.hila
            % ******************************************************
            
        else
            % pakollinen, jos kuvia aiotaan k�ytt��
            % kesken ja ei toimi
            
            p3 = p1(p1>0.01);
            p4 = p2(p2>0.01);
    
            px3 = px1-round((length(p3)-1)/2);
            px4 = px2-round((length(p4)-1)/2);
            %px3a = px1-(length(p3)-1)/2;
            %px4a = px2-(length(p4)-1)/2;
            
            if px3+length(p3)>p.hilakoko || px4+length(p4)>p.hilakoko || px3<1 || px4<1
                % ylilaidan
            else
                p2d = p3'*p4.*aania;
                b.hila(px3:px3+length(p3)-1,px4:px4+length(p4)-1) = b.hila(px3:px3+length(p3)-1,px4:px4+length(p4)-1) + p2d;
                b.puoluehila(px3:px3+length(p3)-1,px4:px4+length(p4)-1,b.puolue(i)) = b.puoluehila(px3:px3+length(p3)-1,px4:px4+length(p4)-1,b.puolue(i)) + p2d;
            end
            
        end
        
    else
        aania = 0;   % miks tää täällä 30.5.2011? Selvitä:
    end
    
end

% ******************************************************
% halutaanko skaalata? ('-' vasttauksilla) 110424 - skaalaa tasaisesti)
% Hetkinen, ainakin ruotsi ja homokysss�riss� kaikka vastanneet, oliko
% raja, ett� kaikki on vastannut kaikkiin kyss�reihin. Edelleen osa
% ehdokkaista j�� huomioimatta, joten siit� virhett�.

% mieti onko tarpeen?!
%b.vastausmaarahila = b.vastausmaarahila/(sum(sum(b.vastausmaarahila)))*length(b.binaaritaulu);
% ******************************************************

% ******************************************************
% mika puolue suurin tietyssa hilapisteessa?
% ******************************************************
b.varihila = zeros(p.hilakoko,p.hilakoko);
b.hilanvari = zeros(p.hilakoko,p.hilakoko);
if p.pelkkahallitus == 1
    b.hilapistepuolueelle = zeros(length(p.hallitusneuvotteluissa),1);
    for i=1:length(p.hallitusneuvotteluissa)
        b.hilahallituspuolue(i) = find(strcmp(m.puolueet,p.hallitusneuvotteluissa{i}));
    end
end
    
% end tässä end?
for i=1:p.hilakoko
    for j=1:p.hilakoko
        
        if p.pelkkahallitus == 0
            b.varihila(i,j) = max(b.puoluehila(i,j,:));
            %b.hilanvari(i,j) = find(b.puoluehila(i,j,:)==b.varihila(i,j));
            apu = find(b.puoluehila(i,j,:)==b.varihila(i,j));
        else
            
            for i2=1:length(p.hallitusneuvotteluissa)
                b.hilapistepuolueelle(i2) = b.puoluehila(i,j,b.hilahallituspuolue(i2));
            end
            
            b.varihila(i,j) = max(b.hilapistepuolueelle);
            
            % *****************
            % IDA idea
            % 110506 TÄNNE vois koodata puolueraja, jos toka suurin lähes
            % yhtä suuri, niin valkoiseksi b.varihila siinä kohtaa
            
            apu = b.hilahallituspuolue(b.hilapistepuolueelle==b.varihila(i,j));
            
        end
        
        b.hilanvari(i,j) = apu(1);  % voi l�yty� montakohtaa
        
        
        if p.kansalaistarkkailu == 0 && p.kuvaantyhjiaruutuja == 1
            if b.varihila(i,j)<0.00024   %1
                b.varihila(i,j) = 0;  % TEKEE valkoista missä ei ehdokasta
            end
        end
        
        
        
    end
end
% ******************************************************
disp('Ehdokkaat ja hilojen v�ri asetettu');

disp('approksimoidaan HS:n vaalikoneeseen vastaamatta j�tt�neet (tai esim veltto - tietokanta)');
% ******************************************************
%%
% LOGO
%p.puolueetL = {'-','IPU','KOK','KA','M2011','PS','PIR','KESK','SKP','KD','SDP','KTP','RKP','VP','VAS','VIHR'};
%p.puolueetL = {'MUUT','IPU','KOK','KA','M2011','PS','PIR','KESK','SKP','KD','SDP','KTP','RKP','VP','VAS','VIHR'};

% kaydaan kaikki lapi
% VAALIKONEESEEN
% SSP Suomen senioripuolue
testi.paaohjelmassanimialisatty = 0;

b.all.n = length(b.binaaritaulu);

for i=1:A.n
    %if sum(strcmp(hs.Nimi(b.b2hs),A.nimi{i}))>0
    if sum(b.ylenimi == i)>0 
        % *****************************
        % on merkattu hilaan
        % paitsi virheet
        % *****************************   
        testi.paaohjelmassanimialisatty = testi.paaohjelmassanimialisatty+1;
        % 1814 - eli kuusi nimeä eri tavalla Ylellä.
        %
        % HUOM! TESTATAAN SUMMALLA, samannimiset voi olla ongelma..
        % 110519
        
        
      
    else
        
        % Visualisointeja varten 28.9.2011
        b.all.n = b.all.n + 1;
        
        % **************************************************
        % 110506 tarkistetaanko onko puuttuva vaalipiirissä
        % **************************************************
        if p.kokosuomi == 0
            p.analyysiinmukaan = 0;
        else
            p.analyysiinmukaan = 1;
        end
        % **************************************************
        if p.analyysiinmukaan ==1
            muokataan = find(strcmp(m.puolueet,A.puolue{i}));
         
            lisataanko = 0;
            if length(muokataan) == 0 || muokataan == 30
            
                %disp(muokataan)
                
                %if strcmp(A.puolue{i},'M11')
                %    muokataan=    find(strcmp(m.puolueet,'M2011'));
                %    disp('TUUUUUR_HA');
            
                %elseif strcmp(A.puolue{i},'ITSP')
                    %    muokataan=    find(strcmp(m.puolueet,'IPU'));
                    %elseif strcmp(A.puolue{i},'KÖY')
                    %    muokataan=    find(strcmp(m.puolueet,'KA'));   
                    %elseif strcmp(A.puolue{i},'K�Y')
                    %    muokataan=    find(strcmp(m.puolueet,'KA'));  
                    %elseif strcmp(A.puolue{i},'SSP')
                        %    muokataan=    find(strcmp(m.puolueet,'SEN'));  
                if strcmp(A.nimi{i}(1:3),'Nau')
                    %disp('NAUK');
                
                    muokataan=    find(strcmp(m.puolueet,'RKP'));  
                    disp('naucler');
                    %sum(sum(b.puoluehila(:,:,muokataan)))
                    %muokataan = 18;
                    % naucler RKP:hen
                elseif strcmp(A.nimi{i}(1:8),'Viitanen') || strcmp(A.nimi{i}(1:8),'Rantanen')
                    muokataan = 29;
                    disp(strcat(A.nimi{i},'->Yl-sit'));
                %elseif strcmp(A.nimi{i}(1:3),'Kor') && strcmp(A.puolue{i},'MUUT')
                %    muokataan=    find(strcmp(m.puolueet,'Korhonen'));  
                %elseif strcmp(A.nimi{i}(1:3),'Val') && strcmp(A.puolue{i},'MUUT')
                %    muokataan=    find(strcmp(m.puolueet,'Valjakka'));   
                else
                    valilyonti = findstr(A.nimi{i},' ')-1;
                    sukunimipuolue = find(strcmp(m.puolueet,A.nimi{i}(1:valilyonti)));
                    if length(sukunimipuolue)>0 && strcmp(A.puolue{i},'MUUT')
                        muokataan = sukunimipuolue;
                        %disp('*');
                        disp(strcat(A.nimi{i},'->',num2str(muokataan)));
                    else
                        disp(strcat('30:',A.nimi{i}));
                    end
                    %disp('KOR');
                end
                %disp('=');
                %disp(num2str(muokataan))
                
                kokohilalisaysnyt = 0;
                if length(muokataan) == 0
                    kokohilalisaysnyt = 1;
                else
                    if sum(sum(b.puoluehila(:,:,muokataan)))==0
                        kokohilalisaysnyt = 1;
                    end
                end
                if kokohilalisaysnyt == 2    
                
                    if p.kansalaistarkkailu == 1
                        % A.nimi(b.ylenimi(487)) esim. timosoini
                        aania = A.aania(i);
                    else
                        % pelkastaan valitut (muille = 0)..
                        aania = A.valitaan(i)==2;
                    end
                
                    
                    rahaa = A.rahaa(i);
                    
                    
                    % vois lisata hilaan tasaisesti.
                    lisataanhilaan = b.hila*(aania/sum(sum(b.hila))); %
                    %TARKISTA ONKO OIKEIN
                    b.hila = b.hila + lisataanhilaan;  % pelk�st��n t�h�n
                    
                    b.rahahila = b.rahahila + lisataanhilaan/aania*rahaa;
                    
                    disp(strcat('Vain Bhila:',A.nimi{i},',',A.puolue{i}))
                    %disp('ei lis�tty puoluehilaan, vain kokohilaan');
                    lisataanko = 1;
                    
                    disp('MOOO');
                    
                    
                
                    
                    %jiiS.vasenoikea((i)) = 100; % VALKOINEN JS
                    %jiiS.libkon((i)) = 100; % VALKOINEN JS
                    
                    
                    
                    
                end
            end
            
            if lisataanko == 0
                %lisataanhilaan = b.hila*(A.aania(i)/sum(sum(b.hila))); % TARKISTA ONKO OIKEIN
                % Puolueen alueella
                if p.kansalaistarkkailu == 1
                    % A.nimi(b.ylenimi(487)) esim. timosoini
                    aania = A.aania(i);
                else
                     % pelkastaan valitut (muille = 0)..
                     if p.kaikkiehdokkaat == 1
                        aania = 1;
                     else
                        % pelkastaan valitut (muille = 0)..
                        %aania = A.valitaan(b.ylenimi(i))==2;  % VALITTU, muuten 0
                        aania = A.valitaan(i)==2;
                     end
                end
            
            
                rahaa = A.rahaa(i);
                
                if sum(sum(b.puoluehila(:,:,muokataan)))>0
                    lisataanhilaan = b.puoluehila(:,:,muokataan)*(aania/sum(sum(b.puoluehila(:,:,muokataan)))); % TARKISTA ONKO OIKEIN
                    b.hila = b.hila + lisataanhilaan;
                    b.puoluehila(:,:,muokataan) = b.puoluehila(:,:,muokataan) + lisataanhilaan;
                
                    
                    b.rahahila = b.rahahila + lisataanhilaan/aania*rahaa;
                    
                    %disp('JEEE');
                    
                    jiiS.vasenoikea((i)) = find(max(b.puoluehila(:,:,muokataan)') == max(max(b.puoluehila(:,:,muokataan))));
                    jiiS.libkon((i)) = find(max(b.puoluehila(:,:,muokataan)) == max(max(b.puoluehila(:,:,muokataan))));
                    
                    
                    b.all.X(b.all.n) = find(max(b.puoluehila(:,:,muokataan)') == max(max(b.puoluehila(:,:,muokataan))));
                    b.all.Y(b.all.n) = find(max(b.puoluehila(:,:,muokataan)) == max(max(b.puoluehila(:,:,muokataan))));
                    b.all.nimi{b.all.n} = strcat('*',A.nimi{i});
                    b.all.aania(b.all.n) = A.aania(i);
                    b.all.puolue{b.all.n} = A.puolue{i};
                    
                    if A.valitaan(i)==2
                        %if aania == 1 && p.kansalaistarkkailu == 0
                        % otetaan ylös ehdokkaan koordinaatit
                        kansislaskuri = kansislaskuri + 1;
                        b.kansanedustaja.X(kansislaskuri) = find(max(b.puoluehila(:,:,muokataan)') == max(max(b.puoluehila(:,:,muokataan))));
                        b.kansanedustaja.Y(kansislaskuri) = find(max(b.puoluehila(:,:,muokataan)) == max(max(b.puoluehila(:,:,muokataan))));
                        b.kansanedustaja.nimi{kansislaskuri} = A.nimi{i};%hs.Nimi{b.b2hs(i)};
                        b.kansanedustaja.approksimoitu(kansislaskuri) = 1;
                        
                        b.kansanedustaja.puolue{kansislaskuri} = A.puolue{i};
                        % nämä lähinnä tarkistusta varten
                        b.kansanedustaja.aania(kansislaskuri) = A.aania(i);
                        b.kansanedustaja.uusi(kansislaskuri) = A.uusi(i);
                        
                        b.kansanedustaja.puhis(kansislaskuri) = 0;
                        
                        b.kansanedustaja.puolueN(kansislaskuri) = muokataan;
                        
                        for ppp = 1:length(p.puhislista) 
                            if strcmp(A.nimi{i},p.puhislista{ppp})
                                % turhaa periaatteessa
                                puhislaskuri = puhislaskuri + 1;
                                b.puheenjohtaja.X(puhislaskuri) = find(max(b.puoluehila(:,:,muokataan)') == max(max(b.puoluehila(:,:,muokataan))));
                                b.puheenjohtaja.Y(puhislaskuri) = find(max(b.puoluehila(:,:,muokataan)) == max(max(b.puoluehila(:,:,muokataan))));
                                b.puheenjohtaja.nimi{puhislaskuri} = A.nimi{i};%hs.Nimi{b.b2hs(i)};
                        
                                b.kansanedustaja.puhis(kansislaskuri) = 1;
                                
                                b.puheenjohtaja.puolue(puhislaskuri) = muokataan;
                            end
                        end
                        
                        
                    end
                else
                    
                    
                     if p.kaikkiehdokkaat == 1 || (A.valitaan(i)==2)
                  
                    
                    disp(muokataan)
                    % ei taida t�� toimia?!?!? 110428 turha muutenkin
                    disp('POOOOOhelminen'); % ei vastauksia - jakauma tasasesti nolla ilman t�t�
                    lisataanhilaan = b.hila(:,:)*(aania/sum(sum(b.hila(:,:)))); % TARKISTA ONKO OIKEIN
                    b.hila = b.hila + lisataanhilaan;
                    b.puoluehila(:,:,muokataan) = b.puoluehila(:,:,muokataan) + lisataanhilaan;
                    
                    
                    b.all.X(b.all.n) = find(max(b.hila(:,:)') == max(max(b.hila(:,:))));
                    b.all.Y(b.all.n) = find(max(b.hila(:,:)) == max(max(b.hila(:,:))));
                    b.all.nimi{b.all.n} = strcat('**',A.nimi{i});
                    b.all.aania(b.all.n) = A.aania(i);
                    b.all.puolue{b.all.n} = A.puolue{i};
                    
                    
                     end
       
                end
            end
        end
    end
end


%%

for i=1:length(b.kansanedustaja.nimi)
    apu =  find(strcmp(hs.Nimi,b.kansanedustaja.nimi{i}));
    if length(apu) == 1
        b.kansanedustaja.hsN(i) = apu;
    else
        if length(apu)==0
            apu = nan; % naucler ei löydy hs tietokannasta ollenkaan
        else
            %disp(length(apu))
            apu = apu(find(strcmp(hs.Puolue(apu),b.kansanedustaja.puolue{i})));
            disp(apu);
        end
        
    end
    b.kansanedustaja.hsN(i) = apu;
    
end

for i=1:length(b.kansanedustaja.nimi)
    if isnan(b.kansanedustaja.hsN(i))
        % naucler
        b.kansanedustaja.ammatti{i} = '?';
        b.kansanedustaja.koulutustaso{i} = '?';
        % ... HS:stä
    else
        b.kansanedustaja.ammatti{i} = hs.Ammatti{b.kansanedustaja.hsN(i)};
        b.kansanedustaja.koulutustaso{i} = hs.Koulutustaso{b.kansanedustaja.hsN(i)};
        % ... HS:stä
    end
end




%%
% laitetaan varit kuosiin 
Zv = zeros(p.hilakoko,p.hilakoko,3);


%p.puoluereunat = 1;
if p.puoluereunat == 1
   for i=1:length(m.puolueet)
        if p.vakiopuoluereuna == 0
            b.puoluereuna(i) = mean(mean(b.puoluehila(:,:,i)))*1.8;
        else
            b.puoluereuna(i) = p.vakiopuoluereuna;
        end
        
        %b.puoluereuna(i) = median(median(b.puoluehila(:,:,i)));
   end
end
    


for i=1:p.hilakoko
    for j=1:p.hilakoko
        %apu = find(strcmp(p.puolueetL, m.puolueet{b.hilanvari(i,j)}));
        apu = b.hilanvari(i,j); %110426 new logo
        %apu=6;
        if length(apu)>0 && p.logotkarttaan == 1
            % JUMAZUIKA - rgb siis 0..1, ja kuva0..255
            % voisko videotallennuskin johtua siit�? no selvit� joskus
            % tuli vain valkoista kuvaa 255,255,255 matriisia..
            
            Zv(i,j,:) = (double(logot{apu}(mod((80-i)*1,79)+1,mod((j)*1,79)+1,:)))/255;
            %for k=1:3
            %    Zv(i,j,4-k) = double(logot{apu}(mod((80-i)*1,79)+1,mod((j)*1,79)+1,k));
            %end
            
            
            % *************** tummenntaa valkoisia kohtia logossa ***
            % ei tarpeen omien logojen j�lkeen
            %if sum(Zv(i,j,:))>2.5
                % vois tehd� kans ehdon, ett� jos ollaan logon laidalla,
                % tai sitten omat logot photarilla
                % esim. V�ri ja teksti - tai palattia.
                % esim. nyt valkonen pois SDP:st�
            %    Zv(i,j,:) = p.vari{b.hilanvari(i,j)};
            %end
            % *****************************************
            
           
            
        else
            
            % ilman varjostusta
            %Zv(i,j,:) = p.vari{b.hilanvari(i,j)};
            
            
            Zv(i,j,:) = p.vari{b.hilanvari(i,j)};%+varinvarjostus*2;
          
        end
        
        if p.varinvarjostus == 1
            
              % Google matlab transperency example - tai opacity
            varinvarjostus = 1 - b.varihila(i,j)/sum(b.puoluehila(i,j,:));
            Zv(i,j,:) = Zv(i,j,:) + varinvarjostus*(3-sum(Zv(i,j,:)))/2; % vaaleita puolueita varjostetaan v�hemm�n
            for izz = 1:3
                if Zv(i,j,izz)>1
                    Zv(i,j,izz) = 1;
                end
            end
            
        end
        
       
        
    end
end


% ***************************** Reunat koodattu vain hallitustan varten
if p.puoluereunat == 1
    if p.pelkkahallitus == 1
        for i=1:p.hilakoko
            for j=1:p.hilakoko
                for apu = b.hilahallituspuolue
                    %apu = b.hilanvari(i,j); %110426 new logo
                    if abs(b.puoluehila(i,j,apu) - b.puoluereuna(apu))<b.puoluereuna(apu).*p.puoluereunanpaksuus
                        Zv(i,j,:) = p.vari{apu}';%*0.5;
                    end
                    
                    % toinen korkeuskäyrä
                    if p.toinenkorkeuskayra>0
                    if abs(b.puoluehila(i,j,apu) - b.puoluereuna(apu)*p.toinenkorkeuskayra)<b.puoluereuna(apu).*p.puoluereunanpaksuus
                        Zv(i,j,:) = p.vari{apu}';%*0.5;
                    end
                    end
                    
                    
                end
            end
        end
    end
end

% on ok kans
%%
% ************************ yhden puolueen piirto
% vois piirt�� laajemmalla, siten muokataan ZV:t� siten, ett� jos v�h�n ees
% kannatusta, niin omalla v�rill�.


% IDA KORJAUKSIA 110512
p.puolueetYLEtulosUnique = unique(A.puolue);

%%
% ******************************************************************
% FIG100: IDA figures 1-3 % report fig 0
% ******************************************************************
cd(p.f);IDAfigure(p,b,A,hs,m,SCORE,Zv,p.idapaperfigure);cd ..
% ******************************************************************
% posteriin 26.10.2011
%%
%cd(p.f);IDAfigure(p,b,A,hs,m,SCORE,Zv,0);cd ..



% **********************************************
% Print posterpic
% käännä käsin sopivaan kulmaan ja kopioi komentoriville
% Print posterpic - ei joka kerta- kun onnistunut kuva, niin 0
p.idaposterprint =0;
% **********************************************
if p.idaposterprint == 1
    
     p.tallennusX = 8;
     p.tallennusY = 6;
    
     view(-123,28)
     grid
     
     
     set(gcf,'PaperUnits','inches','PaperPosition',[0 0 p.tallennusX p.tallennusY]);
            
     % HUOM! ZEKKAA TIEDOSTO
     tallennustiedosto = strcat('posterpics/','k100I1','.png');
     print('-dpng','-r300',tallennustiedosto);% pngWille/testi111015.png -r150
    
     disp('tuhatta ääntä KESK');
     sum(sum(b.puoluehila(:,:,5)))/1000
     
     
end
% **********************************************

%%

% ******************************************************************
% FIG101: Finnish Voters
% ******************************************************************
cd(p.f);Votersfigure(p,b,A,hs,m,SCORE,Zv);cd ..
% ******************************************************************

%%
% **********************************************
% Print posterpic
% käännä käsin sopivaan kulmaan ja kopioi komentoriville
% Print posterpic - ei joka kerta- kun onnistunut kuva, niin 0
p.idaposterprint =0;
% **********************************************
if p.posterprint == 1
    
     p.tallennusX = 8;
     p.tallennusY = 6;
    
     view(-99,64)
     grid
     
     
     set(gcf,'PaperUnits','inches','PaperPosition',[0 0 p.tallennusX p.tallennusY]);
            
     % HUOM! ZEKKAA TIEDOSTO
     tallennustiedosto = strcat('posterpics/','k101','.png');
     print('-dpng','-r300',tallennustiedosto);% pngWille/testi111015.png -r150
    
end
% **********************************************



%%

%%


% KYll� ei.. esim. homo..
figure(444);
Z = b.kysymyshila(:,:,1)./b.vastausmaarahila;
%Z = b.vastausmaarahila;
%surface(-p.hilaskaala:p.hilaskaala,-p.hilaskaala:p.hilaskaala,Z,Zv,'edgecolor',[0.6 0.6 0.6]);
meshc(-p.hilaskaala:p.hilaskaala,-p.hilaskaala:p.hilaskaala,Z);
%surfc(-p.hilaskaala:p.hilaskaala,-p.hilaskaala:p.hilaskaala,Z,Zv,'edgecolor',[0.6 0.6 0.6]);
%title(strcat(hs.k.Kysymys{b.kysymysnumero(6)},':',b.vastaus{6}));

%Z = b.kysymyshila(:,:,2)./b.vastausmaarahila;

%surface(-p.hilaskaala:p.hilaskaala,-p.hilaskaala:p.hilaskaala,Z,Zv,'edgecolor',[0.3 0.3 0.3]);
%title(strcat(hs.k.Kysymys{(p.tarkasteltavakysymys)}));

%%
p.hilavastauspiirretaan = 1;
figure(445);
Z = b.kysymyshila(:,:,p.hilavastauspiirretaan)./b.vastausmaarahila;
%Z = b.vastausmaarahila;
%shading interp
surface(-p.hilaskaala:p.hilaskaala,-p.hilaskaala:p.hilaskaala,Z,Zv,'edgecolor',[0.5 0.5 0.5]);
title(strcat(hs.k.Kysymys{p.tarkasteltavakysymys},':',b.vastaus{p.tarkasteltavat(p.hilavastauspiirretaan)}));

%shading faceted

%%
% tilavuus on 1 = = TOIMII!!!! otettu tosiaan pois painoarvo t�ss�.. '-'
% vastaus on siis suhteutettu..
p.hilavastauspiirretaan = 1;
figure(446);
Z = -b.kysymyshila(:,:,1);
for i=2:2
    Z = Z + b.kysymyshila(:,:,i);
end


% -1..1 hila kantaan
Z = Z./b.vastausmaarahila;

% kerrotaan populalla
Z = Z.*b.hila;
%Z = b.vastausmaarahila;
%shading interp
surface(-p.hilaskaala:p.hilaskaala,-p.hilaskaala:p.hilaskaala,Z,Zv,'edgecolor',[0.5 0.5 0.5]);
title(strcat(hs.k.Kysymys{p.tarkasteltavakysymys},':',b.vastaus{p.hilavastauspiirretaan}));
%%
% TÄHÄN KUVAAN PUOLUIDEN MOODIT (Ilman hlömaaraa, eli mielipide)
p.hilavastauspiirretaan = 1;
figure(448);
%Z = b.kysymyshila(:,:,p.hilavastauspiirretaan)./b.vastausmaarahila.*b.hila;
Z = b.kysymyshila(:,:,p.hilavastauspiirretaan)./b.vastausmaarahila.*b.hila;
%Z = Z - b.kysymyshila(:,:,1)./b.vastausmaarahila;
%Z = b.vastausmaarahila;
%shading interp
%surface(-p.hilaskaala:p.hilaskaala,-p.hilaskaala:p.hilaskaala,Z,Zv,'edgecolor',[0.5 0.5 0.5]);
%surfc(-p.hilaskaala:p.hilaskaala,-p.hilaskaala:p.hilaskaala,Z,Zv,'edgecolor',[0.5 0.5 0.5]);
%meshc(-p.hilaskaala:p.hilaskaala,-p.hilaskaala:p.hilaskaala,Z);
contour(-p.hilaskaala:p.hilaskaala,-p.hilaskaala:p.hilaskaala,Z);
title(strcat(hs.k.Kysymys{p.tarkasteltavakysymys},':',b.vastaus{p.tarkasteltavat(p.hilavastauspiirretaan)}));
disp(strcat(hs.k.Kysymys{p.tarkasteltavakysymys},':',b.vastaus{p.tarkasteltavat(p.hilavastauspiirretaan)}));
disp((b.vastaus{p.tarkasteltavat(p.hilavastauspiirretaan)}));
disp(sum(sum(Z)));


%%
hold off;
xlabel('Arvo');
ylabel('raha');
figure(200);
%subplot(2,2,4);

% mieti viel� kahta eri erottelijaa... konerv. vs. lib

hold on;
for i=1:length(m.puolueEhdokkaita)
    if m.puolueEhdokkaita(i)>2
        
        
    edust = find(hs.Kansanedustaja(m.puolueN{i}));
    if length(edust)>0
        for j=1:length(edust)  % HUOM! tyhji� vastauksia poistettu
            plot(SCORE(edust(j),1),SCORE(edust(j),2),'.','color',p.vari{i}); 
        end
    end
    % ehdokkaat
    kx = mean(SCORE(m.puolueN{i},1));
    ky = mean(SCORE(m.puolueN{i},2));
    plot(kx,ky,'o','color',p.vari{i}); 
    text(kx,ky,m.puolueet{i});
    end
end
hold off;
xlabel('Arvo');
ylabel('raha');

%%


% testailuu:
p.pcakys = [1:5];
b.testi = [];
for i=1:length(p.pcakys)
    b.testi = [b.testi find(b.kysymysnumero==p.pcakys(i))];
end
[CO,SC,la,ts] = princomp(b.binaaritaulu(:,b.testi));

SC2 = SC(:,1) - min(SC(:,1));
SC2 = SC2/max(SC2);

figure(300);
subplot(2,2,1);
plot(SC(:,1),SC(:,2),'.');
subplot(2,2,2);
bar(CO);
%disp(hs.k.Kysymys{p.pcakys});
%disp(b.v{p.pcakys}.v');
subplot(2,2,3);
bar(la/sum(la)*100);
subplot(2,2,4);
plot(ts);

piirra = find(m.puolueEhdokkaita>60); % isot puolueet 12kpl


figure(301);
for i=1:12
    subplot(4,3,i);
    hist(SC2(m.puolueN{piirra(i)}'),0:1:1);
    title(m.puolueet{piirra(i)});
    axis tight;
end




%%
% *******************************************
% t�st� alkaa histogrammit
% *******************************************
figure(500);
kyssari = 23;
clear vastaukset
clear apu
clear e2011
edustajat2007 = find(hs.Kansanedustaja==1);

e2007 = hs.k.Vastaus(edustajat2007,kyssari);

e2012 = hs.k.Vastaus(:,kyssari);
edustajat2011 = 1:2306;
vastaukset = unique(hs.k.Vastaus(:,kyssari));
jakauma = zeros(length(vastaukset),1);
jakauma2 = zeros(length(vastaukset),1);
jakauma3 = zeros(length(vastaukset),1);
jakauma2007hallitus = zeros(length(vastaukset),1);
jakauma2007kok = zeros(length(vastaukset),1);
jakauma2007muut = zeros(length(vastaukset),1);
koklaskuri = 0;
for i=1:length(vastaukset)
    jakauma(i) = sum(strcmp(e2007,vastaukset{i}));
    for j=1:length(edustajat2007)
        if strcmp(hs.Puolue{edustajat2007(j)},'KOK') && strcmp(e2007{j},vastaukset{i})
            koklaskuri = koklaskuri + 1;
            jakauma2007kok(i) = jakauma2007kok(i) + 1;
        end
        if strcmp(hs.Puolue{edustajat2007(j)},'KESK') && strcmp(e2007{j},vastaukset{i})
            koklaskuri = koklaskuri + 1;
            jakauma2007muut(i) = jakauma2007muut(i) + 1;
        end
        if strcmp(hs.Puolue{edustajat2007(j)},'VIHR') && strcmp(e2007{j},vastaukset{i})
            koklaskuri = koklaskuri + 1;
            jakauma2007muut(i) = jakauma2007muut(i) + 1;
        end
        if strcmp(hs.Puolue{edustajat2007(j)},'RKP') && strcmp(e2007{j},vastaukset{i})
            koklaskuri = koklaskuri + 1;
            jakauma2007muut(i) = jakauma2007muut(i) + 1;
        end
    end
end
jakauma2007muut = (jakauma2007muut + jakauma2007kok)/sum(jakauma2007muut + jakauma2007kok)*100;

jakauma2011hallitus = zeros(length(vastaukset),1);
jakauma2011kok = zeros(length(vastaukset),1);
jakauma2011muut = zeros(length(vastaukset),1);
koklaskuri = 0;
for i=1:length(vastaukset)
    jakauma(i) = sum(strcmp(e2012,vastaukset{i}));
    for j=1:length(edustajat2011)
        if strcmp(hs.Puolue{edustajat2011(j)},'KOK') && strcmp(e2012{j},vastaukset{i})
            koklaskuri = koklaskuri + 1;
            jakauma2011kok(i) = jakauma2011kok(i) + 1;
        end
        if strcmp(hs.Puolue{edustajat2011(j)},'SDP') && strcmp(e2012{j},vastaukset{i})
            koklaskuri = koklaskuri + 1;
            jakauma2011muut(i) = jakauma2011muut(i) + 1;
        end
        if strcmp(hs.Puolue{edustajat2011(j)},'PS') && strcmp(e2012{j},vastaukset{i})
            koklaskuri = koklaskuri + 1;
            jakauma2011muut(i) = jakauma2011muut(i) + 1;
        end
        %if strcmp(hs.Puolue{edustajat2011(j)},'RKP') && strcmp(e2012{j},vastaukset{i})
        %    koklaskuri = koklaskuri + 1;
        %    jakauma2011muut(i) = jakauma2011muut(i) + 1;
        %end
    end
end
jakauma2011muut = (jakauma2011muut + jakauma2011kok)/koklaskuri*100;


e2011.puolueet = {'KOK','SDP','PS','KESK','VAS','VIHR','RKP','KD'}
e2011.maarat = [20.2 19.1 19.1 15.9 8.2 7.2 4.3 4.1];



% VIRHE KOODISSA..  TEE viel� eri vaihtoehdoista statsit.

for i=1:8
    nn = find(strcmp(m.puolueet,e2011.puolueet{i}))
    apu = hs.k.Vastaus(m.puolueN{nn},kyssari);
    
    for i=1:length(vastaukset)
        jakauma3(i) = sum(strcmp(apu,vastaukset{i}));
        
    end
    jakauma3 = jakauma3 * e2011.maarat(i);
    jakauma2 = jakauma2 + jakauma3;
end


%jakauma2 = jakauma2./sum(jakauma2)*200;
%jakauma = jakauma./sum(jakauma)*200;

jakauma2 = jakauma2./sum(jakauma2)*100;
jakauma = jakauma./sum(jakauma)*100;


barh([jakauma2007muut jakauma jakauma2 jakauma2011muut jakauma2011kok]);

apu = hs.k.Kysymys{kyssari};
if length(apu)>300
    apu = {apu(1:150),apu(151:300),apu(301:end)}
elseif length(apu)>200
    apu = {apu(1:200),apu(61:end)}
end
title(apu);
set(gca,'YTickLabel',vastaukset)
legend('2007 hallitus','2007 eduskunta','2011 eduskunta',' 2011 Kok','2011 KOK, SDP, PS');
xlabel('%');


%%

aa.kysymys = 1;

aa.edustajat{1} = find(hs.Kansanedustaja==1);   % indeksi 1..2306



e2007 = hs.k.Vastaus(edustajat2007,kyssari);

e2012 = hs.k.Vastaus(:,kyssari);
edustajat2011 = 1:2306;
vastaukset = unique(hs.k.Vastaus(:,kyssari));
jakauma = zeros(length(vastaukset),1);
jakauma2 = zeros(length(vastaukset),1);
jakauma3 = zeros(length(vastaukset),1);
jakauma2007hallitus = zeros(length(vastaukset),1);
jakauma2007kok = zeros(length(vastaukset),1);
jakauma2007muut = zeros(length(vastaukset),1);
koklaskuri = 0;

%%
% 110418
clear all;
close all;
load hs.mat;
%%

%fid = fopen('tulos2011e2.csv');
fid = fopen('tulos2011eXP.csv');
C = textscan(fid,'%s%s%s%s','Delimiter','%','EmptyValue',-Inf);
    
fclose(fid);


k=0;
for i=1:200
    edustaja2011 = find(strcmp(hs.Nimi',C{1}(i)));
    %disp(length(edustaja2011))
    if length(edustaja2011)==0
        disp('Puuttuu:');
        disp(C{1}(i));
        apu = C{1}(i);
        if find(strcmp(apu{1}(1),'E'))
            edustaja2011 =  find(strcmp(hs.Sukunimi,'Elomaa')); % Kike
            disp('kike lisatty');
        end
        if find(strcmp(apu{1}(1),'V'))
            edustaja2011 =  find(strcmp(hs.Sukunimi,'Virtanen')); % Kike
            disp('veltto lisatty');
        else
            % naucler puuttuu..
            disp(apu)
        end
        
    end
    if length(edustaja2011)==2
        disp('Tupla:');
        disp(C{1}(i));
        edustaja2011 = edustaja2011(1);  % tuplanimi, kouvolasta pääs
    end
    
    if length(edustaja2011) == 1
        %disp(k)
        k=k+1;
        tt.edustaja2011(k) = edustaja2011;  % tässä hs-rivinumerot
    end
    
end

tt.edustaja2007 = find(hs.Kansanedustaja==1);  % esim. puolue (hs.Puolue(tt.edustaja20xx))


%%

tt.kysymys = 30;
p.jakokoko = 240;
p.kysymysjako = 267;

tt.vastausvaihtoehdot = unique(hs.k.Vastaus(:,tt.kysymys));


tt.vastauspuoluejako = {'KOK', 'SDP' ,'PS' ,'KESK','VAS' ,'VIHR' ,'RKP','KD' ,'MUUT'};  % hallitukset & mahdollinen
tt.puoluecolor = [0 0 1; 1 0 0; 0 0 0.5; 0 0.6 0; 0.7 0 0; 0 1 0;0.7 0.7 0.2;0.6 0.6 1;0.7 0.7 0.7];
tt.vastaus{1} = zeros(length(tt.vastausvaihtoehdot),length(tt.vastauspuoluejako)); %{vuosi}  6 = [KOK SDP PS KESK VIHR RKP MUUT]
for i=1:length(tt.edustaja2011)
    vastaus = find(strcmp(tt.vastausvaihtoehdot,hs.k.Vastaus(tt.edustaja2011(i),tt.kysymys))); % 1..vastausvaihtoehdot
    loyty = 0;
    for j=1:length(tt.vastauspuoluejako)-1
        if strcmp(hs.Puolue{tt.edustaja2011(i)},tt.vastauspuoluejako{j})
            loyty = j;
        end
    end
    if loyty == 0
        loyty = length(tt.vastauspuoluejako);
    end
    tt.vastaus{1}(vastaus,loyty) = tt.vastaus{1}(vastaus,loyty) + 1;
end

% cp
tt.vastaus{2} = zeros(length(tt.vastausvaihtoehdot),length(tt.vastauspuoluejako)); %{vuosi}  6 = [KOK SDP PS KESK VIHR RKP MUUT]
%tt.vastauspuoluejako = {'KOK', 'SDP' ,'PS' ,'KESK' ,'VIHR' ,'RKP' ,'MUUT'};  % hallitukset & mahdollinen
for i=1:length(tt.edustaja2007)
    vastaus = find(strcmp(tt.vastausvaihtoehdot,hs.k.Vastaus(tt.edustaja2007(i),tt.kysymys))); % 1..vastausvaihtoehdot
    loyty = 0;
    for j=1:length(tt.vastauspuoluejako)-1
        if strcmp(hs.Puolue{tt.edustaja2007(i)},tt.vastauspuoluejako{j})
            loyty = j;
        end
    end
    if loyty == 0
        loyty = length(tt.vastauspuoluejako);
    end
    tt.vastaus{2}(vastaus,loyty) = tt.vastaus{2}(vastaus,loyty) + 1;
end


tt.vastaus{1} = tt.vastaus{1}/length(tt.edustaja2011)*100;
tt.vastaus{2} = tt.vastaus{2}/length(tt.edustaja2007)*100;


%tt.puoluecolor = [0 0 1;1 0 0;0 0 0.5;0 0.6 0;0 1 0;0.7 0.7 0;0.7 0.7 0.7];
colormap(tt.puoluecolor);

apu = strcat('2011: ',hs.k.Kysymys{tt.kysymys});
if length(apu)>p.kysymysjako*2
    apu = {apu(1:p.kysymysjako),apu(p.kysymysjako+1:p.kysymysjako*2),apu(p.kysymysjako*2+1:end)}
elseif length(apu)>p.kysymysjako
    apu = {apu(1:p.kysymysjako),apu(p.kysymysjako+1:end)}
end

figure(1000);

subplot(2,1,1);
%set(gca,'position',[0.033 0.02 0.96 0.955])
set(gca,'position',[0.033 0.53 0.94 0.43])  % [x x lev
barh(tt.vastaus{1},'stacked');
title(apu);
set(gca,'YTickLabel',{})
legend(tt.vastauspuoluejako);
%xlabel('%');



for i=1:length(tt.vastausvaihtoehdot)
    apu2 = (sum(tt.vastaus{1}'));
    apu = tt.vastausvaihtoehdot{i};
    
    if length(apu)>p.jakokoko*2
        apu = {apu(1:p.jakokoko),apu(p.jakokoko+1:p.jakokoko*2),apu(p.jakokoko*2+1:end)}
    elseif length(apu)>p.jakokoko
        apu = {apu(1:p.jakokoko),apu(p.jakokoko+1:end)}
    end
    
    text(apu2(i)+1,i,apu);
end


apu = strcat('2007: ',hs.k.Kysymys{tt.kysymys});
if length(apu)>p.kysymysjako*2
    apu = {apu(1:p.kysymysjako),apu(p.kysymysjako+1:p.kysymysjako*2),apu(p.kysymysjako*2+1:end)}
elseif length(apu)>p.kysymysjako
    apu = {apu(1:p.kysymysjako),apu(p.kysymysjako+1:end)}
end
axis([0 100 0.5 length(tt.vastausvaihtoehdot)+0.5])
subplot(2,1,2);

set(gca,'position',[0.033 0.04 0.94 0.43]) 
barh(tt.vastaus{2},'stacked');
title(apu);
%set(gca,'YTickLabel',tt.vastausvaihtoehdot)
set(gca,'YTickLabel',{})
legend(tt.vastauspuoluejako);
axis([0 100 0.5 length(tt.vastausvaihtoehdot)+0.5])
xlabel('%');
for i=1:length(tt.vastausvaihtoehdot)
    apu2 = (sum(tt.vastaus{2}'));
    
    apu = tt.vastausvaihtoehdot{i};
    
    if length(apu)>p.jakokoko*2
        apu = {apu(1:p.jakokoko),apu(p.jakokoko+1:p.jakokoko*2),apu(p.jakokoko*2+1:end)}
    elseif length(apu)>p.jakokoko
        apu = {apu(1:p.jakokoko),apu(p.jakokoko+1:end)}
    end
    
    text(apu2(i)+1,i,apu);
end

%%
% HALLITuS yll� alla OPPOSITIO
tt.kysymys = 29;
p.jakokoko = 190;
p.kysymysjako = 260;

tt.vastausvaihtoehdot = unique(hs.k.Vastaus(:,tt.kysymys));

tt.hallitus = [1 2 5 7 9];
tt.oppositio = [3 4 6 8 9];
tt.opanalyysi = 1;

tt.hallitus = [1 1 0 0 1 0 1 1 0];
tt.hallitus = [1 1 1 1 1 1 1 1 1];
tt.oppositio = -(tt.hallitus-1);

tt.vastauspuoluejako = {'KOK', 'SDP' ,'PS' ,'KESK','VAS' ,'VIHR' ,'RKP','KD' ,'MUUT'};  % hallitukset & mahdollinen
tt.puoluecolor = [0 0 1; 1 0 0; 0 0 0.5; 0 0.6 0; 0.7 0 0; 0 1 0;0.7 0.7 0.2;0.6 0.6 1;0.7 0.7 0.7];
tt.vastaus{1} = zeros(length(tt.vastausvaihtoehdot),length(tt.vastauspuoluejako)); %{vuosi}  6 = [KOK SDP PS KESK VIHR RKP MUUT]
tt.vastaus{3} = zeros(length(tt.vastausvaihtoehdot),length(tt.vastauspuoluejako)); %{vuosi}  6 = [KOK SDP PS KESK VIHR RKP MUUT]

tt.yhteensaaania = 0;
for i=1:length(tt.edustaja2011)
    vastaus = find(strcmp(tt.vastausvaihtoehdot,hs.k.Vastaus(tt.edustaja2011(i),tt.kysymys))); % 1..vastausvaihtoehdot
    loyty = 0;
    for j=1:length(tt.vastauspuoluejako)-1
        if strcmp(hs.Puolue{tt.edustaja2011(i)},tt.vastauspuoluejako{j})
            loyty = j;
        end
    end
    if loyty == 0
        loyty = length(tt.vastauspuoluejako);
    end
    tt.vastaus{1}(vastaus,loyty) = tt.vastaus{1}(vastaus,loyty) + str2num(C{4}{i});
    tt.yhteensaaania = tt.yhteensaaania + str2num(C{4}{i});
    tt.vastaus{3}(vastaus,loyty) = tt.vastaus{3}(vastaus,loyty) + 1;
end

% cp
tt.vastaus{2} = zeros(length(tt.vastausvaihtoehdot),length(tt.vastauspuoluejako)); %{vuosi}  6 = [KOK SDP PS KESK VIHR RKP MUUT]
%tt.vastauspuoluejako = {'KOK', 'SDP' ,'PS' ,'KESK' ,'VIHR' ,'RKP' ,'MUUT'};  % hallitukset & mahdollinen

if tt.opanalyysi == 0
for i=1:length(tt.edustaja2007)
    vastaus = find(strcmp(tt.vastausvaihtoehdot,hs.k.Vastaus(tt.edustaja2007(i),tt.kysymys))); % 1..vastausvaihtoehdot
    loyty = 0;
    for j=1:length(tt.vastauspuoluejako)-1
        if strcmp(hs.Puolue{tt.edustaja2007(i)},tt.vastauspuoluejako{j})
            loyty = j;
        end
    end
    if loyty == 0
        loyty = length(tt.vastauspuoluejako);
    end
    tt.vastaus{2}(vastaus,loyty) = tt.vastaus{2}(vastaus,loyty) + 1;
end
end



tt.vastaus{1} = tt.vastaus{1}/tt.yhteensaaania*100;
tt.vastaus{3} = tt.vastaus{3}/length(tt.edustaja2011)*100;


%tt.puoluecolor = [0 0 1;1 0 0;0 0 0.5;0 0.6 0;0 1 0;0.7 0.7 0;0.7 0.7 0.7];
colormap(tt.puoluecolor);

apu = strcat('Kansanedustajan ��nim��r�ll� painotettu vastaus:');
if length(apu)>p.kysymysjako*2
    apu = {apu(1:p.kysymysjako),apu(p.kysymysjako+1:p.kysymysjako*2),apu(p.kysymysjako*2+1:end)}
elseif length(apu)>p.kysymysjako
    apu = {apu(1:p.kysymysjako),apu(p.kysymysjako+1:end)}
end

figure(1000);

subplot(2,1,1);
%set(gca,'position',[0.033 0.02 0.96 0.955])
set(gca,'position',[0.033 0.53 0.94 0.43])  % [x x lev

tt.vastaus{2} = tt.vastaus{1};  % OTETAAN TALTAAN ennen hallituksen muodostamista

%tt.vastaus{1} = [tt.vastaus{1}(:,tt.hallitus) tt.vastaus{1}(:,tt.oppositio)*0];

for i=1:9
    for j=1:length(tt.vastausvaihtoehdot)
        if tt.hallitus(i) == 0
            tt.vastaus{1}(j,i) = 0;
        end
    end
end

barh(tt.vastaus{1},'stacked');
title(apu);
set(gca,'YTickLabel',{})
legend(tt.vastauspuoluejako);
%xlabel('%');



for i=1:length(tt.vastausvaihtoehdot)
    apu2 = (sum(tt.vastaus{1}'));
    apu = tt.vastausvaihtoehdot{i};
    
    if length(apu)>p.jakokoko*2
        apu = {apu(1:p.jakokoko),apu(p.jakokoko+1:p.jakokoko*2),apu(p.jakokoko*2+1:end)}
    elseif length(apu)>p.jakokoko
        apu = {apu(1:p.jakokoko),apu(p.jakokoko+1:end)}
    end
    
    text(apu2(i)+1,i,apu);
end



apu = strcat('Hallitus: ',hs.k.Kysymys{tt.kysymys});
if length(apu)>p.kysymysjako*2
    apu = {apu(1:p.kysymysjako),apu(p.kysymysjako+1:p.kysymysjako*2),apu(p.kysymysjako*2+1:end)}
elseif length(apu)>p.kysymysjako
    apu = {apu(1:p.kysymysjako),apu(p.kysymysjako+1:end)}
end
axis([0 100 0.5 length(tt.vastausvaihtoehdot)+0.5])
subplot(2,1,2);

%tt.vastaus{2} = [tt.vastaus{2}(:,tt.hallitus)*0 tt.vastaus{2}(:,tt.oppositio)];
%tt.vastaus{2}.*(tt.oppositio'*ones(1,9))
for i=1:9
    for j=1:length(tt.vastausvaihtoehdot)
        if tt.oppositio(i) == 0
            tt.vastaus{2}(j,i) = 0;
        end
    end
end


tt.vastaus{2} = tt.vastaus{3};

set(gca,'position',[0.033 0.04 0.94 0.43]) 
barh(tt.vastaus{2},'stacked');
title(apu);
%set(gca,'YTickLabel',tt.vastausvaihtoehdot)
set(gca,'YTickLabel',{})
%legend(tt.vastauspuoluejako);
axis([0 100 0.5 length(tt.vastausvaihtoehdot)+0.5])
xlabel('%');
for i=1:length(tt.vastausvaihtoehdot)
    apu2 = (sum(tt.vastaus{2}'));
    
    apu = tt.vastausvaihtoehdot{i};
    
    if length(apu)>p.jakokoko*2
        apu = {apu(1:p.jakokoko),apu(p.jakokoko+1:p.jakokoko*2),apu(p.jakokoko*2+1:end)}
    elseif length(apu)>p.jakokoko
        apu = {apu(1:p.jakokoko),apu(p.jakokoko+1:end)}
    end
    
    text(apu2(i)+1,i,apu);
end

%%




