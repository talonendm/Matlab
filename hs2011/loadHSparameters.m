function p = loadHSparameters(p)
% *************************************************************************
% HS parameters set
% 110516
% talonen
% *************************************************************************
% parameter setting: not depending the data
% *************************************************************************

% *************************************************************************
% kaikki kyssärit  - 110803
% *************************************************************************
p.mitkakyssaritmukana = 1:31;
%p.mitkakyssaritmukana = [1 4:31];
%p.mitkakyssaritmukana = 6:15;
%p.mitkakyssaritmukana = 9:11;
%p.mitkakyssaritmukana = 1:20;
%p.mitkakyssaritmukana = 21;
% *************************************************************************
if length(p.mitkakyssaritmukana) == 31
    % IDA ym.
    p.moneenkovastannuettamukaan = 20;  %35 yleisin.. 
    %p.moneenkovastannuettamukaan = 0;  %35 yleisin.. % JS:ään vähän
    %vastanneista kukaan ei päässyt hallitykseen (1-20 kertaa vastanneet)
else
    % vähintään puoliin kyssäreistä pitää olla vastaus.. 
    p.moneenkovastannuettamukaan = floor(length(p.mitkakyssaritmukana)/2);
end
% *************************************************************************

p.logotkarttaan = 1;
% *************************************************************************
if p.idapaperfigure == 0
    p.kysymyksenpainotusA = 0.5;%0.25;
    p.kysymyksenpainotusB = 0.5;%2;
    % *************************************************************************
    p.vastausmaarallapainottaminen = 0;   % kaikki 0..1 dataa, hs paino vain
    p.vastausmaaraneliojuuressa = 0;
    % *************************************************************************
    p.puoluekohtainenpainotus = 0;%17;%25; %17;%6; % =0 ei tehdä, muuten numero puolueelle 6=KOK
    % *************************************************************************
    p.pelkkahallitus = 0;  % 1;
    p.puoluereunat = 0; p.puoluereunanpaksuus = 0.13; % hattuvakio
    p.toinenkorkeuskayra = 3; % =0 jos ei ole, muuten kerroin ekaan kayraan
    p.vakiopuoluereuna = 0.002;
    p.kaikkiehdokkaat = 1; % = 0 jos pelkästään 200 mukaan
    % *************************************************************************
    p.hallitusneuvotteluissa = {'KOK','KD','PS','SDP','RKP','VAS'};
    p.puheenjohtajatkuvaan = 1;
    p.kansanedustajatkuvaan = 0; % menee sotkuseks
    if p.kansanedustajatkuvaan == 1
        p.kansanedustajapiirtoaaniraja = 10000;
    else
        p.kansanedustajapiirtoaaniraja = 0;
    end
    % *************************************************************************
    p.bw = 1;
    % *************************************************************************
    % p.puhislista = {'Katainen Jyrki','Wallin Stefan','Arhinmäki Paavo','Räsänen Päivi','Soini Timo','Urpilainen Jutta','Sinnemäki Anni','Kiviniemi Mari'};
    % Suomen Senioripuolue SSP Heikki Silván pitäiskö toi a koodata?
    p.puhislista = {'Katainen Jyrki','Wallin Stefan','Arhinmäki Paavo','Räsänen Päivi','Soini Timo','Urpilainen Jutta','Sinnemäki Anni','Kiviniemi Mari','Pesonen Antti','Harju Hannu','Savola Terttu','Keronen Jiri','Palmulehto Pasi','Hakanen Yrjö','Silvan Heikki','Tanski Juhani','Helo Kalevi'};
    % **********
    p.kansalaistarkkailu = 1;  % 1= kerrotaan äänet saaduilla äänillä. nyt 1 ehdokkaalle
    % **********
    p.kompo2maaritys = 0;  % 0 = kaikki kyssarit
    p.vertailekompoissakahta = 0;
    p.regrPCA = 1; % käytetäänkö hyötyniemen toolia
    % *************************************************************************
    p.varinvarjostus = 1;  %1  % BW:ssa aika huono
    p.kuvaantyhjiaruutuja = 0; %1..
    % *************************************************************************
    p.rotaatio = 1;
    p.rotaatiopuolue1 = {'VAS'};
    p.rotaatiopuolue2 = {'KOK'};
    p.rotaatio = 0;
    p.rotaatiopuolue1 = {'-'};
    p.rotaatiopuolue2 = {'-'};
    % *************************************************************************
    p.kokosuomi = 1;
    if p.kokosuomi == 0
        %p.vaalipiirivalinta = 'Helsingin vaalipiiri';
        p.vaalipiirivalinta ='Uudenmaan vaalipiiri';
        %p.vaalipiirivalinta ='Vaasan vaalipiiri'
    end
    % *************************************************************************
    if p.kompo2maaritys ==1
        % HALLITUSKYSYMYS
        p.X1 = [1 6 7 8 9 12 13 27 28 30 31];
        p.X2 = [10 11 16 17 18 20 21];
        p.katainen = [p.X1 p.X2];     % kysymys
        p.katainenpaino = [1.6 1.0 0.5 0.5 1.5 0.6 0.4 0.05 1 0.3 0.2 1.6 1 0.5 0.7 0.5 0.3 0.2]; % massa   
        %p.katainenpaino = [1.2 1.0 1 1 1.5 0.8 0.6 0.05 1 0.6 0.5 1.3 1 0.8 0.7 0.5 0.4 0.05]; % massa
        p.labelX = 'Talouspolitiikan linja';
        p.labelY = 'Suomen EU-, ulko- ja turvallisuuspolitiikan linja';
        % 28.4.2011
        %katainen
        % 2. Talouspolitiikan linja
        %a) Valtiovarainministeriö on arvioinut julkisen talouden kestävyyden turvaamiseksi tarvittavan ylijäämän olevan vuonna 2015 4 % bruttokansantuotteesta. Onko eduskuntaryhmänne valmis kirjaamaan hallitusohjelmaan päätöksen alkavan hallituskauden aikana toteutettavasta toimintasuunnitelmasta, jolla turvataan julkisen talouden kestävyys ja hyvinvointiyhteiskunnan rahoituspohja tulevaisuudessa?
        %b) Mikäli olette tähän valmiita, missä määrin kullakin seuraavista vaihtoehdoista ja millä uusilla konkreettisilla keinoilla vaikutetaan kestävän rahoitusaseman saavuttamiseen (% / bkt tai mrd.)?
        %– työurien pidentämiseen vaikuttavat uudet toimet
        %– kuntarakenteen ja palvelurakenteen uudistaminen
        %– talouskasvua ja työllisyyttä lisäävät uudet toimet esimerkiksi verotusta uudistamalla
        %–valtion budjettitaloutta kiristävät toimet (verotuksen kiristäminen ja menojen leikkaaminen)
        %c) Millaiseksi valtion menotaso tulisi asettaa ja mihin aikatauluun hallituksen tulee sitoutua valtiontalouden velkaantumisen pysäyttämiseksi?
        %3. Suomen EU-, ulko- ja turvallisuuspolitiikan linja
        %a) Mitkä ovat mielestänne Suomen EU-, ulko- ja turvallisuuspolitiikan painopisteet?
        %b) Miten turvaisitte myös jatkossa Suomen kansainvälisen aseman ja vaikutusvallan meihinkin vaikuttavissa kysymyksissä?
        %c) Oletteko valmiit hyväksymään Suomen jo tekemien sitoumuksien loppuunsaattamisen euroalueen sekä Suomen vakauden, talouskasvun ja työllisyyden turvaamiseksi?
    end
    % *************************************************************************
else
    % p.idapaperfigure 1-3
    p.kysymyksenpainotusA = 0.5;%0.25;
    p.kysymyksenpainotusB = 0.5;%2;
    % *************************************************************************
    p.vastausmaarallapainottaminen = 1;   % kaikki 0..1 dataa, hs paino vain
    p.vastausmaaraneliojuuressa = 1;
    % *************************************************************************
    p.puoluekohtainenpainotus = 0;%25; %17;%6; % =0 ei tehdä, muuten numero puolueelle 6=KOK
    % *************************************************************************
    if p.idapaperfigure == 3
        p.pelkkahallitus = 0;  % 1;
        p.kaikkiehdokkaat = 0; % = 0 jos pelkästään 200 mukaan
        p.kansalaistarkkailu = 0;  % 1= kerrotaan äänet saaduilla äänillä. nyt 1 ehdokkaalle
    else
        p.pelkkahallitus = 0;
        p.kaikkiehdokkaat = 1; % = 0 jos pelkästään 200 mukaan
        if p.idapaperfigure == 2
            p.kansalaistarkkailu = 1;  % 1= kerrotaan äänet saaduilla äänillä. nyt 1 ehdokkaalle
        else
            p.kansalaistarkkailu = 0;  % 1= kerrotaan äänet saaduilla äänillä. nyt 1 ehdokkaalle
        end
    end
    p.puoluereunat = 0; p.puoluereunanpaksuus = 0.13; % hattuvakio
    p.toinenkorkeuskayra = 3; % =0 jos ei ole, muuten kerroin ekaan kayraan
    p.vakiopuoluereuna = 0.002;
    
    
    % *************************************************************************
    p.hallitusneuvotteluissa = {'KOK','KD','PS','SDP','RKP','VAS'};
    p.puheenjohtajatkuvaan = 1;
    p.kansanedustajatkuvaan = 0; % menee sotkuseks
    if p.kansanedustajatkuvaan == 1
        p.kansanedustajapiirtoaaniraja = 1000;  % 1000 tulee kaikki
    else
        p.kansanedustajapiirtoaaniraja = 0;
    end
    % *************************************************************************
    p.bw = 1;
    % *************************************************************************
    % p.puhislista = {'Katainen Jyrki','Wallin Stefan','Arhinmäki Paavo','Räsänen Päivi','Soini Timo','Urpilainen Jutta','Sinnemäki Anni','Kiviniemi Mari'};
    % Suomen Senioripuolue SSP Heikki Silván pitäiskö toi a koodata?
    p.puhislista = {'Katainen Jyrki','Wallin Stefan','Arhinmäki Paavo','Räsänen Päivi','Soini Timo','Urpilainen Jutta','Sinnemäki Anni','Kiviniemi Mari','Pesonen Antti','Harju Hannu','Savola Terttu','Keronen Jiri','Palmulehto Pasi','Hakanen Yrjö','Silvan Heikki','Tanski Juhani','Helo Kalevi'};
    % **********
    
    % **********
    p.kompo2maaritys = 0;  % 0 = kaikki kyssarit
    p.vertailekompoissakahta = 0;
    p.regrPCA = 1; % käytetäänkö hyötyniemen toolia
    % *************************************************************************
    p.varinvarjostus = 1;  %1  % BW:ssa aika huono
    p.kuvaantyhjiaruutuja = 0; %1..
    % *************************************************************************
    p.kokosuomi = 1;
    if p.kokosuomi == 0
        p.vaalipiirivalinta = 'Helsingin vaalipiiri';
        %p.vaalipiirivalinta ='Uudenmaan vaalipiiri';
        %p.vaalipiirivalinta ='Vaasan vaalipiiri'
        
    end
    % *************************************************************************
    if p.kompo2maaritys ==1
        % HALLITUSKYSYMYS
        p.X1 = [1 6 7 8 9 12 13 27 28 30 31];
        p.X2 = [10 11 16 17 18 20 21];
        p.katainen = [p.X1 p.X2];     % kysymys
        p.katainenpaino = [1.6 1.0 0.5 0.5 1.5 0.6 0.4 0.05 1 0.3 0.2 1.6 1 0.5 0.7 0.5 0.3 0.2]; % massa   
        %p.katainenpaino = [1.2 1.0 1 1 1.5 0.8 0.6 0.05 1 0.6 0.5 1.3 1 0.8 0.7 0.5 0.4 0.05]; % massa
        p.labelX = 'Talouspolitiikan linja';
        p.labelY = 'Suomen EU-, ulko- ja turvallisuuspolitiikan linja';
        
        
        p.X1 = [1:31];
        p.X2 = [9:11];
        p.katainen = [p.X1 p.X2];     % kysymys
        p.katainenpaino =  ones(length([p.X1 p.X2]),1)'; % massa   
        %p.katainenpaino = [1.2 1.0 1 1 1.5 0.8 0.6 0.05 1 0.6 0.5 1.3 1 0.8 0.7 0.5 0.4 0.05]; % massa
        p.labelX = 'XXXTalouspolitiikan linja';
        p.labelY = 'XXXSuomen EU-, ulko- ja turvallisuuspolitiikan linja';
        
        
        
    end
    % *************************************************************************
      
    % *************************************************************************
    % JOS CR
    p.rotaatio = 0;
    p.rotaatiopuolue1 = {'-'};
    p.rotaatiopuolue2 = {'-'};
    % *************************************************************************
    %p.rotaatio = 1;
    %p.rotaatiopuolue1 = {'VAS'};
    %p.rotaatiopuolue2 = {'KOK'};
    % *************************************************************************
    
end


% TEE solut tästä
% (yleiset) Kysymykset 1-5
% Eläkkeet 6-8
% Talous 9-11
% Verot 12-15
% Puolustus 16-17
% Ulkomaat 18-21
% Kotimaa 22-26
% Kunnat 27-30
% Hallituspohja 31



% *************************************************************************
end
