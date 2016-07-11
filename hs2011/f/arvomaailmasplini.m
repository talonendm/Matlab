function arvomaailmasplini(hs,m,p,num)

n = length(hs.Ika);

figure(num);
%%
% muokataan vastaukset:
unique(hs.k.Vastaus(:,1))  % eri vastausvaihtoehdot
%%
% *************************************************************************
% MIETI KERTOIMET TARKKAAN!
hs.k.v = zeros(n,31);
% 1 Tuloerot ovat kasvaneet Suomessa 1990-luvun puolivälin jälkeen nopeasti. Miten siihen pitäisi suhtautua?
kysymys = 1;
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'Tuloeroja on kavennettava huomattavasti.')*-100;
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'Tuloeroja on kavennettava lievästi.')*-50;
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'Tuloerot saavat edelleen kasvaa hillitysti.')*50;
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'Tuloerot saavat kasvaa vapaasti.')*100;
%hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'-')*0;
%hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'Tuloerot ovat nyt sopivalla tasolla.')*0;

% kyllä / ei
kysymykset = [2 3 15];
for i=1:length(kysymykset)
    kysymys = kysymykset(i);
    hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'Ei.')*1;%*-100;
    hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'Kyllä.')*0;%*100;
    %hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'-')*0;
end

% 4 Lapsilisää maksetaan jokaisesta Suomessa asuvasta lapsesta 17 ikävuoteen asti riippumatta vanhempien tulotasosta. Mitä lapsilisille pitäisi tehdä?
kysymys = 4;
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'Lapsilisän ei pidä vaikuttaa toimeentulotulen määrään, vaan toimeentulotuen varassa elävän tulee saada sekä tuki että lapsilisä.')*-100;
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'Lapsilisän suuruus pitää porrastaa vanhempien tulojen mukaan.')*-50;
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'Lapsilisät tulee poistaa hyvätuloisilta.')*50;
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'Säilytetään nykytilanne.')*100;
%hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'-')*0;

% 5 Sosiaali- ja terveysministeriössä valmistellaan parhaillaan niin sanottua ikälakia, johon mahdollisesti kirjataan yhtenäisiä ikäihmisten hoidon laatuvaatimuksia ja palvelutakuu. Pitäisikö ikäihmisten hoivatakuu kirjata lakiin ja velvoittaa kunnat tarjoamaan vanhuksille subjektiivinen oikeus hoivaan?
kysymys = 5;
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'Ei. Kuntien velvoitteita ei pidä lisätä.')*-100;
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'Ei. Riittää, että ikäihmisten hoidon taso määritellään suosituksilla.')*-50;
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'Kyllä. Ikäihmisillä pitää olla subjektiivinen oikeus hyvään hoitoon.')*100;
%hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'-')*0;

% 6 Vanhuuseläkkeelle siirtymisen ikä on nyt 63–68 vuotta. Eläkeiän alarajaa pitäisi mielestäni 
kysymys = 6;
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'laskea.')*-100;
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'nostaa usealla vuodella.')*100;
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'nostaa vuodella.')*50;
%hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'pitää nykyisellään.')*0;
%hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'-')*0;

% 7 Suomalaiset jäävät nykyään eläkkeelle keskimäärin alle 60-vuotiaina. Keskeisin syy ennenaikaiseen eläkkeeseen on työkyvyttömyys. Työssäjaksamisen ja työhyvinvoinnin kehittämisen lisäksi eläkeiän nostamiseksi on esitetty myös muita keinoja. Mitä seuraavista käyttäisit ensisijaisesti työurien pidentämiseksi? 
kysymys = 7;
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'Nostaisin vanhuuseläkeiän alarajaa 63 vuodesta.')*100;
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'Ottaisin käyttöön nämä kaikki.')*50;
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys), 'Poistaisin mahdollisuuden osa-aikaeläkkeeseen.')*-50;
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys), 'Poistaisin työttömyysputken eli oikeuden työttömyyspäivärahan lisäpäiviin ennen eläkeikää.')*-100;
%hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'En käyttäisi mitään näistä.')*0;
%hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'-')*0;

% 8 Eläkkeitä korotetaan niin sanotulla taitetulla indeksillä, jossa kuluttajahintojen painoarvo on 80 prosenttia ja palkansaajien yleisen palkkakehityksen 20 prosenttia. 1990-luvun puolivälissä käyttöön otetun taitetun indeksin vuoksi eläkeläiset ovat jääneet jälkeen yleisestä tulokehityksestä. Mitä taitetulle indeksille pitäisi tehdä?
kysymys = 8;
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'Indeksi tulisi muuttaa entiselleen, jolloin puolet siitä määräytyisi palkkojen ja puolet hintojen mukaan.')*100;
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'Indeksin tulisi seurata palkkojen nousua hieman nykyistä suuremmalla painolla.')*50;
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys), 'Kuluttajahintojen pitäisi vaikuttaa indeksiin enemmän kuin palkkojen.')*-50;
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys), 'Nykyistä indeksiä ei ole tarpeen muuttaa.')*-100;
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'Palkkakehityksen tulisi vaikuttaa indeksiin enemmän kuin kuluttajahintojen.')*0;
%hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'-')*0;

% Laaja kaikkia miehiä koskeva asevelvollisuus on Euroopan maista vain Kreikassa, Kyproksessa ja Suomessa. Miten järjestäisit Suomen asevelvollisuuden?
% hs.k.Kysymys{16}
kysymys = 16;
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys), 'Myös naisten tulisi olla asevelvollisia.');
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys), 'Nykykäytäntö on hyvä.');
    

% 20 Suomi on sitoutunut YK:n asettamaan tavoitteeseen, jonka mukaan kehitysyhteistyön määrärahat tulee vuoteen 2015 mennessä nostaa 0,7 prosenttiin bruttokansantuotteesta. Viime vuonna Suomi antoi kehitysyhteistyöhön 965,6 miljoonaa euroa, mikä on 0,55 prosenttia bruttokansantuotteesta. Mitä kehitysavulle pitäisi tehdä? 
% *
kysymys = 20; %OK
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'Kehitysyhteistyöhön ei pidä antaa valtion varoja.');
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys), 'Kehitysyhteistyövaroja ei pidä lisätä.');
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys), 'Kehitysyhteistyövaroja pitää leikata.');

% 22 Aselakia kiristettiin syksyllä 2010 mm. nostamalla käsiaseluvan ikäraja 20 vuoteen. Mitä aselaille pitäisi uudessa eduskunnassa tehdä? 
kysymys = 22; %OK
% *    
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'Aselaki on nyt hyvä.');
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'Aselakia pitäisi lieventää.');
% *    
kysymys = 24; %OK
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys),'Kyllä, se kuuluu suomalaiseen kevätjuhlaperinteeseen.');

kysymys = 25; %OK
hs.k.v(:,kysymys) = hs.k.v(:,kysymys) + strcmp(hs.k.Vastaus(:,kysymys), 'Liian löysää.');



%%

X = [hs.k.v(:,2) hs.k.v(:,16) hs.k.v(:,25) hs.k.v(:,20) hs.k.v(:,22) hs.k.v(:,24)];

[COEFF,SCORE,latent,tsquare] = princomp(X);



p.arvomuuttujia = size(X,2);
p.arvomuuttujia = 3;

hsARVO = sum(X')'/p.arvomuuttujia*10;

ARVO = SCORE(:,1) - min(SCORE(:,1));
ARVO = ARVO./max(ARVO)*10;

figure(4);
hist(ARVO);



for i=1:m.n.puolueet
    r.hsarvo(i) = mean(hsARVO(m.puolueN{i}))
end




figure(7);



%varitys = find(m.puolueEhdokkaita>144==1); %HESARIN mukaan
varitys = find(m.puolueEhdokkaita>100==1); %HESARIN mukaan
%varitys = find(m.puolueEhdokkaita<144==1);
varitys = [varitys 18]; % RKP;

%varitys = find((m.puolueEhdokkaita<144).*( m.puolueEhdokkaita>10==1));
vertailussa = m.puolueet(varitys);


p.ikaotsikko = 1;
p.ikal=0;

p.ika1=18;
p.ika2=30;
p.ikaaskel = 10;

p.ika1=p.ika1+p.ikaaskel*p.ikal;
p.ika2=p.ika2+p.ikaaskel*p.ikal;
clear apu;
for i=1:length(vertailussa)
    % hsARVO / ARVO
    %data{i}= ARVO(m.puolueN{strcmp(m.puolueet,vertailussa{i})==1});
    
    
    apu = m.puolueN{strcmp(m.puolueet,vertailussa{i})==1};
    apu = apu(hs.Ika(apu)>=p.ika1);
    apu = apu(hs.Ika(apu)<p.ika2);
    data{i} = hsARVO(apu);
end

%vari{1} = [0 1 0];
%vari{2} = [0.7 0 0];
%vari{3} = [0.7 0.7 0];
%vari{4} = [1 0 0];
%vari{5} = [0 0.6 0];
%vari{6} = [0 0 1];
%vari{7} = [0 0 0.7];
%vari{8} = [0 0 0];
%data2 = ARVO(200:800);
p.normaalisovitus = 0;
%m.puolueN(strcmp(m.puolueet,'KESK')==1)
figure(10);
hold on;
for i=1:length(vertailussa)
keskiarvo(i) = mean(data{i});
varianssi = std(data{i});  % huom - keskihajonta kaavassa

%p1 = pdf('Normal',0:0.1:10,keskiarvo,varianssi)
p1 = normpdf(0:0.1:10,keskiarvo(i),varianssi)*100;

p1 = p1/sum(p1*0.1)*100;

maxp1(i) = max(p1);


%pisteet = hist(data{i},(p.arvomuuttujia+1))./(m.puolueEhdokkaita(varitys(i)))*100;  % VIRHEELLINEN
pisteet = hist(data{i},0:10/p.arvomuuttujia:10)./(m.puolueEhdokkaita(varitys(i)))*100;

x = 0:(10/p.arvomuuttujia):10;
y = pisteet;
cs = spline(x,y);
xx = 0:0.1:10;

if p.normaalisovitus == 1
    plot(0:0.1:10,p1,'color',p.vari{varitys(i)},'linewidth',3);
else
    
    
    
    
    % tulee negatiivisia tästä: ppval(cs,xx)
    yy = ppval(cs,xx);
    yy = (yy>0).*yy;
    kerroin = sum(yy/10);
    yy = yy/kerroin*100;%/10;
    y = y/kerroin*100;%/10;
    disp(sum(yy/10))
    plot(x,y,'o',xx,yy,'-','color',p.vari{varitys(i)},'linewidth',2);
end


maxspline(i) = max(y);
apu = (find(max(y)==y==1)-1)*10/p.arvomuuttujia;
maxsplinex(i) = apu(1);

minspline(i) = min(y(find((y>1)==1)));
apu = (find(minspline(i)==y==1)-1)*10/p.arvomuuttujia;
minsplinex(i) = apu(1);

%plot(0:2:10,pisteet,'color',p.vari{varitys(i)});
%m.puolueEhdokkaita
end
p.tekstiY = 0;
for i=1:length(vertailussa)
    %setBold(text, 'true')
    if p.normaalisovitus == 1
        text(keskiarvo(i),maxp1(i)+p.tekstiY,['\bf' vertailussa{i}],'horizontalalignment','center','verticalalignment','bottom');
        text(keskiarvo(i),maxp1(i)+p.tekstiY, vertailussa{i},'color',p.vari{varitys(i)},'horizontalalignment','center','verticalalignment','bottom');
    else
        if maxsplinex(i) == 0
            sijainti = 'left';
        elseif maxsplinex(i) == 10
            sijainti = 'right';
        else
            sijainti = 'center';
        end
        %text(maxsplinex(i),maxspline(i)+p.tekstiY,['\bf' vertailussa{i}],'horizontalalignment',sijainti,'verticalalignment','bottom');
        text(maxsplinex(i),maxspline(i)+p.tekstiY, vertailussa{i},'color',p.vari{varitys(i)},'horizontalalignment',sijainti,'verticalalignment','bottom');
        
        if minsplinex(i) == 0
            sijainti = 'left';
        elseif minsplinex(i) == 10
            sijainti = 'right';
        else
            sijainti = 'center';
        end
        
        text(minsplinex(i),minspline(i)-p.tekstiY, vertailussa{i},'color',p.vari{varitys(i)},'horizontalalignment',sijainti,'verticalalignment','top');
   
    end
end
%axis([0 10 0 65])
xlabel('<--- liberaali    konservatiivinen --->');
ylabel('%');
%legend('VIHR','VIHR','VIHR','VIHR','VIHR','VIHR','VIHR','VIHR');

if p.ikaotsikko == 1
    title(strcat('Puolueiden arvokonservatiivisuus vaalikonevastausten perusteella (',num2str(p.ika1),'-',num2str(p.ika2),' vuotiaat)'));
else
    title('Puolueiden arvokonservatiivisuus vaalikonevastausten perusteella');
end
    
    
hold off;
%%
figure(6);
hsjarjestys = [4 1 3 2 5 8 6 7];
for i=1:8
    subplot(8,1,i);
    j = hsjarjestys(i);
    boxplot(data{j},'orientation','horizontal');
    title(vertailussa{j});
    axis([0 10 0.8 1.2])
end
%%
figure(9);
hsjarjestys = [4 1 3 2 5 8 6 7];
for i=1:8
    subplot(4,2,i);
    j = hsjarjestys(i);
    hist(data{j},0:(10/p.arvomuuttujia):10,'color',p.vari{varitys(j)});
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',p.vari{varitys(j)},'EdgeColor',[0 0 0])
    title(vertailussa{j});
    axis([-(10/p.arvomuuttujia)/2 10+(10/p.arvomuuttujia)/2 0 max(hist(data{j}))+3])
    ylabel('Ehdokkaita');
    xlabel('arvopisteet');
end

%%

figure(11);
bar(COEFF(:,1:2)')
legend(hs.k.Kysymys{2} ,hs.k.Kysymys{16},  hs.k.Kysymys{25},  hs.k.Kysymys{20},  hs.k.Kysymys{22},  hs.k.Kysymys{24});

title('Pääkomponenttianalyysi');



end