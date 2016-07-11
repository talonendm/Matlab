function Votersfigure(p,b,A,hs,m,SCORE,Zv)
% *************************************************************************
% Finnish Voters
% 110516
% *************************************************************************
figure(101);
%colormap(p.omacolormap);

%Z2 = b.varihila(:,:);  % PUOLUEMAX
%Z2 = b.hila(:,:)/(2931817/4159857);%0.71;      % KOKOKANSA


%Z2 = b.hila(:,:);

%Z2 = b.rahahila(:,:);
Z3 = (b.kysymyshila(:,:,2)./b.vastausmaarahila(:,:)); % esim. ydinvoimassa sum 1787, 1820 on max

p.suhteellinenraha = 0;
if p.suhteellinenraha == 1
    Z2 = b.rahahila(:,:)./b.hila(:,:);
    
    
    % rajotus: 
    % KOODAA tähän et jos Z2 yli 10e eli suhteellinen vaalirahoitus, niin on
    % 10.1.. sitten saa skaalan paremmaksi!!! TEE! 30.5.2011
    p.suhteellisenrajanmaara = 15; % euroa
    
    Z2 = Z2 - (Z2>p.suhteellisenrajanmaara).*Z2 + (Z2>p.suhteellisenrajanmaara).*(p.suhteellisenrajanmaara);
    
    axis([-p.hilaskaala p.hilaskaala -p.hilaskaala p.hilaskaala 0 p.suhteellisenrajanmaara]);
    
else
    Z2 = b.rahahila(:,:);
end
    

% Pakotus... poista jos halut rahan
Z2 = b.hila;



% 


if p.logotkarttaan == 1
    p.lsk = 1;
    
    % OIKEE!
    %surface(-p.hilaskaala*p.lsk:p.lsk:p.hilaskaala*p.lsk,-p.hilaskaala*p.lsk:p.lsk:p.hilaskaala*p.lsk,Z2,(Zv),'FaceColor','texturemap','EdgeColor','none', 'CDataMapping','direct');
    
    
    %surface(-p.hilaskaala*p.lsk:p.lsk:p.hilaskaala*p.lsk,-p.hilaskaala*p.lsk:p.lsk:p.hilaskaala*p.lsk,Z2,(Zv),'edgecolor',[0.7 0.7 0.7]);

    
    %Zopacity =  min(1,b.hila/max(max(b.hila))*100);
    Zopacity =  min(1,b.hila/mean(mean(b.hila))*10);

    % kokeillaan varjostusta
    surface(-p.hilaskaala*p.lsk:p.lsk:p.hilaskaala*p.lsk,-p.hilaskaala*p.lsk:p.lsk:p.hilaskaala*p.lsk,Z2,Zv,'FaceAlpha','flat',...
                'AlphaDataMapping','scaled',...
                'AlphaData',(Zopacity));



    shading interp
    
    %axis off;        
            
            
else
    p.lsk = 1;   % lis�� n�m� teksteihin
    %surface(-p.hilaskaala:p.hilaskaala,-p.hilaskaala:p.hilaskaala,Z2,(Zv));
    
    % oikee
    %surface(-p.hilaskaala*p.lsk:p.lsk:p.hilaskaala*p.lsk,-p.hilaskaala*p.lsk:p.lsk:p.hilaskaala*p.lsk,Z2,(Zv),'edgecolor',[0.5 0.5 0.5]);
    
    % harmaa väritys
    colormap(gray)
    %colormap(cool)
    
    
    % OIKEE
    surface(-p.hilaskaala*p.lsk:p.lsk:p.hilaskaala*p.lsk,-p.hilaskaala*p.lsk:p.lsk:p.hilaskaala*p.lsk,Z2);
    
    
    
    
    
    %shading faceted
%shading flat
    shading interp
    
    axis off;
    
end
    %surface(-p.hilaskaala:p.hilaskaala,-p.hilaskaala:p.hilaskaala,Z2,Zv,'FaceColor','texturemap','EdgeColor','none', 'CDataMapping','direct');
%colormap(bone)
%get(gca)
%set(gca,'ZColor',[0.5 0.5 0.5]);

if p.kompo2maaritys ==1
   ylabel(p.labelX);   % nää sekasin
   xlabel(p.labelY);
end

axis tight;


if p.kansanedustajatkuvaan == 1

    for i=1:length(b.kansanedustaja.nimi)
        
        
        tekstikorkeus =  Z2(b.kansanedustaja.X(i),b.kansanedustaja.Y(i));
        if p.kansalaistarkkailu == 1
            tekstikorkeusviiva =  Z2(b.kansanedustaja.X(i),b.kansanedustaja.Y(i)) + max(max(Z2))/10 + sqrt(b.kansanedustaja.aania(i))/10;%max(max(Z2))/5;
        else
            tekstikorkeusviiva =  Z2(b.kansanedustaja.X(i),b.kansanedustaja.Y(i)) + max(max(Z2))/6;%5; % 2D pieniluku
        end
        line([(b.kansanedustaja.Y(i)-1)-p.hilaskaala (b.kansanedustaja.Y(i)-1)-p.hilaskaala],[(b.kansanedustaja.X(i)-1)-p.hilaskaala (b.kansanedustaja.X(i)-1)-p.hilaskaala],[tekstikorkeus tekstikorkeusviiva],'color',p.vari{b.kansanedustaja.puolueN(i)},'linewidth',1);
        tekstikoko = 10;
        
        
        %kirjoitetaantekstiin = strcat(b.kansanedustaja.nimi{i},' (',b.kansanedustaja.puolue{i},'):',num2str(b.kansanedustaja.aania(i)));
        kirjoitetaantekstiin = strcat(b.kansanedustaja.nimi{i},': ',num2str(b.kansanedustaja.aania(i)));
        
        % HS:n ammatit
        % ehd2011.xls mukaan uusi tietokanta, tee ID!
        %kirjoitetaantekstiin = strcat(b.kansanedustaja.ammatti{i},' (',b.kansanedustaja.puolue{i},'):',num2str(b.kansanedustaja.aania(i)));
      
        
        if b.kansanedustaja.puhis(i) == 1
            tekstikoko = 14;
            varjataan = p.vari{b.kansanedustaja.puolueN(i)};%varjataan = [0 0 0];
        else
            
            varjataan = p.vari{b.kansanedustaja.puolueN(i)};%varjataan = [0.2 0.2 0.2]; 
        end
        
        
        % tähän lisä ehdoilla
        %if (strcmp(b.kansanedustaja.nimi{i},'Soininvaara Osmo')) || b.kansanedustaja.puhis(i) == 1
        if b.kansanedustaja.aania(i)>p.kansanedustajapiirtoaaniraja || b.kansanedustaja.puhis(i) == 1   
        
        
        if b.kansanedustaja.approksimoitu(i) == 0
            
            text((b.kansanedustaja.Y(i)-1)-p.hilaskaala,b.kansanedustaja.X(i)-1-p.hilaskaala,tekstikorkeusviiva,kirjoitetaantekstiin,'color',varjataan,'fontsize',tekstikoko,'horizontalalignment','center','verticalalignment','middle');
        else
            
            if b.kansanedustaja.puhis(i) == 1
                tekstikoko = 12;
                varjataan = p.vari{b.kansanedustaja.puolueN(i)};%[0.7 0.7 0.7];
            else
                tekstikoko = 8;
                varjataan = p.vari{b.kansanedustaja.puolueN(i)};%[0.8 0.8 0.8]; 
            end
            
            
            kirjoitetaantekstiin = strcat('*',kirjoitetaantekstiin);
            
            ehdokkaansatunnaisuus = 16;
            text((b.kansanedustaja.Y(i)-1-p.hilaskaala + rand*ehdokkaansatunnaisuus-ehdokkaansatunnaisuus/2),b.kansanedustaja.X(i)-1-p.hilaskaala + (rand*ehdokkaansatunnaisuus-ehdokkaansatunnaisuus/2),tekstikorkeusviiva,kirjoitetaantekstiin,'color',varjataan,'fontsize',tekstikoko,'horizontalalignment','center','verticalalignment','middle');
        end
        
        end
        
    end
end


if p.puheenjohtajatkuvaan == 1 && p.kansanedustajatkuvaan == 0

    for i=1:length(b.puheenjohtaja.nimi)
        tekstikorkeus =  Z2(b.puheenjohtaja.X(i),b.puheenjohtaja.Y(i));
        %tekstikorkeusviiva =  Z2(b.puheenjohtaja.X(i),b.puheenjohtaja.Y(i)) + max(max(Z2))/5;
        tekstikorkeusviiva =  Z2(b.puheenjohtaja.X(i),b.puheenjohtaja.Y(i)) + max(max(Z2))/6;
        line([(b.puheenjohtaja.Y(i)-1)-p.hilaskaala (b.puheenjohtaja.Y(i)-1)-p.hilaskaala],[(b.puheenjohtaja.X(i)-1)-p.hilaskaala (b.puheenjohtaja.X(i)-1)-p.hilaskaala],[tekstikorkeus tekstikorkeusviiva],'color',p.vari{b.puheenjohtaja.puolue(i)},'linewidth',1);
        text((b.puheenjohtaja.Y(i)-1)-p.hilaskaala,b.puheenjohtaja.X(i)-1-p.hilaskaala,tekstikorkeusviiva,b.puheenjohtaja.nimi{i},'color',p.vari{b.puheenjohtaja.puolue(i)},'fontsize',12,'horizontalalignment','center','verticalalignment','bottom');
    end
end


if p.logotkarttaan == 0 || p.logotkarttaan == 1
for i=1:length(m.puolueet)
    keta = find(b.puolue==i);
    if length(keta)>2
        kxx = max(max(b.puoluehila(:,:,i)));
        kx = find(max(b.puoluehila(:,:,i)')==kxx);
        
        kyy = max(max(b.puoluehila(:,:,i)));
        ky = find(max(b.puoluehila(:,:,i))==kyy);
        
        
        varjataan = p.vari{i};%[0 0 0];
        tekstikorkeus =  Z2(kx(1),ky(1));%*1.001;
        
        tekstikoko= 14;
        if b.puoluehila(kx(1),ky(1),i)<max(b.varihila(kx(1),ky(1),:))
            %kx(1)
            %max(b.puoluehila(kx(1),ky(1),:))
            varjataan = p.vari{i};%[0.3 0.3 0.3];
            %tekstikorkeus = max(b.varihila(kx(1),ky(1),:))*1.03;
            tekstikoko= 12;
        end
            
        
        %tekstikorkeusviiva =  Z2(kx(1),ky(1))*1.2;
        tekstikorkeusviiva =  Z2(kx(1),ky(1)) + max(max(Z2))/11;
        
        %kx = round(mean(SCORE(keta,1)));
        %ky = round(mean(SCORE(keta,2)));
        % pinnassa koordinaatisto toistenp�in!
        %text((ky(1)-1)-p.hilaskaala,(kx(1)-1)-p.hilaskaala,Z2(kx(1),ky(1))+0.01,m.puolueet{i},'color',[0 0 0],'fontsize',12,'horizontalalignment','center','verticalalignment','middle');
       %text((ky(1)-1)-p.hilaskaala,(kx(1)-1)-p.hilaskaala,Z2(kx(1),ky(1))+0.01,m.puolueet{i},'color',[0 0 0],'BackgroundColor',p.vari{i},'fontsize',14,'horizontalalignment','center','verticalalignment','middle');
       
       % esim. ida3-kuvassa - ei pikkupuoluieita
       % purkka if
       if (ky(1)-1)-p.hilaskaala == -100 || sum(sum(b.puoluehila(:,:,i)))==0
       else
           
            line([(ky(1)-1)-p.hilaskaala (ky(1)-1)-p.hilaskaala],[(kx(1)-1)-p.hilaskaala (kx(1)-1)-p.hilaskaala],[tekstikorkeus tekstikorkeusviiva],'color',p.vari{i},'linewidth',2);
      
            text((ky(1)-1)-p.hilaskaala,(kx(1)-1)-p.hilaskaala,tekstikorkeusviiva,m.puolueet{i},'color',varjataan,'fontsize',tekstikoko,'horizontalalignment','center','verticalalignment','bottom');
       end
       
       
        
       
       
        %text(ky,kx,Z2(kx+101,ky+101)+0.1,m.puolueet{i},'color',[0 0 0],'BackgroundColor',p.vari{i},'fontsize',12,'horizontalalignment','center','verticalalignment','middle');
    end
    
        
end
end
%view(-70.5,50); % vois tehd� videoon py�riv�n

p.aitakuvaan = 1; 
p.aidankorkeus = max(max(Z2))/25;
p.infotekstinkorkeus = 2
if p.aitakuvaan == 1
    
    p.aitamaareet = [0.25 0.5 0.75];
    for ii=1:length(p.aitamaareet)
        % aita kuvaan
        
        if ii==1
            p.aidanvari = [1 0.8 0.8];
        elseif ii==2
            p.aidanvari = [1 0.9 0.7];
        else
            p.aidanvari = [0.8 1 0.8];
        end
        
        for i=1:p.hilakoko
            apu0 = abs(Z3(:,i)-p.aitamaareet(ii));
            apu = min(apu0);
        
        
            % HUOM! muuallakin (ehdokkaan sijainnit ovat -1 -1 väärässä
            % paikassa 30.5.2011 - KorjattU:? xx.xx.xxxx
        
            aitaX(i) = find(apu0==apu(1))-1;
            aitaY(i) = i-1;
            %aitaZ(i) = Z3(aitaX(i),aitaY(i));
            aitaZ(i) = Z2(aitaX(i)+1,aitaY(i)+1);   
        end
    
        for i=1:p.hilakoko
            line([aitaY(i)-p.hilaskaala aitaY(i)-p.hilaskaala],[aitaX(i)-p.hilaskaala aitaX(i)-p.hilaskaala],[aitaZ(i) aitaZ(i)+p.aidankorkeus],'color',[1 1 1]);
            if i>1
                line([aitaY(i-1)-p.hilaskaala aitaY(i)-p.hilaskaala],[aitaX(i-1)-p.hilaskaala aitaX(i)-p.hilaskaala],[aitaZ(i-1)+p.aidankorkeus aitaZ(i)+p.aidankorkeus],'color',p.aidanvari);
                line([aitaY(i-1)-p.hilaskaala aitaY(i)-p.hilaskaala],[aitaX(i-1)-p.hilaskaala aitaX(i)-p.hilaskaala],[aitaZ(i-1) aitaZ(i)+p.aidankorkeus],'color',p.aidanvari);
                line([aitaY(i-1)-p.hilaskaala aitaY(i)-p.hilaskaala],[aitaX(i-1)-p.hilaskaala aitaX(i)-p.hilaskaala],[aitaZ(i-1)+p.aidankorkeus aitaZ(i)],'color',p.aidanvari);
            end
            
            if i==50 && ii==1
                infoviivaX = round((aitaX(i-1) + p.hilakoko)/2);
                infoviivaY = aitaY(i-1);
                infoviivaZ = Z2(infoviivaX+1, infoviivaY+1);
                line([infoviivaY-p.hilaskaala infoviivaY-p.hilaskaala],[infoviivaX-p.hilaskaala infoviivaX-p.hilaskaala],[infoviivaZ infoviivaZ+p.infotekstinkorkeus],'color',[1 1 1]);
                text(infoviivaY-p.hilaskaala,infoviivaX-p.hilaskaala,infoviivaZ+p.infotekstinkorkeus,'<25%','color',[1 1 1],'horizontalalignment','center','verticalalignment','bottom');
            end
            if i==50 && ii==3
                infoviivaX = round(aitaX(i-1)/2);
                infoviivaY = aitaY(i-1);
                infoviivaZ = Z2(infoviivaX+1, infoviivaY+1);
                line([infoviivaY-p.hilaskaala infoviivaY-p.hilaskaala],[infoviivaX-p.hilaskaala infoviivaX-p.hilaskaala],[infoviivaZ infoviivaZ+p.infotekstinkorkeus],'color',[1 1 1]);
                text(infoviivaY-p.hilaskaala,infoviivaX-p.hilaskaala,infoviivaZ+p.infotekstinkorkeus,'>75%','color',[1 1 1],'horizontalalignment','center','verticalalignment','bottom');
            end
            
            
            
        end
    end
end

if p.rotaatio == 0
     view(-90,90)
else
    ylabel(strcat(p.rotaatiopuolue1,'-',p.rotaatiopuolue2));
end


%view(-65,57)
%shading interp
%%
%title('Arvio Suomen ��nioikeutettujen arvoista ja mielipiteist� sek� todenn�k�isin puolue tietylle arvolle.');
%zlabel('��nest�ji� kpl');
%text(-100,110,4,'(c)2011 jaakko.talonen@aalto.fi','color',[0.8 0.8 0.8]);
%text(-100,70,5,'(c)2011 talonen@cis.hut.fi','color',[0.8 0.8 0.8]);
%
%shading faceted
%shading flat
%shading interp
%%
%axis off



end
