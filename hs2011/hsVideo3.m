% video
% ERILLINEN LÄPPÄ .. aja eka HSstats..
% 110422 keskenn
%close(111);
%fig=figure;


% ZEKKAA: help axis
% AXIS VIS3D  freezes aspect ratio properties to enable rotation of
 %       3-D objects and overrides stretch-to-fill.
 
        
        


%set(fig,'DoubleBuffer','on');
%set(gca,'xlim',[-80 80],'ylim',[-80 80],...
 %     'NextPlot','replace','Visible','off')
 %mov = avifile('vaalitJT7.avi','fps',5,'compression','Cinepak') 
 %mov = avifile('vaalitJT3.avi','fps',5) 
        
            % mov = avifile('mymovie.avi','fps',videonopeus, 'compression', 'Cinepak');
%mov.Quality = 60;

disp('Asetetaan ehdokkaat hilaan');
% EDITOI! 110421. jokin virhe maksimi pisteen haussa!
p.logotkarttaan = 1;
if p.logotkarttaan == 1
    p.hilakoko=261;
    p.varianssi = 8;  %esim 3..5
else
    p.hilakoko=101;
    p.varianssi = 3;  %esim 3..5
end
p.varinvarjostus = 1;
p.hilaskaala = (p.hilakoko-1)/2;
%p.varianssi = 3;  %esim 3..5
b.hila = zeros(p.hilakoko,p.hilakoko);
b.puoluehila = zeros(p.hilakoko,p.hilakoko,length(m.puolueet));

viewlaskuri = 79;


for i=1:length(b.b2hs)
    if b.ylenimi(i)>0
        b.aania(i) = A.aania(b.ylenimi(i));
    else
        b.aania(i) = 0;
    end
end
%[so1 so2] = sort(b.aania,'descend');
[so1 so2] = sort(b.aania);

%

b.videoi=1190;%1806;%1806;%28;%1806;

ekaijollaaania=1181;%22;
framelaskuri = 0;
for ii=ekaijollaaania:b.videoi%length(b.binaaritaulu)
    i = so2(ii);
    %b.hila(round(SCORE(i,1)+101),round(SCORE(i,2)+101)) = b.hila(round(SCORE(i,1)+101),round(SCORE(i,2)+101)) + rand;

    figure(100);
    
    
    if b.ylenimi(i)>0
        
        px1 = round((SCORE(i,1)+101)/201*p.hilakoko);
        px2 = round((SCORE(i,2)+101)/201*p.hilakoko);
        
        p1 = pdf('Normal',1:1:p.hilakoko,px1,p.varianssi);  % p.varianssi
        p2 = pdf('Normal',1:1:p.hilakoko,px2,p.varianssi);
   
        
        aania = A.aania(b.ylenimi(i));
        viimesinnimi = A.nimi{b.ylenimi(i)};
        tarkastellaanehdokasta = i;
        %aania = 1;
        
        p.optimointihila = 0;  %110422 ei toimi, en löydä vikaa
        
        if p.optimointihila == 0
        
            p2d = p1'*p2.*aania;
            b.hila = b.hila + p2d;
            b.puoluehila(:,:,b.puolue(i)) = b.puoluehila(:,:,b.puolue(i)) + p2d;
        else
            % pakollinen, jos kuvia aiotaan käyttää
            
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
        aania = 0;
    end
    



b.varihila = zeros(p.hilakoko,p.hilakoko);
b.hilanvari = zeros(p.hilakoko,p.hilakoko);
for i=1:p.hilakoko
    for j=1:p.hilakoko
        b.varihila(i,j) = max(b.puoluehila(i,j,:));
        apu = find(b.puoluehila(i,j,:)==b.varihila(i,j));
        b.hilanvari(i,j) = apu(1);  % voi löytyä montakohtaa
        %if b.varihila(i,j)<0.000001
        %    b.varihila(i,j) = nan;
        %end
        
        
        
    end
end

% laitetaan varit kuosiin 
Zv = zeros(p.hilakoko,p.hilakoko,3);
p.puolueetL = {'-','IPU','KOK','KA','M2011','PS','PIR','KESK','SKP','KD','SDP','KTP','RKP','VP','VAS','VIHR'};

for i=1:p.hilakoko
    for j=1:p.hilakoko
        apu = find(strcmp(p.puolueetL, m.puolueet{b.hilanvari(i,j)}));
        %apu=6;
        if length(apu)>0 && p.logotkarttaan == 1
            % JUMAZUIKA - rgb siis 0..1, ja kuva0..255
            % voisko videotallennuskin johtua siitä? no selvitä joskus
            % tuli vain valkoista kuvaa 255,255,255 matriisia..
            
            Zv(i,j,:) = (double(logot{apu}(mod((80-i)*1,79)+1,mod((j)*1,79)+1,:)))/255;
            %for k=1:3
            %    Zv(i,j,4-k) = double(logot{apu}(mod((80-i)*1,79)+1,mod((j)*1,79)+1,k));
            %end
            
            
            if sum(Zv(i,j,:))>2.5
                % vois tehdä kans ehdon, että jos ollaan logon laidalla,
                % tai sitten omat logot photarilla
                % esim. Väri ja teksti - tai palattia.
                % esim. nyt valkonen pois SDP:stä
                Zv(i,j,:) = p.vari{b.hilanvari(i,j)};
            end
            
            
            
            
        else
            
            % ilman varjostusta
            %Zv(i,j,:) = p.vari{b.hilanvari(i,j)};
            
            
            Zv(i,j,:) = p.vari{b.hilanvari(i,j)};%+varinvarjostus*2;
          
        end
        
        if p.varinvarjostus == 1
            
              % Google matlab transperency example - tai opacity
            varinvarjostus = 1 - b.varihila(i,j)/sum(b.puoluehila(i,j,:));
            Zv(i,j,:) = Zv(i,j,:) + varinvarjostus*(3-sum(Zv(i,j,:)))/2; % vaaleita puolueita varjostetaan vähemmän
            for izz = 1:3
                if Zv(i,j,izz)>1
                    Zv(i,j,izz) = 1;
                end
            end
        end
    end
end



%figure(111);
Z = b.hila;
%surface(-p.hilaskaala:p.hilaskaala,-p.hilaskaala:p.hilaskaala,Z,Zv);

if p.logotkarttaan == 1
    p.lsk = 1;
    surface(-p.hilaskaala*p.lsk:p.lsk:p.hilaskaala*p.lsk,-p.hilaskaala*p.lsk:p.lsk:p.hilaskaala*p.lsk,Z,(Zv),'FaceColor','texturemap','EdgeColor','none', 'CDataMapping','direct');
else
    p.lsk = 1;   % lisää nämä teksteihin
    %surface(-p.hilaskaala:p.hilaskaala,-p.hilaskaala:p.hilaskaala,Z2,(Zv));
    surface(-p.hilaskaala*p.lsk:p.lsk:p.hilaskaala*p.lsk,-p.hilaskaala*p.lsk:p.lsk:p.hilaskaala*p.lsk,Z,(Zv),'edgecolor',[0.5 0.5 0.5]);
end


%title('suomalaiset: paljon kutakin arvopistettä äänestettiin. väritys kertoo sen, mikä oli todennäköisimmin äänestän puolue, ja paljon ääniä annettiin yhteensä arvolle');

axis tight;



if p.logotkarttaan == 0 || p.logotkarttaan == 1
for i=1:length(m.puolueet)
    %keta = find(b.puolue==i);
    keta = find(b.puolue(so2(ekaijollaaania:ii))==i);
    if length(keta)>2
        kxx = max(max(b.puoluehila(:,:,i)));
        kx = find(max(b.puoluehila(:,:,i)')==kxx);
        
        kyy = max(max(b.puoluehila(:,:,i)));
        ky = find(max(b.puoluehila(:,:,i))==kyy);
        
        
        varjataan = [0.3 0.3 0.3];
        tekstikorkeus =  Z(kx(1),ky(1));%*1.001;
        
        tekstikoko= 14;
        if b.puoluehila(kx(1),ky(1),i)<max(b.varihila(kx(1),ky(1),:))
            %kx(1)
            %max(b.puoluehila(kx(1),ky(1),:))
            varjataan = [0.7 0.7 0.7];
            %tekstikorkeus = max(b.varihila(kx(1),ky(1),:))*1.03;
            tekstikoko= 10;
        end
            
        
        %tekstikorkeusviiva =  Z(kx(1),ky(1))*1.2;
        tekstikorkeusviiva =  Z(kx(1),ky(1)) + max(max(Z))/12;
        
        %kx = round(mean(SCORE(keta,1)));
        %ky = round(mean(SCORE(keta,2)));
        % pinnassa koordinaatisto toistenpäin!
        %text((ky(1)-1)-p.hilaskaala,(kx(1)-1)-p.hilaskaala,Z(kx(1),ky(1))+0.01,m.puolueet{i},'color',[0 0 0],'fontsize',12,'horizontalalignment','center','verticalalignment','middle');
       %text((ky(1)-1)-p.hilaskaala,(kx(1)-1)-p.hilaskaala,Z(kx(1),ky(1))+0.01,m.puolueet{i},'color',[0 0 0],'BackgroundColor',p.vari{i},'fontsize',14,'horizontalalignment','center','verticalalignment','middle');
       
       line([(ky(1)-1)-p.hilaskaala (ky(1)-1)-p.hilaskaala],[(kx(1)-1)-p.hilaskaala (kx(1)-1)-p.hilaskaala],[tekstikorkeus tekstikorkeusviiva],'color',p.vari{i},'linewidth',2);
      
       text((ky(1)-1)-p.hilaskaala,(kx(1)-1)-p.hilaskaala,tekstikorkeusviiva,m.puolueet{i},'color',varjataan,'fontsize',tekstikoko,'horizontalalignment','center','verticalalignment','bottom');
       
       
       
   
       
       
        %text(ky,kx,Z(kx+101,ky+101)+0.1,m.puolueet{i},'color',[0 0 0],'BackgroundColor',p.vari{i},'fontsize',12,'horizontalalignment','center','verticalalignment','middle');
    end
    
        
end
end


      
        
        
        varjataan = [0.5 0.5 0.5];
        tekstikorkeus =  Z(px1, px2);
        tekstikorkeusviiva =  Z(px1,px2) + max(max(Z))/8;
        tekstikoko= 14;
        
            
        
        %kx = round(mean(SCORE(keta,1)));
        %ky = round(mean(SCORE(keta,2)));
        % pinnassa koordinaatisto toistenpäin!
        %text((ky(1)-1)-p.hilaskaala,(kx(1)-1)-p.hilaskaala,Z(kx(1),ky(1))+0.01,m.puolueet{i},'color',[0 0 0],'fontsize',12,'horizontalalignment','center','verticalalignment','middle');
       %text((ky(1)-1)-p.hilaskaala,(kx(1)-1)-p.hilaskaala,Z(kx(1),ky(1))+0.01,m.puolueet{i},'color',[0 0 0],'BackgroundColor',p.vari{i},'fontsize',14,'horizontalalignment','center','verticalalignment','middle');
       % pinnassa koordinaatisto toistenpäin!
       line([px2-1-p.hilaskaala px2-1-p.hilaskaala],[px1-1-p.hilaskaala px1-1-p.hilaskaala],[tekstikorkeus tekstikorkeusviiva],'color',p.vari{b.puolue(tarkastellaanehdokasta)},'linewidth',2);
       text(px2-1-p.hilaskaala,px1-1-p.hilaskaala,tekstikorkeusviiva,viimesinnimi,'color',varjataan,'fontsize',tekstikoko,'horizontalalignment','center','verticalalignment','bottom');
   

     
viewlaskuri = viewlaskuri + 0.5;
if viewlaskuri>180
    viewlaskuri = -179.5;
end


view(viewlaskuri,60+sin(10+viewlaskuri/50)*20); % vois tehdä videoon pyörivän
title(strcat('(c)2011 jaakko.talonen@aalto.fi (',num2str(aania),')'));
%axis off
%axis tight
  %axis equal
  %set(gca,'nextplot','replacechildren');
  
  framelaskuri = framelaskuri+1;
   saveas(gcf,strcat('C:\Users\Jaakko\Documents\vaaliJPG3\',num2str(framelaskuri),'.jpg'), 'jpg')
  
   disp(framelaskuri);
   %frame = getframe(gca,[1 1 300 300]); %[left bottom width height]
   %frame.colormap = p.omacolormap
   %image(f.cdata)
   %mov = addframe(mov,frame);
  
   %
   %F(framelaskuri) = getframe;
   
  %F(j) = getframe;
  
  %vikakuva auki
  if b.videoi>ii
      close;
  end
end  % FOR SILMULLA


%mov = close(mov)
%
%shading flat
%shading interp
%%
%




%%



