function IDAfigure(p,b,A,hs,m,SCORE,Zv,kuvanro)
% *************************************************************************
% IDA paper figures 1-3 tai 3d = kuvanro = 0
% 110516
% *************************************************************************

% *************************************************************************
u = figure(100);
% *************************************************************************
if p.idapaperfigure > 0 && p.rotaatio == 1
    venytyspcassa = (-m.pca.min(1) + m.pca.max(1))/(-m.pca.min(2) + m.pca.max(2));
    set(u,'position',[0 0 800 800*venytyspcassa])
else
    venytyspcassa = (-m.pca.min(1) + m.pca.max(1))/(-m.pca.min(2) + m.pca.max(2));
    set(u,'position',[0 0 800*venytyspcassa 800])
end
% *************************************************************************
if kuvanro == 0
    p.puoluekuvatasovaikorkeus = 0;  % 3d report part1
    p.ehdokastarkkailu = 0;
else
    p.puoluekuvatasovaikorkeus = 1;
end
% ida
if kuvanro == 1
    p.ehdokastarkkailu = 2;  % 2 = IDA fig3 plot kansa
    % IDA kuvat 1 2
    piir2=1; % ehdokkaat % PAKOTA TÄLLÄ JOS SUMMA KÄYRÄÄ!
elseif kuvanro == 2
    p.ehdokastarkkailu = 1;  % 2 = IDA fig2 plot ehdokkaat
    % IDA kuvat 1 2
    piir2=1; % ehdokkaat % PAKOTA TÄLLÄ JOS SUMMA KÄYRÄÄ!
elseif kuvanro == 3
    p.ehdokastarkkailu = 0;  % 200 eduskunta
    piir2 = [4 5 6 17 18 19 25 26];
end
% *************************************************************************
if p.puoluekuvatasovaikorkeus == 0
    
    if p.puoluekohtainenpainotus>0
        piir = p.puoluekohtainenpainotus
    else
        piir = 5;%26;   %1 = helminen 1773hs datassa, ei vastannut mihink��n..
    end
    Z = b.puoluehila(:,:,piir);
    surface(-p.hilaskaala:p.hilaskaala,-p.hilaskaala:p.hilaskaala,Z,Zv);
    view(65,57)
    shading interp
    axis tight;
    
    
    % *************************************************************************
    % tarkistus
    puoluepoikkeus =0;
    % *************************************************************************
    puolueHS2YLE = find(strcmp(p.puolueetYLEtulosUnique,m.puolueet{piir}));
    if length(puolueHS2YLE)==0
        if piir ==2
            disp(m.puolueet{piir});
            puolueHS2YLE = find(strcmp(p.puolueetYLEtulosUnique,'ITSP'));
            disp(p.puolueetYLEtulosUnique{puolueHS2YLE});
        elseif piir ==3
            disp(m.puolueet{piir});
            puolueHS2YLE = find(strcmp(p.puolueetYLEtulosUnique,'KÖY'));
            disp(p.puolueetYLEtulosUnique{puolueHS2YLE});
        elseif strcmp(m.puolueet{piir},'M2011')
            puolueHS2YLE = find(strcmp(p.puolueetYLEtulosUnique,'M11'));
    
        elseif strcmp(m.puolueet{piir},'SEN')
            puolueHS2YLE = find(strcmp(p.puolueetYLEtulosUnique,'SSP'));   
        else
            %onkomuut?
            puolueHS2YLE = find(strcmp(p.puolueetYLEtulosUnique,'MUUT'));    
            disp('virhe'); 
            if piir==1 || piir == 8
                % Helminen Jyrki
                puoluepoikkeus =1;
            end
        end
    end
    % *************************************************************************
    if puoluepoikkeus == 1
        if piir ==1
            Zaaniatuli = A.aania(strcmp(A.nimi,'Helminen Jyrki'));
        end
        if piir == 8
            Zaaniatuli = A.aania(strcmp(A.nimi,'Korhonen Juho')); 
        end
            
    else
        Zaaniatuli = sum(A.aania(strcmp(A.puolue,p.puolueetYLEtulosUnique{puolueHS2YLE})));  % HUOM! VP T�SS�
    end
    Ztilavuus = sum(sum(Z));
    disp(strcat('Total votes [YLE]:',num2str(Zaaniatuli)));
    disp(strcat('Volume:',num2str(Ztilavuus)));
    %disp(Ztilavuus);
    %disp('ero johtuu siit�, ett� joku puuttuu.. yli menev� tilavuus on skaalattu: (+) liikaa tilavuutta ');
    disp(strcat('Yle votes - Volume V =',num2str(round(Zaaniatuli - Ztilavuus))));
    title(m.puolueet{piir});
    
    % latexiin
    disp('\begin{figure}');
    disp('\centering');
    disp(strcat('\includegraphics[width=150mm]{part1/',num2str(piir),'.eps}'))
    disp(strrep(strcat('\caption{Value grid for Finnish cizizen, who gave vote for the�',m.puolueet{piir},'-party. Total votes:�',num2str(num2str(Zaaniatuli))    ,'}'),'�',' '));
    disp(strcat('\label{fig:',num2str(piir),'}'));
    disp('\end{figure}');


    % latexiin 2 kuvaa p��lkekk�in
    disp('....');
    disp('\begin{figure}[h]');
    disp('\begin{center}');
    disp('\begin{tabular}{c}');
    disp(strcat('\includegraphics[width=0.9\columnwidth]{part1/',num2str(piir),'.eps}\\'))
    disp(strrep(strcat('Value grid for Finnish cizizen, who gave vote for the�',m.puolueet{piir},'-party. Total votes:�',num2str(num2str(Zaaniatuli))    ,'\\'),'�',' '));
    disp(strcat('\includegraphics[width=0.9\columnwidth]{part1/',num2str(piir),'.eps}\\'))
    disp(strrep(strcat('Value grid for Finnish cizizen, who gave vote for the�',m.puolueet{piir},'-party. Total votes:�',num2str(num2str(Zaaniatuli))    ,''),'�',' '));
    disp('\end{tabular}');
    disp(strcat('\label{fig:',num2str(piir),'}'));
    disp('\end{center}');
    disp('\end{figure}');
    
    %shading faceted
    %shading flat
    %shading interp

    % *************************************************************************
    p.tallennuspuoluekuva = 0;
    if p.tallennuspuoluekuva == 1
        cd(p.fig); % nopeuttaa tallentamista
        saveas(gcf,strcat('',num2str(piir),'.eps'), 'psc2'); % save as eps with colors - google
        close;
        cd ..;
    end
% ******************************************************************
    
else
    
    % Colormap toisten päin, ja valeempia tummenetaan
    gray2 = max(0,(flipud(gray))-0.3);
    
    hold on;
    
    % ********************************************************************
    % IDA kuva 3 monta tasoa
    m.puolueet(piir2)
    % ********************************************************************
    for i=1:length(piir2)
        piir = piir2(i);

        if piir2 ==1
            % kuva 1-2
            % MIETI ONKO OK!!!
            Z = b.hila(:,:)/200*2315; % eiku oikein, summa on.. VÄÄRIN.. tässä vain maksimi
            %Z = b.puoluehila(:,:,1);
            %for ii = 2:29
            %    Z(:,:) = Z(:,:) + b.puoluehila(:,:,ii);
            %    disp(ii)
            %end
        else
            % kuva 3
            Z = b.puoluehila(:,:,piir);   
        end
        
        %Z = fliplr(b.puoluehila(:,:,piir)'); % käännellään kuvaa varten
        
        
        %[Co h] = contour(-p.hilaskaala:p.hilaskaala,-p.hilaskaala:p.hilaskaala,Z);
        if p.kansalaistarkkailu == 1
            
            kansakuvapiirretaan = 3;
            if kansakuvapiirretaan == 1
                Z = Z * (1 /  max(pdf('Normal',1:1:p.hilakoko,50,p.varianssi)))^2; % skaalataan siten, että huippu 1 kansalainen.
                [Co h] = contour(-p.hilaskaala:p.hilaskaala,-p.hilaskaala:p.hilaskaala,Z,[10000 20000 50000]);
                
            elseif kansakuvapiirretaan == 2
                % tuhansissa
                Z = (Z * (1 /  max(pdf('Normal',1:1:p.hilakoko,50,p.varianssi)))^2)./1000; % skaalataan siten, että huippu 1 kansalainen.
                [Co h] = contour(-p.hilaskaala:p.hilaskaala,-p.hilaskaala:p.hilaskaala,Z,[20 40 60]);
            elseif kansakuvapiirretaan == 3
                % kaikki edustajat kansalaistarkkailu hilasuma
                Z = (Z * (1 /  max(pdf('Normal',1:1:p.hilakoko,50,p.varianssi)))^2); % skaalataan siten, että huippu 1 kansalainen.
                Z = Z/1000;
                    disp('poo')
                    [Co h] = contour(-p.hilaskaala:p.hilaskaala,-p.hilaskaala:p.hilaskaala,Z,[1 10 100 500 750 1000]);
                
                %[Co h] = contour(-p.hilaskaala:p.hilaskaala,-p.hilaskaala:p.hilaskaala,Z,[10 40 60]);    
            end
        else
            if length(piir2)>1
                % varianssit ym. vaikuttaa tähän
                Z = Z * (1 /  max(pdf('Normal',1:1:p.hilakoko,50,p.varianssi)))^2; % skaalataan siten, että huippu 1.
                [Co h] = contour(-p.hilaskaala:p.hilaskaala,-p.hilaskaala:p.hilaskaala,Z,[2 4 6 10]);
                disp('mmm')
            else
                Z = Z * (1 /  max(pdf('Normal',1:1:p.hilakoko,50,p.varianssi)))^2; % skaalataan siten, että huippu 1.
                if p.kansalaistarkkailu == 0
                    [Co h] = contour(-p.hilaskaala:p.hilaskaala,-p.hilaskaala:p.hilaskaala,Z,[1 100 300 500 700]);
                    disp('IDA2 fig2. ehdokkaat');
                else
                    Z = Z/1000;
                    disp('poo')
                    [Co h] = contour(-p.hilaskaala:p.hilaskaala,-p.hilaskaala:p.hilaskaala,Z,[1 5 10 100]);
                end
                disp('Kaikkien ehdokkaiden arvot')
            end
                
                
        end
        text_handle = clabel(Co,h);
        %text_handle2 =round(text_handle);
        
        if p.bw == 0
            set(text_handle,'color',p.vari{piir});
        else
            set(text_handle,'color',[0 0 0]);
        end
        
        % taustavarilla
        %set(text_handle,'BackgroundColor',sqrt(p.vari{piir}),...
        %'Edgecolor',p.vari{piir})
    end
    
     
    coKy = b.kansanedustaja.Y-1-p.hilaskaala;
    coKx = b.kansanedustaja.X-1-p.hilaskaala;
    %coKx = fliplr(coKx);
    if p.kansalaistarkkailu == 0 && p.ehdokastarkkailu == 0
    for i=1:length(b.kansanedustaja.puolue)
        apu = find(strcmp(m.puolueet,b.kansanedustaja.puolue{i}));
        if length(apu) == 0
            % naucler
            apu = 18;
            %b.kansanedustaja.puolue{i} = 'RKP(a)';
            b.kansanedustaja.puolue{i} = 'RKP';
        end
        if b.kansanedustaja.approksimoitu(i) == 1
            tyyli = 'normal';
            teksti = strcat('(',b.kansanedustaja.puolue{i},')');
        else
            tyyli = 'normal';
            teksti = b.kansanedustaja.puolue{i};
        end
        if b.kansanedustaja.puhis(i) == 1
            koko = 16;
        else
            koko = 8;
        end
        
        if p.bw == 1
            text(coKy(i),coKx(i),teksti,'fontsize',koko,'color',[0 0 0],'fontangle',tyyli);
        else
            text(coKy(i),coKx(i),teksti,'fontsize',koko,'color',p.vari{apu},'fontangle',tyyli);
        end
        
    end
    end
    hold off;
    
    
    axis([-100 100 -100 100]);
    %set(h,'ShowText','on','TextStep',(get(h,'LevelStep')*2))
    colormap(gray2)
    if p.rotaatio == 0
        view(-90,90)
    else
        %ylabel(strcat(p.rotaatiopuolue1,'-',p.rotaatiopuolue2));
    end
    
    
    
end
if p.rotaatio == 0
    ylabel('First component');
    xlabel('Second component');
else
    xlabel(strcat(p.rotaatiopuolue1,'-',p.rotaatiopuolue2));
end

if p.ehdokastarkkailu == 2
    % KESKIARVOT EHDOKKAIDEN
    for i=1:length(m.puolueet)
        keta = find(b.puolue==i);
        if length(keta)>2
            kx = median(SCORE(keta,1));
            ky = median(SCORE(keta,2));
            %text(ky,kx,m.puolueet{i},'color',[0 0 0],'BackgroundColor',p.vari{i},'fontsize',8,'horizontalalignment','center','verticalalignment','middle');

            text(ky,kx,m.puolueet{i},'color',[0 0 0],'fontsize',16,'horizontalalignment','center','verticalalignment','middle');
        end    
    end
end
if p.ehdokastarkkailu == 1
    for i=1:length(m.puolueet)
        keta = find(b.puolue==i);
        if length(keta)>2
            kxx = max(max(b.puoluehila(:,:,i)));
            kx = find(max(b.puoluehila(:,:,i)')==kxx);
        
            kyy = max(max(b.puoluehila(:,:,i)));
            ky = find(max(b.puoluehila(:,:,i))==kyy);
        
        
            varjataan = [0 0 0];
            tekstikorkeus =  Z(kx(1),ky(1));%*1.001;
        
            tekstikoko= 16;
            if b.puoluehila(kx(1),ky(1),i)<max(b.varihila(kx(1),ky(1),:))
                varjataan = [0 0 0];
                tekstikoko= 10;
            end
       
            tekstikorkeusviiva =  Z(kx(1),ky(1)) + max(max(Z))/6;
   
            % esim. ida3-kuvassa - ei pikkupuoluieita
            if (ky(1)-1)-p.hilaskaala == -100 || sum(sum(b.puoluehila(:,:,i)))==0
            else
                text((ky(1)-1)-p.hilaskaala,(kx(1)-1)-p.hilaskaala,m.puolueet{i},'color',varjataan,'fontsize',tekstikoko,'horizontalalignment','center','verticalalignment','bottom');
            end
        end
    end
end





% *************************************************************************
end