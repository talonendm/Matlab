

jiiS.vasenoikea = jiiS.vasenoikea - 100;
jiiS.libkon = jiiS.libkon - 100;
%%

jsydin = zeros(2315,1);




for i=1:length(vastannutkyssareihin)
    jsydin(b.ylenimi(i)) = -b.vastaustaulu(vastannutkyssareihin(i),8)*(b.binaarimassa(vastannutkyssareihin(i),8)+2) + b.vastaustaulu(vastannutkyssareihin(i),9)*(b.binaarimassa(vastannutkyssareihin(i),9)+2);
    
end

%%

ppp = 30;
pp.muutpallo = 0;
pp.kaikki = 1;
pp.aanet =2;  % VAS OIK = 0, 1= aanet, 2= indeksi
m.puolueetYLE = m.puolueet;
m.puolueetYLE{2} = 'IPU';
m.puolueetYLE{3} = 'KA';
m.puolueetYLE{12} = 'M2011';
m.puolueetYLE{20} = 'SEN';

teksti = '';
teksti2 = '';
teksti3 = '';


tekstiK = '';
tekstiKG = '';


%for i=1:200
    %kokoran(i) = round(rand*1000)+1;
 %   
%end



for i=1:ppp
    j=0+i;
    CI{i}.puolue = m.puolueetYLE{j};
    CI{i}.ehdokkaat = find(strcmp(A.puolue,CI{i}.puolue));
    
end


yhteensaaania = zeros(30,1);
for i=1:30
    if length(CI{i}.ehdokkaat)>0
        yhteensaaania(i) = sum(A.aania(CI{i}.ehdokkaat));
    end
end

[yyy aaa] = sort(yhteensaaania,'descend');






% ehdokkaat - sen jälkeen hallitus

disp('*********************');

disp('var parties = {');


k3 = 0;%0;
for i3 = 1:ppp
    i2 = aaa(i3);
    if length(CI{i2}.ehdokkaat)>0 
    disp(strrep(strcat(CI{i2}.puolue,'¤',num2str((sum(A.aania(CI{i2}.ehdokkaat)))),': {'),'¤','_'));
    
    summa = 0;
        for i=1:length(CI{i2}.ehdokkaat)
            ehdokas = A.nimi{CI{i2}.ehdokkaat(i)};
            ehdokas = strrep(ehdokas,'ä','a');
            ehdokas = strrep(ehdokas,'ö','o');
            ehdokas = strrep(ehdokas,'Ä','A');
            ehdokas = strrep(ehdokas,'Ö','O');
            ehdokas = strrep(ehdokas,'Å','A');
            ehdokas = strrep(ehdokas,'å','a');
            ehdokas = strrep(ehdokas,' ','_');
            ehdokas = strrep(ehdokas,'é','e');
            ehdokas = strrep(ehdokas,'-','');
            ehdokas = strrep(ehdokas,'(','');
            ehdokas = strrep(ehdokas,')','');
            ehdokas = strrep(ehdokas,'.','');
            ehdokas = strrep(ehdokas,'ó','o');
            ehdokas = strrep(ehdokas,'á','a');
            
            
        
            
                  
                     if i<length(CI{i2}.ehdokkaat)
                        k3 = k3 +1;
                        jsarvo = k3; % TÄNNE vas-oik arvo
                        teksti = strcat(teksti,',',num2str(jiiS.vasenoikea(CI{i2}.ehdokkaat(i))));
                        teksti2 = strcat(teksti2,',',num2str(jiiS.libkon(CI{i2}.ehdokkaat(i))));
                        
                        % AANET
                        teksti3 = strcat(teksti3,',',num2str(A.aania(CI{i2}.ehdokkaat(i)))); 
                        
                        % kysymys
                        tekstiK = strcat(tekstiK,',',num2str(jsydin(CI{i2}.ehdokkaat(i))));
                        
                        
                        %teksti2
                        disp(strrep(strcat(ehdokas,':¤',num2str(jsarvo),','),'¤',' '));
                     else
                         k3 = k3 +1;
                        jsarvo = k3; % TÄNNE vas-oik arvo
                        teksti = strcat(teksti,',',num2str(jiiS.vasenoikea(CI{i2}.ehdokkaat(i))));
                        teksti2 = strcat(teksti2,',',num2str(jiiS.libkon(CI{i2}.ehdokkaat(i))));
                        teksti3 = strcat(teksti3,',',num2str(A.aania(CI{i2}.ehdokkaat(i)))); 
                        %teksti2
                        disp(strrep(strcat(ehdokas,':¤',num2str(jsarvo)),'¤',' '));
                        
                         % kysymys
                        tekstiK = strcat(tekstiK,',',num2str(jsydin(CI{i2}.ehdokkaat(i))));
                        
                        
                     end
            
            
            
        end
        
       % if i3<ppp || k3==200 && pp.aanet == 2
       %     if yyy(i3+1) == 0 || k3==200 && pp.aanet == 2
       %         disp('}');   
        %    else
                disp('},');   
        %    end
        %else
            %disp('}');   
        %end
    end
end
%disp('}');   
disp('};');  

%teksti2

%var numero = new Array(1,2,3,4,5); 
teksti = strcat('var libkon = new Array(',teksti(2:end),');');
disp(teksti);

teksti2 = strcat('var vasenoikea = new Array(',teksti2(2:end),');');
disp(teksti2);

teksti3 = strcat('var koko = new Array(',teksti3(2:end),');');
disp(teksti3);

tekstiK = strcat('var kysymys3 = new Array(',tekstiK(2:end),');');
disp(tekstiK);


%%
% SIIRRÄ TOHON YLÄKOODIIN 20.5.2011

% hs kysymyksey JS:ään



for i=1:2315
    
    
end



% PAHUS - js.tiedostot on sortattu.. täytyy vääntää noihin vanhoihin
% tiedostoihin.

% LATER !!



%%
