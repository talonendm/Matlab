

jiiS.vasenoikea = jiiS.vasenoikea - 100;
jiiS.libkon = jiiS.libkon - 100;
%%
ppp = 30;
pp.muutpallo = 0;
pp.kaikki = 0;
pp.aanet =2;  % VAS OIK = 0, 1= aanet, 2= indeksi
m.puolueetYLE = m.puolueet;
m.puolueetYLE{2} = 'IPU';
m.puolueetYLE{3} = 'KA';
m.puolueetYLE{12} = 'M2011';
m.puolueetYLE{20} = 'SEN';

teksti = '';
teksti2 = '';
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
    if length(CI{i2}.ehdokkaat)>0 && k3<200 && pp.aanet == 2  % VIRHE TÄSSÄ
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
            
            
        
            
            if pp.aanet == 1
            if i<length(CI{i2}.ehdokkaat) && A.valitaan(CI{i2}.ehdokkaat(i)) == 2
                disp(strrep(strcat(ehdokas,':¤',num2str(A.aania(CI{i2}.ehdokkaat(i))),','),'¤',' '));
            else
                if pp.kaikki == 1
                     if i<length(CI{i2}.ehdokkaat)
                        disp(strrep(strcat(ehdokas,':¤',num2str(A.aania(CI{i2}.ehdokkaat(i))),','),'¤',' '));
                     else
                        disp(strrep(strcat(ehdokas,':¤',num2str(A.aania(CI{i2}.ehdokkaat(i)))),'¤',' '));
                     end
                else
                    summa = summa + A.aania(CI{i2}.ehdokkaat(i));
                    if i == length(CI{i2}.ehdokkaat) && pp.muutpallo == 1
                        disp(strrep(strcat(CI{i2}.puolue,':¤',num2str(summa)),'¤',' '));
                    end
                end
            end
            elseif pp.aanet == 0
                % VARI VAS-OIK
            if i<length(CI{i2}.ehdokkaat) && A.valitaan(CI{i2}.ehdokkaat(i)) == 2
                
                jsarvo = round(rand*200); % TÄNNE vas-oik arvo
                
                disp(strrep(strcat(ehdokas,':¤',num2str(jsarvo),','),'¤',' '));
            else
                if pp.kaikki == 1
                     if i<length(CI{i2}.ehdokkaat)
                        jsarvo = round(rand*200); % TÄNNE vas-oik arvo
                        disp(strrep(strcat(ehdokas,':¤',num2str(jsarvo),','),'¤',' '));
                     else
                        jsarvo = round(rand*200); % TÄNNE vas-oik arvo
                        disp(strrep(strcat(ehdokas,':¤',num2str(jsarvo)),'¤',' '));
                     end
                else
                    % jsarvo  (EI VIELÄ TÄÄLLÄ)
                    summa = summa + A.aania(CI{i2}.ehdokkaat(i));
                    if i == length(CI{i2}.ehdokkaat) && pp.muutpallo == 1
                        disp(strrep(strcat(CI{i2}.puolue,':¤',num2str(summa)),'¤',' '));
                    end
                end
            end    
                
               
            else
                     % VARI #1 ... 2315
            if i<length(CI{i2}.ehdokkaat) && A.valitaan(CI{i2}.ehdokkaat(i)) == 2
                k3 = k3 +1;
                jsarvo = k3;
                
                % TÄNNE vas-oik arvo
                %disp(k3)
                teksti = strcat(teksti,',',num2str(jiiS.vasenoikea(CI{i2}.ehdokkaat(i))));
                teksti2 = strcat(teksti2,',',num2str(jiiS.libkon(CI{i2}.ehdokkaat(i))));
                
               % if A.valitaan(CI{i2}.ehdokkaat(i+1)) < 2
                %    disp(strrep(strcat(ehdokas,':¤',num2str(jsarvo)),'¤',' '));
                %else
                    disp(strrep(strcat(ehdokas,':¤',num2str(jsarvo),','),'¤',' '));
               % end
            else
                if pp.kaikki == 1
                     if i<length(CI{i2}.ehdokkaat)
                        k3 = k3 +1;
                        jsarvo = k3; % TÄNNE vas-oik arvo
                        teksti = strcat(teksti,',',num2str(jiiS.vasenoikea(CI{i2}.ehdokkaat(i))));
                        teksti2 = strcat(teksti2,',',num2str(jiiS.libkon(CI{i2}.ehdokkaat(i))));
                        %teksti2
                        disp(strrep(strcat(ehdokas,':¤',num2str(jsarvo),','),'¤',' '));
                     else
                         k3 = k3 +1;
                        jsarvo = k3; % TÄNNE vas-oik arvo
                        teksti = strcat(teksti,',',num2str(jiiS.vasenoikea(CI{i2}.ehdokkaat(i))));
                        teksti2 = strcat(teksti2,',',num2str(jiiS.libkon(CI{i2}.ehdokkaat(i))));
                        teksti2
                        disp(strrep(strcat(ehdokas,':¤',num2str(jsarvo)),'¤',' '));
                     end
                elseif pp.aanet == 2 && A.valitaan(CI{i2}.ehdokkaat(i)) == 2
                    teksti = strcat(teksti,',',num2str(jiiS.vasenoikea(CI{i2}.ehdokkaat(i))));
                    teksti2 = strcat(teksti2,',',num2str(jiiS.libkon(CI{i2}.ehdokkaat(i))));
                     
                else
                    % jsarvo  (EI VIELÄ TÄÄLLÄ)
                    summa = summa + A.aania(CI{i2}.ehdokkaat(i));
                    if i == length(CI{i2}.ehdokkaat) && pp.muutpallo == 1
                        disp(strrep(strcat(CI{i2}.puolue,':¤',num2str(summa)),'¤',' '));
                        
                    end
                end
                
            end    
                
            
            end
            
            
            
        end
        
        if i3<ppp || k3==200 && pp.aanet == 2
            if yyy(i3+1) == 0 || k3==200 && pp.aanet == 2
                disp('}');   
            else
                disp('},');   
            end
        else
            %disp('}');   
        end
    end
end
%disp('}');   
disp('};');  

teksti2

%var numero = new Array(1,2,3,4,5); 
teksti = strcat('var libkon = new Array(',teksti(2:end),');');
disp(teksti);

teksti2 = strcat('var vasenoikea = new Array(',teksti2(2:end),');');
disp(teksti2);

%%

        disp('KOK: {');
        for i=1:length(CI{2}.ehdokkaat)
            
            ehdokas = A.nimi{CI{2}.ehdokkaat(i)};
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
            
            if i<length(CI{2}.ehdokkaat)
                disp(strrep(strcat(ehdokas,':¤',num2str(A.aania(CI{2}.ehdokkaat(i))),','),'¤',' '));
            else
                disp(strrep(strcat(ehdokas,':¤',num2str(A.aania(CI{2}.ehdokkaat(i)))),'¤',' '));
            
        
            end
        end
        disp('}');   
       
disp('};');    





%%
disp('var flare = {');
    disp('  Government: {');
        disp('KESK: {');
        for i=1:length(CI{1}.ehdokkaat)
            ehdokas = A.nimi{CI{1}.ehdokkaat(i)};
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
            
            if i<length(CI{1}.ehdokkaat)
                disp(strrep(strcat(ehdokas,':¤',num2str(A.aania(CI{1}.ehdokkaat(i))),','),'¤',' '));
            else
                disp(strrep(strcat(ehdokas,':¤',num2str(A.aania(CI{1}.ehdokkaat(i)))),'¤',' '));
            
        
            end
        end
        disp('},');   
        disp('KOK: {');
        for i=1:length(CI{2}.ehdokkaat)
            
            ehdokas = A.nimi{CI{2}.ehdokkaat(i)};
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
            
            if i<length(CI{2}.ehdokkaat)
                disp(strrep(strcat(ehdokas,':¤',num2str(A.aania(CI{2}.ehdokkaat(i))),','),'¤',' '));
            else
                disp(strrep(strcat(ehdokas,':¤',num2str(A.aania(CI{2}.ehdokkaat(i)))),'¤',' '));
            
        
            end
        end
        disp('}');   
    disp('}');     
disp('};');    
   
   
        
     %%   
        
%        AgglomerativeCluster: 3938,
%      CommunityStructure: 3812,
%      HierarchicalCluster: 6714,
%      MergeEdge: 743
    },
    graph: {
      BetweennessCentrality: 3534,
      LinkDistance: 5731,
      MaxFlowMinCut: 7840,
      ShortestPaths: 5914,
      SpanningTree: 3416
    },
    optimization: {
      AspectRatioBanker: 7074
    }
  },
  animate: {
    Easing: 17010,
    FunctionSequence: 5842,
    interpolate: {
      ArrayInterpolator: 1983,
      ColorInterpolator: 2047,
      DateInterpolator: 1375,
      Interpolator: 8746,
      MatrixInterpolator: 2202,
      NumberInterpolator: 1382,
      ObjectInterpolator: 1629,
      PointInterpolator: 1675,
      RectangleInterpolator: 2042
    },
    ISchedulable: 1041,
    Parallel: 5176,
    Pause: 449,
    Scheduler: 5593,
    Sequence: 5534,
    Transition: 9201,
    Transitioner: 19975,
    TransitionEvent: 1116,
    Tween: 6006
  }
  
};