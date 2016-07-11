function [hs A] = KorjaaTietokantojaW7(hs,A)



paikka = find(strcmp(hs.Nimi,'Elomaa Kike'));

hs.Nimi{paikka} = 'Elomaa Kike (Ritva)';
 % sama kuin tuloksissa [YLE]


paikka = find(strcmp(hs.Nimi,'Gutiérrez Sorainen Ana Maria'))

hs.Nimi{paikka} = 'Gutiérrez Sorainen Ana'; % sama kuin tuloksissa [YLE]


paikka = find(strcmp(hs.Nimi,'Czitrom Belá'));

hs.Nimi{paikka} = 'Czitrom Béla'; % sama kuin tuloksissa [YLE]


paikka = find(strcmp(hs.Nimi,'Söderholm Charlotta'));

hs.Nimi{paikka} = 'Söderholm Charlotta (Lotta)'; % sama kuin tuloksissa [YLE]


paikka = find(strcmp(hs.Nimi,'Laasonen Jan-Markus'));

hs.Nimi{paikka} = 'Laasonen Jan-Markus (Jani)'; % sama kuin tuloksissa [YLE]


paikka = find(strcmp(hs.Nimi,'Elmgren Oskar'));

hs.Nimi{paikka} = 'Elmgren Oskar (Jeppe)'; % sama kuin tuloksissa [YLE]


paikka = find(strcmp(hs.Nimi,'Paldán Maarit'));

hs.Nimi{paikka} = 'Paldan Maarit'; % sama kuin tuloksissa [YLE]


paikka = find(strcmp(hs.Nimi,'Jernström Kristian(Kisu)'));

hs.Nimi{paikka} = 'Jernström Kristian (Kisu)'; % sama kuin tuloksissa [YLE]


paikka = find(strcmp(hs.Nimi,'Mavromichalis Nicolas'));

hs.Nimi{paikka} = 'Mavromichalis Nicolas (Black Mike)'; % sama kuin tuloksissa [YLE]



hs.Nimi{1744} = 'Männikkö Tarja';  % koodattu Mannikkö

A.nimi{1837} = 'Pakkasvirta Anneli (Anneli Sauli)';

%paikka = find(strcmp(hs.Nimi,'Pakkasvirta Anneli (Anneli Sauli)'));

%hs.Nimi{paikka} = 'Pakkasvirta Anneli (Anneli Sauli'; % sama kuin tuloksissa [YLE]



%"Jokila Kirsi ""Killi"""

A.nimi{956} = 'Jokila Kirsi (Killi)';

paikka = find(strcmp(hs.Nimi,'Jokila Kirsi "Killi"'));

hs.Nimi{paikka} = 'Jokila Kirsi (Killi)'; % ei vastannut muutenkaan 


A.nimi{74} = 'Virtanen Pertti (Veltto)';

paikka = find(strcmp(hs.Nimi,'Virtanen Pertti "Veltto"'));

hs.Nimi{paikka} = 'Virtanen Pertti (Veltto)'; % ei vastannut muutenkaan hs.
%110516

% muutetaan yksittaiset MUUT - kuten YLElla
%find(strcmp(hs.'Helminen'))



hs.YLEpuolue = hs.Puolue;

hs.YLEpuolue{find(strcmp(hs.Puolue,'Helminen'))} = 'MUUT';

hs.YLEpuolue{find(strcmp(hs.Puolue,'Korhonen'))} = 'MUUT';

hs.YLEpuolue{find(strcmp(hs.Puolue,'Kotajärvi'))} = 'MUUT';

hs.YLEpuolue{find(strcmp(hs.Puolue,'Laajola'))} = 'MUUT';

hs.YLEpuolue{find(strcmp(hs.Puolue,'Lönnroth'))} = 'MUUT';

hs.YLEpuolue{find(strcmp(hs.Puolue,'Marttunen'))} = 'MUUT';

hs.YLEpuolue{find(strcmp(hs.Puolue,'Nikoskinen'))} = 'MUUT';

hs.YLEpuolue{find(strcmp(hs.Puolue,'Niskanen'))} = 'MUUT';

hs.YLEpuolue{find(strcmp(hs.Puolue,'Salminen'))} = 'MUUT';

hs.YLEpuolue{find(strcmp(hs.Puolue,'Simonen'))} = 'MUUT';

hs.YLEpuolue{find(strcmp(hs.Puolue,'Valjakka'))} = 'MUUT';

apu = find(strcmp(hs.Puolue,'Yl-Sit'));

hs.YLEpuolue{apu(1)} = 'MUUT';

hs.YLEpuolue{apu(2)} = 'MUUT';


for i=1:length(hs.YLEpuolue)
    
if strcmp(hs.YLEpuolue{i},'IPU')
    
   hs.YLEpuolue{i} = 'ITSP';
 
   end
   
 if strcmp(hs.YLEpuolue{i},'KA')
   
     hs.YLEpuolue{i} = 'KÖY';
 
   end
 
   if strcmp(hs.YLEpuolue{i},'M2011')
      
  hs.YLEpuolue{i} = 'M11';
  
  end
   
 if strcmp(hs.YLEpuolue{i},'SEN')
    
    hs.YLEpuolue{i} = 'SSP';
  
  end
  
    
    
end



% SUORAAN muokataan YLEe 1105016

for i=1:length(A.nimi)
  
  if strcmp(A.puolue{i},'ITSP')
    
    A.puolue{i} = 'IPU';
  
  end
  
  if strcmp(A.puolue{i},'SSP')
 
       A.puolue{i} = 'SEN';
  
  end
   
 if strcmp(A.puolue{i},'M11')
     
   A.puolue{i} = 'M2011';
  
  end
 
   if strcmp(A.puolue{i},'KÖY')
  
      A.puolue{i} = 'KA';
  
  end

end
