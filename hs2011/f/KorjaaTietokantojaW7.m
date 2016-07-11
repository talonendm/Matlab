function [hs A] = KorjaaTietokantojaW7(hs,A)



paikka = find(strcmp(hs.Nimi,'Elomaa Kike'));

hs.Nimi{paikka} = 'Elomaa Kike (Ritva)';
 % sama kuin tuloksissa [YLE]


paikka = find(strcmp(hs.Nimi,'Guti�rrez Sorainen Ana Maria'))

hs.Nimi{paikka} = 'Guti�rrez Sorainen Ana'; % sama kuin tuloksissa [YLE]


paikka = find(strcmp(hs.Nimi,'Czitrom Bel�'));

hs.Nimi{paikka} = 'Czitrom B�la'; % sama kuin tuloksissa [YLE]


paikka = find(strcmp(hs.Nimi,'S�derholm Charlotta'));

hs.Nimi{paikka} = 'S�derholm Charlotta (Lotta)'; % sama kuin tuloksissa [YLE]


paikka = find(strcmp(hs.Nimi,'Laasonen Jan-Markus'));

hs.Nimi{paikka} = 'Laasonen Jan-Markus (Jani)'; % sama kuin tuloksissa [YLE]


paikka = find(strcmp(hs.Nimi,'Elmgren Oskar'));

hs.Nimi{paikka} = 'Elmgren Oskar (Jeppe)'; % sama kuin tuloksissa [YLE]


paikka = find(strcmp(hs.Nimi,'Pald�n Maarit'));

hs.Nimi{paikka} = 'Paldan Maarit'; % sama kuin tuloksissa [YLE]


paikka = find(strcmp(hs.Nimi,'Jernstr�m Kristian(Kisu)'));

hs.Nimi{paikka} = 'Jernstr�m Kristian (Kisu)'; % sama kuin tuloksissa [YLE]


paikka = find(strcmp(hs.Nimi,'Mavromichalis Nicolas'));

hs.Nimi{paikka} = 'Mavromichalis Nicolas (Black Mike)'; % sama kuin tuloksissa [YLE]



hs.Nimi{1744} = 'M�nnikk� Tarja';  % koodattu Mannikk�

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

hs.YLEpuolue{find(strcmp(hs.Puolue,'Kotaj�rvi'))} = 'MUUT';

hs.YLEpuolue{find(strcmp(hs.Puolue,'Laajola'))} = 'MUUT';

hs.YLEpuolue{find(strcmp(hs.Puolue,'L�nnroth'))} = 'MUUT';

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
   
     hs.YLEpuolue{i} = 'K�Y';
 
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
 
   if strcmp(A.puolue{i},'K�Y')
  
      A.puolue{i} = 'KA';
  
  end

end
