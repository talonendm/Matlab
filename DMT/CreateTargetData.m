function [d p] = CreateTargetData(d,p);

if length(p.variables)==0 
    d.target.variables = SelectAlmostAll(p.muuttujassamuutoksia, d.ori.g.data(p.aloitus:p.tiheys:p.lopetus,:));
else
    d.target.variables = p.variables;
end
d.target.nimet = deal(d.ori.g.names(d.target.variables)); 
d.target.selitys = deal(d.ori.g.exps(d.target.variables)); 
d.target.t = d.ori.g.time(p.aloitus:p.tiheys:p.lopetus);
d.target.yksikko = deal(d.ori.g.units(d.target.variables));

% *************************************************************************
%d.target.esikasittely = 0; % 0/1
% *************************************************************************
cd(p.tyoskentelyhakemisto);
d.target.aX = deal(d.ori.g.data(p.aloitus:p.tiheys:p.lopetus,d.target.variables));

% Luo target.(local/global) eaX:n - skaalatun datan

p.skaalaustapa
%p.skaalaustapa = 1   % jossain asettuuu 2:ksi.

[d.target p] = ScaleDataV2(d.target,p);  % skaalaustapa 1 = normvarianssi 2=tietokanta
% alkuper�inen recursive std
%d.target
%d.target.kanta

d.target.sX = StdData(d.target.aX,p.window.sX,0);  % 100531 Tässä aina ollu 1 eikä 0, eli ennustusta 1 ruutu.. menee jälkeen
d.target.dX = DifferenceData(d.target.aX,p.window.sX,1);

d.target.fX = SmoothData(d.target.aX,p.window.X,1);

% 081107
% tehd��n Mat-2.3128 mukaisesti harkka8
d.target.kausivaihtelu = d.target.aX - d.target.fX;
d.target.dkausivaihtelu = SmoothData(d.target.kausivaihtelu,p.window.X,1);
d.target.kohina = d.target.kausivaihtelu - d.target.dkausivaihtelu;
% 081107

%30.5.2007
d.target.sdX = StdData(d.target.dX,p.window.sX,0);
d.target.sdefX = StdData(d.target.defX,p.window.sX,0);
% Mirror
d.target.mfX = MirrorDataV2(d.target.efX,p);
%d.target.mX = MirrorDataV2(d.target.eaX,p);  % t�m� hieman turha
%d.target.mdefX = MirrorDataV2(d.target.defX,p);
d.target.skew = SkewData(d.target.efX,p.window.X,0);

% *************************************************************************
% Siirretty SCALEDATAAN
% skaalattu recursive std
%d.target.seX = StdData(d.target.eaX,p.window.sX,1);
% filter�ity ja skaalattu X
%d.target.efX = SmoothData(d.target.eaX,p.window.X,2);
% kiinnostava: d(filter�ity ja skaalattu X)
%d.target.defX = DifferenceData(d.target.efX,p.window.dX,1);
% *************************************************************************


%d.target.aX = SmoothData(d.target.aX,d.target.X_valinpituus,1);
%d.target.dX = DifferenceData(d.target.X,p.window.dX,1);
%d.target.daX = DifferenceData(d.target.aX,p.window.dX,1);
%d.target.edX = DifferenceData(d.target.eX,p.window.dX,1);
% *************************************************************************

% 24.4.2007 - LEIKATAAN DATASTA LOPPU JA ALKU POIS!
d.target.aX = d.target.aX(p.window.X:end,:);
d.target.fX = d.target.fX(p.window.X:end,:);
d.target.eaX = d.target.eaX(p.window.X:end,:);
d.target.seX = d.target.seX(p.window.X:end,:);

d.target.sdX = d.target.sdX(p.window.X:end,:);
d.target.sdefX = d.target.sdefX(p.window.X:end,:);
d.target.efX = d.target.efX(p.window.X:end,:);
d.target.defX = d.target.defX(p.window.X:end,:);
d.target.dfX = d.target.dfX(p.window.X:end,:);
d.target.sX = d.target.sX(p.window.X:end,:);
d.target.dX = d.target.dX(p.window.X:end,:);

d.target.kausivaihtelu = d.target.kausivaihtelu(p.window.X:end,:);
d.target.dkausivaihtelu = d.target.dkausivaihtelu(p.window.X:end,:);
d.target.kohina = d.target.kohina(p.window.X:end,:);

d.target.deafX = d.target.deafX(p.window.X:end,:);

d.target.mfX = d.target.mfX(p.window.X:end,:);
%d.target.mX = d.target.mX(p.window.X:end,:);
%d.target.mdefX = d.target.mdefX(p.window.X:end,:);
d.target.skew = d.target.skew(p.window.X:end,:);


d.target.t = d.target.t(p.window.X:end);

% TESTI 30.11.2007: Saa helposti PCA analyysiin dX:n
%d.target.efX = d.target.defX;

%30.5.2007 - v�liaikainen - PCA kuitenkin dX datalle...
%d.target.efX = d.target.defX;
%d.target.mX = d.target.mdefX;


% 6.8.2007 - Smoothing index - kuvaa tasoituksen m��r�� DI
%figure(66)
%j = abs(d.target.efX(:,:)-d.target.eaX(:,:));
%bar(sum(j)/length(j));
%title('Smoothing index for each variable');

%figure(67);
%subplot(1,2,1);
%hold on;
%k=1;
%apu = SmoothData(d.target.aX(:,k),p.window.sX,1);
%plot(d.target.t,apu);
%plot(d.target.t,[apu + d.target.sX(:,k)*2 apu - d.target.sX(:,k)*2],'color',[0.6 0.6 0.6])

%plot(d.target.t, d.target.aX(:,1),'color',[0 0 0]);
%plot(d.target.t,apu,'color',[0.6 0.6 0.6]);
%hold off;

%xlabel('Time [s]');
%ylabel('[kg/s]');
%title('Coolant flow rate');
%subplot(1,2,2);
%hold on;
%plot(d.target.t, d.target.eaX(:,1),'color',[0.6 0.6 0.6]);
%plot(d.target.t, d.target.efX(:,1),'color',[0 0 0]);
%xlabel('Time [s]');
%title('Preprocessed coolant flow rate');
%hold off;



%subplot(2,2,4);
%plot(d.target.t, d.target.dfX(:,1));

%figure(68);
%subplot(1,2,1);
%plot(d.target.t, d.target.dX(:,1),'color',[0 0 0]);
%xlabel('Time [s]');
%ylabel('[d(kg/s)]');
%title('Gas flow');
%subplot(1,2,2);
%hold on;
%plot(d.target.t, d.target.eaX(:,1),'color',[0.6 0.6 0.6]);
%plot(d.target.t, d.target.sX(:,1),'color',[0 0 0]);
%xlabel('Time [s]');
%title('Preprocessed gas flow');
%hold off;

%figure(69);
%subplot(1,2,1);
%plot(d.target.t, d.target.dX(:,1),'color',[0 0 0]);
%xlabel('Time [s]');
%ylabel('[d(kg/s)]');
%title('Gas flow');
%subplot(1,2,2);
%hold on;
%plot(d.target.t, d.target.eaX(:,1),'color',[0.6 0.6 0.6]);
%plot(d.target.t, sum(d.target.sX(:,:)')','color',[0 0 0]);
%xlabel('Time [s]');
%title('Volatility');
%hold off;

%figure(71);
%subplot(1,2,1);
%plot(d.target.t, sum(d.target.seX(:,:)')','color',[0 0 0]);
%xlabel('Time [s]');
%title('Volatility');
%subplot(1,2,2);
%plot(d.target.t, sum(d.target.sdefX(:,:)')','color',[0 0 0]);
%xlabel('Time [s]');
%title('Volatility');

disp('Target data is ready');
%DataInformationV2
% *************************************************************************
% 28.8.2007
%englanniksi osa muuttujista:
if p.DI.english == 1
if  p.DI.englanniksiMuutettu == 0
    p.DI.englanniksiMuutettu = 1;
    
    disp('OLE TARKKANA, meneek� englanniksi oikeisiin kohtiin! Onko?');
    
    d.target.selitys{3} = 'sum of steam flows';
    
    disp(d.target.selitys{4});
    d.target.selitys{4} = 'steam flow 1';
    disp(d.target.selitys{4});
    d.target.selitys{5} = 'steam flow 2';
    d.target.selitys{6} = 'steam flow 3';
    d.target.selitys{7} = 'steam flow 4';
    d.target.selitys{8} = 'feed water flow';
    d.target.selitys{9} = 'HP-inlet pressure';
    %d.target.selitys{10} = 'piston position of high pressure feed water';
    d.target.selitys{10} = 'HPFW piston position';
    %d.target.selitys{11} = '413V503 piston position';
    %d.target.selitys{12} = 'piston position of high pressure feed water'; 
    d.target.selitys{12} = 'HPFW piston position';
    d.target.selitys{15} = 'idle power';
    
    if p.segmentation.start == 0
        d.target.selitys{44} = 'output voltage P5';
        d.target.selitys{37} = '400 kV voltage';
        d.target.selitys{32} = 'output of controller 535';
        
        disp(d.target.selitys{36});
        d.target.selitys{36} = 'output of level controller 537';
        disp(d.target.selitys{36});
    else
        d.target.selitys{32} = '400 kV voltage';
        d.target.selitys{30} = 'output of controller 535';
        
        disp(d.target.selitys{34});
        d.target.selitys{34} = 'output of level controller 537';
        disp(d.target.selitys{34});
    end
      
    %d.target.selitys{29} = 'APRM';
end
end

%\textbf{output of level controller 537}&98.70&101.19&85.95&102.98\\\hline
%\textbf{output voltage P5}&0.77&0.78&0.65&0.82\\\hline
%\textbf{output of controller 535}&5219.96&5265.84&4912.11&5345.05\\\hline
%\textbf{piston position of high pressure feed water}&81.98&87.60&55.83&89.54\\\hline
%\textbf{steam flow 4}&314.03&320.78&275.06&327.77\\\hline
%\textbf{400 kV voltage}&62.39&63.48&49.94&69.70\\\hline
%\textbf{idle power}&-0.47&-0.46&-1.09&-0.20\\\hline
%\textbf{steam flow 1}&1266.28&1289.61&1114.81&1303.59\\\hline


