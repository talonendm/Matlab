function varargout = DMT2(varargin)
% 
% **************************************************************
% Data management tool v2.0
% **************************************************************
% 
% talonen 
% 28.9.2007
% 25.6.2008: ajobugi huomattu: yksik�t muuttujissa vaihtuu ainakin
% 1.setiss� esik�sittelyn j�lkeen: ks. steam flow 4. PS. onko t�m� oikea
% DMT2-versio?
% 
% **************************************************************

% **************************************************************
% DMT2 M-file for DMT2.fig
%      DMT2, by itself, creates a new DMT2 or raises the existing
%      singleton*.
%
%      H = DMT2 returns the handle to a new DMT2 or the handle to
%      the existing singleton*.
%
%      DMT2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DMT2.M with the given input arguments.
%
%      DMT2('Property','Value',...) creates a new DMT2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DMT2_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DMT2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DMT2

% Last Modified by GUIDE v2.5 28-Sep-2007 14:50:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DMT2_OpeningFcn, ...
                   'gui_OutputFcn',  @DMT2_OutputFcn, ...
                   'gui_LayoutFcn',  @DMT2_LayoutFcn, ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before DMT2 is made visible.
function DMT2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DMT2 (see VARARGIN)

% Choose default command line output for DMT2
handles.output = hObject;


%set(handles.hpop2,'String', {'hsv','hot','kdk','gray'});


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DMT2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% **************************************************************
% Callbackit:
% **************************************************************
set(handles.infotext2, 'String','Real high frequency data sets: 1-11 & 12-39 (rXX), OL design based events 1-28 (dXX), OL simulator data 1-4 (sX)');

% **************************************************************
% --- Executes on button press in closePush.
function closePush_Callback(hObject, eventdata, handles)
% hObject    handle to closePush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.infotext, 'String'),'Please start data mining')
    joo = 1;
    % Ei edes a�loitettu dataminingia
else
p = handles.p;

% PDF report
if p.tool.report > 1
    d = handles.d;
    taulu = d.o.taulu;
    % Create PDF
    DMT2createPDF(p,d,taulu); %turhaan toi taulu tonne menee erikseen
    
    
    disp('END');
    %!latex reports.tex

end
end
% THE END
close all;
% **************************************************************




% **************************************************************
% **************************************************************
% Tekstieditointi kent�t
% **************************************************************
% **************************************************************

% **************************************************************
% Tiedoston nimi: editointi
% **************************************************************
function inEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% ***********
% Tiedoston nimi: Callback
% ***********
function inEdit_Callback(hObject, eventdata, handles)
% hObject    handle to inEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inEdit as text
%        str2double(get(hObject,'String')) returns contents of inEdit as a double
% **************************************************************
function ss_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
function ss_Callback(hObject, eventdata, handles)
% **************************************************************
function fsX_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
function fsX_Callback(hObject, eventdata, handles)
% **************************************************************
function fsdX_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
function fsdX_Callback(hObject, eventdata, handles)
% **************************************************************
function fssX_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
function fssX_Callback(hObject, eventdata, handles)
% **************************************************************


% **************************************************************
function variance_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
function variance_Callback(hObject, eventdata, handles)
% **************************************************************
function alpha_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
function alpha_Callback(hObject, eventdata, handles)
% **************************************************************
function max_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
function max_Callback(hObject, eventdata, handles)
% **************************************************************
function min_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
function min_Callback(hObject, eventdata, handles)
% **************************************************************
% **************************************************************
function reportfile_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
function reportfile_Callback(hObject, eventdata, handles)
% **************************************************************

% **************************************************************
function leakagespeed_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
function leakagespeed_Callback(hObject, eventdata, handles)
% **************************************************************
% **************************************************************
function leakageacceleration_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
function leakageacceleration_Callback(hObject, eventdata, handles)
% **************************************************************
% **************************************************************
function leakagestarttime_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
function leakagestarttime_Callback(hObject, eventdata, handles)
% **************************************************************







% **************************************************************
% Tiedosto (TEXTBOX) 4xx: Variable fast selection
% **************************************************************
function vfs_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
function vfs_Callback(hObject, eventdata, handles)
% **************************************************************





% **************************************************************
% Lataa datasetti: Callback
% **************************************************************
% --- Executes on button press in loadPush.
function loadPush_Callback(hObject, eventdata, handles)
% hObject    handle to loadPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dataset_file = get(handles.inEdit,'String');
if ~isempty(dataset_file)
  
% ******************************************
p.tyoskentelyhakemisto = '/share/dsnpp/doc/DMT2demo/DM';
%%
p.tyoskentelyhakemisto = '/m/fs/project1/dsnpp/data/olkiluoto-realdata/data_processed/Jaakko';
%%
p.tyoskentelyhakemisto = cd;
p.tyoskentelyhakemisto = strcat(p.tyoskentelyhakemisto,'/DM')

%%

cd(p.tyoskentelyhakemisto);
p = LoadParameters();

dataset_file = get(handles.inEdit,'String');

apu = str2num(dataset_file(2:end));
datasetfile = 1;
if apu>0
    if apu<38
        p.tilanne = apu;
        if apu>11
            datasetfile = 2
        end
    end
end

if strcmp(dataset_file(1),'r')
    p.datasettityyppi = 1;
elseif strcmp(dataset_file(1),'d')
    p.datasettityyppi = 2;
elseif strcmp(dataset_file(1),'s')
    p.datasettityyppi = 3;  % s = simulaattori
else
    p.datasettityyppi = 4;
end

disp('Loading ...');
% kentt� p�ivittyy vasta kun tulee funktiosta ulos
% set(handles.infotext, 'String','Loading data set ...');
[d.ori p] = LoadDataSet(p,datasetfile);
%p.kayttis = '/m/fs/home/talonen/Matlab';
cd(p.kayttis);
% **************************************************************
p.tool.phase = 1;    % 1=ladattu data

handles.p = p;
handles.d = d;
set(handles.createPush, 'Enable','on');
set(handles.applyPush, 'Enable','on');
set(handles.selectPush, 'Enable','off');
set(handles.selectmodelPush, 'Enable','off');
set(handles.modelPush, 'Enable','off');

set(handles.ss,'enable','on');
set(handles.fsX,'enable','on');
set(handles.fssX,'enable','on');
set(handles.fsdX,'enable','on');

set(handles.reportfile,'enable','on');

% **************************************************************
% KUVA PREVIEW
set(handles.Ima,'HandleVisibility','ON');
axes(handles.Ima);
plot(d.ori.g.time,d.ori.g.data(:,1));
title({d.ori.g.names{1},d.ori.g.exps{1}});
ylabel(d.ori.g.units{1});
xlabel('Time [s]')
%axis equal;
axis tight;
set(handles.Ima,'HandleVisibility','OFF');
% **************************************************************
% SLIDER2 property:
set(handles.intSlider2,'Value',1);
apu = (size(d.ori.g.data,2));
set(handles.intSlider2,'Max',apu);
set(handles.intSlider2,'SliderStep',[1/(apu-1) 1/(apu-1)]);
% **************************************************************
infoo = strcat(num2str(p.tilanne),': ',d.ori.tiedosto,', ',d.ori.aika, d.ori.kuvaus,', size of data set [n,m]: [',num2str(size(d.ori.g.data)),']')
set(handles.infotext2, 'String',infoo);
% **************************************************************
set(handles.infotext, 'String','Data set loaded, please start DM');
guidata(gcbo, handles);
end
% **************************************************************

% **************************************************************
% Leakage: Check button
% **************************************************************
% --- Executes on button press in transfCheck.
function transfCheck_Callback(hObject, eventdata, handles)
% hObject    handle to transfCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of transfCheck
% **************************************************************


% **************************************************************
% Duunaa datasetti: Callback
% **************************************************************
% --- Executes on button press in loadPush.
function createPush_Callback(hObject, eventdata, handles)
% hObject    handle to loadPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% **************************************************************
p = handles.p;
d = handles.d;

p.segmentation.start = str2num(get(handles.valText3,'String'));
p.segmentation.end = str2num(get(handles.valText4,'String'));
p.lopetus = p.datalength;
p.aloitus = round(p.lopetus*p.segmentation.start)+1;  %0.25 %kerro nollalla, niin alusta asti
p.lopetus = round(p.lopetus*p.segmentation.end);  % 0.6 hyv�
p.muuttujassamuutoksia = str2num(get(handles.valText,'String'));


apuHzP = findstr(d.ori.aika,'Hz');
if length(apuHzP) == 0
    p.naytettasek = 1;
    disp(strcat('Naytteita sekunnissa: (arvaus)',num2str(p.naytettasek)));
else
    %disp(apuHzP)
    disp(d.ori.aika(apuHzP-13:apuHzP-1))
    if p.datasettityyppi == 4
        p.naytettasek = str2num(d.ori.aika(apuHzP-13:apuHzP-1));
    else
        p.naytettasek = str2num(d.ori.aika(apuHzP-3:apuHzP-1));
    end
    disp(strcat('Naytteita sekunnissa:',num2str(p.naytettasek)));
end

p.window.dX = str2num(get(handles.fsdX,'String'));
p.window.X = str2num(get(handles.fsX,'String'));
p.window.sX = str2num(get(handles.fssX,'String'));
p.tiheys = str2num(get(handles.ss,'String'));

set(handles.ss,'enable','off');
set(handles.fsX,'enable','off');
set(handles.fssX,'enable','off');
set(handles.fsdX,'enable','off');

cd(p.tyoskentelyhakemisto);
[d p] = CreateTargetData(d,p);
cd(p.kayttis);
p.tool.phase = 2;    % 2 = target data luotu;



% POP UP
set(handles.hpop2,'String', d.target.selitys);
% **************************************************************

handles.d = d;
handles.p = p;
set(handles.selectPush, 'Enable','on');
set(handles.storevariablePush, 'Enable','on');
%set(handles.applyPush, 'Enable','on');
set(handles.createPush, 'Enable','off');

% **************************************************************
% SLIDER2 property:
set(handles.intSlider2,'Value',1);
apu = (size(d.target.aX,2));
set(handles.intSlider2,'Max',apu);
set(handles.intSlider2,'SliderStep',[1/(apu-1) 1/(apu-1)]);
% **************************************************************

set(handles.infotext, 'String','Target data set was created');
guidata(gcbo, handles);
% **************************************************************




% **************************************************************
% Valitse muuttujia datasetist�: Callback
% **************************************************************
% --- Executes on button press in loadPush.
function selectPush_Callback(hObject, eventdata, handles)
% hObject    handle to loadPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% **************************************************************
muuttujat = str2num(get(handles.vfs,'String'));
tarketti = get(handles.u0Text, 'String');
if length((muuttujat)) == 0 && strcmp(tarketti,'-')
    set(handles.infotext, 'String','Choose target variable or model variables');
else
% **************************************************************
p = handles.p;
d = handles.d;
% *********
p.DI.nopeamuuttujavalinta = 0;
p.WRLS.summahoyry = 0;
p.DI.nopeatmuuttujat = [];
% *********


if length((muuttujat))>0
    % MUUTTUJA FAST VALINTA
    %muuttujat = str2num(get(handles.vfs,'String'));
    
        p.DI.nopeamuuttujavalinta = 1;
        p.DI.nopeatmuuttujat = (muuttujat)';
        for i=1:length(muuttujat)
            if strcmp('CH311K035',d.target.nimet{muuttujat(i)})
                p.WRLS.summahoyry = 1;  % jos 3:n mukana, muuten 0
            end
        end
        p.PCA.target = muuttujat(1);
        p.PCA.target0 = muuttujat(1);
        set(handles.u0Text,'String',d.target.selitys{muuttujat(1)});
else
    
    uu = get(handles.u0Text,'String');
    p.PCA.target = find(strcmp(d.target.selitys,uu),1);
    p.PCA.target0 = find(strcmp(d.target.selitys,uu),1);
end



% *********

% *********
cd(p.tyoskentelyhakemisto);
d = SelectFeatures(d,p);
[d.features p] = FeatureSelection(d.features,p);
cd(p.kayttis);
% *********
p.tool.reportfile = get(handles.reportfile,'String');
if length(p.tool.reportfile) == 0 || strcmp(p.tool.reportfile,'report')
    p.tool.report = 0;
else
    p.tool.report = 1;
end
% *********
if p.tool.report > 0
    p.tool.figuuri = 301;
    p = DMT2report(p);
end

p.tool.phase = 3;
handles.p = p;   % my�s p:t� muutetaan
handles.d = d;
set(handles.selectmodelPush, 'Enable','on');
set(handles.reportfile,'enable','off');
%set(handles.selectPush, 'Enable','off');
set(handles.infotext, 'String','Variable selection was done');
guidata(gcbo, handles);
% **************************************************************
end
% **************************************************************






% **************************************************************
% Valitse muuttujia datasetist�: Callback
% **************************************************************
% --- Executes on button press in loadPush.
function selectmodelPush_Callback(hObject, eventdata, handles)
% hObject    handle to loadPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% **************************************************************
p = handles.p;
d = handles.d;
cd(p.tyoskentelyhakemisto);

% ********* Ylim. skeida **** kuva ulos
% 31.10.2007
% muualle muutokset muistiin (dd).. haetaan vain kuva
dd = SelectFeatures(d,p);

% ******************
% 29.1.2008 - reportti seteist� varten
if p.PCA.selection.type == 4
figure(321);  % h�t�sesti koodattu
for i=1:p.PCA.clusters
    subplot(3,2,i);
    plot(d.target.t, d.target.aX(:,dd.features.variables(i)));
    title(dd.features.selitys{i});
    ylabel(dd.features.yksikko{i});
end
end
if p.PCA.selection.type == 0
figure(321);  % h�t�sesti koodattu
for i=1:length(dd.features.variables)
    subplot(4,3,i);
    
    disp('dd.features.variables(i) (selectiontype =0)');
    disp(dd.features.variables(i));
    plot(d.target.t, d.target.aX(:,dd.features.variables(i)));
    title(strcat('virhe',dd.features.selitys{i}));
    ylabel(dd.features.yksikko{i});
end
end

%cd(p.tyoskentelyhakemisto);
% ******************
[dd.features pp] = FeatureSelection(dd.features,p);

disp('dd.features.nimet');
dd.features.nimet

cd(p.kayttis);
if p.tool.report > 0
    if p.PCA.selection.type == 4
        p.tool.figuuri = 321;
        p = DMT2report(p);
    end
end
if p.tool.report > 0
    p.tool.figuuri = 301;
    p = DMT2report(p);
end

cd(p.tyoskentelyhakemisto);




% ****************** tunnarit jota tPDF taulukkoon 4.12.2007
d.o.taulu.taulukko.rowLabels = dd.features.selitys;
d.o.taulu.taulukko.columnLabels = {'Tunnus','Yksikk�','Min','Max','Mean','Smin','Smax'};

for i=1:length(d.o.taulu.taulukko.rowLabels)
    %for j=1:length(d.o.taulu.taulukko.columnLabels)
        
        % AIVAN TYPER� SILMUKKA
        %if j==1 
            data{i,1} = (dd.features.nimet{i});
            data{i,1} = strrep(data{i,1},'_','\_');
            %data2{i,1} = d.target.nimet{dd.features.variables(i)};
       
            data{i,2} = (dd.features.yksikko{i});
            data{i,2} = strrep(data{i,2},'%','\%');
            
        %elseif j==2
            %data{i,2} = num2str(min(dd.features.data(:,i)));
            data{i,3} = num2str(min(d.target.aX(:,dd.features.variables(i))));
       
        %elseif j==3
            %data{i,3} = num2str(max(dd.features.data(:,i)));
            data{i,4} = num2str(max(d.target.aX(:,dd.features.variables(i))));
       
        %else
            %data{i,4} = num2str(mean(dd.features.data(:,i)));
        %end
        
            data{i,5} = num2str(mean(d.target.aX(:,dd.features.variables(i))));
        
            data{i,6} = num2str(min(dd.features.data(:,i)));
            data{i,7} = num2str(max(dd.features.data(:,i)));
            
            
            
    %end
end

d.o.taulu.taulukko.taulu = data;
%d.o.taulu.taulukko.taulu2 = data2;
% ******************

% ********* Ylim. skeida loppuu **** kuva ulos

[d p] = ModelSelection(d,p);
cd(p.kayttis);
handles.p = p;
handles.d = d;
set(handles.modelPush, 'Enable','on');
set(handles.selectmodelPush, 'Enable','off');
set(handles.selectPush, 'Enable','off');
set(handles.infotext, 'String','Model features was selected');
guidata(gcbo, handles);
% **************************************************************


% **************************************************************
% Valitse muuttujia datasetist�: Callback
% **************************************************************
% --- Executes on button press in loadPush.
function modelPush_Callback(hObject, eventdata, handles)
% hObject    handle to loadPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% **************************************************************
p = handles.p;
d = handles.d;

p.WRLS.alpha = str2num(get(handles.alpha,'String'));
p.WRLS.normalvariance = str2num(get(handles.variance,'String'));
p.WRLS.lambda.min = str2num(get(handles.min,'String'));
p.WRLS.lambda.max = str2num(get(handles.max,'String')); 
p.WRLS.lambda.start = p.WRLS.lambda.min;

if get(handles.transfCheck,'Value')
    p.WRLS.vuotosimulointi = 2;
    p.WRLS.leak.size2 = str2num(get(handles.leakageacceleration,'String')); 
    p.WRLS.leak.size = str2num(get(handles.leakagespeed,'String')); 
    p.WRLS.leak.starttime = str2num(get(handles.leakagestarttime,'String'));
else
    p.WRLS.vuotosimulointi = 0;
end

%disp(p.WRLS.leak.size2)
%disp(p.WRLS.leak.size)

%p.tool.reportfile = get(handles.reportfile,'String');
%if length(p.tool.reportfile) == 0 || strcmp(p.tool.reportfile,'report')
%    p.tool.report = 0;
%end


cd(p.tyoskentelyhakemisto);
d.model = WRLSData(d.features,p);
PlotModel(d,p);  
cd(p.kayttis);

if p.tool.report > 0
    p.tool.figuuri = 40;
    p = DMT2report(p);
end
if p.tool.report > 0
    p.tool.figuuri = 343;
    p = DMT2report(p);
end
if p.tool.report > 0
    p.tool.figuuri = 701;
    p = DMT2report(p);
end
if p.tool.report > 0
    p.tool.figuuri = 503;
    p = DMT2report(p);
end

handles.p = p;
handles.d = d;
set(handles.infotext, 'String','Model was created');
guidata(gcbo, handles);
% **************************************************************

% **************************************************************
function applyPush_Callback(hObject, eventdata, handles)
p = handles.p;
d = handles.d;

t = get(handles.intSlider2,'value');
t2 = floor(t);

p.tool.val = get(handles.hpop2,'Value');
val = p.tool.val;

figure(124);
if p.tool.phase < 2
  
    plot(d.ori.g.time,d.ori.g.data(:,t2));
    title({d.ori.g.names{t2},d.ori.g.exps{t2}});
    ylabel(d.ori.g.units{t2});
    xlabel('Time [s]');
    %axis equal;
    axis tight;
    koodi = strcat('v',d.ori.g.names{t2}(3:5),'.png');
    onkuva = 0;
    for i=1:length(p.fig.virtauskuvakoodit)
        if str2num(koodi(2:4)) == p.fig.virtauskuvakoodit(i)
            onkuva =1;
        end
    end
    if onkuva == 0
        koodi = 'Paavirtaus.png';
    end
    figure(125);
    cd(p.virtauskaaviot);
    I = imread(koodi);
    image(I);
    axis off;
    title({d.ori.g.names{t2},d.ori.g.exps{t2}});
else

    plot(d.target.t,d.target.aX(:,t2));
    title({d.target.nimet{t2},d.target.selitys{t2}});
    ylabel(d.target.yksikko{t2});
    xlabel('Time [s]');
    %axis equal;
    axis tight;
    koodi = strcat('v',d.target.nimet{t2}(3:5),'.png');
    onkuva = 0;
    for i=1:length(p.fig.virtauskuvakoodit)
        if str2num(koodi(2:4)) == p.fig.virtauskuvakoodit(i)
            onkuva =1;
        end
    end
    if onkuva == 0
        koodi = 'Paavirtaus.png';
    end
    figure(125);
    cd(p.virtauskaaviot);
    I = imread(koodi);
    image(I);
    axis off;
    title({d.target.nimet{t2},d.target.selitys{t2}});
    % 081107 - kausivaihtelun tarkastelu
    figure(126);
    subplot(2,3,1);
    plot(d.target.t,d.target.aX(:,t2));
    title('aX');
    subplot(2,3,2);
    plot(d.target.t,d.target.fX(:,t2));
    title('fX');
    subplot(2,3,3);
    plot(d.target.t,d.target.kausivaihtelu(:,t2));
    title('kausivaihtelu');
    subplot(2,3,4);
    plot(d.target.t,d.target.dkausivaihtelu(:,t2));
    title('dkausivaihtelu');
    subplot(2,3,5);
    plot(d.target.t,d.target.kohina(:,t2));
    title('kohina');
    subplot(2,3,6);
    plot(d.target.t,d.target.dX(:,t2));
    title('dX');
    
    figure(127);
    plot(d.target.kohina(:,t2),d.target.kohina(:,val),'x')
        
    X = (d.target.kohina(:,t2));
    Y = (d.target.kohina(:,val));
    B = X \ Y;
    B = [ones(length(X),1) X] \ Y;

    %for j=2:length(vp2{i});
    %    line([d.k(vp2{i}(j)) d.k(vp2{i}(j-1))],[d.h(vp2{i}(j)) d.h(vp2{i}(j-1))],'color',[i/length(vp2),0,0]);
    %end
    
    line([min(X) max(X)],[B(1)+B(2)*min(X) B(1)+B(2)*max(X)],'color',[0.3,0.1,0.3]);
   
    correlaatio = corr(d.target.kohina(:,t2),d.target.kohina(:,val));
    title(strcat('Kohinakorrelaatio',d.target.nimet{t2},' vs. ',d.target.nimet{val},' corr: ',num2str(correlaatio)));
    
end 



cd(p.kayttis);


% **************************************************************


% **************************************************************
function storevariablePush_Callback(hObject, eventdata, handles)
p = handles.p;
d = handles.d;

%save d_p.mat d p -mat
save d_p100531.mat d p -mat
dtarget = d.target.efX;

%dtargetefX = d.target.efX;
save dtarget.mat dtarget -mat
%save dtargetefx dtargetefX -csv

disp('variables d & p stored');
% **************************************************************





% **************************************************************



% **************************************************************
% Slider
% **************************************************************
% --- Executes during object creation, after setting all properties.
function intSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to intSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% --- Executes on slider movement.
% ***********
% Slider: Callback
% ***********
function intSlider_Callback(hObject, eventdata, handles)
% hObject    handle to intSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
t = get(handles.intSlider,'value');
t2 = floor(t);
% VAIN  KOKONAISLUKUJA
set(handles.intSlider,'value',t2);
set(handles.valText,'String',num2str(t2));
% **************************************************************


% **************************************************************
% Slider2
% **************************************************************
function intSlider2_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% ***********
% Slider: Callback
% ***********
function intSlider2_Callback(hObject, eventdata, handles)
t = get(handles.intSlider2,'value');
t2 = floor(t);

set(handles.intSlider2,'value',t2);
set(handles.valText2,'String',num2str(t2));
% **************************************************************
% KUVA PREVIEW
d = handles.d;
p = handles.p;
if p.tool.phase < 2
    set(handles.Ima,'HandleVisibility','ON');
    axes(handles.Ima);
    plot(d.ori.g.time,d.ori.g.data(:,t2));
    title({d.ori.g.names{t2},d.ori.g.exps{t2}});
    ylabel(d.ori.g.units{t2});
    xlabel('Time [s]');
    %axis equal;
    axis tight;
    set(handles.Ima,'HandleVisibility','OFF');
else
    set(handles.Ima,'HandleVisibility','ON');
    axes(handles.Ima);
    plot(d.target.t,d.target.aX(:,t2));
    title({d.target.nimet{t2},d.target.selitys{t2}});
    ylabel(d.target.yksikko{t2});
    xlabel('Time [s]');
    %axis equal;
    axis tight;
    set(handles.Ima,'HandleVisibility','OFF');
end 
% **************************************************************

% **************************************************************


% **************************************************************
% Slider3 - segmentation start
% **************************************************************
function intSlider3_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% ***********
function intSlider3_Callback(hObject, eventdata, handles)
t = get(handles.intSlider3,'value');
%t2 = floor(t);
t2=t;
set(handles.intSlider3,'value',t2);
set(handles.valText3,'String',num2str(t2));
% **************************************************************

% **************************************************************
% Slider4 - segmentation start
% **************************************************************
function intSlider4_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% ***********
function intSlider4_Callback(hObject, eventdata, handles)
t = get(handles.intSlider4,'value');
t2 = floor(t);
t2=t;
set(handles.intSlider4,'value',t2);
set(handles.valText4,'String',num2str(t2));
% **************************************************************







% **************************************************************
% Radio buttons
% **************************************************************
function u0_Callback(hObject, eventdata, handles)
set(handles.u1,'Value',0);
set(handles.u0,'Value',1);
set(handles.u2,'Value',0);
function u1_Callback(hObject, eventdata, handles)
set(handles.u1,'Value',1);
set(handles.u0,'Value',0);
set(handles.u2,'Value',0);
function u2_Callback(hObject, eventdata, handles)
set(handles.u2,'Value',1);
set(handles.u1,'Value',0);
set(handles.u0,'Value',0);
% **************************************************************



% **************************************************************
% POP UP 1
% **************************************************************
function hpop_Callback(hObject, eventdata, handles)
% **************************************************************
% KUVA PREVIEW
d = handles.d;
p = handles.p;
p.tool.val = get(handles.hpop2,'Value');
p.tool.tyyli = get(handles.hpop,'Value');


if get(handles.u2,'Value') == 1
     kuva = 2;
elseif get(handles.u1,'Value') == 1
     kuva = 1;
else
    kuva = 0;
end

if kuva > 0
if p.tool.phase > 1
    t = get(handles.intSlider2,'value');
    t2 = floor(t);
    p.ttt2 = t2;
    
    % SOM:n muuttujat textboxista 081126 ei kai toimi
    muuttujat = str2num(get(handles.vfs,'String'))
    
    if length(muuttujat)>0
        p.alarm.som.variables = muuttujat;
    end
    
    if kuva == 1
        set(handles.Ima2,'HandleVisibility','ON');   
        axes(handles.Ima2);   
        DMT2fig2(d,p);
        axes(handles.Ima2);
    else
        set(handles.Ima3,'HandleVisibility','ON');   
        axes(handles.Ima3);   
        DMT2fig2(d,p);
        axes(handles.Ima3);   
    end
    %axis equal;
    axis tight;
    if p.tool.tyyli == 3
       axis([d.target.t(1),d.target.t(end),0,1]);
    end
    if kuva == 1
        set(handles.Ima2,'HandleVisibility','OFF');
    else
        set(handles.Ima3,'HandleVisibility','OFF');
    end
end
end


% **************************************************************
% POP UP 2: Kuvaaja tyyli (eax, fx ym.)
% **************************************************************

function hpop2_Callback(hObject, eventdata, handles)
% **************************************************************
% KUVA PREVIEW
d = handles.d;
p = handles.p;
p.tool.val = get(handles.hpop2,'Value');
p.tool.tyyli = get(handles.hpop,'Value');


if get(handles.u2,'Value') == 1
    set(handles.u2Text,'String',d.target.selitys{get(handles.hpop2,'Value')});
    kuva = 2;
    
    
elseif get(handles.u1,'Value') == 1
    set(handles.u1Text,'String',d.target.selitys{get(handles.hpop2,'Value')});
    kuva = 1;
else
    set(handles.u0Text,'String',d.target.selitys{get(handles.hpop2,'Value')});
    kuva = 0;
end

if kuva > 0
if p.tool.phase > 1
    t = get(handles.intSlider2,'value');
    t2 = floor(t);
    p.ttt2 = t2;
    if kuva == 1
        set(handles.Ima2,'HandleVisibility','ON');   
        axes(handles.Ima2);   
        DMT2fig2(d,p);
        axes(handles.Ima2);
        
    else
        set(handles.Ima3,'HandleVisibility','ON');   
        axes(handles.Ima3);   
        DMT2fig2(d,p);
        axes(handles.Ima3);   
        
    end
    %axis equal;
    axis tight;
    if p.tool.tyyli == 3
       axis([d.target.t(1),d.target.t(end),0,1]);
    end
    if kuva == 1
        set(handles.Ima2,'HandleVisibility','OFF');
    else
        set(handles.Ima3,'HandleVisibility','OFF');
    end
end

end


% **************************************************************
% POP UP 3
% **************************************************************
function hpop3_Callback(hObject, eventdata, handles)
kuvaaja = get(handles.hpop3,'Value');
d = handles.d;
p = handles.p;

muuttujia = length(d.target.variables);

if p.tool.phase > 1
    
    set(handles.Ima4,'HandleVisibility','ON');
    axes(handles.Ima4); 
    
    if kuvaaja == 1
        % volatility
        if length(p.PCA.variables)>0
            plot(d.target.t, sum(d.target.seX(:,p.PCA.variables)')'/length(p.PCA.variables),'color',[0 0 0]);
        else
            plot(d.target.t, sum(d.target.seX(:,:)')'/muuttujia,'color',[0 0 0]);
        end
        xlabel('Time [s]');
    elseif kuvaaja == 2
        % smoothing indices
        if length(p.PCA.variables)>0
            j = abs(d.target.efX(:,p.PCA.variables)-d.target.eaX(:,p.PCA.variables));
            bar(sum(j)/length(j));
        else
            j = abs(d.target.efX(:,:)-d.target.eaX(:,:));
            bar(sum(j)/length(j));
        end
        xlabel('Variable number');
        
    elseif kuvaaja ==3
        % margin = max-min: bar plot
        if length(p.PCA.variables)>0
            j = max(d.target.efX(:,p.PCA.variables)) - min(d.target.efX(:,p.PCA.variables));
            bar(j);
        else
            j = max(d.target.efX(:,:)) - min(d.target.efX(:,:));
            bar(j);
        end
        xlabel('Variable number');
        ylabel('0..1');
    elseif kuvaaja == 4
        if length(p.PCA.variables)>0
            boxplot(d.target.eaX(:,p.PCA.variables));
            
        else
            boxplot(d.target.eaX);
        end
        title('Boxplot');
    elseif kuvaaja == 5
        if length(get(handles.u1Text,'String'))>1
            if length(get(handles.u2Text,'String'))>1  
                
                uu = get(handles.u1Text,'String');
                yy = get(handles.u2Text,'String');
                uuu = find(strcmp(d.target.selitys,uu),1);
                yyy = find(strcmp(d.target.selitys,yy),1);
                
                kohde = [uuu yyy];
                viive = DMT2delayestimation(d,p,kohde);
                % 1 ja 4 autokorrelaatiot
                plot(viive.tau, viive.risti(:,2),'color',[0 0 0]);
                xlabel('\tau')
                %title('Cross-correlation: from u to y');
                
            end
        end
    elseif kuvaaja == 6
        %Segmentation 12.8.2008
        %setti = 8;
        klustereita = 2;
        cd(p.tyoskentelyhakemisto);

        %figure(1000);

        %k=0;
        %[a b] = k_means(vika{setti}.g.data,klustereita);
        %[a b] = k_means(d.target.aX,klustereita);
        %for i=1:5
        %    for j=1:4
        %        k=k+1;
        %        subplot(5,4,k);
        %        if k==20
        %            plot(a,'x')
        %        else
        %            %plot(vika{setti}.g.data(:,k));
        %            plot(d.target.aX(:,k));
        %        end
        %    end
        %end

       

        

        k=0;
        %[a b] = k_means(vika{setti}.g.data,klustereita);
        [a b] = k_means(abs(d.target.defX),klustereita);
        [a2 b2] = k_means(abs(d.target.efX),klustereita);
        
        % pikkuikkuna
        plot(d.target.t, a,'x','color',[0 0 0]);
        xlabel('Segments')
        
        % segmentointi - oma idea 0808 & 081127
        figure(1501);
        subplot(2,1,1);
        plot(d.target.t, a,'x','color',[0 0 0]);
        xlabel('Segments')
        subplot(2,1,2);
        plot(d.target.t, a2,'x','color',[0 0 0]);
        xlabel('Segments')
        
        
        figure(1502);
        subplot(2,1,1);
        plot(b(1,:),b(2,:),'.');
        for i=1:length(b)
            text(b(1,i),b(2,i),d.target.selitys{i});
        end
        title('cluster centroid locations - using defX');
        xlabel('b1')
        ylabel('b2')
        
        subplot(2,1,2);
        plot(b2(1,:),b2(2,:),'.');
        for i=1:length(b2)
            text(b2(1,i),b2(2,i),d.target.selitys{i});
        end
            %size(b)
        title('cluster centroid locations - using efX');
        xlabel('b1')
        ylabel('b2')
        
        disp('Ehdotus segmentoinniksi:');
        segmenttipituus = 0;
        segmentteja = 1;
        virhesegmentteja = 0;
        
        viimeisinvirhe = 0;
        
        segmenttia{1} = 1;
        segmenttil{1} = 0;
        
        for i=2:length(a)
            if a(i) == a(i-1)
                segmenttipituus = segmenttipituus + 1;
            else
                if segmenttipituus>p.window.X
                    segmenttil{segmentteja} = i - p.window.X;
                    segmentteja = segmentteja + 1;
                    segmenttia{segmentteja} = i - p.window.X + 1;
                    disp(segmentteja);
                    disp(d.target.t(i));
                else
                    virhesegmentteja = virhesegmentteja + 1;
                    viimeisinvirhe = i;
                    
                end
                segmenttipituus = 0;
                
            end
        end
        if virhesegmentteja>0
            if viimeisinvirhe>segmenttia{segmentteja-1}
                segmenttil{segmentteja} = viimeisinvirhe - p.window.X;
            else
                segmenttil{segmentteja} = i;
            end
        else
            segmenttil{segmentteja} = i;
        end
        
        disp('virhesegmentteja');
        disp(virhesegmentteja);
        
        if virhesegmentteja>0
        if viimeisinvirhe + p.window.X +1<segmenttil{segmentteja}
            segmentteja = segmentteja + 1;
            segmenttia{segmentteja} = viimeisinvirhe - p.window.X + 1;
            segmenttil{segmentteja} = viimeisinvirhe + p.window.X;
            
            segmentteja = segmentteja + 1;
            segmenttia{segmentteja} = viimeisinvirhe + p.window.X + 1;
            segmenttil{segmentteja} = length(a);
            
        end
        end
        
        for i=1:segmentteja
            disp(strcat('segment:',num2str(i)));
            disp(d.target.t(segmenttia{i}));
            disp(d.target.t(segmenttil{i}));
        end
        
        
        
        %figure(1002);
        %for i=1:5
        %    for j=1:4
        %        k=k+1;
         %       subplot(5,4,k);
        %        if k==20
         %           plot(a,'x')
         %       else
                    %plot(vika{setti}.g.data(:,k));
        %            plot(d.target.aX(:,k));
        %        end
        %    end
       % end
        
        cd(p.kayttis);
        
    end
    
    axis tight;
    if kuvaaja>1
        axis tight;
    end
    set(handles.Ima4,'HandleVisibility','OFF');
end
% **************************************************************















% **************************************************************
% T�h�n p��ttyi k�ytt�liittym�n koodi
% **************************************************************
% Alla oleva koodi m��rittelee K�ytt�liittym�n ulkoasun
% **************************************************************
% Exported FIG
% **************************************************************
% --- Outputs from this function are returned to the command line.
function varargout = DMT2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Creates and returns a handle to the GUI figure. 
function h1 = DMT2_LayoutFcn(policy)
% policy - create a new figure or use a singleton. 'new' or 'reuse'.

persistent hsingleton;
if strcmpi(policy, 'reuse') & ishandle(hsingleton)
    h1 = hsingleton;
    return;
end

appdata = [];
appdata.GUIDEOptions = struct(...
    'active_h', [], ...
    'taginfo', struct(...
    'figure', 2), ...
    'override', 0, ...
    'release', 13, ...
    'resize', 'none', ...
    'accessibility', 'callback', ...
    'mfile', 1, ...
    'callbacks', 1, ...
    'singleton', 1, ...
    'syscolorfig', 1, ...
    'blocking', 0, ...
    'lastSavedFile', '/m/fs/home/talonen/Material_Matlab/gui/DMT2.m');
appdata.lastValidTag = 'figure1';
appdata.GUIDELayoutEditor = [];

h1 = figure(...
'Units','characters',...
'Color',[0.6 0.6 0.6],...
'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
'IntegerHandle','off',...
'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
'MenuBar','none',...
'Name','DMT2',...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'Position',[10 30 120 62],...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1',...
'UserData',[],...
'Behavior',get(0,'defaultfigureBehavior'),...
'Visible','on',...
'CreateFcn', {@local_CreateFcn, '', appdata} );

% **************************************************************
% Buttons here
% **************************************************************

% **************************************************************
% Close
% **************************************************************
appdata = [];
appdata.lastValidTag = 'closePush';

h9 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','DMT2(''closePush_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[5 5 25 2],...
'String','Close',...
'Tag','closePush',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

% **************************************************************
% Tiedosto NUM
% **************************************************************
appdata = [];
appdata.lastValidTag = 'inEdit';

h3 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','DMT2(''inEdit_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[5 55 5 2],...
'String','r8',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'DMT2(''inEdit_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','inEdit',...
'Behavior',get(0,'defaultuicontrolBehavior'));

% **************************************************************
% Tiedosto avaa
% **************************************************************
appdata = [];
appdata.lastValidTag = 'loadPush';

h4 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','DMT2(''loadPush_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[10 55 20 2],...
'String','Load data set',...
'Tag','loadPush',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************

% **************************************************************
% Leakage or not?!
% **************************************************************
appdata = [];
appdata.lastValidTag = 'transfCheck';

h5 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','DMT2(''transfCheck_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[5 9.8 25 2],...
'String','Artificial leakage',...
'Style','checkbox',...
'Tag','transfCheck',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );







% **************************************************************
% Create target data
% **************************************************************
appdata = [];
appdata.lastValidTag = 'createPush';

h101 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','DMT2(''createPush_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[5 52 25 2],...
'String','Create target data',...
'Tag','createPush',...
'Enable','off',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************


% **************************************************************
% PCA to target data
% **************************************************************
appdata = [];
appdata.lastValidTag = 'selectPush';

h102 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','DMT2(''selectPush_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[5 49 25 2],...
'String','Select variables',...
'Tag','selectPush',...
'Enable','off',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************

% **************************************************************
% Model selection
% **************************************************************
%appdata = [];
%appdata.lastValidTag = 'selectmodelPush';
h103 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','DMT2(''selectmodelPush_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[5 46 25 2],...
'String','Select model features',...
'Tag','selectmodelPush',...
'Enable','off',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************


% **************************************************************
% Model selection
% **************************************************************
%appdata = [];
%appdata.lastValidTag = 'modelPush';
h104 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','DMT2(''modelPush_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[5 43 25 2],...
'String','Modeling',...
'Tag','modelPush',...
'Enable','off',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************


% **************************************************************
% APPLY preprocess
% **************************************************************
%appdata = [];
%appdata.lastValidTag = 'applyPush';
h105 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','DMT2(''applyPush_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[5 25 15 2],...
'String','Draw #',...
'Tag','applyPush',...
'Enable','off',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************


% **************************************************************
% STOREVARIABLES 
% **************************************************************
h106 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','DMT2(''storevariablePush_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[20 25 10 2],...
'String','Save d&p',...
'Tag','storevariablePush',...
'Enable','off',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************


% **************************************************************
% Slider (huom. pos Y-akseli sijainti alhaalta p�in..)
% **************************************************************
%appdata = [];
%appdata.lastValidTag = 'intSlider';
h6 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback','DMT2(''intSlider_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Max',100,...
'Min',1,...
'Position',[32 54.3 20 0.7],...
'String',{  '' },...
'Style','slider',...
'Value',100,...
'SliderStep',[1/99 1/99],...
'CreateFcn', {@local_CreateFcn, 'DMT2(''intSlider_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','intSlider',...
'Behavior',get(0,'defaultuicontrolBehavior'));
% **************************************************************
% **************************************************************
% Slider text
% **************************************************************
%appdata = [];
%appdata.lastValidTag = 'valText';
h7 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Position',[32 53 20 1.2],...
'String','100',...
'Style','text',...
'Tag','valText',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************


% **************************************************************
% Slider2 (huom. pos Y-akseli sijainti alhaalta p�in..)
% **************************************************************
h201 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback','DMT2(''intSlider2_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Max',2,...
'Min',1,...
'Position',[32 56.3 20 0.7],...
'String',{  '' },...
'Style','slider',...
'Value',1,...
'SliderStep',[0.1 0.1],...
'CreateFcn', {@local_CreateFcn, 'DMT2(''intSlider2_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','intSlider2',...
'Behavior',get(0,'defaultuicontrolBehavior'));
% **************************************************************
% **************************************************************
% Slider2 text
% **************************************************************
h202 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Position',[32 55 20 1.2],...
'String','1',...
'Style','text',...
'Tag','valText2',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************

% **************************************************************
% Slider - SegmentationStart
% **************************************************************
h203 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback','DMT2(''intSlider3_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Max',1,...
'Min',0,...
'Position',[32 52.3 10 0.7],...
'String',{  '' },...
'Style','slider',...
'Value',0,...
'SliderStep',[0.01 0.02],...
'CreateFcn', {@local_CreateFcn, 'DMT2(''intSlider3_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','intSlider3',...
'Behavior',get(0,'defaultuicontrolBehavior'));
% **************************************************************
% **************************************************************
% Slider text - SegmentationStart
% **************************************************************
h204 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Position',[32 51 10 1.2],...
'String','0.00',...
'Style','text',...
'Tag','valText3',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************

% **************************************************************
% Slider - SegmentationStop
% **************************************************************
h205 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback','DMT2(''intSlider4_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Max',1,...
'Min',0,...
'Position',[42 52.3 10 0.7],...
'String',{  '' },...
'Style','slider',...
'Value',0.92,...
'SliderStep',[0.02 0.04],...
'CreateFcn', {@local_CreateFcn, 'DMT2(''intSlider4_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','intSlider4',...
'Behavior',get(0,'defaultuicontrolBehavior'));
% **************************************************************
% **************************************************************
% Slider text - SegmentationStop
% **************************************************************
h206 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Position',[42 51 10 1.2],...
'String','0.92',...
'Style','text',...
'Tag','valText4',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************







% **************************************************************
% Figure (3XX): Original variables
% **************************************************************
h301 = axes(...
'Parent',h1,...
'Units','characters',...
'Position',[62 45.5 55 10.5],...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'XColor',get(0,'defaultaxesXColor'),...
'YColor',get(0,'defaultaxesYColor'),...
'ZColor',get(0,'defaultaxesZColor'),...
'Tag','Ima',...
'Behavior',get(0,'defaultaxesBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
%h15 = axes(...
%'Parent',h1,...
%'Units','characters',...
%'Position',[60 12 30 15],...
%'CameraPosition',[0.5 0.5 9.16025403784439],...
%'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
%'Color',get(0,'defaultaxesColor'),...
%'ColorOrder',get(0,'defaultaxesColorOrder'),...
%'LooseInset',[18 4 13 3],...
%'XColor',get(0,'defaultaxesXColor'),...
%'YColor',get(0,'defaultaxesYColor'),...
%'ZColor',get(0,'defaultaxesZColor'),...
%'Tag','Ima',...
%'Behavior',get(0,'defaultaxesBehavior'),...
%'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************

% **************************************************************
% Figure: Time series
% **************************************************************
h302 = axes(...
'Parent',h1,...
'Units','characters',...
'Position',[62 31 55 10.5],...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'XColor',get(0,'defaultaxesXColor'),...
'YColor',get(0,'defaultaxesYColor'),...
'ZColor',get(0,'defaultaxesZColor'),...
'Tag','Ima2',...
'Behavior',get(0,'defaultaxesBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************

% **************************************************************
% Figure: Time series
% **************************************************************
h303 = axes(...
'Parent',h1,...
'Units','characters',...
'Position',[62 16.5 55 10.5],...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'XColor',get(0,'defaultaxesXColor'),...
'YColor',get(0,'defaultaxesYColor'),...
'ZColor',get(0,'defaultaxesZColor'),...
'Tag','Ima3',...
'Behavior',get(0,'defaultaxesBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************

% **************************************************************
% Figure: Time series
% **************************************************************
h304 = axes(...
'Parent',h1,...
'Units','characters',...
'Position',[62 2.5 55 9],...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'XColor',get(0,'defaultaxesXColor'),...
'YColor',get(0,'defaultaxesYColor'),...
'ZColor',get(0,'defaultaxesZColor'),...
'Tag','Ima4',...
'Behavior',get(0,'defaultaxesBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************


% **************************************************************
% **************************************************************
% TEXTBOXES
% **************************************************************
% **************************************************************

% **************************************************************
% Tiedosto (TEXTBOX) (4xx): Sample size: tiheys
% **************************************************************
h400 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','DMT2(''ss_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[5 21 5 2],...
'String','1',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'DMT2(''ss_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','ss',...
'Behavior',get(0,'defaultuicontrolBehavior'));
% **************************************************************

% **************************************************************
% Tiedosto (TEXTBOX) 4xx: Frame size: X
% **************************************************************
h401 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','DMT2(''fsX_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[5 18 5 2],...
'String','5',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'DMT2(''fsX_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','fsX',...
'Behavior',get(0,'defaultuicontrolBehavior'));
% **************************************************************

% **************************************************************
% Tiedosto (TEXTBOX) 4xx: Frame size: dX
% **************************************************************
h402 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','DMT2(''fsdX_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[5 15 5 2],...
'String','1',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'DMT2(''fsdX_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','fsdX',...
'Behavior',get(0,'defaultuicontrolBehavior'));
% **************************************************************

% **************************************************************
% Tiedosto (TEXTBOX) 4xx: Frame size: sX
% **************************************************************
h403 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','DMT2(''fssX_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[5 12 5 2],...
'String','20',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'DMT2(''fssX_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','fssX',...
'Behavior',get(0,'defaultuicontrolBehavior'));
% **************************************************************

% **************************************************************
% Tiedosto (TEXTBOX) 4xx: Variable fast selection
% **************************************************************
h404 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','DMT2(''vfs_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[32 49 20 2],...
'String','',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'DMT2(''vfs_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','vfs',...
'Behavior',get(0,'defaultuicontrolBehavior'));
% **************************************************************
% [7 3 9 10] esim..



% **************************************************************
% Tiedosto (TEXTBOX) (4xx): normalvariance
% **************************************************************
h405 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','DMT2(''variance_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[28 21 12 2],...
'String','0.000008',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'DMT2(''variance_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','variance',...
'Behavior',get(0,'defaultuicontrolBehavior'));
% **************************************************************

% **************************************************************
% Tiedosto (TEXTBOX) 4xx: alpha
% **************************************************************
h406 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','DMT2(''alpha_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[28 18 12 2],...
'String','1000',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'DMT2(''alpha_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','alpha',...
'Behavior',get(0,'defaultuicontrolBehavior'));
% **************************************************************

% **************************************************************
% Tiedosto (TEXTBOX) 4xx: max
% **************************************************************
h407 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','DMT2(''max_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[28 15 12 2],...
'String','0.999',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'DMT2(''max_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','max',...
'Behavior',get(0,'defaultuicontrolBehavior'));
% **************************************************************

% **************************************************************
% Tiedosto (TEXTBOX) 4xx: min
% **************************************************************
h408 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','DMT2(''min_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[28 12 12 2],...
'String','0.98',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'DMT2(''min_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','min',...
'Behavior',get(0,'defaultuicontrolBehavior'));
% **************************************************************

% **************************************************************
% Tiedosto (TEXTBOX) 4xx: min
% **************************************************************
h409 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','DMT2(''reportfile_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[5 2 25 2],...
'String','report',...
'Style','edit',...
'Visible','on',...
'CreateFcn', {@local_CreateFcn, 'DMT2(''reportfile_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','reportfile',...
'Behavior',get(0,'defaultuicontrolBehavior'));
% **************************************************************

% **************************************************************
% vuoto nopeus (TEXTBOX) 4xx: 
% **************************************************************
h410 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','DMT2(''leakagespeed_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[5 8 15 2],...
'String','0.0003',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'DMT2(''leakagespeed_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','leakagespeed',...
'Behavior',get(0,'defaultuicontrolBehavior'));
% **************************************************************

% **************************************************************
% vuoto kiihtyvyys (TEXTBOX) 4xx: 
% **************************************************************
h411 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','DMT2(''leakageacceleration_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[28 8 10 2],...
'String','0.00003',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'DMT2(''leakageacceleration_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','leakageacceleration',...
'Behavior',get(0,'defaultuicontrolBehavior'));
% **************************************************************

% **************************************************************
% vuoto alkaminen (TEXTBOX) 4xx: 
% **************************************************************
h412 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','DMT2(''leakagestarttime_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[40 8 5 2],...
'String','0.75',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'DMT2(''leakagestarttime_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','leakagestarttime',...
'Behavior',get(0,'defaultuicontrolBehavior'));
% **************************************************************







% **************************************************************
% **************************************************************
% TEXTS
% **************************************************************
% **************************************************************

% **************************************************************
% Sample Size text
% **************************************************************
h500 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Position',[10 20.5 18 2],...
'String','Sample size',...
'Style','text',...
'Tag','ssText',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************

% **************************************************************
% Frame Size X text
% **************************************************************
h501 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Position',[10 17.5 18 2],...
'String','Frame size X',...
'Style','text',...
'Tag','fsXText',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************

% **************************************************************
% Frame Size dX text
% **************************************************************
h502 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Position',[10 14.5 18 2],...
'String','Frame size dX',...
'Style','text',...
'Tag','fsdXText',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************

% **************************************************************
% Frame Size sX text
% **************************************************************
h503 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Position',[10 11.5 18 2],...
'String','Frame size sX',...
'Style','text',...
'Tag','fssXText',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************

% **************************************************************
% Frame Size u0 text: Target
% **************************************************************
h504 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Position',[22 32.4 30 1.5],...
'String','-',...
'Style','text',...
'HorizontalAlignment','left',...
'Tag','u0Text',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************

% **************************************************************
% Frame Size u1 text: FIG1
% **************************************************************
h505 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Position',[22 29.9 30 1.5],...
'String','-',...
'Style','text',...
'HorizontalAlignment','left',...
'Tag','u1Text',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************

% **************************************************************
% Frame Size u2 text: FIG2
% **************************************************************
h506 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Position',[22 27.4 30 1.5],...
'String','-',...
'Style','text',...
'HorizontalAlignment','left',...
'Tag','u2Text',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************

% **************************************************************
% Info-box: loading, calcuting ... etc. 
% **************************************************************
h507 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Position',[5 57 50 2],...
'String','Please start data mining',...
'Style','text',...
'Tag','infotext',...
'HorizontalAlignment','left',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************
% Lis�� noita m��reit�: doc uicontrol
% lis�t��n tonne v�leihin vaan...
% **************************************************************

% **************************************************************
% Info-box2: DATASET 
% **************************************************************
h508 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Position',[5 59.5 110 2],...
'String',' ',...
'Style','text',...
'Tag','infotext2',...
'HorizontalAlignment','left',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************

% **************************************************************
% Normalvarince text
% **************************************************************
h509 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Position',[40 20.5 15 2],...
'String','\sigma',...
'Style','text',...
'Tag','varianceText',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************

% **************************************************************
% alpha text
% **************************************************************
h510 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Position',[40 17.5 15 2],...
'String','\alpha',...
'Style','text',...
'Tag','alphaText',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************

% **************************************************************
% lambda max text
% **************************************************************
h511 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Position',[40 14.5 15 2],...
'String','\lambda max',...
'Style','text',...
'Tag','maxText',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************

% **************************************************************
% lambda min text
% **************************************************************
h512 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Position',[40 11.5 15 2],...
'String','\lambda min',...
'Style','text',...
'Tag','minText',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );
% **************************************************************







% **************************************************************
% Radiobuttons *3
% **************************************************************
h600 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Callback','DMT2(''u0_Callback'',gcbo,[],guidata(gcbo))',...
'Position',[5 32.5 16 2],...
'String','Target:',...
'Style','Radio',...
'Tag','u0');
% **************************************************************
h601 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Value',1,...
'Callback','DMT2(''u1_Callback'',gcbo,[],guidata(gcbo))',...
'Position',[5 30 16 2],...
'String','Figure 1:',...
'Style','Radio',...
'Tag','u1');
% **************************************************************
h602 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Callback','DMT2(''u2_Callback'',gcbo,[],guidata(gcbo))',...
'Position',[5 27.5 16 2],...
'String','Figure 2:',...
'Style','Radio',...
'Tag','u2');
% **************************************************************

% **************************************************************
%appdata = [];
%appdata.lastValidTag = 'hpop';

h701 = uicontrol('Style', 'popup',...
       'Parent',h1,...
       'String', 'Original (tight)|Scaled (tight)|Scaled|Difference|Scaled difference|MSDV|Original|Correlation|Alarms|SOM',...
       'Position', [30 515 270 10],...
       'Callback','DMT2(''hpop_Callback'',gcbo,[],guidata(gcbo))',...
       'Tag','hpop');

h702 = uicontrol('Style', 'popup',...
       'Parent',h1,...
       'String', '-',...
       'Position', [30 555 270 10],...
       'Callback','DMT2(''hpop2_Callback'',gcbo,[],guidata(gcbo))',...
       'Tag','hpop2');

h703 = uicontrol('Style', 'popup',...
       'Parent',h1,...
       'String', 'Volatility|Smoothing indeces|Margins|Boxplot|Delay estimation|Segmentation',...
       'Position', [389 175 300 10],...
       'Callback','DMT2(''hpop3_Callback'',gcbo,[],guidata(gcbo))',...
       'Tag','hpop3');  
   
% Margin = max-min values for variables   

%set(h601,'SelectionChangeFcn',@selcbk);
%set(h601,'SelectedObject',[]);  % No selection
%set(h601,'Visible','on');


%h601 = uibuttongroup('visible','off','Position',[0 0 .2 1]);
% Create three radio buttons in the button group.
%u0 = uicontrol('Style','Radio','String','Option 1',...
%    'pos',[10 350 100 30],'parent',h1,'HandleVisibility','off');
%u1 = uicontrol('Style','Radio','String','Option 2',...
%    'pos',[10 250 100 30],'parent',h1,'HandleVisibility','off');
%u2 = uicontrol('Style','Radio','String','Option 3',...
%    'pos',[10 150 100 30],'parent',h1,'HandleVisibility','off');
% Initialize some button group properties. 
%set(h,'SelectionChangeFcn',@selcbk);
%set(h,'SelectedObject',[]);  % No selection
%set(h,'Visible','on');






% **************************************************************

hsingleton = h1;


% --- Set application data first then calling the CreateFcn. 
function local_CreateFcn(hObject, eventdata, createfcn, appdata)

if ~isempty(appdata)
   names = fieldnames(appdata);
   for i=1:length(names)
       name = char(names(i));
       setappdata(hObject, name, getfield(appdata,name));
   end
end

if ~isempty(createfcn)
   eval(createfcn);
end


% --- Handles default GUIDE GUI creation and callback dispatch
function varargout = gui_mainfcn(gui_State, varargin)

gui_StateFields =  {'gui_Name'
                    'gui_Singleton'
                    'gui_OpeningFcn'
                    'gui_OutputFcn'
                    'gui_LayoutFcn'
                    'gui_Callback'};
gui_Mfile = '';
for i=1:length(gui_StateFields)
    if ~isfield(gui_State, gui_StateFields{i})
        error('Could not find field %s in the gui_State struct in GUI M-file %s', gui_StateFields{i}, gui_Mfile);        
    elseif isequal(gui_StateFields{i}, 'gui_Name')
        gui_Mfile = [gui_State.(gui_StateFields{i}), '.m'];
    end
end

numargin = length(varargin);

if numargin == 0
    % DMT2
    % create the GUI
    gui_Create = 1;
elseif isequal(ishandle(varargin{1}), 1) && ispc && iscom(varargin{1}) && isequal(varargin{1},gcbo)
    % DMT2(ACTIVEX,...)    
    vin{1} = gui_State.gui_Name;
    vin{2} = [get(varargin{1}.Peer, 'Tag'), '_', varargin{end}];
    vin{3} = varargin{1};
    vin{4} = varargin{end-1};
    vin{5} = guidata(varargin{1}.Peer);
    feval(vin{:});
    return;
elseif ischar(varargin{1}) && numargin>1 && isequal(ishandle(varargin{2}), 1)
    % DMT2('CALLBACK',hObject,eventData,handles,...)
    gui_Create = 0;
else
    % DMT2(...)
    % create the GUI and hand varargin to the openingfcn
    gui_Create = 1;
end

if gui_Create == 0
    varargin{1} = gui_State.gui_Callback;
    if nargout
        [varargout{1:nargout}] = feval(varargin{:});
    else
        feval(varargin{:});
    end
else
    if gui_State.gui_Singleton
        gui_SingletonOpt = 'reuse';
    else
        gui_SingletonOpt = 'new';
    end
    
    % Open fig file with stored settings.  Note: This executes all component
    % specific CreateFunctions with an empty HANDLES structure.
    
    % Do feval on layout code in m-file if it exists
    if ~isempty(gui_State.gui_LayoutFcn)
        gui_hFigure = feval(gui_State.gui_LayoutFcn, gui_SingletonOpt);
        % openfig (called by local_openfig below) does this for guis without
        % the LayoutFcn. Be sure to do it here so guis show up on screen.
        movegui(gui_hFigure,'onscreen')
    else
        gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt);            
        % If the figure has InGUIInitialization it was not completely created
        % on the last pass.  Delete this handle and try again.
        if isappdata(gui_hFigure, 'InGUIInitialization')
            delete(gui_hFigure);
            gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt);            
        end
    end
    
    % Set flag to indicate starting GUI initialization
    setappdata(gui_hFigure,'InGUIInitialization',1);

    % Fetch GUIDE Application options
    gui_Options = getappdata(gui_hFigure,'GUIDEOptions');
    
    if ~isappdata(gui_hFigure,'GUIOnScreen')
        % Adjust background color
        if gui_Options.syscolorfig 
            set(gui_hFigure,'Color', get(0,'DefaultUicontrolBackgroundColor'));
        end

        % Generate HANDLES structure and store with GUIDATA. If there is
        % user set GUI data already, keep that also.
        data = guidata(gui_hFigure);
        handles = guihandles(gui_hFigure);
        if ~isempty(handles)
            if isempty(data)
                data = handles;
            else
                names = fieldnames(handles);
                for k=1:length(names)
                    data.(char(names(k)))=handles.(char(names(k)));
                end
            end
        end
        guidata(gui_hFigure, data);
    end
    
    % If user specified 'Visible','off' in p/v pairs, don't make the figure
    % visible.
    gui_MakeVisible = 1;
    for ind=1:2:length(varargin)
        if length(varargin) == ind
            break;
        end
        len1 = min(length('visible'),length(varargin{ind}));
        len2 = min(length('off'),length(varargin{ind+1}));
        if ischar(varargin{ind}) && ischar(varargin{ind+1}) && ...
                strncmpi(varargin{ind},'visible',len1) && len2 > 1
            if strncmpi(varargin{ind+1},'off',len2)
                gui_MakeVisible = 0;
            elseif strncmpi(varargin{ind+1},'on',len2)
                gui_MakeVisible = 1;
            end
        end
    end
    
    % Check for figure param value pairs
    for index=1:2:length(varargin)
        if length(varargin) == index || ~ischar(varargin{index})
            break;
        end
        try set(gui_hFigure, varargin{index}, varargin{index+1}), catch break, end
    end

    % If handle visibility is set to 'callback', turn it on until finished
    % with OpeningFcn
    gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
    if strcmp(gui_HandleVisibility, 'callback')
        set(gui_hFigure,'HandleVisibility', 'on');
    end
    
    feval(gui_State.gui_OpeningFcn, gui_hFigure, [], guidata(gui_hFigure), varargin{:});
    
    if ishandle(gui_hFigure)
        % Update handle visibility
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);
        
        % Make figure visible
        if gui_MakeVisible
            set(gui_hFigure, 'Visible', 'on')
            if gui_Options.singleton 
                setappdata(gui_hFigure,'GUIOnScreen', 1);
            end
        end

        % Done with GUI initialization
        rmappdata(gui_hFigure,'InGUIInitialization');
    end
    
    % If handle visibility is set to 'callback', turn it on until finished with
    % OutputFcn
    if ishandle(gui_hFigure)
        gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
        if strcmp(gui_HandleVisibility, 'callback')
            set(gui_hFigure,'HandleVisibility', 'on');
        end
        gui_Handles = guidata(gui_hFigure);
    else
        gui_Handles = [];
    end
    
    if nargout
        [varargout{1:nargout}] = feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    else
        feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    end
    
    if ishandle(gui_hFigure)
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);
    end
end    

function gui_hFigure = local_openfig(name, singleton)

% this application data is used to indicate the running mode of a GUIDE
% GUI to distinguish it from the design mode of the GUI in GUIDE.
setappdata(0,'OpenGuiWhenRunning',1);

% openfig with three arguments was new from R13. Try to call that first, if
% failed, try the old openfig.
try 
    gui_hFigure = openfig(name, singleton, 'auto');
catch
    % OPENFIG did not accept 3rd input argument until R13,
    % toggle default figure visible to prevent the figure
    % from showing up too soon.
    gui_OldDefaultVisible = get(0,'defaultFigureVisible');
    set(0,'defaultFigureVisible','off');
    gui_hFigure = openfig(name, singleton);
    set(0,'defaultFigureVisible',gui_OldDefaultVisible);
end
rmappdata(0,'OpenGuiWhenRunning');

