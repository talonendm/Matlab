function p = pvarit(m,p)

for i=1:m.n.puolueet
    p.vari{i} = [0.9 0.9 0.9];
end

p.vari{2} = [0.8 0.3 0.8];
p.vari{3} = [0.8 0.3 0.4];
p.vari{4} = [0.6 0.6 0.9];
p.vari{5} = [0 0.7 0];
p.vari{6} = [0.1 0.1 1];
p.vari{7} = [0.45 0 0];

p.vari{12} = [0.3 0.2 0.1];
p.vari{16} = [0.3 0 0];
p.vari{17} = [0.2 0.2 0.5];
p.vari{18} = [1 0.7 0.2];
p.vari{19} = [1 0 0];
%p.vari{20} = [0.3 0.3 0.5];
p.vari{20} = [0.65 0 0];
p.vari{21} = [0.6 0 0];
p.vari{22} = [0.7 0 0];

p.vari{25} = [0.70 0 0];
p.vari{26} = [0 0.9 0];
p.vari{27} = [0.4 0.9 0.5];


% 'Helminen'
 %   'IPU'
  %  'KA'
   % 'KD'
%    'KESK'
%    'KOK'
%    'KTP'
 %   'Korhonen'
%    'Kotajärvi'
 %   'Laajola'
%    'Lönnroth'
%    'M2011'
%    'Marttunen'
%    'Nikoskinen'
 %   'Niskanen'
%    'PIR'
%    'PS'
%    'RKP'
%    'SDP'
%    'SEN'
%    'SKP'
%    'STP'
%    'Salminen'
 %   'Simonen'
 %   'VAS'
%    'VIHR'
%    'VP'
%    'Valjakka'
%    'Yl-Sit'
end