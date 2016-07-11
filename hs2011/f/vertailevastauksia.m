function vertailevastauksia(hs,p,h1,h2)
h1 = find(strcmp(hs.Nimi,h1));
h2 = find(strcmp(hs.Nimi,h2));
disp('************* ');
j=0;
%disp(p.labelX);
for i=p.X1
    j=j+1;
    disp('-----');
    disp(p.katainenpaino(j));
    disp(hs.k.Kysymys{i});
    disp(hs.k.Vastaus{h1,i});
    if strcmp(hs.k.Vastaus{h1,i},hs.k.Vastaus{h2,i})
        disp('SAMA vastaus');
    else
        disp(hs.k.Vastaus{h2,i});
    end
end
disp('************ ');
%disp(p.labelY);
for i=p.X2
    j=j+1;
    disp('-----');
    disp(p.katainenpaino(j));
    disp(hs.k.Kysymys{i});
    disp(hs.k.Vastaus{h1,i});
    if strcmp(hs.k.Vastaus{h1,i},hs.k.Vastaus{h2,i})
        disp('SAMA vastaus');
    else
        disp(hs.k.Vastaus{h2,i});
    end
end
disp('************ ');
end