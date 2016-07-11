function b=kysymyksetjavastaukset(hs,m,p)
% *************************************************************************
%tyhjia = 0;
for i=1:31
    disp(' ');
    disp(num2str(i))
    disp(hs.k.Kysymys{i});
    apu = unique(hs.k.Vastaus(:,i));
    k=0;
    for j=1:length(apu)
        disp(apu{j});
        if strcmp(apu{j},'-')
            %tyhjia = tyhjia + 1;  % 31
        else
            k=k+1;
            b.v{i}.v{k} = apu{j};
        end
    end
    
end
% *************************************************************************
% 21 ja 31 monivalinta kysymyksia. Pilkotaan osiin.
% *************************************************************************


end