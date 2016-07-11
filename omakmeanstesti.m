


setti = 8;
klustereita = 3;


figure(1);

k=0;
%[a b] = k_means(vika{setti}.g.data,klustereita);
[a b] = k_means(d.target.aX,klustereita);
for i=1:5
    for j=1:4
        k=k+1;
        subplot(5,4,k);
        if k==20
            plot(a,'x')
        else
            %plot(vika{setti}.g.data(:,k));
            plot(d.target.aX(:,k));
        end
    end
end


figure(2);

k=0;
%[a b] = k_means(vika{setti}.g.data,klustereita);
[a b] = k_means(d.target.defX,klustereita);
for i=1:5
    for j=1:4
        k=k+1;
        subplot(5,4,k);
        if k==20
            plot(a,'x')
        else
            %plot(vika{setti}.g.data(:,k));
            plot(d.target.aX(:,k));
        end
    end
end
