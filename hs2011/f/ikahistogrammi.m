function ikahistogrammi(hs,m,p,num)
figure(num);
for i=1:m.n.vaalipiirit
    subplot(7,2,i);
    hist(hs.Ika(m.ehdokkaat{i}),p.ikahistX);
    axis tight;
    title(m.vaalipiirit{i});
    %xlabel('ikä');
    
end
end