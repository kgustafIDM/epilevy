features_Levy = [pdf_slopes; disp_slopes; fp_binary];

kk = 0;
clearcases = [];    
for ii = 1:6
    for jj = [1:3,9:11]
        
        kk = kk + 1;
        clearcases(kk) = jj+11*(ii-1);
    end
end

clearsub = [];
clearsup = [];
kk = 0;

for ii = 1:6
    for jj = 1:3
        
        kk = kk + 1;
        clearsub(kk) = jj+11*(ii-1);
    end
end

kk = 0;
for ii = 1:6
    for jj = 9:11
        
        kk = kk + 1;
        clearsup(kk) = jj+11*(ii-1);
    end
end

pcc = pdf_slopes(clearcases,:);
dcc = disp_slopes(clearcases,:);

features_Levy = [pdf_slopes(clearcases,:) disp_slopes(clearcases,:) fp_binary(clearcases,:)];