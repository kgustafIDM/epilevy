% plot of gamma(mu,nu)

mu_space = 0:0.1:4;
nu_space = 0:0.1:2;

gamma = zeros(length(mu_space),length(nu_space));

for ii = 1:length(mu_space)
    for jj = 1:length(nu_space)
        mu = mu_space(ii);
        nu = nu_space(jj);
        muprime = mu - 1 + 1/nu;
        muprimenu = muprime*nu;
        smpn = muprimenu - 2*nu;
        if muprimenu >= 2
            if smpn >= 1
                gamma(ii,jj) = 1;
            elseif smpn < 1 && smpn > 0
                gamma(ii,jj) = 2 - smpn;
            end
        elseif muprimenu < 2 && muprimenu > 1
            if smpn >= 1
                gamma(ii,jj) = muprimenu - 1;
            elseif smpn < 1 && smpn > 0
                gamma(ii,jj) = 2*nu;
            end
        end
    end
end

figure; surf(mu_space,nu_space,gamma'); view([0 0 90]);
%set(gca,'Xlim',[2 4],'Ylim',[0 2]);
xlabel('\mu'); ylabel('\nu');
