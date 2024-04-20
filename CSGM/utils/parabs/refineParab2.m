function parab = refineParab2(parab, params)
%% Check if pixels are invalid

x = -parab.b./(2*parab.a);
idx_min = (parab.a>-1e-32);
parab.a(idx_min) = -1e-32;
parab.b(idx_min) = x(idx_min).*(-2.*parab.a(idx_min)); % b - is adjusted so x is the minimum x = b/-2a -> b = -2ax
nanparab = isnan(parab.b);
parab.b(nanparab) = 0;
parab.a(nanparab) = -1e-102;

end


