function C = C_ijkl(E, DF, D, nsd, aux_size)
DF = permute(DF, [2 1 3]);
DF_minus_T = geopdes_inv__(DF);
C = zeros(nsd,nsd,nsd,nsd,size(DF,3));
for i=1:nsd
    for j=1:nsd
        for k=1:nsd
            for ell=1:nsd
                Eik = Ehat(E,DF_minus_T,i,k);
                Eik = permute(Eik,[2 1 3]);
                tmp = multiprod(Eik,D);
                Ejl = Ehat(E,DF_minus_T,j,ell);
                C(i,j,k,ell,:) = multiprod(tmp,Ejl);
            end
        end
    end
end
DF_det = geopdes_det__(DF);
C = reshape(C,[nsd^4,size(DF,3)]);
C = C.*DF_det';
C = reshape(C,[nsd,nsd,nsd,nsd,size(DF,3)]);
Cee = cell(nsd,nsd);
for i=1:nsd
    for j=1:nsd
        tmp = squeeze(C(i,j,:,:,:));
        tmp = reshape(tmp,[nsd, nsd, aux_size]);
        Cee{i,j} = permute(tmp,[3:2+nsd 1 2]);
    end
end
C = Cee;
end