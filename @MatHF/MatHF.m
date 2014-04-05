classdef MatHF < handle
    properties
        matpsi;
        BasisProp;
        Int;
        HFpars;
    end
    methods
        
        % constructor
        function res = MatHF(geomstr, basisname)
            res.matpsi = MatPsi({geomstr, basisname});
            nbasis = res.matpsi.nbasis();
            res.BasisProp.nbasis = res.matpsi.nbasis();
            res.Int.S = res.matpsi.overlap();
            res.Int.H1 = res.matpsi.kinetic() + res.matpsi.potential();
            res.Int.H2 = zeros(nbasis, nbasis, nbasis, nbasis);
            h2uniq = res.matpsi.tei_alluniq();
            for i = 1:nbasis
                for j = 1:nbasis
                    for k = 1:nbasis
                        for l = 1:nbasis
                            res.Int.H2(i,j,k,l) = h2uniq(MatHF.ij2I(MatHF.ij2I(i, j), MatHF.ij2I(k, l)));
                        end
                    end
                end
            end
            res.HFpars.eps = 1.0e-8;
            res.HFpars.maxIter = 100;
        end
    end
        
    methods (Static)
        
        function res = ij2I(i, j)
            if(i<j)
                tmp = i;
                i = j;
                j = tmp;
            end
            res = i.*(i-1)./2 + j;
        end
        
    end
end