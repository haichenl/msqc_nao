classdef FastJK < handle
    
    properties
        
        % bigN by bigN matrices
        teiMatJ; % J and K in total occupy ~ n^4 / 2 * 8 bytes of memory 
        teiMatK; % e.g. nbasis = 316 => 40 gb; but the initialization 
                 % requires ~ n^4 / 1.5 * 8 bytes, as nbasis = 294 => 40 gb
        
        % frequently used numbers 
        nbasis; % number of basis functions 
        bigN; % nchoosek(nbasis, 2) + nbasis 
        uniqN; % nchoosek(bigN, 2) + bigN 
        
        % an array holding information for getting rows 
        arr;
        
    end
    
    methods
        
        function res = InCoreJK(matpsi)
            
            [teiVecJ, teiVecK] = matpsi.tei_alluniqJK();
            res.nbasis = matpsi.nbasis();
            res.bigN = res.nbasis * ( res.nbasis + 1 ) / 2;
            res.uniqN = ( res.nbasis * ...
                        ( res.nbasis + 1 ) * ...
                        ( res.nbasis * res.nbasis + res.nbasis + 2 ) )...
                        / 8;
                    
            res.teiMatJ = zeros(res.bigN);
            res.teiMatJ(triu(true(res.bigN))) = teiVecJ;
            teiVecJ = []; % frees memory 
            res.teiMatK = zeros(res.bigN);
            res.teiMatK(triu(true(res.bigN))) = teiVecK;
            teiVecK = []; % frees memory 
            % Hermitianize J and K 
            for i = 1:res.bigN
                res.teiMatJ(i+1:res.bigN, i) = res.teiMatJ(i, i+1:res.bigN);
                res.teiMatK(i+1:res.bigN, i) = res.teiMatK(i, i+1:res.bigN);
            end

            res.arr = 1 : res.bigN;
            res.arr = res.arr .* ( res.arr - 1 ) ./ 2 + 1;

        end
        
        function jmat = ComputeJ(obj, densMat)
            % in-core algorithm
            
            % scale density vector 
            densVec = 2 * densMat(triu(true(obj.nbasis)));
            diagind = 1 : obj.nbasis;
            diagind = diagind .* ( diagind + 1 ) / 2;
            densVec(diagind) = 0.5 * densVec(diagind);
            
            % form the uniq J vector 
            juniqvec = res.teiMatJ * densVec;
            
            % form the full J matrix, i.e. Hermitianize 
            jmat = zeros(obj.nbasis);
            jmat(triu(true(obj.nbasis))) = juniqvec;
            jmat = jmat + jmat' - diag(diag(jmat));
            
        end
        
        function kmat = Kmat(obj, densMat)
            % in-core algorithm
            
            densVec = -0.5 * densMat(triu(true(obj.nbasis)));
            diagind = 1 : obj.nbasis;
            diagind = diagind .* ( diagind + 1 ) / 2;
            densVec(diagind) = 0.5 * densVec(diagind);
            
            kuniqvec = res.teiMatK * densVec;
            
            kmat = zeros(obj.nbasis);
            kmat(triu(true(obj.nbasis))) = kuniqvec;
            kmat = kmat + kmat' - diag(diag(kmat));
            
        end
        
    end
    
end