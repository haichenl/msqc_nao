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
        
        function obj = FastJK(matpsi)
            
            [teiVecJ, teiVecK] = matpsi.tei_alluniqJK();
            obj.nbasis = matpsi.nbasis();
            obj.bigN = obj.nbasis * ( obj.nbasis + 1 ) / 2;
            obj.uniqN = ( obj.nbasis * ...
                        ( obj.nbasis + 1 ) * ...
                        ( obj.nbasis * obj.nbasis + obj.nbasis + 2 ) )...
                        / 8;
                    
            obj.teiMatJ = zeros(obj.bigN);
            obj.teiMatJ(triu(true(obj.bigN))) = teiVecJ;
            teiVecJ = []; % frees memory 
            obj.teiMatK = zeros(obj.bigN);
            obj.teiMatK(triu(true(obj.bigN))) = teiVecK;
            teiVecK = []; % frees memory 
            % Hermitianize J and K 
            for i = 1:obj.bigN
                obj.teiMatJ(i+1:obj.bigN, i) = obj.teiMatJ(i, i+1:obj.bigN);
                obj.teiMatK(i+1:obj.bigN, i) = obj.teiMatK(i, i+1:obj.bigN);
            end

            obj.arr = 1 : obj.bigN;
            obj.arr = obj.arr .* ( obj.arr - 1 ) ./ 2 + 1;

        end
        
        function jmat = Jmat(obj, densMat)
            % in-core algorithm
            
            % scale density vector 
            densVec = 2 * densMat(triu(true(obj.nbasis)));
            diagind = 1 : obj.nbasis;
            diagind = diagind .* ( diagind + 1 ) / 2;
            densVec(diagind) = 0.5 * densVec(diagind);
            
            % form the uniq J vector 
            juniqvec = obj.teiMatJ * densVec;
            
            % form the full J matrix, i.e. Hermitianize 
            jmat = zeros(obj.nbasis);
            jmat(triu(true(obj.nbasis))) = juniqvec;
            jmat = jmat + jmat' - diag(diag(jmat));
            
        end
        
        function kmat = Kmat(obj, densMat)
            % in-core algorithm
            
            densVec = densMat(triu(true(obj.nbasis)));
            diagind = 1 : obj.nbasis;
            diagind = diagind .* ( diagind + 1 ) / 2;
            densVec(diagind) = 0.5 * densVec(diagind);
            
            kuniqvec = obj.teiMatK * densVec;
            
            kmat = zeros(obj.nbasis);
            kmat(triu(true(obj.nbasis))) = kuniqvec;
            kmat = kmat + kmat' - diag(diag(kmat));
            
        end
        
    end
    
end