classdef InCoreJK < handle
    
    properties
        
        % 2*uniqN length vector, the first uniqN length holding a "normal"
        % two-electron integral unique vector, the last uniqN length
        % holding an arranged one for the forming of K 
        teiVecJ; % occupy ~ n^4 / 4 * 8 bytes of memory. nbasis = 375 => 40 gb
        teiVecK; % occupy ~ n^4 / 4 * 8 bytes of memory. nbasis = 375 => 40 gb
        
        % frequently used numbers 
        nbasis; % number of basis functions 
        bigN; % nchoosek(nbasis, 2) + nbasis 
        uniqN; % nchoosek(bigN, 2) + bigN 
        
        % an array holding information for getting rows 
        arr;
        
    end
    
    methods
        
        function res = InCoreJK(matpsi)
            [res.teiVecJ, res.teiVecK] = matpsi.tei_alluniqJK();
            res.nbasis = matpsi.nbasis();
            res.bigN = res.nbasis * ( res.nbasis + 1 ) / 2;
            res.uniqN = ( res.nbasis * ...
                        ( res.nbasis + 1 ) * ...
                        ( res.nbasis * res.nbasis + res.nbasis + 2 ) )...
                        / 8;
            res.arr = 1 : res.bigN;
            res.arr = res.arr .* ( res.arr - 1 ) ./ 2 + 1;

        end
        
        function jmat = ComputeJ(obj, densMat)
            % in-core algorithm
            
            % scale density vector 
            densVec = 2 * densMat(triu(true(obj.nbasis)));
            for i = 1 : obj.nbasis
                densVec( i*(i+1)/2 ) = 0.5 * densVec( i*(i+1)/2 );
            end
            
            % allocate J vector 
            juniqvec = zeros(1, obj.bigN);
            
            % extract every row and multiply it with the scaled density vector 
            for i = 1 : obj.bigN
                juniqvec(i) = obj.JRow( i ) * densVec;
            end
            
            % forming the full J matrix 
            jmat = zeros(obj.nbasis);
            jmat(triu(true(obj.nbasis))) = juniqvec;
            jmat = jmat + jmat' - diag(diag(jmat));
            
        end
        
        function kmat = ComputeK(obj, densMat)
            % in-core algorithm
            
            densVec = -0.5 * densMat(triu(true(obj.nbasis)));
            for i = 1 : obj.nbasis
                densVec( i*(i+1)/2 ) = 0.5 * densVec( i*(i+1)/2 );
            end
            kuniqvec = zeros(1, obj.bigN);
            for i = 1 : obj.bigN
                kuniqvec(i) = obj.KRow( i ) * densVec;
            end
            kmat = zeros(obj.nbasis);
            kmat(triu(true(obj.nbasis))) = kuniqvec;
            kmat = kmat + kmat' - diag(diag(kmat));
            
        end
        
        function jrow = JRow(obj, ind)
            if(ind == obj.bigN)
                append = [];
            else
                tmp = obj.teiVecJ( obj.arr + ind - 1);
                append = tmp(end-(obj.bigN-ind)+1 : end);
            end
            jrow = [ (obj.teiVecJ( (ind-1)*ind/2+1 : ind*(ind+1)/2 ))', append'];
        end
        
        function krow = KRow(obj, ind)
            if(ind == obj.bigN)
                append = [];
            else
                tmp = obj.teiVecK(obj.arr + ind - 1);
                append = tmp( end-(obj.bigN-ind)+1 : end );
            end
            krow = [ (obj.teiVecK( (ind-1)*ind/2+1 : ind*(ind+1)/2 ))', append'];
        end
        
    end
    
end