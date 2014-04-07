classdef MatHF < handle
    
    properties
        matpsi;
        incorejk;
        
        MolProp;
        BasisProp;
        OEI; % one-electron integrals 
        HFpars;
        
        HFresults;
        
    end
    
    methods
        
        % constructor
        function res = MatHF(geomstr, basisname)
            res.matpsi = MatPsi({geomstr, basisname});
            
            res.MolProp.Enuc = res.matpsi.Enuc();
            res.MolProp.natom = res.matpsi.natom();
            res.MolProp.nelec = res.matpsi.nelec();
            
            res.BasisProp.nbasis = res.matpsi.nbasis();
            res.BasisProp.basisAtom = res.matpsi.func2center();
            
            res.OEI.S = res.matpsi.overlap();
            res.OEI.H1 = res.matpsi.kinetic() + res.matpsi.potential();
            
            res.HFpars.eps = 1.0e-8;
            res.HFpars.maxIter = 100;
            res.HFpars.minIter = 5;
            
            res.HFresults.density = [];
            res.HFresults.Ehf = [];
            res.HFresults.Eorb = [];
            res.HFresults.orb = [];
        end
        
        % incorejk object initializer 
        function UseInCoreJK(obj)
            obj.incorejk = InCoreJK(obj.matpsi);
        end
        
        function doHF(obj)
            filled = 1:(obj.MolProp.nelec/2);
            H1 = obj.OEI.H1;
            % Step 3 -- Calculate transformation matrix (eq. 3.167).
            X = inv(sqrtm(obj.OEI.S));
            
            % Step 4 -- Guess at density matrix -- all zeros right now.
            if ( isempty(obj.HFresults.density) )
                Pn = zeros(obj.BasisProp.nbasis);
            else
                Pn = obj.HFresults.density;
            end
            
            iter = 0;
            finished = false;
            
            Ehfsave = 0;
            % Begin iteration through.
            while (~finished)  % Step 11 -- Test convergence
                
                P = Pn;
                
                % Step 6 -- Obtain F (fock matrix). In-core right now
                F = H1 + obj.incorejk.ComputeJ(P) + obj.incorejk.ComputeK(P);
                
                % Step 7 -- Calculate the transformed F matrix.
                Ft = X'*F*X;
                
                % Step 8 -- Find e and the transformed expansion coefficient matrices.
                [Ct1, e1] = eig(Ft);
                e2 = diag(e1);
                [e, i1] = sort(e2);
                Ct = Ct1(:, i1);
                
                % Step 9 -- Transform Ct back to C.
                C = X*Ct;
                
                % Step 10 -- Calculate the new density matrix.
                
                Pn = 2* C(:,filled)*( C(:,filled)');
                iter = iter + 1;
                
                changeInDensity = max(max(abs(P - Pn)));
                Ehftemp = sum(sum(Pn.*(H1+F)));
                changeInEnergy = abs(Ehftemp - Ehfsave);
                if (iter > obj.HFpars.maxIter)
                    finished = true;
                elseif (iter > obj.HFpars.minIter)
                    if (changeInDensity < obj.HFpars.eps...
                            || changeInEnergy < obj.HFpars.eps)
                        finished = true;
                    end
                end
                Ehfsave = Ehftemp;
            end
            % End of iteration of steps 5-11.
            obj.HFresults.density = Pn;
            
            % Step 12: Output.
            obj.HFresults.Ehf = Ehfsave/2 + obj.MolProp.Enuc;
            
            % Orbital energies.
            obj.HFresults.Eorb = e;
            
            % Molecular orbital components.
            obj.HFresults.orb = C;
            
            
            
            if (iter+1 > obj.HFpars.maxIter)
                disp('You are living on the edge.. hartree fock didn''t converge');
            end
            
        end
        
    end
        
end