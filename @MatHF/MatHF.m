classdef MatHF < handle
    
    properties
        matpsi;
        incorejk;
        fastjk;
        
        % Molecule properties 
        Enuc;
        natom;
        nelec;
        
        % Basis set properties 
        nbasis;
        basisAtom;
        
        % one-electron integrals 
        S;
        H1
        
        % Hartree-Fock settings 
        eps;
        maxIter;
        minIter;
        
        % Hartree-Fock results
        density;
        Ehf;
        Eorb;
        orb;
        
    end
    
    methods
        
        % constructor
        function res = MatHF(geomstr, basisname)
            res.matpsi = MatPsi({geomstr, basisname});
            
            res.Enuc = res.matpsi.Enuc();
            res.natom = res.matpsi.natom();
            res.nelec = res.matpsi.nelec();
            
            res.nbasis = res.matpsi.nbasis();
            res.basisAtom = res.matpsi.func2center();
            
            res.S = res.matpsi.overlap();
            res.H1 = res.matpsi.kinetic() + res.matpsi.potential();
            
            res.eps = 1.0e-8;
            res.maxIter = 100;
            res.minIter = 5;
            
            res.density = [];
            res.Ehf = [];
            res.Eorb = [];
            res.orb = [];
        end
        
        % incorejk object initializer 
        function UseInCoreJK(obj)
            obj.incorejk = InCoreJK(obj.matpsi);
        end
        
        function UseFastJK(obj)
            obj.fastjk = FastJK(obj.matpsi);
        end
        
        function doHF(obj)
            
            obj.UseFastJK();
            
            filled = 1:(obj.nelec/2);
            % Step 3 -- Calculate transformation matrix (eq. 3.167).
            X = inv(sqrtm(obj.S));
            
            % Step 4 -- Guess at density matrix -- all zeros right now.
            if ( isempty(obj.density) )
                Pn = zeros(obj.nbasis);
            else
                Pn = obj.density;
            end
            
            iter = 0;
            finished = false;
            
            Ehfsave = 0;
            % Begin iteration through.
            while (~finished)  % Step 11 -- Test convergence
                
                P = Pn;
                
                % Step 6 -- Obtain F (fock matrix). In-core right now
                F = obj.H1 + obj.fastjk.ComputeJ(P) + obj.fastjk.ComputeK(P);
                
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
                Ehftemp = sum(sum(Pn.*(obj.H1+F)));
                changeInEnergy = abs(Ehftemp - Ehfsave);
                if (iter > obj.maxIter)
                    finished = true;
                elseif (iter > obj.minIter)
                    if (changeInDensity < obj.eps...
                            || changeInEnergy < obj.eps)
                        finished = true;
                    end
                end
                Ehfsave = Ehftemp;
            end
            % End of iteration of steps 5-11.
            obj.density = Pn;
            
            % Step 12: Output.
            obj.Ehf = Ehfsave/2 + obj.Enuc;
            
            % Orbital energies.
            obj.Eorb = e;
            
            % Molecular orbital components.
            obj.orb = C;
            
            
            
            if (iter+1 > obj.maxIter)
                disp('You are living on the edge.. hartree fock didn''t converge');
            end
            
        end
        
    end
        
end