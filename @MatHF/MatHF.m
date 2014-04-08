classdef MatHF < handle
    
    properties
        zmat;
        matpsi;
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
        
        % Hartree-Fock results
        density;
        Ehf;
        Eorb;
        orb;
        
    end
    
    methods
        
        % constructor
        function obj = MatHF(zmat_, basisname)
            obj.zmat = zmat_;
            obj.matpsi = MatPsi({obj.zmat.build_molstr(), basisname});
            
            obj.Enuc = obj.matpsi.Enuc();
            obj.natom = obj.matpsi.natom();
            obj.nelec = obj.matpsi.nelec();
            
            obj.nbasis = obj.matpsi.nbasis();
            obj.basisAtom = obj.matpsi.func2center();
            
            obj.S = obj.matpsi.overlap();
            obj.H1 = obj.matpsi.kinetic() + obj.matpsi.potential();
            
            obj.eps = 1.0e-8;
            obj.maxIter = 100;
            
            obj.density = [];
            obj.Ehf = [];
            obj.Eorb = [];
            obj.orb = [];
            
            obj.UseFastJK();
            
        end
        
        % fastjk object initializer 
        function UseFastJK(obj)
            obj.fastjk = FastJK(obj.matpsi);
        end
        
        function Ehf_ = doHF(obj)
            
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
                
                % Step 6 -- Obtain F (fock matrix); Fast algorithm 
                F = obj.H1 + 2 .* obj.fastjk.Jmat(P) - obj.fastjk.Kmat(P);
                
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
                
                Pn = C(:,filled)*( C(:,filled)');
                iter = iter + 1;
                
                DeltaDens = reshape(P - Pn, 1, obj.nbasis * obj.nbasis);
                rmsDeltaDens = DeltaDens * DeltaDens';
                Ehftemp = sum(sum(Pn.*(obj.H1+F)));
                DeltaE = abs(Ehftemp - Ehfsave);
                if (iter > obj.maxIter)
                    finished = true;
                elseif (rmsDeltaDens < obj.eps || DeltaE < obj.eps)
                    finished = true;
                end
                Ehfsave = Ehftemp;
            end
            % End of iteration of steps 5-11.
            obj.density = Pn;
            
            % Step 12: Output.
            Ehf_ = Ehfsave + obj.Enuc;
            obj.Ehf = Ehf_;
            
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