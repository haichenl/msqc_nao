classdef ZMatrix < handle
    %ZMatrix is a list of atoms from the ZAtom class
    %Pars can be a class that is constructed as needed.
    %Pars should be a list of bond lengths, angles and dihedral values
    %For more info on pars class, go to Pars_template class
    
    properties
        atoms = {};
        pars = {};      %Pars is in the following form, which does have to
                        %be it's own class. This could change if the 
                        %implementation is improved for auto-generating
                        %config files. Here is an example of what the
                        %class needs to have implemented for the 
                        %functions here to work:
                        %
                        %   pars.bond_pars = {}  List of bond vals in order
                        %   pars.ang_pars = {}   List of ang vals in order
                        %   pars.di_pars = {}    List of di vals in order
                        %
                        %Any other stuff the class has would be for the
                        %user, as the project needs.
                        
        file_name_base;
    end
    
    methods
        function zmat = ZMatrix()
            %Constructor
            %
            %file_name_base is automatically initialized to a random unique
            %name but can be changed as needed

            zmat.atoms = {};
            zmat.pars.bond_pars = [];
            zmat.pars.ang_pars = [];
            zmat.pars.di_pars = [];
        end
           
        function make_atom( zmat, z, bond_ref, ang_ref, di_ref )
            %z should be an int
            %All inputs after z should be nonnegative integers
            %None of those inputs should be larger than the length of atoms
            %Has error messages for if conditions are not met (for range
            %and not having overlapping references, not for type checks)
            
            % Num of args and error handling
            atom_num = length( zmat.atoms ) + 1;
            if nargin > 4
                if di_ref ~= 0
                    if di_ref == ang_ref
                        error('di_ref cannot equal ang_ref');
                    end
                    if di_ref == bond_ref
                        error('di_ref cannot equal ang_ref');
                    end
                    if di_ref > atom_num
                        error('di_ref must reference an existing atom in zmat.atoms');
                    end
                end
            else
                di_ref = 0;
            end
            
            if nargin > 3
                if ang_ref ~= 0
                    if ang_ref == bond_ref
                        error('ang_ref cannot equal bond_ref');
                    end
                    if ang_ref > atom_num
                        error('ang_ref must reference an existing atom in zmat.atoms');
                    end
                end
            else
                ang_ref = 0;
            end
            
            if nargin > 2
                if bond_ref > atom_num
                    error('bond_ref must reference an existing atom in zmat.atoms');
                end
            else
                bond_ref = 0;
            end
            
            if nargin == 1
                z = 0;
            end
            %
            atom = ZAtom( z, atom_num, bond_ref, ang_ref, di_ref );
            zmat.atoms{ atom_num } = atom;
            
            if bond_ref ~= 0
                zmat.atoms{ bond_ref }.up_bond_total()
            end
        end
        
        function text = build_atoms( zmat, bq )
            %Calls print_atom for each atom in the zmat
            newLine = char(10);
            if nargin < 2
                bq = 0;
            end
            
            text = '';
            for i = 1:length( zmat.atoms )
                atom = zmat.atoms{i};
                if i == bq || bq == 0
                    text = [text, atom.atom_text( 0)];
                else
                    text = [text, atom.atom_text( 1 )];
                end
            end
            text = [text, newLine];
        end
        
        
        
        function text = build_molstr( zmat, bq )
            
            if nargin < 2
                bq = 0;
            end
            text = zmat.build_atoms( bq );
            text = [text, zmat.build_molstr_pars()];
        end
        
        function text = build_molstr_pars( zmat )
            %Prints B, A, and D vals with the PAR instead of numbers
            newLine = char(10);
            formSpec = '%f';
            
            num_atoms = length(zmat.atoms);
            num_bonds = num_atoms - 1;
            num_angs = num_atoms - 2;
            num_dis = num_atoms - 3;
            text = '';
            for i = 1:num_bonds
                num_str = num2str( real(zmat.pars.bond_pars(i)), formSpec );
                i_str = num2str(i);
                temp = ['   B', i_str, '     =        ', num_str, newLine];
                text = [text, temp];
            end
            for i = 1:num_angs
                num_str = num2str( real(zmat.pars.ang_pars(i)), formSpec );
                i_str = num2str(i);
                temp = ['   A', i_str, '     =        ', num_str, newLine];
                text = [text, temp];
            end
            for i = 1:num_dis
                num_str = num2str( real(zmat.pars.di_pars(i)), formSpec );
                i_str = num2str(i);
                temp = ['   D', i_str, '     =        ', num_str, newLine];
                text = [text, temp];
            end
            text = [text, newLine];
        end
        
    end
    
end

