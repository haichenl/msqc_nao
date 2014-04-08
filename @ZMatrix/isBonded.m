function res = isBonded(obj,iatom,jatom)
if (nargin == 1)
   natoms = length(obj.atoms);
   res = zeros(natoms,natoms);
   for i = 1:natoms
      for j = (i+1):natoms
         res(i,j) = checkBond(obj,i,j);
         res(j,i) = res(i,j);
      end
   end
else
   res = checkBond(obj,iatom,jatom);
end
   
end

function bonding = checkBond(obj,iatom,jatom)
iatom_bond_ref = obj.atoms{iatom}.bond_ref;
jatom_bond_ref = obj.atoms{jatom}.bond_ref;
if (jatom_bond_ref == iatom || iatom_bond_ref == jatom)
   bonding = 1;
else
   bonding = 0;
end
end

