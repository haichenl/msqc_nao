% zmat = ZMatrix();
% zmat.file_name_base = 'met';
% zmat.make_atom(6,0,0,0);
% zmat.make_atom(1,1,0,0);
% zmat.make_atom(1,1,2,0);
% zmat.make_atom(1,1,2,3);
% zmat.make_atom(1,1,2,3);
% zmat.pars.bond_pars = [1.09 1.09 1.09 1.09];
% zmat.pars.ang_pars = [109.4712 109.4712 109.4712];
% zmat.pars.di_pars = [120.0 -120];

zmat.file_name_base = 'ethane';
zmat.make_atom(6,0,0,0); 
zmat.make_atom(6,1,0,0); 
zmat.make_atom(1,1,2,0); 
zmat.make_atom(1,1,2,3); 
zmat.make_atom(1,1,2,3); 
zmat.make_atom(1,2,1,3); 
zmat.make_atom(1,2,1,3); 
zmat.make_atom(1,2,1,3); 
zmat.pars.bond_pars = [1.54, 1.07, 1.07, 1.07, 1.07, 1.07, 1.07];
zmat.pars.ang_pars = [109.5, 109.5, 109.5, 109.5, 109.5, 109.5];
zmat.pars.di_pars = [120, -120, 180, 60, -60];

basisname = 'sto-3g';

mathf = MatHF(zmat, basisname);
