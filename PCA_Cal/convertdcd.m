
function convertdcd(filename, file_no)

% xyz = readdcd(filename, indices)
% reads an dcd and puts the x,y,z coordinates corresponding to indices 
% in the rows of x,y,z

h = read_dcdheader(filename)
nsets = h.NSET;
natoms = h.N;
numind = length(natoms);

x = zeros(natoms,1);
y = zeros(natoms,1);
z = zeros(natoms,1);

if nsets == 0
  xyz = zeros(1,3*numind, 'single');
  nsets = 99999;
else
  xyz = zeros(nsets, 3*natoms, 'single');
end

for i=1:nsets
  pos = ftell(h.fid);
  if pos == h.endoffile 
    break;
  end
  [x,y,z] = read_dcdstep(h);
  xyz(i,1:3:3*natoms) = x(1:natoms)';
  xyz(i,2:3:3*natoms) = y(1:natoms)';
  xyz(i,3:3:3*natoms) = z(1:natoms)';
end

TRJ=xyz';
ss=['TRJ' num2str(file_no) '= TRJ;'];
eval(ss);
save(['TRJ' num2str(file_no) '.mat'],['TRJ' num2str(file_no)],'-v7.3');

end
