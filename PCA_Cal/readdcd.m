
function xyz = readdcd(filename, ind)

% xyz = readdcd(filename, indices)
% reads an dcd and puts the x,y,z coordinates corresponding to indices 
% in the rows of x,y,z

h = read_dcdheader(filename)
nsets = h.NSET;
natoms = h.N;
numind = length(ind);

x = zeros(natoms,1);
y = zeros(natoms,1);
z = zeros(natoms,1);

if nsets == 0
  xyz = zeros(1,3*numind);
  nsets = 99999;
else
  xyz = zeros(nsets, 3*2615);
end

for i=1:nsets
  pos = ftell(h.fid);
  if pos == h.endoffile 
    break;
  end
  [x,y,z] = read_dcdstep(h);
  xyz(i,1:3:3*numind) = x(1:2615)';
  xyz(i,2:3:3*numind) = y(1:2615)';
  xyz(i,3:3:3*numind) = z(1:2615)';
end

fclose(h);

