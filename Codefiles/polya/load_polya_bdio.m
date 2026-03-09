% load_polya_bdio.m
%
% reads the run.polya.bdio file
%
% p = load_polya_bdio(fname)
%
% fname: path to the run.polya.bdio file
%
% p: structure containing the data from the file
%     p.L:         vector with lattice dimensions
%     p.bc:        cell-array with boundary conditions
%     p.N:         number of configurations
%     p.nc:        vector with configuration numbers
%    
%     p.C:         N x 4 matrix with average polyakov lines in each of the
%                  four directions on every configuration
%
%

% Tomasz Korzec 2022
function p = load_polya_bdio(fname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% constants
%
BDIO_PATH = '~/bdio/matlab/';

% open bdio file
finfo = dir(fname);
if isempty(finfo)
   P.N=-1;
   return;
end
addpath(BDIO_PATH);
bdio_define_constants;

fid = bdio_open(fname,'r','Generic Correlator Format 1.0');
if fid==-1
   p.N=-1;
   return;
end

% read lattice description record
i=bdio_seek_record(fid);
if i==-1
   p.N=-1;
   fprintf('No records in old data file found\n');
   return;
end
if bdio_get_ruinfo(fid)~=0
   p.N=-1;
   fprintf('Data file does not start with a lattice record\n');
   return;
end
if bdio_get_rfmt(fid)~=BDIO_ASC_GENERIC
   p.N=-1;
   fprintf('Wrong record type');
   return;
end
str=char(bdio_read(bdio_get_rlen(fid), fid)');
if length(str)~=bdio_get_rlen(fid)
   p.N=-1;
   fprintf('Could not read from data file\n');
   return;
end

p.ensemble=find_val(str,'ENSEMBLE=');
p.L(1) = find_val_num(str,'L0=');
p.L(2) = find_val_num(str,'L1=');
p.L(3) = find_val_num(str,'L2=');
p.L(4) = find_val_num(str,'L3=');
p.bc{1}= find_val(str,'BC0=');
p.bc{2}= find_val(str,'BC1=');
p.bc{3}= find_val(str,'BC2=');
p.bc{4}= find_val(str,'BC3=');

% read the correlator record
j=bdio_seek_record(fid);
if j~=0
   return;
end
if bdio_get_ruinfo(fid)~=1
   return;
end
% found a correlator description record. read it
if bdio_get_rfmt(fid)~=BDIO_ASC_GENERIC
   p.N=-1;
   fprintf('Correlator description record has wrong format\n');
   return;
end
str=char(bdio_read(bdio_get_rlen(fid), fid)');
if length(str)~=bdio_get_rlen(fid)
   p.N=-1;
   fprintf('Could not read from data file\n');
   return;
end
corrID = find_val_num(str,'CORR_ID=');
if corrID ~=0 
   p.N=-1;
   printf('Unexpected corr_ID\n');
   return;
end
   
if find_val_num(str,'NDIM=') ~= 1
   p.N=-1;
   printf('Unexpected NDIM\n');
   return;
end
nmu = find_val_num(str,'D0=');
if nmu~=4
   p.N=-1;
   fprintf('correlator has wrong dimension\n');
   return;
end
if ~strcmp(find_val(str,'DATATYPE='),'complex')
   p.N=-1;
   fprintf('Data should be complex\n');
   return;
end

p.N=0;
ic=0;
j=bdio_seek_record(fid);
while j==0
   % read next configuration
   if bdio_get_rfmt(fid)~=BDIO_ASC_GENERIC || bdio_get_ruinfo(fid)~=4
       return;
   end
   str=char(bdio_read(bdio_get_rlen(fid), fid)');
   if length(str) ~= bdio_get_rlen(fid)
      fprintf('could not read past cnfg entry no %d\n',p.N);
      return;
   end
   ic = ic+1;
   p.N=ic;
   p.nc(ic) = find_val_num(str,'NC=');
   j=bdio_seek_record(fid);
   if bdio_get_rfmt(fid)~=BDIO_BIN_F64LE || bdio_get_ruinfo(fid)~=5
      p.N = p.N-1;
      p.nc= p.nc(1:end-1);
      p.C = p.C(1:p.N,:);
      fprintf('Unexpected record encountered. Last valid cnfg entry: %d\n',p.N);
      return;
   end
   dat = bdio_read_f64(bdio_get_rlen(fid),fid);
   if length(dat)*8 ~= bdio_get_rlen(fid)
      p.N = p.N-1;
      p.nc= p.nc(1:end-1);
      for ii=1:ncorr
         p.C = p.C(1:p.N,:);
      end
      fprintf('Incomplete record. Last valid cnfg entry: %d\n',p.N);
      return;
   end
   if dat(1)~=p.nc(end) || dat(2) ~= 0
      p.N = p.N-1;
      p.nc= p.nc(1:end-1,:);
      p.C = p.C(1:p.N,:);
      fprintf('Unexpected record order. Last valid cnfg: %d\n',p.N);
      return;
   end
   dat = dat(3:end);
   p.C(ic,:) = dat(1:2:end-1)+1i*dat(2:2:end);

   j=bdio_seek_record(fid);
end
bdio_close(fid);
end

function valstr = find_val(str, tag)
   % scan string for the occurence of tag, return pointer to substring
   % between tag and the next '\0'
   t = regexp(str,sprintf('%s([^\\x0]*)',tag),'tokens','once');
   if length(t) ~= 1
      valstr='';
   else
      valstr=t{1};
   end
   return;
end

function val = find_val_num(str, tag)
   % scan string for the occurence of tag, return substring
   % between tag and the next '\0'
   t = regexp(str,sprintf('%s([^\\x0]*)',tag),'tokens','once');
   if length(t) ~= 1
      val=-1;
   else
      val=str2num(t{1});
   end
   return;
end
