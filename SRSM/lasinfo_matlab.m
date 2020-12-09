% Get lasinfo with directory + filename + arguments (see below for options)
function lasinfo_matlab(dirData, filename, args)
% Example: lasinfo_matlab("D:\workingData\LiDAR2017\", "17_4784_dc.las", "");

if nargin<3
    args="";
    if nargin<2
        filename="";
    end
end

current_dir=pwd;
cd C:\OSGeo4W64
[status,cmdout]=system(strcat("OSGeo4W.bat lasinfo ", dirData, filename, args),"-echo");
cd(current_dir)
end

% ----------------------------------------------------------
% lasinfo options:
%   -h [ --help ]         produce help message
%   -i [ --input ] arg    input LAS file
%   -v [ --verbose ]      Verbose message output
%   --no-vlrs             Don't show VLRs
%   --no-schema           Don't show schema
%   --no-check            Don't scan points
%   --xml                 Output as XML
%   -p [ --point ] arg    Display a point with a given id.  --point 44
%   --locale              Use the environment's locale for output
% 
% Filtering options:
%   -e [ --extent ] arg    Extent window that points must fall within to keep.
%                          Use a comma-separated or quoted, space-separated list,
%                          for example,
%                           -e minx, miny, maxx, maxy
%                           or
%                           -e minx, miny, minz, maxx, maxy, maxz
%                           -e "minx miny minz maxx maxy maxz"
%   --minx arg             Extent must be greater than or equal to minx to be
%                          kept.
%                           --minx 1234.0
%   --miny arg             Extent must be greater than or equal to miny to be
%                          kept.
%                           --miny 5678.0
%   --minz arg             Extent must be greater than or equal to minz to be
%                          kept. If maxx and maxy are set but not minz *and maxz,
%                          all z values are kept.
%                           --minz 0.0
%   --maxx arg             Extent must be less than or equal to maxx to be kept.
%                           --maxx 1234.0
%   --maxy arg             Extent must be less than or equal to maxy to be kept.
%                           --maxy 5678.0
%   --maxz arg             Extent must be less than or equal to maxz to be kept.
%                          If maxx and maxy are set but not maxz *and minz, all z
%                          values are kept.
%                           --maxz 10.0
%   -t [ --thin ] arg (=0) Simple decimation-style thinning.
%                          Thin the file by removing every t'th point from the
%                          file.
%   --last-return-only     Keep last returns (cannot be used with
%                          --first-return-only)
%   --first-return-only    Keep first returns (cannot be used with
%                          --last-return-only
%   --keep-returns arg     A list of return numbers to keep in the output file:
%                          --keep-returns 1 2 3
%   --drop-returns arg     Return numbers to drop.
%                          For example, --drop-returns 2 3 4 5
%   --valid_only           Keep only valid points
%   --keep-classes arg     A list of classifications to keep:
%                          --keep-classes 2 4 12
%                          --keep-classes 2
%   --drop-classes arg     A list of classifications to drop:
%                          --drop-classes 1 7 8
%                          --drop-classes 2
%   --keep-intensity arg   Range in which to keep intensity.
%                          The following expression types are supported:
%                          --keep-intensity 0-100
%                          --keep-intensity <200
%                          --keep-intensity >400
%                          --keep-intensity >=200
%   --drop-intensity arg   Range in which to drop intensity.
%                          The following expression types are supported:
%                          --drop-intensity <200
%                          --drop-intensity >400
%                          --drop-intensity >=200
%   --keep-time arg        Range in which to keep time.
%                          The following expression types are supported:
%                          --keep-time 413665.2336-414092.8462
%                          --keep-time <414094.8462
%                          --keep-time >413665.2336
%                          --keep-time >=413665.2336
%   --drop-time arg        Range in which to drop time.
%                          The following expression types are supported:
%                          --drop-time <413666.2336
%                          --drop-time >413665.2336
%                          --drop-time >=413665.2336
%   --keep-scan-angle arg  Range in which to keep scan angle.
%                          The following expression types are supported:
%                          --keep-scan-angle 0-100
%                          --keep-scan-angle <100
%                          --keep-scan-angle <=100
%   --drop-scan-angle arg  Range in which to drop scan angle.
%                          The following expression types are supported:
%                          --drop-scan-angle <30
%                          --drop-scan-angle >100
%                          --drop-scan-angle >=100
%   --keep-color arg       Range in which to keep colors.
%                          Define colors as two 3-tuples (R,G,B-R,G,B):
%                          --keep-color '0,0,0-125,125,125'
%   --drop-color arg       Range in which to drop colors.
%                          Define colors as two 3-tuples (R,G,B-R,G,B):
%                          --drop-color '255,255,255-65536,65536,65536'
% 
% For more information, see the full documentation for lasinfo at:
%  http://liblas.org/utilities/lasinfo.html
% ----------------------------------------------------------