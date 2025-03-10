% U matrix from Euler angles phi1, PHI, phi2.
% The formalism follows the ID11-3DXRD specs
%
%   U = euler2u(phi1, PHI, phi2)
%
% INPUT: phi, PHI, and phi2 in radians
% OUTPUT [U11 U12 U13; U21 U22 U23; U31 U32 U33]
%
%  Henning Poulsen, Risoe 15/6 2002.

% Changed input angles to be in radians instead of degrees
% Henning Osholm S�rensen, Ris� National Laboratory, June 23, 2006.
%

function U = euler2u(phi1, PHI, phi2)

U11 =  cos(phi1)*cos(phi2)-sin(phi1)*sin(phi2)*cos(PHI);
U21 =  sin(phi1)*cos(phi2)+cos(phi1)*sin(phi2)*cos(PHI);
U31 =  sin(phi2)*sin(PHI);
U12 =  -cos(phi1)*sin(phi2)-sin(phi1)*cos(phi2)*cos(PHI);
U22 =  -sin(phi1)*sin(phi2)+cos(phi1)*cos(phi2)*cos(PHI);
U32 =  cos(phi2)*sin(PHI);
U13 =  sin(phi1)*sin(PHI);   
U23 =  -cos(phi1)*sin(PHI);
U33 =  cos(PHI);
U = [U11 U12 U13; U21 U22 U23; U31 U32 U33];
					
