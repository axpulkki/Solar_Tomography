function [G_T,G_P,G_R,G_tot] = G_Thomson(x_obs,y_obs,z_obs,x,y,z,u)
%G_THOMSON Computation of the Thomson coefficients for the LOS integration.
%   
% [G_T,G_P,G_R,G_tot] = G_Thomson(x_obs,y_obs,z_obs,x,y,z,u) ... (x,y,z)
% the location of the plasma packet.
%

% Notation from Howard and Tappin (Space Sci Rev (2009) 147: 31?54 DOI 10.1007/s11214-009-9542-5) used.
% Antti Pulkkinen, March 2017.

% Solar radius.
Rs = 695700e3; % m.
% Thomson cross section.
sigma_e = 7.95e-30; % m^2/sr.

% Half-angular dimater of the Sun seen at (x,y,z). https://en.wikipedia.org/wiki/Angular_diameter.
omega = asin(Rs./(2*sqrt(x.^2 + y.^2 + z.^2)));

% Angles between the observer-(x,y,x)-Sun lines.
r_obs = [x_obs ; y_obs ; z_obs]; chi = zeros(size(x));
for ii = 1:length(x),
    
    r_plasma_packet = [x(ii) ; y(ii) ; z(ii)];
    r_packet_to_obs =  r_plasma_packet - r_obs;
    
    % Use dot product to determine the angle.
    chi(ii) = acos( r_plasma_packet.'*r_packet_to_obs/(sqrt(r_plasma_packet.'*r_plasma_packet)*sqrt(r_packet_to_obs.'*r_packet_to_obs)) );
    
end;

% Van de Hulst coefficients.
A = cos(omega).*sin(omega).^2;
B = -(1/8)*( 1 - 3*sin(omega).^2 - (cos(omega).^2./sin(omega)).*(1 + 3*sin(omega).^2).*log((1 + sin(omega))./cos(omega)));
C = 4/3 - cos(omega) - (cos(omega).^3)/3;
D = (1/8)*( 5 + sin(omega).^2 - (cos(omega).^2./sin(omega)).*(5 - sin(omega).^2).*log((1 + sin(omega))./cos(omega)));

G_T = (pi*sigma_e/2)*( (1 - u)*C + u*D );
G_P = (pi*sigma_e/2)*sin(chi).^2.*( (1 - u)*A + u*B );
G_R = G_T - G_P;
G_tot = G_T + G_R;
