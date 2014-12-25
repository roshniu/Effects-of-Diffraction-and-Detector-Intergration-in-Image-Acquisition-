function psf=modal_psf(fnum,lambda,a,b,deltax,deltay,M,N)

%
% psf = model_psf( fnum, lambda, a, b, dx, dy, M, N ) 
%
% Impulse invariant PSF model for diffraction limited optical system with
% detector integration.
%
% REQUIRED INPUTS
%
% fnum - F-number of optics (focal length / aperture)
% lambda - Wavelength of light considered
% a - Detector active area width 
% b - Detector active area height
% dx - Horizontal detector pitch
% dy - Vertical detector pitch
% M - Output PSF width (must be odd)
% N - Output PSF height (must be odd)
%
% 
% OUTPUTS
%
% psf - Discrete impulse invariant PSF
%
% Author:Roshni Uppala
% Date: 4/4/2013
%
% COPYRIGHT © 2013 ------------. ALL RIGHTS RESERVED.
% This is UNPUBLISHED PROPRIETARY SOURCE CODE of -------------; the
% contents of this file may not be disclosed to third parties, copied or
% duplicated in any form, in whole or in part, without the prior written
% permission of --------------.

% size of the filter

K=15;
L=15;

K=15;
L=15;
cutoff_freq=1/(lambda*fnum); % Cut off frequency
f1=([1:K]-(floor(K/2)+1))/K;% the value of f1 must be |f1|<0.5
f2=([1:L]-(floor(L/2)+1))/L; %% the value of f2 must be |f2|<0.5
u=f1./deltax;
v=f1./deltay;
[U,V]=meshgrid(u,v); % grid the axis
row_2=sqrt(U.^2+V.^2);

%COMPUTING FREQUENCY RESPONSE
% Obtaining the diffraction OTF
H_dif= (2/pi)*acos(row_2./cutoff_freq)-2*row_2./(pi*cutoff_freq).*sqrt(1-(row_2./cutoff_freq).^2);
H_dif(row_2>cutoff_freq)=0;% give the value of zero beyond cut_off frequency

% a and b are the physical dimensions of a single detector in FPA in mm. in
% the x and y dimensions respectively. Given that the detector dimensions
% correspond to 100% fill factor detectors . Therefore, the dimensions of a
% and b are equal to the respective pitches. 
deltax=1/(2*cutoff_freq);
deltay=1/(2*cutoff_freq);
a=deltax;
b=deltay;
H_det= [(sin(pi*a.*u))/(pi*a.*u)].*[(sin(pi*b.*v))/(pi*b.*v)];%Detector integration OTF for rectangular detectors 

H_dsft= (H_dif)*(H_det); % Obtaining the Hdsft( f1, f2) substituting the u and v values 
H_dft=ifftshift(H_dsft);%ifftshift to put DFT in the first quadrant 

h2=ifft2(H_dft); %compute the ifft2 
h_zeroc=fftshift(h2);% use fftshift to shift impilse response to zero centered

psf=h_zeroc; %Discrete impulse invariant PSF zero centered as the output to this funtion