%%ECE563 Image Processing Project 2
%Simulating the Effects of Diffraction Limited Optics Using an Impulse Invariant System
% Roshni Uppala


%% reading the file ( image )and  loading the parameters

cl % clearing th past data
load kett                   % loading the input image kett.m
kett=double(kett);          % Converting image file to double
[M,N]=size(kett);           % Getting the size information of the input image file
lambda= 0.5*10^-3;          % The wavelength in mm.
fnum=3;                     % focal number 
cutoff_freq=1/(lambda*fnum) %cut off frequency

%% Computing the continuous frequency response OTF and generating its high
%%resolution ( 25 x 256 or more )

u=linspace(-1000,1000,4*M);  % Indices in X axis 
v=linspace(-1000,1000,4*N);  % Indices in Y axix
[U,V]= meshgrid(u,v); % Grid the axis 

row=sqrt(U.^2+V.^2); % Computing value of row in cycles/mm
H_dif= (2/pi)*acos(row./cutoff_freq)-2*row./(pi*cutoff_freq).*sqrt(1-(row./cutoff_freq).^2); % Obtaining the diffraction OTF
H_dif(row>cutoff_freq)=0; % Beyond cutoff frequency assigning value zero to the diffraction OTF

%Respective pitches in x and y directions 
%Define deltax & deltay as a function of sampling rate(2*cut_off frequency)
deltax=1/(2*cutoff_freq);
deltay=1/(2*cutoff_freq);
% a and b are the physical dimensions of a single detector in FPA in mm. in
% the x and y dimensions respectively. Given that the detector dimensions
% correspond to 100% fill factor detectors . Therefore, the dimensions of a
% and b are equal to the respective pitches. 
a=deltax ;
b=deltay;
H_det= ((sin(pi*a*U))./(pi*a*U)).*((sin(pi*b*V))./(pi*b*V));%Detector integration OTF for rectangular detectors 

H_conti= (H_dif).*(H_det);% The continuous OTF with the two components

 max(H_conti(:))% Magnitude of thE CONTI FREQ OTF

%display

surfl(u,v,H_conti);
shading interp;
colormap(gray(256));
xlabel('u (cycles/mm)');
ylabel('v (cycles/mm)');
zlabel('|H_conti(u,v)|');
title( ' Continuous OTF ');
axis tight;

%% High resolution surface plot of the desired impulse invariant discrete
%%frequency response Hd( f1,f2)

u=linspace(-cutoff_freq,cutoff_freq,4*M);  % Indices in X axis 
v=linspace(-cutoff_freq,cutoff_freq,4*N);  % Indices in Y axix

[U,V]= meshgrid(u,v); % Grid the axis 

row=sqrt(U.^2+V.^2); % Computing value of row in cycles/mm
H_dif= (2/pi)*acos(row./cutoff_freq)-2*row./(pi*cutoff_freq).*sqrt(1-(row./cutoff_freq).^2); % Obtaining the diffraction OTF
H_dif(row>cutoff_freq)=0; % Beyond cutoff frequency assigning value zero to the diffraction OTF

%Respective pitches in x and y directions 
%Define deltax & deltay as a function of sampling rate(2*cut_off frequency)
deltax=1/(2*cutoff_freq);
deltay=1/(2*cutoff_freq);
% a and b are the physical dimensions of a single detector in FPA in mm. in
% the x and y dimensions respectively. Given that the detector dimensions
% correspond to 100% fill factor detectors . Therefore, the dimensions of a
% and b are equal to the respective pitches. 
a=deltax ;
b=deltay;
H_det= ((sin(pi*a*U))./(pi*a*U)).*((sin(pi*b*V))./(pi*b*V));%Detector integration OTF for rectangular detectors 

H_conti= (H_dif).*(H_det);% The continuous OTF with the two components

%display

f1=u.*deltax;
f2=v.*deltay;
[F1,F2]=meshgrid(f1,f2);% Getting the grid discrete axis

u=F1/deltax; v=F2/deltay; % u and v in cycles/samples

%display
figure;
surfl(f1,f2,H_conti);
shading interp;
colormap(gray(256));
% Plotting from -1/2 to 1/2 cycles/ sample in each dimension
xlabel('f1 (cyc/sample)');
ylabel('f2 (cyc/sample)');
zlabel('|H_dsft(f1,f2)|');
title('Impulse invariant discrete frequency response');


%% Executing the function modal_psf() and plotting the designed 15 x 15 discrete impulse response which is also zero centered.
%Calling function
h_zeroc=modal_psf(fnum,lambda,a,b,deltax,deltay,M,N); % Obtaining the zero centered discrete impulse reponse from the function modal_psf()

[M,N]=size(h_zeroc); %get the size info
m=[-(M-1)/2:(M-1)/2];
n=[-(N-1)/2:(N-1)/2];
figure;
image2=mesh(m,n,h_zeroc);
axis tight
xlabel('m');
ylabel('n');
zlabel('h(m,n)');
title('Impulse response of 15x15 filter(zero centered)');


%% High resolution DFT of 15 X 15 discrete filter displayed as if it were a
% DSFT( using surfl with units of cycles /sample). also showing the magnitude
% surface.




H_dft_zeroc=fft2(h_zeroc,512,512); % computing the  DFT of the impulse response
[P Q]=size(H_dft_zeroc); %Getting the size information
w1=linspace(-pi,pi,P);  %Indices in x direction
w2=linspace(-pi,pi,Q);  %Indices in y direction
%Display
figure;
surfl(w1,w2,fftshift(abs(H_dft_zeroc)))
shading interp
colormap(gray(256));
xlabel('\omega_1 (rads/sample)')
ylabel('\omega_2 (rads/sample)')
zlabel('|H_{DFT}(\omega_1,\omega_2)|');
title('High resolution DFT of the 15x15 filter');
axis tight


%% Input image and output image processed with the 15 X15 filter.
% Filter size K and L
K=15;
L=15;
% Using padding with valid to produce same sizes of the input and output
% image and using convolution for the input image with the zero centered
% impulse response
g=conv2(padarray(kett,[floor(K/2) floor(L/2)],'symmetric','both'),h_zeroc,'valid');
im(kett) %Display the input image
title('Input image')
im(g)  %Display the output-processedimage
title('Output image(processed image)')
