%	template_scancon.m
%	template script for converting from r-sin(theta) data to x-y image
%	must load r-sin(theta) data: rsdata

% --> QUESTION l. <--
% Scan convert the r-sin(theta) buffer to produce a sector scan image.
% Use bilinear interpolation to compute the image values on the
% sector scan image grid.  Matlab's "interp2" function will help you
% do bilinear interpolation.

% compute values needed for interp2
x = linspace(-35,35,512);
z = linspace(0,70,512);
[XI,ZI] = meshgrid(x,z);
[theta_i,r_i]=meshgrid(sin_theta,r);
r_ii = sqrt(XI.^2+ZI.^2);
theta_ii = XI./r_ii;

% Create image w/ bilinear interpolation

im = interp2(theta_i, r_i, abs(rsdata), theta_ii, r_ii, 'bilinear');
tt = find(isnan(im));
im(tt) = zeros(size(tt));

im_i = interp2(theta_i, r_i, abs(rsdata_i), theta_ii, r_ii, 'bilinear');
tt_i = find(isnan(im_i));
im(tt_i) = zeros(size(tt_i));

im_j = interp2(theta_i, r_i, abs(rsdata_j), theta_ii, r_ii, 'bilinear');
tt_j = find(isnan(im_j));
im(tt_j) = zeros(size(tt_j));

im_k = interp2(theta_i, r_i, abs(rsdata_i), theta_ii, r_ii, 'bilinear');
tt_k = find(isnan(im_k));
im(tt_k) = zeros(size(tt_k));

% do similarly all r-sin(theta) buffers

% --> QUESTION l. <--
% Use two images on a logarithmic scale to answer this question:
% one on a 40dB scale, the other on a 20dB scale
figure(10); showimage3(im, 1, 20,70/512,70/512); axis('image')		% Display 20 dB scale image
figure(11); showimage3(im, 1, 40,70/512,70/512); axis('image')		% Display 40 dB scale image

figure(12); showimage3(im_i, 1, 20,70/512,70/512); axis('image')		% Display 20 dB scale image
figure(13); showimage3(im_i, 1, 40,70/512,70/512); axis('image')		% Display 40 dB scale image

figure(14); showimage3(im_j, 1, 20,70/512,70/512); axis('image')		% Display 20 dB scale image
figure(15); showimage3(im_j, 1, 40,70/512,70/512); axis('image')		% Display 40 dB scale image

figure(16); showimage3(im_k, 1, 20,70/512,70/512); axis('image')		% Display 20 dB scale image
figure(17); showimage3(im_k, 1, 40,70/512,70/512); axis('image')		% Display 40 dB scale image