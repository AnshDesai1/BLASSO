function y = linear_dir_vec (N, d, theta, plotFlag)
% assumes electronic steering of zero degrees (boresight)
% INPUTS
% N = number of elements
% d = element spacing (wavelengths)
% theta = angle (degrees) (this can be a vector)
% plotFlag = boolean to determine whether to plot pattern

% OUTPUTS
% y = radiation pattern vector(ne x len(theta)) (complex voltage)

if nargin < 4 || ~exist('plotFlag','var') || isempty(plotFlag)
  plotFlag = 0;
end

theta = theta * pi /180.;                        % convert degrees to radians
theta=reshape(theta,1,length(theta));            % reshape vector to be a row vector
y=exp(1j.*(0:N-1)'.*2*pi*d*sin(theta))/sqrt(N);  % create steering vector

if (length(theta) > 2) && plotFlag          % automatically plot pattern if given more than one angle of evaluation
   eps=1e-6;                                % eliminate the log of zero
   theta=0;                                 %project beam at desired angle
   y1=exp(1j.*(0:N-1)'.*2*pi*d*sin(theta)); %boresight beam is zero degrees, (same as ones, or sum channel)
   yy=y1'*y*1/sqrt(N);                      %beam projected onto desired beam
   figure
   plot(theta*180/pi,20*log10(abs(yy+eps)),'LineWidth',2)
   grid on, zoom on
   title(['antenna pattern for linear array of ',num2str(N),...
         'elements, ',num2str(d),'\lambda spacing']);
   xlabel('degrees')
   ylabel('antenna response (dB)')
   axis([min(theta*180/pi) max(theta*180/pi) -60 1])
end


