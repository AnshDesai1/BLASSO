function x = signal_gen (N, d, theta, SNR, K)
    %function [x] = signal_gen (ne, d, theta, snr, Nsamp)
    % N = number of elements
    % d = element spacing (wavelengths)
    % theta = signal angle (degrees) (this can be a vector)
    % SNR = signal to noise ratio (dB) (SNR should also be a vector that
    %       corresponds to length(theta)
    % K = samples (snapshots) of signal
    % x = radiation pattern vector for K snapshots(ne x K) (complex voltage)

    s=[];
    for ii=1:length (theta) % For each source
       A(:,ii)=linear_dir_vec(N,d,theta(ii)); % Determine steering vector that corresponds to source DoA
       % sample realizations
       if ((length(theta)>1) &&(ii==1)) % coherent source
           %s1=ones(1,K) * 10 ^(SNR(ii)/20);   % coherent source
           s1=(randn(1,K)+1i*randn(1,K))/sqrt(2) * 10 ^(SNR(ii)/20);  %noise source
       else      
           %s1=ones(1,Nsamp) * 10 ^(snr(ii)/20);   % coherent source
           s1=(randn(1,K)+1i*randn(1,K))/sqrt(2) * 10 ^(SNR(ii)/20);  %noise source
       end
       s=[s;s1]; % Complex Signal Coefficient vector
    end
    for kk= 1:N  % thermal noise per element
       n(kk,:)=(randn(1,K)+1i*randn(1,K))/sqrt(2);
    end
    
    x=A*s + n;
end

