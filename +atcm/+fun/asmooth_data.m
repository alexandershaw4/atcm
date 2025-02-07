function smoothed_data = asmooth_data(data, sigma)
    % smooth_data: Smooth a vector of data using a Gaussian filter with scale parameter sigma.
    %
    % Inputs:
    %   data  - Input data vector (1D array).
    %   sigma - Standard deviation of the Gaussian kernel (scale of smoothing).
    %
    % Output:
    %   smoothed_data - Smoothed data vector.

    % Ensure the input is a column vector
    data = data(:);
    
    % Define the kernel size based on sigma
    kernel_radius = ceil(3 * sigma); % 3 standard deviations as a rule of thumb
    x = -kernel_radius:kernel_radius;
    
    % Create a Gaussian kernel
    gaussian_kernel = exp(-x.^2 / (2 * sigma^2));
    gaussian_kernel = gaussian_kernel / sum(gaussian_kernel); % Normalize
    
    % Apply convolution to smooth the data
    smoothed_data = conv(data, gaussian_kernel, 'same'); % 'same' keeps the output size same as input
end
