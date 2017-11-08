function sigmas = VestBMS_sensoryNoise(model,srange,sigmazero,w)
%VESTBMS_SENSORYNOISE Compute standard deviation of sensory noise.
%
%   SIGMAS = VESTBMS_SENSORYNOISE(MODEL,SRANGE,SIGMAZERO,W) returns
%   the standard deviation SIGMAS of sensory noise model MODEL for stimuli
%   SRANGE and parameters SIGMAZERO and SZERO.
%
%   The first character of char array MODEL specifies the noise model:
%     'q'   quadratic noise model, additive sigmazero
%     'c'   cosine noise model, additive sigmazero
%     'a'   absolute-valued sinusoidal noise model, additive sigmazero
%     'Q'   quadratic noise model, multiplicative sigmazero
%     'C'   cosine noise model, multiplicative sigmazero
%     'A'   absolute-valued sinusoidal noise model, multiplicative sigmazero

if w == 0
    sigmas = sigmazero;
else
    switch model(1)
        case 'Q'; sigmas = sigmazero .* sqrt(1 + (w.*srange).^2);    
        case 'C'; sigmas = sigmazero .* sqrt(1 + 2*(90/pi)^2*(1 - cos(srange*pi/90)).*w.^2);
        case 'A'; sigmas = sigmazero .* (1 + (90/pi)*abs(sin(srange*pi/90)).*w);
        case 'D'; sigmas = sigmazero .* (1 + (180/pi)*abs(sin(srange*pi/180)).*w);
        case 'q'; sigmas = sqrt(sigmazero.^2 + (w.*srange).^2);    
        case 'c'; sigmas = sqrt(sigmazero.^2 + 2*(90/pi)^2*(1 - cos(srange*pi/90)).*w.^2);
        case 'a'; sigmas = sigmazero + (90/pi)*abs(sin(srange*pi/90)).*w;
    end
end

end