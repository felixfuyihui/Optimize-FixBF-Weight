clear;
clc;
for angle = 0:20:360
    micloc = [-0.05, 0.05];
    Fs = 16000;
    Nf = 257;
    %angle = 180;
    gamma = 0;

    N     = 180;
    M     = 2;%length(micloc);
    d     = angle / 180 * pi;
    theta = linspace(0, 2*pi, N);
    lobe  = 20;

    wng_cst = 10 ^ (gamma / 20);
    pattern = zeros(N, Nf);
    weights = zeros(M, Nf);

    for f = 1: Nf
        fprintf('Processing %d/%d...\n', f, Nf);
        freq = (f - 1) / (2 * Nf - 2)  * Fs;

        % steer vector from all doa
        Gs = zeros(N, M);
        % target beampattern on current frequency
        b  = zeros(N, 1);

        for n = 1: N
            steer_vector = transfer(micloc, theta(n), freq);
            Gs(n, :) = steer_vector;
            if abs(theta(n) - d) < lobe / 180 * pi
                b(n) = 1;
            end
        end

        % expect doa
        % 1 x M
        in_transfer = transfer(micloc, d, freq);

        % x: M x 1 beam weights for current frequency
        cvx_begin quiet
            variable x(M) complex
            minimize(norm(Gs * conj(x) - b, 2))
            subject to
                in_transfer * conj(x) == 1;
                x' * x <= 1 / wng_cst;
        cvx_end

        weights(:, f) = x;
        pattern(:, f) = abs(Gs * conj(x));

    end

    imagesc(pattern);

    freq_ticks = linspace(0, Nf - 1, 9);
    doa_ticks  = linspace(0, N, 7);

    xticks(freq_ticks);
    xticklabels(freq_ticks / (Nf - 1) * 8);
    yticks(doa_ticks);
    yticklabels(doa_ticks * 2);

    xlabel('Frequency(kHz)');
    ylabel('DOA(\circ)');

    save(['2mic_20lobe/doa' num2str(angle)], 'weights');
    save(['2mic_20lobe/pattern' num2str(angle)], 'pattern');

end
function steer_vector = transfer(micloc, phi, freq) 
	delays = cos(phi) * micloc / 340;
	steer_vector = exp(-2j * pi * delays * freq);
end
% function steer_vector = transfer(micloc, phi, freq)
%     delays = [cos(phi), sin(phi), -cos(phi), -sin(phi)] * 0.05 / 340;
% 	steer_vector = exp(-2j * pi * delays * freq);
% end
% function steer_vector = transfer(micloc, phi, freq)
%     R = 0.0425;
%     delays = [0,-cos(phi),-cos(pi/3-phi),-cos(2*pi/3-phi),cos(phi),cos(pi/3-phi),cos(2*pi/3-phi)] * R / 340;
% 	steer_vector = exp(-2j * pi * delays * freq);
% end