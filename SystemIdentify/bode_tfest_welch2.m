function [sys, coh, f_w] = bode_tfest_welch2(r, u, y, Fs, fmin, fmax)
    assert (length(u) == length(y));
    assert (length(u) > 1);

    assert (length(u) > 100*Fs);

    % segments = floor(length(u) / Fs);
    % window_length = floor(length(u) / segments);

    window_length = 2*Fs;
    welch_win = hann(window_length);
    overlap = floor(window_length*49/50);

    Nfft = window_length;
    [puu_welch, f_w] = pwelch(u, welch_win, overlap, Nfft, Fs);
    [pyy_welch, ~] = pwelch(y, welch_win, overlap, Nfft, Fs);
    [prr_welch, ~] = pwelch(r, welch_win, overlap, Nfft, Fs);
    [puy_welch, ~] = cpsd(y, u, welch_win, overlap, Nfft, Fs);  % 互功率谱
    [pry_welch, ~] = cpsd(y, r, welch_win, overlap, Nfft, Fs);  % 互功率谱
    [pru_welch, ~] = cpsd(u, r, welch_win, overlap, Nfft, Fs);  % 互功率谱

    sense = pry_welch ./ prr_welch;
    process = pru_welch ./ prr_welch;

    sys = sense ./ process;

    coh1 = abs(pru_welch).^2 ./ (prr_welch .* puu_welch);
    coh2 = abs(pry_welch).^2 ./ (prr_welch .* pyy_welch);
    coh3 = abs(puy_welch).^2 ./ (puu_welch .* pyy_welch);
    % coh = coh1 .* coh2;
    coh = coh2;
    coh = max(min(coh, 1), 0); % 限制在 [0,1] 范围内（数值误差修正）

    sys = sys((fmin <= f_w) & (f_w <= fmax));
    coh = coh((fmin <= f_w) & (f_w <= fmax));
    f_w = f_w((fmin <= f_w) & (f_w <= fmax));
end

