[x1, Fs1] = wavread('mic1.wav'); 
% Change the path to your first wav file.
[x2, Fs2] = wavread('mic2.wav'); 
% Change the path to your second wav file.

m = size(x1,1);
n = 2;
A = randn(n, n); 
% A random matrix that linearly combines the sound  signals.
x = A*[x1';x2'];


% x is the input. If you already have a linear combination, like something from two microphones, feed it here.
c = cov(x');
sq = inv(sqrtm(c));
mx = mean(x, 2)';
xx = x - mx'*ones(1, size(x, 2));
xx = sq*xx;

w1 = randn(n, 1);
w1 = w1/norm(w1,2);
w0 = randn(n, 1);
w0 = w0/norm(w0, 2);

while abs(abs(w0'*w1)-1) > 0.001
    w0 = w1;
    w1 = xx*G(w1'*xx)'/m - mean(DG(w1'*xx), 2)*w1;
    w1 = w1/norm(w1, 2);
end

w2 = randn(n, 1);
w2 = w2/norm(w2,2);
w0 = randn(n, 1);
w0 = w0/norm(w0, 2);

while abs(abs(w0'*w2)-1) > 0.001
    w0 = w2;
    w2 = xx*G(w2'*xx)'/m - mean(DG(w2'*xx), 2)*w2;
    w2 = w2 - w2'*w1*w1;
    w2 = w2/norm(w2, 2);
end


w = [w1 w2];
s = w*x;
s1 = s(1,:);
s2 = s(2,:);

% Writes out the extracted sound signals into two different wav files.
wavwrite( s1', 'out1.wav' );
wavwrite( s2', 'out2.wav' );
