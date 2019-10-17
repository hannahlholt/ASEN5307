



scales = 1:120;
coefs = cwt(t_series.Den, scales, 'amor');
f = scale2frq(scales, 'amor', dt);
contour(t,f,abs(coefs), 'Linestyle', 'none', 'Linecolor', [0 0 0], 'Fill', 'on');
