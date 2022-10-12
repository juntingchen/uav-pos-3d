function f = evaluavpos_sf(PosUE, PosUAV, PosBS, Maps, meterPerPixel, ...
                           map_x0, K, funs_sf, SFcor3, Xrange, Yrange, Zrange)
% All positions in 3D


los = IsLosK_discrete(PosUE, PosUAV, Maps, meterPerPixel, map_x0);
ks = round((1 - los) * (K - 1) + 1);
du2u = norm(PosUAV - PosUE);
du2b = norm(PosUAV - PosBS);

[~, i] = min(abs(Xrange - PosUAV(1)));
[~, j] = min(abs(Yrange - PosUAV(2)));
[~, k] = min(abs(Zrange - PosUAV(3)));
sf = SFcor3(i, j, k);

f = funs_sf{ks}(du2u, du2b, sf);
