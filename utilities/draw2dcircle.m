function [h,out] = draw2dcircle(centre, rad, varargin)
%DRAW2DCIRCLE draw a 2d circle with center and radius
	opt.n = 50;

    [opt,arglist] = tb_optparse(opt, varargin);

    % compute points on circumference
	th = [0:opt.n-1]'/ opt.n*2*pi;
    x = rad*cos(th) + centre(1);
    y = rad*sin(th) + centre(2);

    % add extra row if z-coordinate is specified, but circle is always in xy plane
    if length(centre) > 2
        z = ones(size(x))*centre(3);
        p = [x y z]';
    else
        p = [x y]';
    end

    out = p;
    % else plot the circle
    p = [p p(:,1)]; % make it close
    if numrows(p) > 2
        h = plot3(p(1,:), p(2,:), p(3,:), arglist{:});
    else
        h = plot(p(1,:), p(2,:), arglist{:});
    end
end

