function PlotMirror(savedir,idx,phi_net_all_full,xpup,ypup,xi,CM_PD_full,IM_PD_full)

%% Show IF Wavefronts
mirrorposidx=[17 25 33 41 10 18 26 34 42 50 3 11 19 27 35 43 51 59 4 12 20 28 36 44 52 60 5 13 21 29 37 45 53 61 6 14 22 30 38 46 54 62 15 23 31 39 47 55 24 32 40 48 ];
paddedmatCMD=nan(9);
paddedmatslope=nan(9);
minx = min(xpup(:));
maxx = max(xpup(:));
miny = min(ypup(:));
maxy = max(ypup(:));

hFigIF=figure('WindowStyle','docked');
haxIF=tiledlayout(8,8);

hSlope=figure('WindowStyle','docked');
haxSlope=tiledlayout(5,4);

hCMD=figure('WindowStyle','docked');
haxCMD=tiledlayout(5,4);
Sx=size(phi_net_all_full,1);
Sy=Sx;


for zern=1:20

    mirrormatCMD=nan(8);
    mirrormatslope=nan(8);


    for a=1:52

        mirrormatCMD(mirrorposidx(a))=CM_PD_full(a,zern);


        mirrormatslope(mirrorposidx(a))=IM_PD_full(a,zern);


        if zern==1
            figure(hFigIF)
            nexttile(mirrorposidx(a))

            showWavefront = nan(Sx,Sy);
            wftemp=phi_net_all_full(:,:,a);
            showWavefront(idx) = wftemp(idx);
            pcolor(xi,xi,showWavefront);
            colormap turbo
            shading interp;
            axis square;
            axis off
            xlim([minx-5 maxx+5])
            ylim([miny-5 maxy+5])
        end


    end

    if zern==1
        figure(hFigIF)
        drawnow
        haxIF.TileSpacing='loose';
        haxIF.Padding='tight';
        haxIF.TileSpacing='none';
        box on
    end


    %% Show Slopes and Commands
    figure(hCMD)
    nexttile(zern)
    paddedmatCMD(1:8,1:8)=mirrormatCMD;
    pcolor(1:9,1:9,paddedmatCMD)
    colormap parula
    axis image
    axis off
    shading flat
    % clim(voltlims(zern,:))
    colorbar
    title(['Zernike ' num2str(zern)])
    drawnow


    figure(hSlope)
    nexttile(zern)
    paddedmatslope(1:8,1:8)=mirrormatslope;
    pcolor(1:9,1:9,paddedmatslope)
    colormap parula
    axis image
    axis off
    shading flat
    % clim(slopelims(zern+2,:))
    % clim([-1 1])
    colorbar
    title(['Zernike ' num2str(zern)])
    drawnow

end

haxCMD.TileSpacing='loose';
box on
title(haxCMD,'Command Voltage for Zernike Modes (pseudoinverted slopes) (V/µm)')



haxSlope.TileSpacing='loose';
box on
title(haxSlope,'Influence Slopes for Zernike Modes (µm/V)')

exportgraphics(hFigIF,[savedir '/InfluenceFunctionMirror.png'])
saveas(hFigIF,[savedir '/InfluenceFunctionMirror.fig'])
close(hFigIF)

exportgraphics(hSlope,[savedir '/InfluenceZernSlopesMirror.png'])
saveas(hSlope,[savedir '/InfluenceZernSlopesMirror.fig'])
close(hSlope)

exportgraphics(hCMD,[savedir '/ZernVoltageMirror.png'])
saveas(hCMD,[savedir '/ZernVoltageMirror.fig'])
close(hCMD)
end