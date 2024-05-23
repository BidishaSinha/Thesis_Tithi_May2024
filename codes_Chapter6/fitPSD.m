function [FfbrAc, FfbrC, FfbrRS, FfbrTen, FfbrV]=fitPSD(yallsB, xx, rf)
%%
w=1./(xx);
for i=1:rf
    yy=yallsB(:,i);
    
        [xData, yData, weights] = prepareCurveData( xx, yy, w );
        g=fittype('psdsimfunc6(AT,cytvis,ff,sig,x)');
        topts = fitoptions( g );
       % opts.Display = 'Off';
        topts.Lower = [1.0,1.0,0.0,0.0];%changed here
        topts.StartPoint = [1.0,20000.0,0.0,1.0];
        topts.Upper = [10.0,100000000.0,1.0,1000.0];
        topts.Weights = weights;
%         opt.Method='NonlinearLeastSquares';
        
        [f,g]=fit(xx, yy, g, topts);            
        FfbrAc(1,i)=f.AT;FfbrV(1,i)=0.01*(f.cytvis);FfbrC(1,i)=f.ff;FfbrTen(1,i)=50*(f.sig);
        FfbrRS(1,i)=g.rsquare;
end