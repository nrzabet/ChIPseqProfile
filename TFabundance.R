
bcd.min=10000;
bcd.max=30000;

cad.min=(bcd.min*(900/730));
cad.max=(bcd.max*(900/730));

gt.min=(bcd.min*(1700/730));
gt.max=(bcd.max*(1700/730));

hb.min=(bcd.min*(820/730));
hb.max=(bcd.max*(820/730));

kr.min=(bcd.min*(2700/730));
kr.max=(bcd.max*(2700/730));


bcd.min.prop=2000/bcd.max;
bcd.max.prop=2000/bcd.min;

cad.min.prop=10000/cad.max;
cad.max.prop=10000/cad.min;

gt.min.prop=1000/gt.max;
gt.max.prop=1000/gt.min;

hb.min.prop=2000/hb.max;
hb.max.prop=2000/hb.min;

kr.min.prop=20000/kr.max;
kr.max.prop=20000/kr.min;


TFabundance=matrix(0, nrow=5, ncol=4);
rownames(TFabundance)=c("BCD", "CAD","GT","HB","KR");
colnames(TFabundance)=c("min abundance","max abundance","min specifically bound","max specifically bound");

TFabundance[1,1]=bcd.min;TFabundance[1,2]=bcd.max;TFabundance[1,3]=bcd.min.prop;TFabundance[1,4]=bcd.max.prop;
TFabundance[2,1]=cad.min;TFabundance[2,2]=cad.max;TFabundance[2,3]=cad.min.prop;TFabundance[2,4]=cad.max.prop;
TFabundance[3,1]=gt.min;TFabundance[3,2]=gt.max;TFabundance[3,3]=gt.min.prop;TFabundance[3,4]=gt.max.prop;
TFabundance[4,1]=hb.min;TFabundance[4,2]=hb.max;TFabundance[4,3]=hb.min.prop;TFabundance[4,4]=hb.max.prop;
TFabundance[5,1]=kr.min;TFabundance[5,2]=kr.max;TFabundance[5,3]=kr.min.prop;TFabundance[5,4]=kr.max.prop;

signif(TFabundance,2)

cbind(paste(signif(TFabundance[,1],2)," - ",signif(TFabundance[,2],2)," molecules",sep=""), 
paste(round(TFabundance[,3]*100),"% - ",round(TFabundance[,4]*100),"%",sep=""))


TFabundancedf=data.frame("min abundance"=format(signif(TFabundance[,1],2),digits=0,scientific=F),"max abundance"=format(signif(TFabundance[,2],2),digits=0,scientific=F), 
                         "min specifically bound"=paste(round(TFabundance[,3]*100),"\\%",sep=""), "max specifically bound"=paste(round(TFabundance[,4]*100),"\\%",sep="") )
colnames(TFabundancedf)=c("min abundance","max abundance","min specifically bound","max specifically bound");


tableCaption="\\emph{The TF abundance in the nucleus and the percentage of specifically bound TF}. In the first two columns, we list the number of molecules that are in the nucleus for the five TFs. For our estimates for Bicoid nuclear abundance, see the Discussion section in the main text. In \\cite{zamparo_2009}, the authors reanalysed the FlyEx data \\cite{pisarev_2009} and proposed a lower limit for the nuclear abundance of the five TFs, but the proposed values are underestimates of the real values. For the last four TFs (Caudal, Giant, Hunchback and Kruppel), we considered the nuclear abundances of the four TFs relative to Bicoid nuclear abundance, as estimated in \\cite{zamparo_2009} using the Poisson method, and then we multiplied these relative abundances with our estimates for the Bicoid nuclear abundance. In the last two columns we list the percentage of specifically bound TFs, based on the estimations of our method (Table 2 in the main text) and the values in the first two columns.";  
tableLatex=xtable(TFabundancedf, label ="tab:TFabundanceEstimates", caption =tableCaption,align="|l|r|r|r|r|");
print.xtable(tableLatex, sanitize.text.function = function(x) x,file="tables/TFabundanceEstimates.tex",
             hline.after=-1:nrow(TFabundancedf),floating.environment="table", include.rownames=TRUE)



