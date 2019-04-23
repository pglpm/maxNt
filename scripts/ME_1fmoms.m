(* Homog maxent with 4 moments - corresponds to constraining variances *)

Clear[ME1fmoms,ME1fmoms0]; 

ME1fmoms[eemoms_, nn_, wp_:MachinePrecision,pg_:4] := 
  Block[{ll=Length[eemoms],js=Table[Unique["js"],Length[eemoms]], freeenm, soluzm, prob, bvalues,fmoms,emoms},
	emoms=1-eemoms;
	fmoms[s_]=Table[1-Binomial[s,r]/Binomial[nn,r], {r,ll}];
	freeenm = Log@Sum[(1-ss+nn)/((nn+1)*(nn+2))*
			 Exp[js.(fmoms[ss]-emoms)],{ss, 0, nn}];

  soluzm = 
   NMinimize[freeenm, js, MaxIterations -> 10000, 
    PrecisionGoal -> pg, AccuracyGoal -> Infinity,WorkingPrecision->wp];

  prob = Table[
    (1-ss+nn)/((nn+1)*(nn+2))*
    Exp[(js /. soluzm[[2]]).fmoms[ss]] ,{ss, 0, nn}];
  totalprob = Total[prob];prob=prob/totalprob;
  
  bvalues=(Total[Table[fmoms[ss],{ss,0,nn}]*prob]/emoms-1);
  Print["rel. discrepancies %: ", bvalues*100//N];
  
 {js /. soluzm[[2]],prob,totalprob}]


ME1fmoms0[eemoms_, nn_, wp_:MachinePrecision,pg_:4] := 
  Block[{ll=Length[eemoms],js=Table[Unique["js"],Length[eemoms]], freeenm, soluzm, prob, bvalues,fmoms,emoms},
	emoms=1-eemoms;
	fmoms[s_]=Table[1-Binomial[s,r]/Binomial[nn,r], {r,ll}];
	freeenm = Log@Sum[Binomial[nn,ss]*
			 Exp[js.(fmoms[ss]-emoms)],{ss, 0, nn}]/2^nn;

  soluzm = 
   NMinimize[freeenm, js, MaxIterations -> 10000, 
    PrecisionGoal -> pg, AccuracyGoal -> Infinity,WorkingPrecision->wp];

  prob = Table[
    Binomial[nn,ss]*
    Exp[(js /. soluzm[[2]]).fmoms[ss]] ,{ss, 0, nn}]/2^nn;
  totalprob = Total[prob];prob=prob/totalprob;
  
  bvalues=(Total[Table[fmoms[ss],{ss,0,nn}]*prob]/emoms-1);
  Print["rel. discrepancies %: ", bvalues*100//N];
  
 {js /. soluzm[[2]],prob,totalprob}]
