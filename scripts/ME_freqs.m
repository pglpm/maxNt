(* Homog maxent with 4 moments - corresponds to constraining variances *)

Clear[MEfreqs]; 

MEfreqs[freqs_, n_,nn_, wp_:MachinePrecision,pg_:4] := 
  Block[{r=Length[freqs]-1,js=Table[Unique["js"],Length[freqs]], freeenm, varm, soluzm, prob, bvalues},
	
       freeenm = Log@Sum[(1-ss+nn)/((nn+1)*(nn+2))*
			 Exp[js.(Table[Binomial[n,s]*Binomial[nn-nn,ss-s]/Binomial[nn,ss],{s,0,r}]-freqs)]
       ,{ss, 0, nn}];

  soluzm = 
   NMinimize[freeenm, js, MaxIterations -> 10000, 
    PrecisionGoal -> pg, AccuracyGoal -> Infinity,WorkingPrecision->wp];

  prob = Table[
    (1-ss+nn)/((nn+1)*(nn+2))*
    Exp[js.(Table[Binomial[n,s]*Binomial[nn-nn,ss-s]/Binomial[nn,ss],{s,r}]-freqs)]
	 /. soluzm[[2]],{ss, 0, nn}];
  totalprob = Total[prob];prob=prob/totalprob;
  bvalues=(Table[
    Table[Binomial[n,s]*Binomial[nn-nn,ss-s]/Binomial[nn,ss],{ss,0,nn}].prob
	   ,{s,0,r}]/freqs-1);
 Print["rel. discrepancies %: ", bvalues*100//N];
 {js /. soluzm[[2]],prob,totalprob}]
