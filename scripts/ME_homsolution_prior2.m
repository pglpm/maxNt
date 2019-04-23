(* Homog maxent *)

Clear[MEhomsolution,MEhomsolutiontruncated,MEhomprob]; 

MEhomsolution::usage = "Input: pop-avgd activity [0,1], pop-avgd couple activity, pop size, machine precision. Output: Lagr multipl, Lagr multipl, distribution, Z"
(* it uses normalized values in the optimization; usually gives more precise results than MEhomsolution2, although it reports difficulties in finding the optimum *)
MEhomsolution[meanmm_, meancorr_, nn_, wp_:MachinePrecision,pg_:6] := 
  Block[{totcorr = Rationalize[meancorr,0],
	 totmm = Rationalize[meanmm,0],totalprob,
   zzzm, js2, js1, freeenm, varm, soluzm, prob,ran=Range[0,nn],bvalues},
       zzzm = Sum[(1-m+nn)/((nn+1)*(nn+2))*Exp[js2*(m^2-m)/(nn^2-nn) + js1*m/nn], {m, 0, nn}];
  freeenm = Log[zzzm] - js2*totcorr - js1*totmm;
  soluzm = 
   NMinimize[freeenm, {js1, js2}, MaxIterations -> 10000, 
    PrecisionGoal -> pg,Method->Automatic, AccuracyGoal -> Infinity,WorkingPrecision->wp];
  prob = Table[
    (1-m+nn)/((nn+1)*(nn+2))*Exp[js2*(m^2-m)/(nn^2-nn) + js1*m/nn] /. soluzm[[2]], {m, 
     0, nn}];
  totalprob = Total[prob];prob=prob/totalprob;
  bvalues={(ran.prob/nn)/meanmm-1,((ran^2-ran).prob/(nn^2-nn))/meancorr-1};
 Print["rel. discrepancies %: ", bvalues*100//N];
  {js1 /. soluzm[[2]], js2/. soluzm[[2]], prob,totalprob}]

MEhomsolutiont[meanmm_, meancorr_, nn_, wp_:MachinePrecision] := 
 Block[{totcorr = meancorr,
   totmm = meanmm,totalprob,
   js2, js1, freeenm, varm, soluzm, prob,ran=Range[0,nn],bvalues},
       freeenm = Log@Sum[(1-m+nn)/((nn+1)*(nn+2))*Exp[js2*((m^2-m)/(nn^2-nn)-totcorr) + js1*(m/nn-totmm)], {m, 0, nn}];
  soluzm = 
   NMinimize[freeenm, {js1, js2}, MaxIterations -> 10000, 
    PrecisionGoal -> 3,Method->Automatic, AccuracyGoal -> Infinity,WorkingPrecision->wp];
  prob = Table[
    (1-m+nn)/((nn+1)*(nn+2))*Exp[js2*(m^2-m)/(nn^2-nn) + js1*m/nn] /. soluzm[[2]], {m, 
     0, nn}];
  totalprob = Total[prob];prob=prob/totalprob;
  bvalues={(ran.prob/nn)/meanmm-1,((ran^2-ran).prob/(nn^2-nn))/meancorr-1};
 Print["rel. discrepancies %: ", bvalues*100//N];
  {js1 /. soluzm[[2]], js2/. soluzm[[2]], prob,totalprob}]

(* This approximates the sum with an integral, but doesn't give good results *)
MEhomsolutioni[meanmm_, meancorr_, nn_, wp_:MachinePrecision] := 
 Block[{totcorr = meancorr,
   totmm = meanmm,totalprob,
   js2, js1, freeenm, varm, soluzm, prob,bvalues},
       freeenm[js1_,js2_] :=Log@NIntegrate[Exp[nn*(
	 -x*Log[x]-(1-x)*Log[1-x]+ js1*(x-totmm)+js2*(x^2-totcorr))], {x,0,1},
					   PrecisionGoal -> 6, AccuracyGoal -> Infinity,WorkingPrecision->wp];
  soluzm = 
   NMinimize[freeenm[js1,js2], {js1, js2}, MaxIterations -> 10000, 
	     PrecisionGoal -> 6,Method->Automatic, AccuracyGoal -> Infinity,WorkingPrecision->wp];
  
  {jss1,jss2}={js1,js2}/.soluzm[[2]];

  totalprob=NIntegrate[Exp[nn*(
	 -x*Log[x]-(1-x)*Log[1-x]+ jss1*x+jss2*x^2)], {x,0,1},
					   PrecisionGoal -> 3, AccuracyGoal -> Infinity,WorkingPrecision->wp];
  bvalues=NIntegrate[{x,x^2}*Exp[nn*(
	 -x*Log[x]-(1-x)*Log[1-x]+ jss1*x+jss2*x^2)]/totalprob, {x,0,1},
					   PrecisionGoal -> 3, AccuracyGoal -> Infinity,WorkingPrecision->wp]/{meanmm,meancorr}-1;

  Print["rel. discrepancies %: ", bvalues*100//N];
  {jss1, jss2,totalprob}]

 MEhomsolutiontruncated[meanmm_, meancorr_, nn_, trunc_, wp_:MachinePrecision] := 
 Block[{totcorr = meancorr,
   totmm = meanmm,totalprob,
   zzzm, js2, js1, freeenm, varm, soluzm, prob,ran=Range[0,nn],bvalues},
  zzzm = Sum[(1-m+nn)/((nn+1)*(nn+2))*Exp[js2*(m^2-m)/(nn^2-nn) + js1*m/nn], {m, 0, trunc}];
  freeenm = Log[zzzm] - js2*totcorr - js1*totmm;
  soluzm = 
   NMinimize[freeenm, {js1, js2}, MaxIterations -> 10000, 
    PrecisionGoal -> 6, AccuracyGoal -> Infinity,WorkingPrecision->wp];
  prob =PadRight[Table[
    (1-m+nn)/((nn+1)*(nn+2))*Exp[js2*(m^2-m)/(nn^2-nn) + js1*m/nn] /. soluzm[[2]], {m, 
     0, trunc}],nn+1];
 probtrue =Table[
    (1-m+nn)/((nn+1)*(nn+2))*Exp[js2*(m^2-m)/(nn^2-nn) + js1*m/nn] /. soluzm[[2]], {m, 
     0, nn}];
 totalprob = Total[prob];prob=prob/totalprob;probtrue=probtrue/Total[probtrue];
 bvalues={(ran.prob/nn)/meanmm-1,((ran^2-ran).prob/(nn^2-nn))/meancorr-1};
 Print["rel. discrepancies: ", bvalues//N];
  {js1/nn /. soluzm[[2]], 2*js2/(nn^2-nn) /. soluzm[[2]], prob,probtrue,totalprob}]

 MEhomprobn[js1v_,js2v_,nn_] :=
  Block[{js1=js1v*nn, js2=js2v*(nn^2-nn)/2, prob},
  prob = Table[
    (1-m+nn)/((nn+1)*(nn+2))*Exp[js2*(m^2-m)/(nn^2-nn) + js1*m/nn], {m, 
     0, nn}];
  prob/Total[prob]];

 MEhomprob[js1_,js2_,nn_] :=
  Block[{prob},
  prob = Table[
    (1-m+nn)/((nn+1)*(nn+2))*Exp[js2*(m^2-m)/(nn^2-nn) + js1*m/nn], {m, 
     0, nn}];
  prob/Total[prob]];
