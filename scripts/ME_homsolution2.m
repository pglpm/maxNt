(* Homog maxent *)

Clear[MEhomsolution2,MEhomsolutiontruncated,MEhomprob]; 

MEhomsolution2::usage = "Input: pop-avgd activity [0,1], pop-avgd couple activity, pop size, machine precision. Output: Lagr multipl, Lagr multipl, distribution, Z"
(* it uses absolute values in the optimization; usually gives much less precise results than MEhomsolutio2, but it doesn't complain about difficulties in finding the optimum *)
MEhomsolution2[meanmm_, meancorr_, nn_, wp_:MachinePrecision] := 
  Block[{totcorr = meancorr*(nn^2-nn)/2,(* it uses absolute values in the optimization *)
   totmm = meanmm*nn,totalprob,
   zzzm, js2, js1, freeenm, varm, soluzm, prob,ran=Range[0,nn],bvalues},
  zzzm = Sum[Binomial[nn, m]*Exp[js2*(m^2-m)/2 + js1*m], {m, 0, nn}]/2^nn;
  freeenm = Log[zzzm] - js2*totcorr - js1*totmm;
  soluzm = 
   NMinimize[freeenm, {js1, js2}, MaxIterations -> 10000, 
    PrecisionGoal -> 6,Method->Automatic, AccuracyGoal -> Infinity,WorkingPrecision->wp];
  prob = Table[
    Binomial[nn, m]*Exp[js2*(m^2-m)/2 + js1*m]/2^nn /. soluzm[[2]], {m, 
     0, nn}];
  totalprob = Total[prob];prob=prob/totalprob;
  bvalues={(ran.prob/nn)/meanmm-1,((ran^2-ran).prob/(nn^2-nn))/meancorr-1};
 Print["rel. discrepancies %: ", bvalues*100//N];
  {js1*nn /. soluzm[[2]], js2*(nn^2-nn)/2 /. soluzm[[2]], prob,totalprob}]

 MEhomsolutiontruncated[meanmm_, meancorr_, nn_, trunc_, wp_:MachinePrecision] := 
 Block[{totcorr = meancorr,
   totmm = meanmm,totalprob,
   zzzm, js2, js1, freeenm, varm, soluzm, prob,ran=Range[0,nn],bvalues},
  zzzm = Sum[Binomial[nn, m]*Exp[js2*(m^2-m)/(nn^2-nn) + js1*m/nn], {m, 0, trunc}];
  freeenm = Log[zzzm] - js2*totcorr - js1*totmm;
  soluzm = 
   NMinimize[freeenm, {js1, js2}, MaxIterations -> 10000, 
    PrecisionGoal -> 6, AccuracyGoal -> Infinity,WorkingPrecision->wp];
  prob =PadRight[Table[
    Binomial[nn, m]*Exp[js2*(m^2-m)/(nn^2-nn) + js1*m/nn] /. soluzm[[2]], {m, 
     0, trunc}],nn+1];
 probtrue =Table[
    Binomial[nn, m]*Exp[js2*(m^2-m)/(nn^2-nn) + js1*m/nn] /. soluzm[[2]], {m, 
     0, nn}];
 totalprob = Total[prob];prob=prob/totalprob;probtrue=probtrue/Total[probtrue];
 bvalues={(ran.prob/nn)/meanmm-1,((ran^2-ran).prob/(nn^2-nn))/meancorr-1};
 Print["rel. discrepancies: ", bvalues//N];
  {js1/nn /. soluzm[[2]], 2*js2/(nn^2-nn) /. soluzm[[2]], prob,probtrue,totalprob}]

 MEhomprob[js1v_,js2v_,nn_] :=
  Block[{js1=js1v*nn, js2=js2v*(nn^2-nn)/2, prob},
  prob = Table[
    Binomial[nn, m]*Exp[js2*(m^2-m)/(nn^2-nn) + js1*m/nn], {m, 
     0, nn}];
  prob/Total[prob]];
