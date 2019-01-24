(* Homog maxent with 4 moments - corresponds to constraining variances *)

Clear[MEhom4moments,MEhom4thmoment,MEhomvmoments,MEhom4momentsprob]; 

MEhom4moments[m1_, m2_, m3_, m4_, nn_, wp_:MachinePrecision] := 
 Block[{totalprob,
	zzzm, js1, js2, js3, js4, freeenm, varm, soluzm, prob,ran=Range[0,nn],
	nn2=Binomial[nn,2], nn3=Binomial[nn,3], nn4=Binomial[nn,4], 
	bvalues},
       freeenm = Log@Sum[Binomial[nn, m]*
			 Exp[js1*(m/nn-m1) + js2*(Binomial[m,2]/nn2-m2) +
			     js3*(Binomial[m,3]/nn3-m3) +
			     js4*(Binomial[m,4]/nn4-m4)],
       {m, 0, nn}];
  soluzm = 
   NMinimize[freeenm, {js1, js2, js3, js4}, MaxIterations -> 10000, 
    PrecisionGoal -> 4, AccuracyGoal -> Infinity,WorkingPrecision->wp];
  prob = Table[
	  Binomial[nn, m]*Exp[js1*m/nn + js2*Binomial[m,2]/nn2 +
		      js3*Binomial[m,3]/nn3 +
		      js4*Binomial[m,4]/nn4]/. soluzm[[2]],
	  {m, 0, nn}];
  totalprob = Total[prob];prob=prob/totalprob;
    bvalues={(ran.prob/nn)/m1-1,
	     ((Binomial[ran,2]/nn2).prob)/m2-1,
	     ((Binomial[ran,3]/nn3).prob)/m3-1,
	     ((Binomial[ran,4]/nn4).prob)/m4-1};
 Print["rel. discrepancies %: ", bvalues*100//N];
 {js1 /. soluzm[[2]], js2 /. soluzm[[2]],
  js3 /. soluzm[[2]], js4 /. soluzm[[2]],prob,totalprob}]

MEhom4thmoment[m1_, m2_, m4_, nn_, wp_:MachinePrecision] := 
 Block[{totcorr = meancorr,
   totmm = meanmm,totalprob,
	zzzm, js1, js2, js4, freeenm, varm, soluzm, prob,ran=Range[0,nn],
	nn2=nn*(nn-1), nn3=nn*(nn-1)*(nn-2),nn4=nn*(nn-1)*(nn-2)*(nn-3),
	bvalues},
       zzzm = Sum[Binomial[nn, m]*
		  Exp[js1*m/nn + js2*m*(m-1)/nn2 +
		      js4*m*(m-1)*(m-2)*(m-3)/nn4],
       {m, 0, nn}];
  freeenm = Log[zzzm] - js1*m1 - js2*m2 -js4*m4;
  soluzm = 
   NMinimize[freeenm, {js1, js2, js4}, MaxIterations -> 10000, 
    PrecisionGoal -> 6, AccuracyGoal -> Infinity,WorkingPrecision->wp];
  prob = Table[
	  Binomial[nn, m]*
		  Exp[js1*m/nn + js2*m*(m-1)/nn2 +
		      js4*m*(m-1)*(m-2)*(m-3)/nn4]/. soluzm[[2]],
	  {m, 0, nn}];
  totalprob = Total[prob];prob=prob/totalprob;
  bvalues={(ran.prob/nn)/m1-1,((ran*(ran-1)).prob/nn2)/m2-1,
	   ((ran*(ran-1)*(ran-2)*(ran-3)).prob/nn4)/m4-1};
 Print["rel. discrepancies: ", bvalues//N];
 {js1 /. soluzm[[2]], js2 /. soluzm[[2]],
  js4 /. soluzm[[2]],prob,totalprob}]

MEhomvmoments[m1_, m2_, mv_, nn_, wp_:MachinePrecision] := 
 Block[{totcorr = meancorr,
   totmm = meanmm,totalprob,
	zzzm, js1, js2, jsv, freeenm, varm, soluzm, prob,ran=Range[0,nn],
	nn2=Binomial[nn,2], bvalues},
       zzzm = Sum[Binomial[nn, m]*
		  Exp[js1*m/nn + js2*Binomial[m,2]/nn2 +
		      jsv*(Binomial[m,2]/nn2 - (Binomial[m,2]/nn2)^2)],
       {m, 0, nn}];
  freeenm = Log[zzzm] - js1*m1 - js2*m2 -jsv*mv;
  soluzm = 
   NMinimize[freeenm, {js1, js2, jsv}, MaxIterations -> 10000, 
    PrecisionGoal -> 6, AccuracyGoal -> Infinity,WorkingPrecision->wp];
  prob = Table[
	  Binomial[nn, m]*
		  Exp[js1*m/nn + js2*Binomial[m,2]/nn2 +
		      jsv*(Binomial[m,2]/nn2 - (Binomial[m,2]/nn2)^2)]/. soluzm[[2]],
	  {m, 0, nn}];
  totalprob = Total[prob];prob=prob/totalprob;
  bvalues={(ran.prob/nn)/m1-1,
	   (Binomial[ran,2]/nn2).prob/m2-1,
	   ((Binomial[ran,2]/nn2)-(Binomial[ran,2]/nn2)^2).prob/mv-1};
 Print["rel. discrepancies: ", bvalues//N];
 {js1 /. soluzm[[2]], js2 /. soluzm[[2]],
  jsv /. soluzm[[2]],prob,totalprob}]

 MEhom4momentsprob[js1_, js2_,js3_, js4_, nn_] :=
	Block[{nn2=nn*(nn-1), nn3=nn*(nn-1)*(nn-2),nn4=nn*(nn-1)*(nn-2)*(nn-3),
	       prob},
  prob = Table[
    Binomial[nn, m]*Exp[js1*m/nn + js2*m*(m-1)/nn2 +
		      js3*m*(m-1)*(m-2)/nn3 +
		      js4*m*(m-1)*(m-2)*(m-3)/nn4], {m, 0, nn}];
  prob/Total[prob]];
