/*  prederrw.c	*/
/*
**  Input
**    n - no of subjects
**    time - time points from data
**    status - event indicators from data
**    nsurv - no of time points in tsurv
**    ncens - no of time points in tcens
**    nout - no of time points in tout
**    tsurv - time points for survival (first=0)
**    survmat - nsurv x n matrix with survival probabilities
**    tcens - time points for censoring (first=0)
**    censmat - ncens x n matrix with censoring probabilities
**    w - the window width
**    tout - time points at which prediction error is evaluated
**    FUN - 1 for Brier, 2 for Kullback-Leibler
**
** Output
**    err - vector of length nall with prediction errors
**
**  Work
**    work[n] - working vectors/matrices
*/
#include <math.h>
#include <R.h>

void prederrw(
			int    *sn, /* n, no of subjects */
			double *time, /* time points from data */
			int    *status, /* event indicators from data */
			int    *snsurv, /* nsurv, no of time points in tsurv */
			int    *sncens, /* ncens, no of time points in tcens */
			int    *snout, /* nall, no of time points in tout */
			double *tsurv, /* time points for survival */
			double *survmat, /* nsurv x n matrix with survival probabilities */
			double *tcens, /* time points for censoring */
			double *censmat, /* ncens x n matrix with censoring probabilities */
			double *w, /* the window width */
			double *tout, /* time points at which prediction error is evaluated */
			int    *sFUN, /* 1 for Brier, 2 for Kullback-Leibler */
			double *err, /* OUTPUT, containing prediction errors at tout */
			double *work /* WORKING */
			)
{
    int n, nsurv, ncens, nout, FUN;
    int i, iout, j, nrisk, idx, idxsurvt, idxcenst, idxsurvtw, idxcenstw;
    double t, tw, y, score, derrtw, errtw, censtw, survtw;
/*
  	Rprintf("Entering prederrw() ...\n");
  	R_FlushConsole();
*/
    n = *sn; nsurv = *snsurv; ncens = *sncens, nout = *snout, FUN = *sFUN;
/*
  	Rprintf("n = %d, nsurv = %d, ncens = %d, nout = %d, FUN = %d, w = %6.4f\n\n",
              n, nsurv, ncens, nout, FUN, *w);
  	R_FlushConsole();
*/
    for (iout=0; iout<nout; iout++) {
      t = tout[iout]; tw = t+*w;
      /* Find out how many at risk */
      nrisk = n;
      while ((time[nrisk-1]>t) && (nrisk>0)) nrisk--;
      nrisk = n - nrisk;
/*
  	  Rprintf("\niout=%d: t=%6.4f, tw=%6.4f, n=%d, nrisk=%d\n",iout,t,tw,n,nrisk);
      R_FlushConsole();
*/
      /* find indices for censoring and survival at t and tw */
      idxsurvt=0; j=1; while ((j<nsurv) && (tsurv[j]<=t)) {j++; idxsurvt++;}
      idxcenst=0; j=1; while ((j<ncens) && (tcens[j]<=t)) {j++; idxcenst++;}
      idxsurvtw=0; j=1; while ((j<nsurv) && (tsurv[j]<=tw)) {j++; idxsurvtw++;}
      idxcenstw=0; j=1; while ((j<ncens) && (tcens[j]<=tw)) {j++; idxcenstw++;}
/*
      Rprintf("\tidxsurvt=%d, idxcenst=%d, idxsurvtw=%d, idxcenstw=%d\n",idxsurvt,idxcenst,idxsurvtw,idxcenstw);
   	  R_FlushConsole();
*/
      errtw = 0.0;
      if (nrisk==0) err[iout] = errtw/nrisk;
      else {
		  for (i=0; i<nrisk; i++) {
			/* loop over individuals at risk at time t; i=0 is first individual at risk */
			idx = i+n-nrisk; /* idx is used to access "individual i" in data */
			survtw = survmat[idx*nsurv+idxsurvtw]/survmat[idx*nsurv+idxsurvt]; /* conditional survival S(t+w|t,xi) */
/*
			Rprintf("\tidx=%d, time[idx]=%6.4f, status[idx]=%d, survtw=%6.4f/%6.4f=%6.4f\n",
				idx,time[idx],status[idx],survmat[idx*nsurv+idxsurvtw],survmat[idx*nsurv+idxsurvt],survtw);
			R_FlushConsole();
*/
			if ((time[idx]<=tw) && (status[idx]==1)) { /* a case */
			  y = 0;
			  censtw = 1.0; /* censoring needed at ti-, called here censtw nonetheless */
			  j=1; while ((j<ncens) && (tcens[j]<=time[idx])) {censtw = censmat[idx*ncens+(j-1)]; j++;}
			  censtw = censtw/censmat[idx*ncens+idxcenst];
			  if (FUN==1) score = (y-survtw)*(y-survtw);
			  else score = -(y*log(survtw)+(1-y)*log(1-survtw));
			  derrtw = score/censtw;
			  errtw += derrtw;
/*
			  Rprintf("\t\tcase: C(ti)=%6.4f/%6.4f=%6.4f,\n\t\tderr=(%1.0f-%6.4f)*(%1.0f-%6.4f)/%6.4f=%6.4f/%6.4f=%6.4f, errtw=%6.4f\n",
						censtw*censmat[idx*ncens+idxcenst],censmat[idx*ncens+idxcenst],censtw,y,survtw,y,survtw,censtw,score,censtw,derrtw,errtw);
			  R_FlushConsole();
*/
			}
			else if (time[idx]>tw) { /* a control (or not yet case) */
			  y = 1;
			  censtw = censmat[idx*ncens+idxcenstw]/censmat[idx*ncens+idxcenst];
			  if (FUN==1) score = (y-survtw)*(y-survtw);
			  else score = -(y*log(survtw)+(1-y)*log(1-survtw));
			  derrtw = score/censtw;
			  errtw += derrtw;
/*
			  Rprintf("\t\tcontrol: C(ti)=%6.4f/%6.4f=%6.4f, derr=(%1.0f-%6.4f)*(%1.0f-%6.4f)/%6.4f=%6.4f/%6.4f=%6.4f, errtw=%6.4f\n",
						censmat[idx*ncens+idxcenstw],censmat[idx*ncens+idxcenst],censtw,y,survtw,y,survtw,censtw,score,censtw,derrtw,errtw);
			  R_FlushConsole();
*/
			}
		  }
		  err[iout] = errtw/nrisk;
	  }
/*
  	  Rprintf("Prediction error at %6.4f: %6.4f\n",t,err[iout]);
      R_FlushConsole();
*/
    }
}
