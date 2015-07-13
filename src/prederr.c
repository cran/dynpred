/*  prederr.c	*/
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

void prederr(
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
			double *tout, /* time points at which prediction error is evaluated */
			int    *sFUN, /* 1 for Brier, 2 for Kullback-Leibler */
			double *err, /* OUTPUT, containing prediction errors at tout */
			double *work /* WORKING */
			)
{
    int n, nsurv, ncens, nout, FUN;
    int i, iout, j, idxsurvt0, idxcenst0;
    double t0, y, score, derrt0, errt0, censt0, survt0;
/*
  	Rprintf("Entering prederr() ...\n");
  	R_FlushConsole();
*/
    n = *sn; nsurv = *snsurv; ncens = *sncens, nout = *snout, FUN = *sFUN;
/*
  	Rprintf("n = %d, nsurv = %d, ncens = %d, nout = %d, FUN = %d\n\n",
              n, nsurv, ncens, nout, FUN);
  	R_FlushConsole();
*/
    for (iout=0; iout<nout; iout++) {
      t0 = tout[iout];
/*
  	  Rprintf("\niout=%d: t0=%6.4f\n",iout,t0);
      R_FlushConsole();
*/
      errt0 = 0.0;
      /* find indices for censoring and survival at t0 */
      idxsurvt0=0; j=1; while ((j<nsurv) && (tsurv[j]<=t0)) {j++; idxsurvt0++;}
      idxcenst0=0; j=1; while ((j<ncens) && (tcens[j]<=t0)) {j++; idxcenst0++;}
/*
      Rprintf("\tidxsurvt0=%d, idxcenst0=%d\n",idxsurvt0,idxcenst0);
   	  R_FlushConsole();
*/
      for (i=0; i<n; i++) {
        survt0 = survmat[i*nsurv+idxsurvt0];
/*
        Rprintf("\ti=%d, time[i]=%6.4f, status[i]=%d, survt0=%6.4f\n",i,time[i],status[i],survt0);
      	R_FlushConsole();
*/
        if ((time[i]<=t0) && (status[i]==1)) { /* a case */
          y = 0;
          censt0 = 1.0; /* censoring needed at time[i]-, called here censt0 nonetheless */
          j=1; while ((j<ncens) && (tcens[j]<=time[i])) {censt0 = censmat[i*ncens+(j-1)]; j++;}
          if (FUN==1) score = (y-survt0)*(y-survt0);
          else score = -(y*log(survt0)+(1-y)*log(1-survt0));
          derrt0 = score/censt0;
          errt0 += derrt0;
/*
          Rprintf("\t\tcase: C(ti)=%6.4f, derr=(%1.0f-%6.4f)*(%1.0f-%6.4f)/%6.4f=%6.4f/%6.4f=%6.4f, errt0=%6.4f\n",
          			censt0,y,survt0,y,survt0,censt0,score,censt0,derrt0,errt0);
       	  R_FlushConsole();
*/
        }
        else if (time[i]>t0) { /* a control (or not yet case) */
          y = 1;
          censt0 = censmat[i*ncens+idxcenst0];
          if (FUN==1) score = (y-survt0)*(y-survt0);
          else score = -(y*log(survt0)+(1-y)*log(1-survt0));
          derrt0 = score/censt0;
          errt0 += derrt0;
/*
          Rprintf("\t\tcontrol: C(ti)=%6.4f, derr=(%1.0f-%6.4f)*(%1.0f-%6.4f)/%6.4f=%6.4f/%6.4f=%6.4f, errt0=%6.4f\n",
          			censt0,y,survt0,y,survt0,censt0,score,censt0,derrt0,errt0);
       	  R_FlushConsole();
*/
        }
      }
      err[iout] = errt0/n;
/*
  	  Rprintf("Prediction error at %6.4f: %6.4f\n",t0,err[iout]);
      R_FlushConsole();
*/
    }
}
