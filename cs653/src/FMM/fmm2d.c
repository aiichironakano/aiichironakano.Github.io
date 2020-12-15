/*------------------------------------------------------------------------------
Fast multipole method (FMM) algorithm in 2 dimensions.
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "fmm2d.h"

/*----------------------------------------------------------------------------*/
int main(int argc, char **argv) {
/*----------------------------------------------------------------------------*/
  int j;
  double diff,max_diff = 0.0,error;
  clock_t t1,t2,t3;

  initialize();  /* Randomly generate particle positions & charges */
  t1 = clock();
  mp_leaf();  /* Compute multipoles at the leaf cells */
  upward();  /* Upward pass to compute multipoles at all quadtree levels */
  downward();  /* Dowaward pass to compute local-expansion terms */
  nn_direct();  /* Evaluate potentials at all particle positions */
  t2 = clock();
  all_direct();  /* All-pairs direct evaluation of potentials for validation */
  t3 = clock();

  for (j=0; j<Npar; j++) {
    diff = (pot[j]-pot_direct[j])/pot_direct[j];
    diff = diff>0 ? diff : -diff;
    max_diff = diff>max_diff ? diff : max_diff;
  }
  error = (eng-eng_direct)/eng_direct;
  printf("===== Max potential difference = %e =====\n",max_diff);
  printf("===== Total FMM vs. direct energies & error = %e %e %e =====\n",
    eng,eng_direct,error);
  tfmm    = (t2-t1)/(double)CLOCKS_PER_SEC;
  tdirect = (t3-t2)/(double)CLOCKS_PER_SEC;
  printf("===== FMM & direct CPU times = %e %e =====\n",
    tfmm,tdirect);

  return 0;
}

/*----------------------------------------------------------------------------*/
void initialize() {
/*------------------------------------------------------------------------------
  Randomly generates particle positions in the range [0,BOX] & charges in [0,1];
  also initializes basic data structures.
------------------------------------------------------------------------------*/
  int j,b,l;

  /* Particle positions & charges */
  for (j=0; j<Npar; j++) {
    for (b=0; b<2; b++) z[j][b] = rand()/(double)RAND_MAX*BOX;
    q[j] = rand()/(double)RAND_MAX;
  }

  /* Cell-index offsets */
  for (l=0; l<=L; l++) c0[l] = (pow(4,l)-1)/3;
}

/*----------------------------------------------------------------------------*/
void mp_leaf() {
/*------------------------------------------------------------------------------
  Computes multipoles for all cells at the leaf level.
------------------------------------------------------------------------------*/
  int nc,lc,c,a,b,j,cj[2];
  double rc,zjc[2],qz[2];

  nc = pow(4,L);  /* # of quadtree cells at the leaf level */
  lc = pow(2,L);  /* # of cells in each direction */
  rc = BOX/lc;    /* Leaf-cell length */

  /* Clear multipoles */
  for (c=c0[L]; c<c0[L]+nc; c++)
    for (a=0; a<=P; a++)
      cini(0.0,0.0,phi[c][a]);

  /* Scan particles to add their multipoles */
  for (j=0; j<Npar; j++) {
    for(b=0; b<2; b++) cj[b] = z[j][b]/rc;  /* Particle-to-cell mapping */
    c = c0[L]+cj[0]*lc+cj[1];  /* Serial cell index */
    for (b=0; b<2; b++) zjc[b] = z[j][b]-(cj[b]+0.5)*rc;
    cini(q[j],0.0,qz);
    cadd(1.0,phi[c][0],1.0,qz,phi[c][0]);
    for (a=1; a<=P; a++) {
      cmul(qz,zjc,qz);
      cadd(1.0,phi[c][a],-1.0/a,qz,phi[c][a]);
    }
  }
}

/*----------------------------------------------------------------------------*/
void upward() {
/*------------------------------------------------------------------------------
  Upward pass to compute multipoles for all cells at all quadtree levels.
------------------------------------------------------------------------------*/
  int nc,lc,l,vc[2],c,vcd[2],cd,a,b,g;
  double rc,zdm[2],pz[2],zg[2],w[2];

  /* Ascend the quadtree upward */
  for (l=L-1; l>=0; l--) {
    nc = pow(4,l);  /* # of cells at quadtree level l */
    lc = pow(2,l);  /* # of cells in each direction */
    rc = BOX/lc;  /* Cell length */
    /* Loop over mother cells at level l */
    for (c=c0[l]; c<c0[l]+nc; c++) {
      for (a=0; a<=P; a++)
        cini(0.0,0.0,phi[c][a]);
      vc[0] = (c-c0[l])/lc; vc[1] = (c-c0[l])%lc; /* Mother's vector cell ID */
      /* Loop over 4 daughter cells at level l+1 */
      for (vcd[0]=2*vc[0]; vcd[0]<=2*vc[0]+1; (vcd[0])++)
      for (vcd[1]=2*vc[1]; vcd[1]<=2*vc[1]+1; (vcd[1])++) {
        cd = c0[l+1]+vcd[0]*(2*lc)+vcd[1]; /* Daughter's serial cell index */
        for (b=0; b<2; b++)
          zdm[b] = (vcd[b]+0.5)*(rc/2)-(vc[b]+0.5)*rc;  /* Zdaughter-Zmother */
        cadd(1.0,phi[c][0],1.0,phi[cd][0],phi[c][0]);
        smul(phi[cd][0],1.0,pz);
        for (a=1; a<=P; a++) {
          cmul(pz,zdm,pz);
          cadd(1.0,phi[c][a],-1.0/a,pz,phi[c][a]);
          for (g=0; g<=a-1; g++) {
            if (g == 0)
              cini(1.0,0.0,zg);
            else
              cmul(zg,zdm,zg);
            cmul(phi[cd][a-g],zg,w);
            cadd(1.0,phi[c][a],comb(a-1,a-g-1),w,phi[c][a]);
          }
        } /* End for multipole terms a */
      } /* End for daughter cells vcd */
    } /* End for mother cells c */
  } /* End for quadtree levels l */
}

/*----------------------------------------------------------------------------*/
void downward() {
/*------------------------------------------------------------------------------
  Downward pass to compute local-expansion terms for all cells at all levels.
------------------------------------------------------------------------------*/
  int nc,lc,c,a,b,l,vc[2],cm,vcm[2],g,vcb[2],vce[2],vci[2],ci;
  double rc,zdm[2],zg[2],w[2],zdi[2],lz[2],zi[2],zib[2],zia[2],zim[2],w0[2];

  nc = pow(4,1);  /* # of cells at quadtree level 1 */
  for (c=0; c<c0[1]+nc; c++)  /* Clear local expansion terms at levels 0 & 1 */
    for (a=0; a<=P; a++)
      cini(0.0,0.0,psi[c][a]);

  /* Decend the quadtree downward */
  for (l=2; l<=L; l++) {
    nc = pow(4,l);  /* # of cells at quadtree level l */
    lc = pow(2,l);  /* # of cells in each direction */
    rc = BOX/lc;    /* Cell length */

    /* Local-to-local transformation from the mother */
    for (c=c0[l]; c<c0[l]+nc; c++) {  /* Loop over daughter cells */
      vc[0] = (c-c0[l])/lc; vc[1] = (c-c0[l])%lc;  /* Daughter's vector cell ID */
      for (b=0; b<2; b++) vcm[b] = vc[b]/2;  /* Mother's vector cell ID */
      cm = c0[l-1]+vcm[0]*(lc/2)+vcm[1];  /* Mother's serial cell ID */
      for (b=0; b<2; b++)
        zdm[b] = (vc[b]+0.5)*rc-(vcm[b]+0.5)*(2*rc);  /* Zdaughter-Zmother */
      for (a=0; a<=P; a++) {
        cini(0.0,0.0,psi[c][a]);
        for (g=0; g<=P-a; g++) {
          if (g == 0)
            cini(1.0,0.0,zg);
          else
            cmul(zg,zdm,zg);
          cmul(psi[cm][a+g],zg,w);
          cadd(1.0,psi[c][a],comb(a+g,a),w,psi[c][a]);
        }
      }  /* End for local-expansion terms a */
    } /* End for daughter cells c */

    /* Multipole-to-local transfomation from the interactive cells */
    for (c=c0[l]; c<c0[l]+nc; c++) {  /* Loop over cells c*/
      vc[0] = (c-c0[l])/lc; vc[1] = (c-c0[l])%lc;  /* Vector cell ID */
      for (b=0; b<2; b++) {  /* Beginning & ending interactive-cell indices */
        vcb[b] = 2*(vc[b]/2-1);
        vcb[b] = vcb[b]>=0 ? vcb[b] : 0;
        vce[b] = 2*(vc[b]/2+1)+1;
        vce[b] = vce[b]<lc ? vce[b] : lc-1;
      }
      /* Loop over cells vci[] within mother's nearest neighbor */
      for (vci[0]=vcb[0]; vci[0]<=vce[0]; (vci[0])++)
      for (vci[1]=vcb[1]; vci[1]<=vce[1]; (vci[1])++) {
        /* Exclude daughter's nearest neighbors */
        if ( (vci[0]-vc[0]<-1)||(vci[0]-vc[0]>1)||
             (vci[1]-vc[1]<-1)||(vci[1]-vc[1]>1) ) {
          ci = c0[l]+vci[0]*lc+vci[1];  /* Serial ID of interactive cell */
          for (b=0; b<2; b++)
            zdi[b] = (vc[b]-vci[b])*rc;  /* Zdaughter-Zinteractive */
          clgn(zdi,lz);
          cmul(phi[ci][0],lz,w);
          cadd(1.0,psi[c][0],1.0,w,psi[c][0]);
          cinv(zdi,zi);
          cini(1.0,0.0,zib);
          for (b=1; b<=P; b++) {
            cmul(zib,zi,zib);
            cmul(phi[ci][b],zib,w);
            cadd(1.0,psi[c][0],1.0,w,psi[c][0]);
          }
          smul(zi,-1.0,zim);
          cini(1.0,0.0,zia);
          for (a=1; a<=P; a++) {
            cmul(zia,zim,zia);
            cmul(phi[ci][0],zia,w);
            cadd(1.0,psi[c][a],-1.0/a,w,psi[c][a]);
            cini(0.0,0.0,w0);
            cini(1.0,0.0,zib);
            for (b=1; b<=P; b++) {
              cmul(zib,zi,zib);
              cmul(phi[ci][b],zib,w);
              cadd(1.0,w0,comb(a+b-1,b-1),w,w0);
            }
            cmul(zia,w0,w0);
            cadd(1.0,psi[c][a],1.0,w0,psi[c][a]);
          } /* End for local-expansion terms a */
        } /* End if not daughter's nearest neighbor */
      } /* End for interactive cells vci[] */
    } /* End for cells c */

  } /* End for quadtree level l */
}

/*----------------------------------------------------------------------------*/
void nn_direct() {
/*------------------------------------------------------------------------------
  Direct calculation of the electrostatic potentials between the nearest-
  neighbor leaf cells, along with the evaluation of the local expansion.
------------------------------------------------------------------------------*/
  int nc,lc,j,k,a,b,vc[2],c,vc1[2],c1,vcb[2],vce[2];
  double rc,zjc[2],zjk[2],rjk,cpot[2],w[2],za[2];

  /* Make linked lists, lscl & evaluate local expansion-----------------------*/

  nc = pow(4,L);  /* # of quadtree cells at the leaf level */
  lc = pow(2,L);  /* # of cells in each direction */
  rc = BOX/lc;    /* Leaf-cell length */
  for (c=c0[L]; c<c0[L]+nc; c++) head[c] = EMPTY;  /* Reset the headers, head */

  /* Scan particles to construct headers, head, & linked lists, lscl */
  for (j=0; j<Npar; j++) {
    for (b=0; b<2; b++) vc[b] = z[j][b]/rc;  /* Vector cell index */
    c = c0[L]+vc[0]*lc+vc[1];  /* Vector-to-scalar cell-index mapping */
    lscl[j] = head[c];  /* Link to last occupant (or EMPTY if you're 1st) */
    head[c] = j;  /* The last one goes to the header */

    /* Evaluate the local expansion */
    for (b=0; b<2; b++) zjc[b] = z[j][b]-(vc[b]+0.5)*rc;
    cini(0.0,0.0,cpot);
    for (a=0; a<=P; a++) {
      if (a == 0)
        cini(1.0,0.0,za);
      else
        cmul(za,zjc,za);
      cmul(psi[c][a],za,w);
      cadd(1.0,cpot,1.0,w,cpot);
    }
    pot[j] = cpot[0];
  }

  /* Calculate pair interactions----------------------------------------------*/

  /* Scan inner cells */
  for (vc[0]=0; vc[0]<lc; (vc[0])++)
  for (vc[1]=0; vc[1]<lc; (vc[1])++) {
    c = c0[L]+vc[0]*lc+vc[1];  /* Scalar cell index */
    if (head[c] == EMPTY) continue;  /* Skip this cell if empty */

    /* Scan the neighbor cells (including itself) of cell c */
    for (b=0; b<2; b++) {
      vcb[b] = (vc[b]-1 >= 0) ? vc[b]-1 : 0;
      vce[b] = (vc[b]+1 < lc) ? vc[b]+1 : lc-1;
    }
    for (vc1[0]=vcb[0]; vc1[0]<=vce[0]; (vc1[0])++)
    for (vc1[1]=vcb[1]; vc1[1]<=vce[1]; (vc1[1])++) {
      c1 = c0[L]+vc1[0]*lc+vc1[1];  /* Scalar cell index of the neighbor cell */
      if (head[c1] == EMPTY) continue;  /* Skip this neighbor cell if empty */

      /* Scan atom i in cell c */
      j = head[c];
      while (j != EMPTY) {

        /* Scan atom j in cell c1 */
        k = head[c1];
        while (k != EMPTY) {

        /* Avoid double counting of pairs */
        if (j < k) {
          /* Pair vector zjk[] = z[j][]-z[k][] */
          for (rjk=0.0, b=0; b<2; b++) {
            zjk[b] = z[j][b]-z[k][b];
            rjk += zjk[b]*zjk[b];
          }
          rjk = sqrt(rjk);
          pot[j] += q[k]*log(rjk);
          pot[k] += q[j]*log(rjk);
        } /* End if i<j */

        k = lscl[k];
        } /* End while k not empty */

      j = lscl[j];
      } /* End while j not empty */

    } /* End for neighbor cells c1 */
  } /* End for central cells c */

  /* Compute the electrostatic energy */
  eng = 0.0;
  for (j=0; j<Npar; j++)
    eng += q[j]*pot[j];
  eng *= 0.5;
}

/*----------------------------------------------------------------------------*/
void all_direct() {
/*------------------------------------------------------------------------------
  All-pair direct calculation of the electrostatic potential.
------------------------------------------------------------------------------*/
  int j,k,a;
  double zjk[2],rjk,pot_jk;

  /* All-pair calculation of the electrostatic potentials */
  for (j=0; j<Npar; j++)
    pot_direct[j] = 0.0;
  for (j=0; j<Npar-1; j++)
    for (k=j+1; k<Npar; k++) {
      for (rjk=0.0, a=0; a<2; a++) {
        zjk[a] = z[j][a]-z[k][a];
        rjk += zjk[a]*zjk[a];
      }
      rjk = sqrt(rjk);
      pot_direct[j] += q[k]*log(rjk);
      pot_direct[k] += q[j]*log(rjk);
    }

  /* Compute the electrostatic energy */
  eng_direct = 0.0;
  for (j=0; j<Npar; j++)
    eng_direct += q[j]*pot_direct[j];
  eng_direct *= 0.5;
}
