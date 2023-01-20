/* To be put in a header file *************************************************/
#define NMAX 100           /* Max # of atoms */
#define NAUC 8             /* # of atoms per unit cell */
#define LCNS (1.0*5.43)    /* Lattice constant of Si in angstrom */
int nAtom;                 /* # of atoms */
double r[NMAX][3];         /* r[i][0|1|2] is the x|y|z coordinate of atom i */
int InitUcell[3];          /* # of unit cells */

/*----------------------------------------------------------------------------*/
void InitConf() {
/*------------------------------------------------------------------------------
	r are initialized to diamond lattice positions.  
------------------------------------------------------------------------------*/
	double gap[3];      /* Unit cell size */
	double c[3];
	int j,k,nX,nY,nZ;
	/* Atom positions in a unit diamond crystalline unit cell */
	double origAtom[NAUC][3] = {{0.0, 0.0, 0.0 }, {0.0, 0.5, 0.5 },
                                {0.5, 0.0, 0.5 }, {0.5, 0.5, 0.0 },
                                {0.25,0.25,0.25}, {0.25,0.75,0.75},
                                {0.75,0.25,0.75}, {0.75,0.75,0.25}};

	/* Read the # of unit cells in the x, y & z directions */
	scanf("%d%d%d",&InitUcell[0],&InitUcell[1],&InitUcell[2]);

	/* Sets up a diamond lattice */
	for (k=0; k<3; k++) gap[k] = LCNS;
	nAtom = 0;
	for (nZ=0; nZ<InitUcell[2]; nZ++) {
		c[2] = nZ*gap[2];
		for (nY=0; nY<InitUcell[1]; nY++) {
			c[1] = nY*gap[1];
			for (nX=0; nX<InitUcell[0]; nX++) {
				c[0] = nX*gap[0];
				for (j=0; j<NAUC; j++) {
					for (k=0; k<3; k++)
						r[nAtom][k] = c[k] + gap[k]*origAtom[j][k];
					++nAtom;
				}
			}
		}
	}
}

