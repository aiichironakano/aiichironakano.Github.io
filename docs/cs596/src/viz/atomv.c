/***********************************************************************
  Program atomv.c--ball representation of atoms.
  Required files
    atomv.h:   Include file
    md.conf:   MD configuration file containing atomic coordinates
***********************************************************************/
#include "atomv.h"
#include <stdio.h>
#include <math.h>
#include <OpenGL/gl.h>    /* Header file for the OpenGL library */
#include <OpenGL/glu.h>   /* Header file for the GLu library */
#include <GLUT/glut.h>    /* Header file for the GLut library */

GLuint sphereid;          /* display-list id of atom sphere geom */
GLuint atomsid;           /* display-list id of all atoms */
GLdouble fovy, aspect, near_clip, far_clip;  
                          /* parameters for gluPerspective() */
FILE *fp;                 /* pointer to open an MD-configuration file */

/* Function prototypes ************************************************/
void reshape(int, int);
void makeFastNiceSphere(GLuint, double);
void makeAtoms(void);
void makeCurframeGeom(void);
void drawScene(void);
void display(void);
void initView(float *, float *);
void readConf(void);

/**********************************************************************/
void reshape (int w, int h) {
/***********************************************************************
  Callback for glutReshapeFunc()
***********************************************************************/
  /* set the GL viewport to match the full size of the window */
  glViewport(0, 0, (GLsizei)w, (GLsizei)h);
  aspect = w/(float)h;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(fovy,aspect,near_clip,far_clip);
  glMatrixMode(GL_MODELVIEW);
}

/**********************************************************************/
void makeFastNiceSphere(GLuint listid, double radius) {
/***********************************************************************
Called once to generate and compile sphere geometry into the given
display list id.
***********************************************************************/
  int i,j;
  float lon,lat;
  float loninc,latinc;
  float x,y,z;

  loninc = 2*M_PI/nlon;
  latinc = M_PI/nlat;

  glNewList(listid,GL_COMPILE);

    /* South-pole triangular fan */
    glBegin(GL_TRIANGLE_FAN);
      glNormal3f(0,-1,0);
      glVertex3f(0,-radius,0);
      lon = 0;
      lat = -M_PI/2 + latinc;
      y = sin(lat);
      for (i=0; i<=nlon; i++) {
        x = cos(lon)*cos(lat);
        z = -sin(lon)*cos(lat);
        glNormal3f(x,y,z);
        glVertex3f(x*radius,y*radius,z*radius);
        lon += loninc;
      }
    glEnd();

    /* Quadrilateral stripes to cover the sphere */
    for (j=1; j<nlat-1; j++) {
      lon = 0;
      glBegin(GL_QUAD_STRIP);
        for (i=0; i<=nlon; i++) {
          x = cos(lon)*cos(lat);
          y = sin(lat);
          z = -sin(lon)*cos(lat);
          glNormal3f(x,y,z);
          glVertex3f(x*radius,y*radius,z*radius);
          x = cos(lon)*cos(lat+latinc);
          y = sin(lat+latinc);
          z = -sin(lon)*cos(lat+latinc);
          glNormal3f(x,y,z);
          glVertex3f(x*radius,y*radius,z*radius);
          lon += loninc;
        }
      glEnd();
      lat += latinc;
    }

    /* North-pole triangular fan */
    glBegin(GL_TRIANGLE_FAN);
      glNormal3f(0,1,0);
      glVertex3f(0,radius,0);
      y = sin(lat);
      lon = 0;
      for (i=0; i<=nlon; i++) {
        x = cos(lon)*cos(lat);
        z = -sin(lon)*cos(lat);
        glNormal3f(x,y,z);
        glVertex3f(x*radius,y*radius,z*radius);
        lon += loninc;
      }
    glEnd();

  glEndList();
}

/**********************************************************************/
void makeAtoms() {
/***********************************************************************
  Makes display-list of all atoms in the current frame using spheres.
***********************************************************************/
  int i;
  float rval,gval,bval;

  glNewList(atomsid, GL_COMPILE);
  rval = Ratom; gval = Gatom; bval = Batom;  /* RGB color of an atom */
  for (i=0; i < natoms; i++) {
    glPushMatrix();
    glTranslatef(atoms[i].crd[0],atoms[i].crd[1],atoms[i].crd[2]);
    glColor3f(rval,gval,bval);
    glCallList(sphereid);
    glPopMatrix();
  }
  glEndList();
}

/**********************************************************************/
void makeCurframeGeom() {
/***********************************************************************
  Reads the atoms information for the current time frame and makes the
  display-list of all the atoms' geometry.
***********************************************************************/
  makeAtoms();
}

/**********************************************************************/
void drawScene() {
/***********************************************************************
  Called by display() to draw the view of the current scene.
***********************************************************************/
  /* Define viewing transformation */
  gluLookAt(
    (GLdouble)eye[0],(GLdouble)eye[1],(GLdouble)eye[2],
    (GLdouble)center[0],(GLdouble)center[1],(GLdouble)center[2],
    (GLdouble)up[0],(GLdouble)up[1],(GLdouble)up[2]);
  glCallList(atomsid);
}

/**********************************************************************/
void display() {
/***********************************************************************
  Callback for glutDisplayFunc().  It clears the frame and depth 
  buffers and draws the atoms in the current frame.
***********************************************************************/
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  drawScene();
  glutSwapBuffers();
}

/**********************************************************************/
void initView (float *min_ext, float *max_ext) {
/***********************************************************************
  Initializes global viewing, lighting, and projection values.
***********************************************************************/
  GLfloat light_diffuse[]   = {1.0, 1.0, 1.0, 1.0};
  GLfloat light_position1[] = {0.5, 0.5, 1.0, 0.0};
  float dif_ext[3],dis;
  int i;

  /* Define normal light */
  glLightfv(GL_LIGHT0,GL_DIFFUSE,light_diffuse);
  glLightfv(GL_LIGHT0,GL_POSITION,light_position1);

  /* Enable a single OpenGL light */
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  /* Use depth buffering for hidden surface elimination */
  glEnable(GL_DEPTH_TEST);

  /* get diagonal and average distance of extent */
  for (i=0; i<3; i++) dif_ext[i] = max_ext[i] - min_ext[i];
  dis = 0.0;
  for (i=0; i<3; i++) dis += dif_ext[i]*dif_ext[i];
  dis = (float)sqrt((double)dis);

  /* set center in world space */
  for (i=0; i<3; i++) center[i] = min_ext[i] + dif_ext[i]/2.0;

  /* set initial eye & look at location in world space */
  eye[0] = center[0];
  eye[1] = center[1];
  eye[2] = center[2] + dis;
  up[0] = 0.0;
  up[1] = 1.0;
  up[2] = 0.0;

  /* set parameters for gluPerspective() */
  /* Near- & far clip-plane distances */
  near_clip = (GLdouble)( 0.5*(dis-0.5*dif_ext[2]) );
  far_clip  = (GLdouble)( 2.0*(dis+0.5*dif_ext[2]) );
  /* Field of view */
  fovy = (GLdouble)( 0.5*dif_ext[1]/(dis-0.5*dif_ext[2]) );
  fovy = (GLdouble)( 2*atan((double)fovy)/M_PI*180.0 );
  fovy = (GLdouble)(1.2*fovy);

  /* Enable the color material mode */
  glEnable(GL_COLOR_MATERIAL);
}

/**********************************************************************/
void readConf() {
/***********************************************************************
Read atomic coordinates from an MD-configuration file & allocates 
necessary arrays.
***********************************************************************/
  int l, j;

  /* Open an MD-configuration file */
  fp = fopen("md.conf","r");
  /* Read the # of atoms */
  fscanf(fp,"%d",&natoms);
  /* allocate atoms array */
  atoms = (AtomType *) malloc(sizeof(AtomType)*natoms);
  /* Maximum & minimum extent of system in angstroms */
  for (l=0; l<3; l++) fscanf(fp,"%f%f",&min_ext[l],&max_ext[l]);
  /* Atomic coordinates */
  for (j=0; j<natoms; j++)
    fscanf(fp,"%f %f %f",&(atoms[j].crd[0]),&(atoms[j].crd[1]),
                         &(atoms[j].crd[2]));
  fclose(fp);
}

/**********************************************************************/
int main(int argc, char **argv) {
/**********************************************************************/

  glutInit(&argc, argv);

  /* Read atomic coordinates from an MD-configuration file */
  readConf();

  /* Set up an window */
  /* Initialize display mode */
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
  /* Specify window size */
  glutInitWindowSize(winx, winy);
  /* Open window */
  glutCreateWindow("Lennard-Jones Atoms");

  /* Initialize view */
  initView(min_ext, max_ext);

  /* Set a glut callback functions */
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);

  /* generate an OpenGL display list for single sphere */
  sphereid = glGenLists(1);
  makeFastNiceSphere(sphereid,atom_radius);
  
  /* generate an OpenGL display list for the atoms' geometry */
  atomsid = glGenLists(1);
  /* make the geometry of the current frame's atoms */
  makeCurframeGeom();

  /* Start main display loop */
  glutMainLoop();
  
  return 0;
}
/**********************************************************************/
