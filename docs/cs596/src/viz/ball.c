/* ball.c
     A trivial CAVE demo program, demonstrating the most basic CAVE library
     functions. This program just draws a red ball in the front of the
     CAVE. No interaction (outside of moving around), and nothing changes.
*/

#include <cave_ogl.h>
#include <GL/glu.h>

void init_gl(void),draw_ball(void);

main(int argc,char **argv)
{
/* Initialize the CAVE */
 CAVEConfigure(&argc,argv,NULL);
 CAVEInit();
/* Give the library a pointer to the GL initialization function */
 CAVEInitApplication(init_gl,0);
/* Give the library a pointer to the drawing function */
 CAVEDisplay(draw_ball,0);
/* Wait for the escape key to be hit */
 while (!CAVEgetbutton(CAVE_ESCKEY))
	sginap(10);    /* Nap so that this busy loop doesn't waste CPU time */
/* Clean up & exit */
 CAVEExit();
}


static GLUquadricObj *sphereObj;

/* init_gl - GL initialization function. This function will be called
  exactly once by each of the drawing processes, at the beginning of
  the next frame after the pointer to it is passed to CAVEInitApplication.
  It defines and binds the light and material data for the rendering,
  and creates a quadric object to use when drawing the sphere. */
void init_gl(void)
{
 float redMaterial[] = { 1, 0, 0, 1 };
 glEnable(GL_LIGHT0);
 glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, redMaterial);
 sphereObj = gluNewQuadric();
}

/* draw_ball - the display function. This function is called by the
  CAVE library in the rendering processes' display loop. It draws a
  ball 1 foot in radius, 4 feet off the floor, and 1 foot in front
  of the front wall (assuming a 10' CAVE). */
void draw_ball(void)
{
 glClearColor(0., 0., 0., 0.);
 glClear(GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT);
 glEnable(GL_LIGHTING);
 glPushMatrix();
 glTranslatef(0.0, 4.0, -4.0);
 gluSphere(sphereObj, 1.0, 8, 8);
 glPopMatrix();
 glDisable(GL_LIGHTING);
}
