#include <iostream>
#include <GL/freeglut.h>
#include <vector>
#include <chrono>
#include <string>

using std::cout;
using std::cin;
using std::vector;

//type definitions
struct cPoint
{
	double x, y, z;
};

struct cColor
{
	float r, g, b;
};

typedef vector<cPoint> ptArr;
typedef vector<ptArr> ptMat;

//global variables
ptMat oFuncData;
double zRot = 0;
int wx = 800; int wy = 600;
double sx1, sx2, sy1, sy2, minx, miny, minz, maxz;
int useFunc = 0;
std::chrono::steady_clock::time_point lastUpdate = std::chrono::steady_clock::now();

//function forward declarations
double normalizeTo(double val, double minVal, double maxVal, double normVal);
float getDeltaTime();
double func(double x, double y, int funcTp);
void calcFuncData(ptMat &fd, double xs, double xe, double ys, double ye, int n);
void cbReshape(int x, int y);
void setupProjection(double cx, double cy, double cz);
void DrawChartLine(ptArr pts, cColor clr);
void DrawLine(cPoint pt1, cPoint pt2, cColor clr);
void cbDisplay();
void cbIdle();
void prepGLUT(int argc, char **argv);
void quickMin(double xs, double xe, double ys, double ye, double prec, double &minx_r, double &miny_r, double &minz_r);
void coordDescend(double x1, double y1, double x2, double y2, double sprec, double error, double &minx_r, double &miny_r, double &minz_r);

//main
int main(int argc, char **argv)
{
	std::string ans;
	cout << "Manual setup (M) or auto (A)?\n";
	cin >> ans;
	
	double sp = 0, err = 0;
	int meshDensity = 100;

	if (ans == "M" || ans == "m")
	{
		cout << "Enter X1, Y1, X2, Y2\n";
		cin >> sx1 >> sy1 >> sx2 >> sy2;

		double b = 0;
		if (sx2 < sx1)
		{
			b = sx1;
			sx1 = sx2;
			sx2 = b;
		}
		if (sy2 < sy1)
		{
			b = sy1;
			sy1 = sy2;
			sy2 = b;
		}

		cout << "Enter serach precision (i.e. 1000) and error allowance (i.e. 0.001)\n";
		cin >> sp >> err;
		cout << "Enter mesh density (i.e. 100)\n";
		cin >> meshDensity;

		cout << "Enter function type (0/1/2)\n";
		cin >> useFunc;
	}
	else
	{
		sx1 = -10;
		sx2 = 10;
		sy1 = -10;
		sy2 = 10;
		sp = 1000;
		err = 0.001;
		meshDensity = 100;
		useFunc = 1;
	}


	coordDescend(sx1,sy1,sx2,sy2,sp,err,minx,miny,minz);
	printf("Found minimum: x = %f; y = %f; z = %f\n", minx, miny, minz);

	calcFuncData(oFuncData,sx1,sx2,sy1,sy2, meshDensity);

	prepGLUT(argc, argv);

	system("pause");
}

//function implementations

double normalizeTo(double val, double minVal, double maxVal, double normVal)
{
	return normVal * (val / (maxVal - minVal));
}

float getDeltaTime()
{
	auto now = std::chrono::steady_clock::now();
	float deltaTime = std::chrono::duration_cast<std::chrono::microseconds>(now - lastUpdate).count() / 1000000.0f;
	lastUpdate = now;
	return deltaTime;
}

double func(double x, double y, int funcTp)
{

	if (funcTp == 0)
		return pow(x, 2) - x * y + pow(y, 2) + 9 * x - 6 * y + 20;
	if (funcTp == 1)
		return abs(3 * pow(x, 3))*pow(y, 2) + 20 * x*y + 1000 * sin(0.5*(x + y));
	if (funcTp == 2)
		return -(1/(1+pow(x,2))+1/(1+pow(y,2)));
}

void calcFuncData(ptMat &fd, double xs, double xe, double ys, double ye, int n)
{
	fd.clear();

	double stepx = (xe - xs) / (double)n;
	double stepy = (xe - xs) / (double)n;

	double cx = xs;
	double cy = ys;
	double cz = 0;

	for (int i = 0; i < n; i++)
	{

		fd.push_back({});

		cx = xs + (double)i * stepx;

		for (int j = 0; j < n; j++)
		{
			cy = ys + (double)j * stepy;
			cz = func(cx, cy, useFunc);

			if (i == 0) maxz = cz;
			if (maxz < cz) maxz = cz;

			fd[i].push_back({ cx,cy,cz });
		}
	}
}

void cbReshape(int x, int y)
{
	wx = x;
	wy = y;
	glViewport(0, 0, wx, wy);
}

void setupProjection(double cx, double cy, double cz)
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(45, wx / wy, 0, 1000);

	gluLookAt(cx, cy, cz, 0, 0, 0, 0, 0, 1);
}

void DrawChartLine(ptArr pts, cColor clr)
{
	int l = pts.size();
	glColor3d(clr.r, clr.g, clr.b);

	glBegin(GL_LINE_STRIP);
	for (int i = 0; i < l; i++)
	{
		glVertex3d(pts[i].x, pts[i].y,
			normalizeTo(pts[i].z, minz, maxz, 50));
	}
	glEnd();
}


void DrawLine(cPoint pt1, cPoint pt2, cColor clr)
{
	glColor3d(clr.r, clr.g, clr.b);

	glBegin(GL_LINES);
	glVertex3d(pt1.x, pt1.y,
		normalizeTo(pt1.z, minz, maxz, 50));
	glVertex3d(pt2.x, pt2.y,
		normalizeTo(pt2.z, minz, maxz, 50));
	glEnd();
}

void cbDisplay()
{

	float dt = getDeltaTime();

	setupProjection(sx1 * 2, sy1 * 2, 100);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glRotated(zRot, 0, 0, 1);

	zRot += 25 * (double)dt;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glLineWidth(2.0f);

	for (int i = 0; i < oFuncData.size(); i++)
		DrawChartLine(oFuncData[i], { 1.0,1.0,1.0 });

	DrawLine({ minx - 1,miny,minz }, { minx + 1,miny,minz }, { 1,0,0 });
	DrawLine({ minx,miny - 1,minz }, { minx,miny + 1,minz }, { 1,0,0 });
	DrawLine({ minx,miny,minz - 1 }, { minx,miny,minz + 1 }, { 1,0,0 });

	glutSwapBuffers();
}

void cbIdle()
{
	cbDisplay();
}

void prepGLUT(int argc, char **argv)
{
	glutInit(&argc, argv); //initializing GLUT

	glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);

	glutInitWindowSize(wx, wy);

	glutCreateWindow("GLUT WINDOW");

	//OpenGL setup block
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

	//GLUT launch block
	glutIdleFunc(cbIdle);
	glutReshapeFunc(cbReshape);
	glutDisplayFunc(cbDisplay);

	glutMainLoop();
}

void quickMin(double xs, double xe, double ys, double ye, double prec, double &minx_r, double &miny_r, double &minz_r)
{
	//not quick, actually; better use some proper
	//minimize algorithm
	double stepx = abs(xe - xs) / prec;
	double stepy = abs(ye - ys) / prec;
	int nOp = 10; //should be swapped by one of the next
	if (abs(stepx) > 0) nOp = abs(xe - xs) / stepx;
	if (abs(stepy) > 0) nOp = abs(ye - ys) / stepy;

	double dir = 1;
	if (xe < xs) dir = -1;

	double cx = xs;
	double cy = ys;
	double minx = cx;
	double miny = cx;
	double cz = func(cx, cy, useFunc);
	double minz = cz;

	for (int i = 0; i < nOp; i++)
	{
		cx = xs + i * dir*stepx;
		cy = ys + i * dir*stepy;
		cz = func(cx, cy, useFunc);
		if (cz < minz) { minx = cx; miny = cy; minz = cz; }
	}

	minx_r = minx;
	miny_r = miny;
	minz_r = minz;
}

void coordDescend(double x1, double y1, double x2, double y2, double sprec, double error, double &minx_r, double &miny_r, double &minz_r)
{
	//specifically the coord descent method
	double xl = x1, yl = y1;
	double xprev = xl, yprev = yl;
	double dx = 0, dy = 0;

	double minx = xl, miny = yl, minz = 0;

	int iteration = 0;

	do
	{

		xprev = xl; yprev = yl;

		if (iteration == 0) { quickMin(xl, x2, yl, yl, sprec, minx, miny, minz); iteration = 1; xl = minx; }
		else if (iteration == 1) { quickMin(xl, xl, y1, y2, sprec, minx, miny, minz); iteration = 0; yl = miny; }

		dx = abs(xl - xprev);
		dy = abs(yl - yprev);

		printf("Iteration complete;\ndx = %f; dy = %f\n", dx, dy);
		printf("mx = %f; my = %f mz = %f\n", minx, miny, minz);

		if ((dx == 0) && (dy == 0))
		{
			cout << "Breaking descent\n";
			break;
		}

	} while ((dx > error) || (dy > error));

	minx_r = minx;
	miny_r = miny;
	minz_r = minz;

}