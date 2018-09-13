#include<bits/stdc++.h>
#include<GL/gl.h>
#include <GL/glut.h>
///#include<string>

using namespace std;

#define FOR(i, a, b) for(int i = a; i <= b; i++)
#define ROF(i, a, b) for(int i = a; i >= b; i--)
#define F first
#define S second
#define PB push_back
#define PDD pair<double, double>
#define PII pair<int, int>

#define sim_spd 1000
#define fst_spd 25
#define mn_spd 1

#define vec_rep_idx 1
#define scale_rotate_vec_idx 2
#define par_rep_idx 3
#define line_int_idx 4
#define conv_idx 5
#define cpop_idx 6
#define fpop_idx 7
#define maer_idx 8
#define back_idx 9


#include "imagehandler.h"

void Sprint( float x, float y, string st)   {///printing text
    //cout << "y : " << y << "\n";
    ///glColor3f(0.0,1.0,0.7);
    //glDisable(GL_LIGHTING);
    glRasterPos2f( x, y); // location to start printing text
    for(int i = 0; i < (int)st.size(); i++) // loop until i is greater then l
    {
       glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, st[i]);

    }
}

double print_doc(double start_x, double start_y, string ss, int ch_num, int row_st) {///return the last row position
    ///start_x determines the x position of first row
    ///ch_num means how many characters will occupie one row
    ///row_st determines from which row it should show
    int row_no = 0;
    FOR(i, 0, (int)ss.size() - 1)   {
        int j;
        for(j = i; (j < ss.size()) && ((ss[j] == '\n') || (j < i + ch_num)) ; j++)   {
            if(ss[j] == '\n')    {
                j++;
                break;
            }
            if(ss[j] == ' ')    {
                continue;
            }
            int k = j;
            bool val = 1;
            while((k < ss.size()) && (ss[k] != ' ')) {
                ///cout << "ayy\n";
                /***/
                if(ss[k] == '\n')    {
                    break;
                }
                if(k == i + ch_num) {
                    val = 0;
                    break;
                }
                k++;
            }
            if(val)    {
                if((k < ss.size()) && (ss[k] == '\n'))    {
                    j = k - 1;
                }
                else    {
                    j = k;
                }
            }
            else    {
                break;
            }
        }
        string tmpss = "";
        FOR(pp, i, (int)j - 1)   {
            if(ss[pp] == '\n')    {
                continue;
            }
            tmpss += ss[pp];
        }
        Sprint(start_x, start_y, tmpss);
        start_y -= 0.1;
        i = j - 1;
    }
    return 0.0;///return something for now
    ///cout << "hi\n";
}

double fx(double yo)    {
    return (2 * yo / glutGet(GLUT_WINDOW_WIDTH)) - 1;
}

double fy(double yo)    {
    return (2 * yo / glutGet(GLUT_WINDOW_HEIGHT)) - 1;
}


int spd = fst_spd;
void update(int value) {
    ///cout << "timer chole ki?\n";
	glutPostRedisplay(); //Tell GLUT that the display has changed

	//Tell GLUT to call update again in 25 milliseconds
	glutTimerFunc(spd, update, 0);
}

int ent_pressed = 0;
int page_no = 1;

bool take_inp = 0;
///this is used for taking mouse input
vector< pair<double, double> > vv;///this is mouse input
void mouse(int button, int state, int x, int y) {
    if(take_inp == 0)    {
        return;
    }
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) { // Pause/resume
        ///vv.push_back(make_pair(ballx, bally));
        int windowHeight = glutGet(GLUT_WINDOW_HEIGHT);
        int nx = x;
        int ny = windowHeight-y;
        if(fx(nx) >= -0.005)    {
            return;
        }
        vv.push_back(make_pair(fx(nx), fy(ny)));
    }
}

double col(double cc)   {///color decimal to float
    return cc / 255.0;
}

double dis(complex<double> a, complex<double> b)   {
    return abs(a - b);
}

void draw_entire_region() {///drawing half of the screen for instruction
    glColor3f(col(255), col(204), col(204));///green color
    glBegin(GL_POLYGON);
        glVertex3f(-1, 1, 0);
        glVertex3f(1, 1, 0);
        glVertex3f(1, -1, 0);
        glVertex3f(-1, -1, 0);
    glEnd();
    glColor3f(col(0), col(0), col(0));///black color
}

void draw_right_region() {///drawing half of the screen for instruction
    glColor3f(col(255), col(204), col(204));///green color
    glBegin(GL_POLYGON);
        glVertex3f(0, 1, 0);
        glVertex3f(1, 1, 0);
        glVertex3f(1, -1, 0);
        glVertex3f(0, -1, 0);
    glEnd();
    glColor3f(col(0), col(0), col(0));///black color
}

void draw_1st_page()  {
    /**
	glColor3f(0.0, 0.0, 0.0);///black color

    string ss =
    "   Geometry Simulator"
    "\nClick Enter To Continue"
    ;
    print_doc(-0.15, 0, ss, 70, 0);
    ///cout << "hi2\n";

    //glutSwapBuffers();
    */
    drawCover();
}

void draw_highlighter(double up_y, double low_y) {
    glColor3f(col(255), col(204), col(204));///green color
    glBegin(GL_POLYGON);
        glVertex3f(-1, up_y, 0);
        glVertex3f(1, up_y, 0);
        glVertex3f(1, low_y, 0);
        glVertex3f(-1, low_y, 0);
    glEnd();
    glColor3f(col(0), col(0), col(0));///black color
}

int p2_menu_no = 1;
int p2_menu_min = 1;
int p2_menu_max = 4;

void draw_2nd_page()  {
    drawBackground();
    string ss =
    "  Menu"
    "\n  About"
    "\n  Reference"
    "\n  Exit"
    ;
    ///cout << "menu_no : " << p2_menu_no << "\n";
    double st = 0;///this is start y value
    double var = st + 0.1;
    for(int i = 1; i < p2_menu_no; i++)   {
        var -= 0.1;
    }
    double ts = 0.03;///threshold value
    draw_highlighter(var - ts, var - 0.1 - ts);
    glColor3f(0.0, 0.0, 0.0);///black color
    print_doc(-0.1, st, ss, 70, 0);
    ///cout << "hi2\n";

    //glutSwapBuffers();
}

void draw_about_page()  {
    draw_entire_region();
    string ss =
    "This simulator gives a introduction to computational geometry and this whole simulator is based on vector geometry."
    "\nThis simulator is useful for anyone who would like to learn or teach basic computational geometry"
    "\n\n                                                      Click Enter To Go Back"
    ;
    print_doc(-0.6, 0.7, ss, 100, 0);
    ///cout << "hi2\n";
}

void draw_ref_page()  {
    draw_entire_region();
    string ss =
    "[1]. Competitive Programmer’s Handbook by Antti Laaksonen"
    "\n[2]. Programming contest(data structure and algoritm) by Md Mahbubul Hasan"
    "\n[3]. http://codeforces.com/blog/entry/46162"
    "\n[3]. http://www.ahinson.com/algorithms_general/Sections/Geometry/ParametricLineIntersection.pdf"
    "\n[4]. The Rotating Calipers: An Efficient, Multipurpose, Computational Tool by Godfried T. Toussaint"
    "\nCover credit : http://www.deepgreensea.net"
    "\nBackgroud credit : https://www.freepik.com/free-vector/bright-background-with-dots_1252906.htm"
    "\n\n                                                      Click Enter To Go Back"
    ;
    print_doc(-0.6, 0.7, ss, 100, 0);
    ///cout << "hi2\n";
}

int p3_menu_no = 1;
int p3_menu_min = 1;
int p3_menu_max = 9;

void draw_3rd_page()  {
    drawBackground();
    string ss =
    "  Vector Representation"
    "\n  Scaling and rotating a vector"
    "\n  Parametric representation of a line"
    "\n  Intersecting point of two parametric lines"
    "\n  Convex Hull"
    "\n  Closest Pair Of Points"
    "\n  Furthest Pair Of Points"
    "\n  Minimum area enclosing rectangle"
    "\n  Back"
    ;
    ///cout << "menu_no : " << p2_menu_no << "\n";
    double st_x = -0.2;///this is start x value
    double st_y = 0.5;///this is start y value
    double var = st_y + 0.1;
    for(int i = 1; i < p3_menu_no; i++)   {
        var -= 0.1;
    }
    double ts = 0.03;///threshold value
    draw_highlighter(var - ts, var - 0.1 - ts);
    glColor3f(0.0, 0.0, 0.0);///black color
    print_doc(st_x, st_y, ss, 70, 0);
    ///cout << "hi2\n";
}

void draw_axis() {
    glColor3f(0.0, 0.0, 0.0);///black color
    for(double i = -1; i < 0; i += 0.1)   {
        glBegin(GL_LINES);
            glVertex3f(i, -1, 0);
            glVertex3f(i, 1, 0);
        glEnd();
    }
    for(double j = -1; j < 1; j += 0.1)   {
        glBegin(GL_LINES);
            glVertex3f(-1, j, 0);
            glVertex3f(0, j, 0);
        glEnd();
    }
    glColor3f(1.0, 0.0, 0.0);///red color
    glBegin(GL_POINTS);
        glVertex3f(-0.5, 0, 0);
    glEnd();

}

///vector representation part start
void draw_vec_ins_1()  {
    string ss =
    "A vector can be represented with a complex number where a complex number is a number of the form x+yi, where i is the imaginary unit."
    "\nx + yi is rectengular form of a complex number representing a two-dimensional point (x,y) or a vector from the origin to a point (x,y)."
    "\n(r, theta) is polar form of a comples number where r is magnitude of the vector and theta is angle of that vector with x axis"
    "\n\n\n                                   Plot 1 point on white space"
    "\n                                                      Then"
    "\n                                        Click Enter To Continue"
    ;
    print_doc(0.01, 0.9, ss, 70, 0);

    ///cout << "hi2\n";
}

void draw_vec_ins_2(PDD yo)  {

    ostringstream ss1, ss2, r, theta;
    ss1 << yo.F;
    ss2 << yo.S;

    complex<double> cp_yo = {yo.F, yo.S};
    r << abs(cp_yo);
    theta << arg(cp_yo);

    string ss =
    "\n\nPoint is situated at (" + ss1.str() + ", " + ss2.str() + ")."
    "\nComplex number representation:"
    "\nRectangular form: " + ss1.str() + " + " + ss2.str() + "i"
    "\nPolar form: (" + r.str() + ", " + theta.str() + ")"
    "\n\n\n\n                                        Click Enter To Continue"
    ;
    print_doc(0.01, 0.9, ss, 70, 0);

    ///cout << "hi2\n";
}

void draw_vec_rep()  {
    draw_axis();
    draw_right_region();

    glColor3f(0.0, 0.0, 0.0);///black color
    if(vv.size())    {
        glBegin(GL_POINTS);
            glVertex3f(vv[0].F, vv[0].S, 0);
        glEnd();
	}

    if(ent_pressed == 0)    {
        draw_vec_ins_1();
	}
	else if(ent_pressed == 1)    {
        draw_vec_ins_2({(vv[0].F + 0.5) * 10, (vv[0].S) * 10});
        glBegin(GL_LINES);
            glVertex3f(-0.5, 0, 0);
            glVertex3f(vv[0].F, vv[0].S, 0);
        glEnd();
	}


}
///vector representation part end

double pnt_x_axis(double x)  {
    return (x / 10) - 0.5;
}
double pnt_y_axis(double y)  {
    return (y / 10);
}

///scaling and rotating vector part start
double scarot_r_init, scarot_r_inc;
double scarot_theta_init, scarot_theta_inc;
int scarot_div_final, scarot_div_init;

void init_scale_rotate_vec(double scale, double rot, int div) {
    scarot_r_init = 1;///(1, 0) that's why magnitude only in x axis
    scarot_theta_init = 0;

    scarot_r_inc = (scale - scarot_r_init) / div;
    scarot_theta_inc = rot / div;

    /**
    vector<double> yo = {scarot_r_init, scarot_theta_init};
    yo *= polar(scale, rot);
    scarot_r_final = yo.real();
    scarot_theta_final = yo.imag();
    */

    scarot_div_final = div;
    scarot_div_init = 0;
}

void draw_scale_rotate_ins_1()  {
    string ss =
    "A vector can be scaled and rotated with polar form (r, theta)."
    "\nIf we multiply x with r, vector is scaled by x."
    "\nIf we multiply y with theta, vector is rotated by y"
    "\n\nHere we have a (1, 0) vector. Let's scale it by 5 and rotate it by 90 degree."

    "\n                                        Click Enter To Continue"
    ;
    print_doc(0.01, 0.9, ss, 70, 0);
    ///cout << "hi2\n";
}

void draw_scale_rotate_vec() {
    draw_axis();
    draw_right_region();

    if(ent_pressed == 0)    {
        glColor3f(1.0, 0.0, 0.0);///red color
        glBegin(GL_LINES);
            glVertex3f(pnt_x_axis(0), pnt_y_axis(0), 0);
            glVertex3f(pnt_x_axis(1), pnt_y_axis(0), 0);
        glEnd();

        glColor3f(0.0, 0.0, 0.0);///black color
        draw_scale_rotate_ins_1();
    }
    else if(ent_pressed == 1)    {
        if(scarot_div_init == scarot_div_final)    {
            ent_pressed = 2;
            return;
        }

        complex<double> yo = {1.0 / 10.0, 0};
        yo *= polar(scarot_r_init, scarot_theta_init);

        glColor3f(1.0, 0.0, 0.0);///red color
        glBegin(GL_LINES);
            glVertex3f(0 - 0.5, 0, 0);
            glVertex3f(yo.real() - 0.5, yo.imag(), 0);
        glEnd();

        scarot_r_init += scarot_r_inc;
        scarot_theta_init += scarot_theta_inc;
        cout << scarot_r_init << " " << scarot_theta_init << "\n";
        scarot_div_init++;
    }
    else if(ent_pressed == 2)    {
        complex<double> yo = {1.0 / 10.0, 0};
        yo *= polar(scarot_r_init, scarot_theta_init);

        glColor3f(1.0, 0.0, 0.0);///red color
        glBegin(GL_LINES);
            glVertex3f(0 - 0.5, 0, 0);
            glVertex3f(yo.real() - 0.5, yo.imag(), 0);
        glEnd();
    }

}
///scaling and rotating vector part end

///parametric representation of a line start
void draw_par_rep_1()  {
    string ss =
    "If we are given two vectors (can also represent it as points) A and B, we can represent a line which goes through these points in parametric form."
    "\n\nParametric form of Line : A + t(B - A), where t is scaling factor. and giving t different values means giving different points on that line"
    "\nt = 0 means point A"
    "\nt = 1 means point B"
    "\nt = [0 to 1] means point inside vector [A to B]"
    "\nt = (greater than 1) means points beyond B"
    "\nt = [-x] means points in opposite direction"

    "\n\n\n                                   Plot 2 points on white space"
    "\n                                                      Then"
    "\n                                        Click Enter To Continue"
    ;
    print_doc(0.01, 0.9, ss, 70, 0);

    ///cout << "hi2\n";
}
void draw_par_rep_2()  {
    string ss =
    "\n\nA line is shown which goes through two points and it is formed by parametric equation with different values of t"
    "\n\n\n                                        Click Enter To Continue"
    ;
    print_doc(0.01, 0.9, ss, 70, 0);

    ///cout << "hi2\n";
}

void draw_par_rep()  {
    glColor3f(0.0, 0.0, 0.0);///black color
    draw_right_region();

    glColor3f(0.0, 0.0, 0.0);///black color
    if(vv.size() >= 1)    {
        glBegin(GL_POINTS);
            FOR(i, 0, min(1, (int)vv.size() - 1))   {
                glVertex3f(vv[i].F, vv[i].S, 0);
            }
        glEnd();
	}
	if(ent_pressed == 0)    {
	    draw_par_rep_1();
	}
	if(ent_pressed == 1)    {
        complex<double> vec_a = {vv[0].F, vv[0].S};
        complex<double> vec_b = {vv[1].F, vv[1].S};
        complex<double> vec_c, vec_d;///parametric equation
        double t1 = -11;///scaling factor 1
        double t2 = 11;///scaling factor 2
        vec_c = vec_a + (t1 * (vec_b - vec_a));///parametric equation
        vec_d = vec_a + (t2 * (vec_b - vec_a));///parametric equation

        glColor3f(0.0, 0.0, 0.0);///black color
        glBegin(GL_LINES);
            glVertex3f(vec_c.real(), vec_c.imag(), 0);
            glVertex3f(vec_d.real(), vec_d.imag(), 0);
        glEnd();

        glColor3f(1.0, 0.0, 0.0);///red color
        glBegin(GL_POINTS);
            glVertex3f(vv[0].F, vv[0].S, 0);
            glVertex3f(vv[1].F, vv[1].S, 0);
        glEnd();

        glColor3f(0.0, 0.0, 0.0);///black color
        draw_right_region();
	    draw_par_rep_2();
	}
}
///parametric representation of a line end

///intersection point between parametric lines start

complex<double> intersecting_point(complex<double> a, complex<double> b, complex<double> c, complex<double> d)    {
    ///return s and t which is sequentially scaling factor of 1st and 2nd parametric line
    ///1st line line goes throung a and b
    ///2nd line line goes throung a and b
    double s1 = ((d.real() - c.real()) * (c.imag() - a.imag())) - ((c.real() - a.real()) * (d.imag() - c.imag()));
    double s2 = ((d.real() - c.real()) * (b.imag() - a.imag())) - ((b.real() - a.real()) * (d.imag() - c.imag()));

    double s = s1 / s2;

    /**
    ///or,
    double t1 = ((b.real() - a.real()) * (c.imag() - a.imag())) - ((c.real() - a.real()) * (b.imag() - a.imag()));
    double t2 = ((d.real() - c.real()) * (b.imag() - a.imag())) - ((b.real() - a.real()) * (d.imag() - c.imag()));
    double t = t1 / t2;
    */
    return (a + (s * (b - a)));///parametric equation;
}

void draw_line_int_1()  {
    string ss =
    "Let's say vector A and B forms a line and vector C and D forms a line."
    "\nParametric form of 1st line : A + s(B - A)"
    "\nParametric form of 2nd line : C + t(D - C)"
    "\nIt can be solved with A + s(B - A) = C + t(D - C), where s and t are unknown"
    "\nWe can divide 1st line with respect to it's real and imaginary part"
    "\nx axis : Ax + s(Bx - Ax)"
    "\ny axis : Ay + s(By - Ay)"
    "\nWe can divide 2nd line with respect to it's real and imaginary part"
    "\nx axis : Cx + t(Dx - Cx)"
    "\ny axis : Cy + t(Dy - Cy)"
    "\nNow, we get, "
    "\nAx + s(Bx - Ax) = Cx + t(Dx - Cx)"
    "\nAy + s(By - Ay) = Cy + t(Dy - Cy)"
    "\nsolving this equations would get us s and t"
    "\n                                    Plot 2 points to form a line"
    "\n                                                      Then"
    "\n                                        Click Enter To Continue"
    ;
    print_doc(0.01, 0.9, ss, 70, 0);

    ///cout << "hi2\n";
}
void draw_line_int_2()  {
    string ss =
    "After solving that we get, "
    "\n\ns = s1 / s2 .... where"
    "\ns1 = (Dx - Cx)(Cy - Ay) - (Cx - Ax)(Dy - Cy)"
    "\ns2 = (Dx - Cx)(By - Ay) - (Bx - Ax)(Dy - Cy)"

    "\n\nt = t1 / t2 .... where"
    "\nt1 = (Bx - Ax)(Cy - Ay) - (Cx - Ax)(By - Ay)"
    "\nt2 = (Dx - Cx)(By - Ay) - (Bx - Ax)(Dy - Cy)"

    "\n\nPlotting s in parametric equation of 1st line or t in parametric equation of 2nd line line would give intersecting point"
    "\n\n                                    Plot 2 points to form a line"
    "\n                                                      Then"
    "\n                                        Click Enter To Continue"
    ;
    print_doc(0.01, 0.9, ss, 70, 0);

    ///cout << "hi2\n";
}
void draw_line_int_3()  {
    string ss =
    "\n\n                                    Blue point is intersecting point"
    "\n\n\n                                        Click Enter To Continue"
    ;
    print_doc(0.01, 0.9, ss, 70, 0);

    ///cout << "hi2\n";
}

void draw_line_int() {
    draw_right_region();

    if(ent_pressed == 0)    {
        glColor3f(0.0, 0.0, 0.0);///black color
        draw_line_int_1();

        ///drawing input points
        FOR(i, 0, min(1, (int)vv.size() - 1))   {
            glBegin(GL_POINTS);
                glVertex3f(vv[i].F, vv[i].S, 0);
            glEnd();
        }
    }
    else if(ent_pressed == 1)    {
        glColor3f(0.0, 0.0, 0.0);///black color
        draw_line_int_2();

        ///drawing 1st line
        glColor3f(1.0, 0.0, 0.0);///red color
        glBegin(GL_LINES);
            glVertex3f(vv[0].F, vv[0].S, 0);
            glVertex3f(vv[1].F, vv[1].S, 0);
        glEnd();
        ///drawing input points
        glColor3f(0.0, 0.0, 0.0);///black color
        FOR(i, 0, min(3, (int)vv.size() - 1))   {
            glBegin(GL_POINTS);
                glVertex3f(vv[i].F, vv[i].S, 0);
            glEnd();
        }
    }
    else if(ent_pressed == 2)    {
        glColor3f(0.0, 0.0, 0.0);///black color
        draw_line_int_3();

        ///drawing 1st line
        glColor3f(1.0, 0.0, 0.0);///red color
        glBegin(GL_LINES);
            glVertex3f(vv[0].F, vv[0].S, 0);
            glVertex3f(vv[1].F, vv[1].S, 0);
        glEnd();
        ///drawing 2nd line
        glColor3f(0.0, 1.0, 0.0);///green color
        glBegin(GL_LINES);
            glVertex3f(vv[2].F, vv[2].S, 0);
            glVertex3f(vv[3].F, vv[3].S, 0);
        glEnd();
        complex<double> int_pnt = intersecting_point({vv[0].F, vv[0].S}, {vv[1].F, vv[1].S}, {vv[2].F, vv[2].S}, {vv[3].F, vv[3].S});
        glColor3f(0.0, 0.0, 1.0);///blue color
        glBegin(GL_POINTS);
            glVertex3f(int_pnt.real(), int_pnt.imag(), 0);
        glEnd();
    }
}
///intersection point between parametric lines end

///convex hull part start
vector<PDD> up_hull;///this contains the points of upper_hull of convex hull
vector<PDD> down_hull;///this contains the points of lower_hull of convex hull
vector<PDD> hull;///this contains the whole hull
vector<string> conv_sim1, conv_sim2;
vector<PDD> conv_simxy1, conv_simxy2;
vector< pair<double, double> > conv_ara1, conv_ara2;
int conv_pos1 = -1, conv_pos2;
int conv_pos_sim1 = -1, conv_pos_sim2 = -1;

void init_convex_hull()  {
    up_hull.clear();down_hull.clear();hull.clear();
    conv_sim1.clear(), conv_sim2.clear();
    conv_simxy1.clear(); conv_simxy2.clear();
    conv_ara1.clear();conv_ara2.clear();
    conv_pos1 = -1, conv_pos2;
    conv_pos_sim1 = -1, conv_pos_sim2 = -1;
}

void run_convex_hull()   {

    conv_ara1 = vv;
    sort(conv_ara1.begin(), conv_ara1.end());
    conv_ara1.erase(unique(conv_ara1.begin(), conv_ara1.end()), conv_ara1.end());

    up_hull.clear();
    FOR(i, 0, (int)conv_ara1.size() - 1)   {
        while(1) {
            if(up_hull.size() <= 1) {
                break;
            }
            complex<double> p = complex<double>(up_hull[ up_hull.size() - 2 ].F, up_hull[ up_hull.size() - 2 ].S);
            complex<double> n = complex<double>(up_hull[ up_hull.size() - 1 ].F, up_hull[ up_hull.size() - 1 ].S);
            complex<double> x = complex<double>(conv_ara1[i].F, conv_ara1[i].S);
            if((conj(n - p) * (x - p)).imag() > 0)    {///cw when -ve, linear when 0
                conv_sim1.PB("pop");
                up_hull.pop_back();
            }
            else    {
                break;
            }
        }
        conv_sim1.PB("push");
        up_hull.PB(conv_ara1[i]);
    }
    ///sorting x small to big then sorting y big to small
    conv_ara2 = conv_ara1;
    FOR(i, 0, (int)conv_ara2.size() - 1)   {
        int j = i;
        while((j < conv_ara2.size()) && (conv_ara2[j].F == conv_ara2[i].F)) {
            j++;
        }
        reverse(conv_ara2.begin() + i, conv_ara2.begin() + j);
        i = j - 1;
    }
    //**
    down_hull.clear();
    ROF(i, (int)conv_ara2.size() - 1, 0)   {
        while(1) {
            if(down_hull.size() <= 1) {
                break;
            }
            complex<double> p = complex<double>(down_hull[ down_hull.size() - 2 ].F, down_hull[ down_hull.size() - 2 ].S);
            complex<double> n = complex<double>(down_hull[ down_hull.size() - 1 ].F, down_hull[ down_hull.size() - 1 ].S);
            complex<double> x = complex<double>(conv_ara2[i].F, conv_ara2[i].S);
            //int n = up_hull[ up_hull.size() - 1 ];
            if((conj(n - p) * (x - p)).imag() > 0)    {///cw when -ve, linear when 0
                conv_sim2.PB("pop");
                down_hull.pop_back();
            }
            else    {
                break;
            }
        }
        conv_sim2.PB("push");
        down_hull.PB(conv_ara2[i]);
    }
    hull.clear();
    FOR(i, 0, (int)up_hull.size() - 2)   {
        hull.PB(up_hull[i]);
    }
    FOR(i, 0, (int)down_hull.size() - 2)   {
        hull.PB(down_hull[i]);
    }

    FOR(i, 0, (int)vv.size() - 1)   {
        cout << i << " : " << vv[i].F << " " << vv[i].S << "\n";
    }
    /**
    FOR(i, 0, (int)hull.size() - 1)   {
        cout << i << " : " << hull[i].F << " " << hull[i].S << "\n";
    }
    */
    //*/
}
void draw_conv_instruction_1()  {
    string ss =
    "A convex hull is the smallest convex polygon that contains all points of a given set."
    "\n\nAndrew’s algorithm provides an easy way to construct the convex hull for a set of points in O(nlogn) time."
    "It uses sweepline technique"
    "The algorithm first locates the leftmost and rightmost points, and then constructs the convex hull in two parts: "
    "first the upper hull and then the lower hull."

    "\n\n\n                       Plot atleast 3 non-linear points on white space"
    "\n                                                      Then"
    "\n                                        Click Enter To Continue"
    ;
    print_doc(0.01, 0.9, ss, 70, 0);
    ///cout << "hi2\n";
}

void draw_conv_instruction_upper_hull()  {
    string ss =
    "Upper Hull Construction : "
    "First, we sort(ascending order) the points primarily according to x coordinates and secondarily according to y coordinates."
    "After this, we go through the points and add each point to the hull."
    "Always after adding a point to the hull, we make sure that the last line segment in the hull does not turn ccw(counter clockwise)."
    "As long as it turns ccw, we repeatedly remove the second last point from the hull."
    "\n\n\n                                        Click Enter To Continue"
    ;
    print_doc(0.01, 0.9, ss, 70, 0);
    ///cout << "hi2\n";
}

void draw_conv_instruction_lower_hull()  {
    string ss =
    "Upper Hull Complete"
    "\n\nLower Hull Construction : "
    "First, we sort(descending order) the points primarily according to x coordinates and secondarily according to y coordinates."
    "After this, we go through the points and add each point to the hull."
    "Always after adding a point to the hull, we make sure that the last line segment in the hull does not turn ccw(counter clockwise)."
    "As long as it turns ccw, we repeatedly remove the second last point from the hull."
    "\n\n\n                                        Click Enter To Continue"
    ;
    print_doc(0.01, 0.9, ss, 70, 0);
    ///cout << "hi2\n";
}

void draw_conv_valid()  {
    string ss =
    "\n\n\n\n\n\n\n\n"
    "                                         turning clockwise"
    "\n                                         pushing into hull"
    ;
    print_doc(0.01, 0.9, ss, 70, 0);
    ///cout << "hi2\n";
}

void draw_conv_invalid()  {
    string ss =
    "\n\n\n\n\n\n\n\n"
    "                                        turning counterclockwise"
    "\n                                  removing second last point from hull"
    ;
    print_doc(0.01, 0.9, ss, 70, 0);
    ///cout << "hi2\n";
}

void drawConvexHull()    {

	draw_right_region();

	if(ent_pressed == 0)    {
        draw_conv_instruction_1();
	}
	if(ent_pressed == 1)    {
        draw_conv_instruction_upper_hull();
	}

    ///printing upper hull which has already been constructed
    glColor3f(0.0, 0.0, 1.0);
    FOR(i, 0, (int)conv_simxy1.size() - 2)   {
        glBegin(GL_LINES);
            glVertex3f(conv_simxy1[i].F, conv_simxy1[i].S, 0);
            glVertex3f(conv_simxy1[i + 1].F, conv_simxy1[i + 1].S, 0);
        glEnd();

    }
    if(ent_pressed == 2)    {
        ///check if full upper hull has been constructed
        if(conv_pos_sim1 != (int)conv_sim1.size() - 1)    {
            ///if not then consider next line vertex
            conv_pos_sim1++;
            if(conv_sim1[conv_pos_sim1] == "push")    {
                conv_pos1++;
                if(conv_simxy1.size())    {
                    draw_conv_valid();
                    glColor3f(0.0, 1.0, 0.0);///green color

                    glBegin(GL_LINES);
                        glVertex3f(conv_simxy1[conv_simxy1.size() - 1].F, conv_simxy1[conv_simxy1.size() - 1].S, 0);
                        glVertex3f(conv_ara1[conv_pos1].F, conv_ara1[conv_pos1].S, 0);
                    glEnd();
                }
                conv_simxy1.PB(conv_ara1[conv_pos1]);

            }
            else    {
                draw_conv_invalid();
                glColor3f(1.0, 0.0, 0.0);///red color

                glBegin(GL_LINES);
                    glVertex3f(conv_simxy1[conv_simxy1.size() - 1].F, conv_simxy1[conv_simxy1.size() - 1].S, 0);
                    glVertex3f(conv_ara1[conv_pos1 + 1].F, conv_ara1[conv_pos1 + 1].S, 0);
                glEnd();
                glBegin(GL_LINES);
                    glVertex3f(conv_simxy1[conv_simxy1.size() - 1].F, conv_simxy1[conv_simxy1.size() - 1].S, 0);
                    glVertex3f(conv_simxy1[conv_simxy1.size() - 2].F, conv_simxy1[conv_simxy1.size() - 2].S, 0);
                glEnd();

                conv_simxy1.pop_back();
            }
        }
        else    {
            ent_pressed = 3;
            spd = sim_spd;
        }
    }


    else if(ent_pressed == 3) {
        draw_conv_instruction_lower_hull();
    }

    ///printing lower hull which has already been constructed
    glColor3f(0.0, 0.0, 1.0);
    FOR(i, 0, (int)conv_simxy2.size() - 2)   {
        glBegin(GL_LINES);
            glVertex3f(conv_simxy2[i].F, conv_simxy2[i].S, 0);
            glVertex3f(conv_simxy2[i + 1].F, conv_simxy2[i + 1].S, 0);
        glEnd();

    }
    /**down hull**/
    if(ent_pressed == 4)    {
        if((conv_pos_sim1 == (int)conv_sim1.size() - 1) && (conv_pos_sim2 != (int)conv_sim2.size() - 1))    {
            ///spd = 2000;
            conv_pos_sim2++;
            if(conv_sim2[conv_pos_sim2] == "push")    {
                conv_pos2--;
                if(conv_simxy2.size())    {
                    draw_conv_valid();
                    glColor3f(0.0, 1.0, 0.0);
                    glBegin(GL_LINES);
                        glVertex3f(conv_simxy2[conv_simxy2.size() - 1].F, conv_simxy2[conv_simxy2.size() - 1].S, 0);
                        glVertex3f(conv_ara2[conv_pos2].F, conv_ara2[conv_pos2].S, 0);
                    glEnd();
                }
                conv_simxy2.PB(conv_ara2[conv_pos2]);

            }
            else    {
                draw_conv_invalid();
                glColor3f(1.0, 0.0, 0.0);

                glBegin(GL_LINES);
                    glVertex3f(conv_ara2[conv_pos2 - 1].F, conv_ara2[conv_pos2 - 1].S, 0);
                    glVertex3f(conv_simxy2[conv_simxy2.size() - 1].F, conv_simxy2[conv_simxy2.size() - 1].S, 0);

                glEnd();
                glBegin(GL_LINES);
                    glVertex3f(conv_simxy2[conv_simxy2.size() - 2].F, conv_simxy2[conv_simxy2.size() - 2].S, 0);
                    glVertex3f(conv_simxy2[conv_simxy2.size() - 1].F, conv_simxy2[conv_simxy2.size() - 1].S, 0);
                glEnd();

                conv_simxy2.pop_back();
            }
        }
        else    {
            ent_pressed = 5;
            spd = sim_spd;
        }
    }

    glColor3f(0.0, 0.0, 0.0);
    ///printing all points
    glBegin(GL_POINTS);
        for(int i = 0; i < (int)vv.size(); i++)   {
            glVertex3f(vv[i].first, vv[i].second, 0);
        }
    glEnd();

	//glutSwapBuffers();
}
///convex hull part end

///closest pair of points start
double cpop_tmp_dd = 1000000000;///minimum distance between points with bruteforce
vector<PDD> cpop_ara;///copy of input to work with
multiset<PDD> cpop_ara_y;///{y, x}
priority_queue<PDD> cpop_valid;///contains valid points between x - d
double cpop_dd;///minimum distance between points
int cpop_pos = -1;
PDD cpop_ap1, cpop_ap2;///this points represent minimum line

void draw_cpop_instruction_1()  {
    string ss =
    "Given a set of n points, closest pair of problem is to find two points whose Euclidean distance is minimum."
    "\nThis problem can be solved in O(nlogn) time using a sweep line algorithm."
    "We go through the points from left to right and maintain a value d: the minimum distance between two points seen so far."
    "At each point, we find the nearest point to the left. "
    "If the distance is less than d, it is the new minimum distance and we update the value of d."
    "If the current point is (x,y) and there is a point to the left within a distance of less than d, the x coordinate of such a point must be between [x-d,x] and the y coordinate must be between [y-d,y+d]."
    "Thus, it suffices to only consider points that are located in those ranges, which makes the algorithm efficient."
    "\nImplementation details : We can keep a set/bst to hold points inside range [x - d, x] and can check only values from range [y-d,y+d] from set/bst only using binary search."
    "\n\n                             Plot atleast 2 points on white space"
    "\n                                                      Then"
    "\n                                        Click Enter To Continue"
    ;
    ///cout << "hi1\n";
    print_doc(0.01, 0.9, ss, 70, 0);
    ///cout << "hi2\n";
}

void draw_cpop_sim_instruction()  {
    string ss =
    "Yellow Region : Active region(points to consider from, range x=>{x,x-d}, y=>{y-d,y+d})."
    "\nGray Region : Semi-Active region(points inside our set/bst, range x=>{x,x-d})."
    "\n\nGreen Line : Lines representing current minimum distance."
    "\nRed Line : Lines representing past minimum distance which just got removed."
    "\nGray Line : Lines which are inside active region and considered for checking minimum distance."
    "\n\n\n                                        Click Enter To End Simulation."
    ;
    ///cout << "hi1\n";
    print_doc(0.01, 0.9, ss, 70, 0);
    ///cout << "hi2\n";
}

void drawClosestPair()    {


	draw_right_region();
	if(ent_pressed == 0)    {
        glColor3f(0.0, 0.0, 0.0);
        draw_cpop_instruction_1();
	}

    ///drawing all points
    glColor3f(0.0, 0.0, 0.0);
    glBegin(GL_POINTS);
        for(int i = 0; i < (int)vv.size(); i++)   {
            glVertex3f(vv[i].first, vv[i].second, 0);
        }
    glEnd();


	if(ent_pressed == 1 && cpop_pos != cpop_ara.size() - 1)    {
        draw_cpop_sim_instruction();

        cpop_pos++;
        ///drawing (x - d) region
        ///glColor3f(0.6, 0.6, 1);
        ///glColor3f(col(102), col(102), col(153));
        glColor3f(col(128), col(128), col(128));
        glBegin(GL_POLYGON);
            glVertex3f(cpop_ara[cpop_pos].F - cpop_dd, 1, 0);
            glVertex3f(cpop_ara[cpop_pos].F, 1, 0);
            glVertex3f(cpop_ara[cpop_pos].F, -1, 0);
            glVertex3f(cpop_ara[cpop_pos].F - cpop_dd, -1, 0);
        glEnd();
        ///drawing active region (x - d), (y - d) to (y + d)
        glColor3f(1.0, 1.0, 0.6);
        glBegin(GL_POLYGON);
            glVertex3f(cpop_ara[cpop_pos].F - cpop_dd, cpop_ara[cpop_pos].S + cpop_dd, 0);
            glVertex3f(cpop_ara[cpop_pos].F, cpop_ara[cpop_pos].S + cpop_dd, 0);
            glVertex3f(cpop_ara[cpop_pos].F, cpop_ara[cpop_pos].S - cpop_dd, 0);
            glVertex3f(cpop_ara[cpop_pos].F - cpop_dd, cpop_ara[cpop_pos].S - cpop_dd, 0);
        glEnd();

        vector<PDD> act;
        ///drawing all points
        glBegin(GL_POINTS);
            for(int i = 0; i < (int)vv.size(); i++)   {
                if((vv[i].F >= cpop_ara[cpop_pos].F - cpop_dd) && (vv[i].F <= cpop_ara[cpop_pos].F) && (vv[i].S >= cpop_ara[cpop_pos].S - cpop_dd) && (vv[i].S <= cpop_ara[cpop_pos].S + cpop_dd))    {
                    ///if point is inside active region
                    ///glColor3f(0.6, 0.6, 1);
                    ///glColor3f(col(102), col(102), col(153));
                    glColor3f(col(128), col(128), col(128));
                    act.PB(vv[i]);
                }
                else    {
                    glColor3f(0.0, 0.0, 0.0);
                }
                glVertex3f(vv[i].first, vv[i].second, 0);
            }
        glEnd();


        ///drawing lines with current point to points inside active region
        for(int i = 0; i < (int)act.size(); i++)   {
            ///glColor3f(1.0, 0.0, 0.0);
            ///glColor3f(0.6, 0.6, 1);
            ///glColor3f(col(102), col(102), col(153));
            glColor3f(col(128), col(128), col(128));
            glBegin(GL_LINES);
                glVertex3f(act[i].first, act[i].second, 0);
                glVertex3f(cpop_ara[cpop_pos].first, cpop_ara[cpop_pos].second, 0);
            glEnd();
        }
        ///drawing past minimum line
        if(cpop_dd != 1000000000)    {
            glColor3f(1.0, 0.0, 0.0);
            glBegin(GL_LINES);
                glVertex3f(cpop_ap1.F, cpop_ap1.S, 0);
                glVertex3f(cpop_ap2.F, cpop_ap2.S, 0);
            glEnd();
        }


        ///drawing current point
        glColor3f(0.0, 1.0, 0.0);
        glBegin(GL_POINTS);
            glVertex3f(cpop_ara[cpop_pos].F, cpop_ara[cpop_pos].S, 0);
        glEnd();


        ///drawing end, algo part starts
        while((!cpop_valid.empty()) && (cpop_valid.top().F < cpop_ara[cpop_pos].F - cpop_dd)) {
            ///remove this point from set as it is outside x - d region
            cpop_ara_y.erase(cpop_ara_y.find({cpop_valid.top().S, cpop_valid.top().F}));
            cpop_valid.pop();
        }
        ///we have to search points which are in range y - dd to y + dd
        ///x doesn't matter now as we have already handled all outside range x - dd above
        double x = -1000000000, y = cpop_ara[cpop_pos].S - cpop_dd;
        while(1) {
            auto it = cpop_ara_y.upper_bound({y, x});
            if(x == -1000000000)    {
                it = cpop_ara_y.lower_bound({y, x});
            }

            if(it == cpop_ara_y.end())    {
                break;
            }
            x = (*it).S;
            y = (*it).F;
            if(y > cpop_ara[cpop_pos].S + cpop_dd)    {
                break;
            }
            if(dis({x, y}, {cpop_ara[cpop_pos].F, cpop_ara[cpop_pos].S}) < cpop_dd)    {
                cpop_dd = dis({x, y}, {cpop_ara[cpop_pos].F, cpop_ara[cpop_pos].S});
                cpop_ap1 = {x, y};cpop_ap2 = {cpop_ara[cpop_pos].F, cpop_ara[cpop_pos].S};
            }
            ///cpop_dd = min(cpop_dd, dis({x, y}, {cpop_ara[cpop_pos].F, cpop_ara[cpop_pos].S}));
        }
        cpop_valid.push(cpop_ara[cpop_pos]);
        cpop_ara_y.insert({cpop_ara[cpop_pos].S, cpop_ara[cpop_pos].F});
    }
    else if(ent_pressed == 1)   {
        ent_pressed = 2;
    }
    else if(ent_pressed == 2)    {///simulation done
        draw_cpop_sim_instruction();
    }

    ///drawing current minimum line if any
    if(cpop_dd != 1000000000)    {
        glColor3f(0.0, 1.0, 0.0);
        glBegin(GL_LINES);
            glVertex3f(cpop_ap1.F, cpop_ap1.S, 0);
            glVertex3f(cpop_ap2.F, cpop_ap2.S, 0);
        glEnd();
        if(abs(dis({cpop_ap1.F, cpop_ap1.S}, {cpop_ap2.F, cpop_ap2.S}) - cpop_dd) >= 1e-9)    {
            ///cout << "na : " << cpop_dd << " " << dis({cpop_ap1.F, cpop_ap1.S}, {cpop_ap2.F, cpop_ap2.S}) << "\n";
        }
    }
    /**
    if(ent_pressed == 2 && (cpop_pos == cpop_ara.size()-1) && cpop_tmp_dd >= 0)    {
        cout << "dis : " << cpop_dd << " " << cpop_tmp_dd << "\n";
        cout << cpop_ap1.F << " " << cpop_ap1.S << ", " << cpop_ap2.F << " " << cpop_ap2.S << "\n";
        cpop_tmp_dd = -1;
    }
    //*/

}
///closest pair of points end


///furthest pair of points start

PDD fpop_mn_x = {1000000000, -1};///
int fpop_mn_xpos = -1;
PDD fpop_mx_x = {-1000000000, -1};
int fpop_mx_xpos = -1;

PDD fpop_support1 = {0, 1};///this is 90 degree up
PDD fpop_support2 = {0, -1};///this is 90 degree down
int fpop_pos1 = -1;
int fpop_pos2 = -1;

double fpop_brute_ans = 0, fpop_ans = 0;
double fpop_tot_rot = 0;///total rotation
int fpop_ans_pnt1, fpop_ans_pnt2;


void init_fpop() {
    fpop_mn_x = {1000000000, -1};///point whose value in x axis is minimum
    fpop_mn_xpos = -1;///position of point whose value in x axis is minimum
    fpop_mx_x = {-1000000000, -1};///point whose value in x axis is max
    fpop_mx_xpos = -1;///position of point whose value in x axis is max

    fpop_support1 = {0, 1};///this is 90 degree up
    fpop_support2 = {0, -1};///this is 90 degree down
    fpop_pos1 = -1;fpop_pos2 = -1;
    fpop_brute_ans = 0;fpop_ans = 0;
    fpop_tot_rot = 0;///total rotation
    fpop_ans_pnt1 = -1;fpop_ans_pnt1 = -1;
}

double pi = 3.1415;

double get_angle(complex<double> a, complex<double> b) {
    ///get angle between two vector where a goes clockwise to b
    ///this function is written for convex polygon where angle between them is less then 180
    double ang1 = arg(a);
    double ang2 = arg(b);
    if(ang1 < 0)    {
        if(ang2 < 0)    {
            return abs(ang2 - ang1);
        }
        else if(ang2 > 0)    {
            return ((pi) - abs(ang1)) + ((pi) - ang2);
        }
    }
    else if(ang1 >= 0)    {
        if(ang2 < 0)    {
            return abs(ang2) + ang1;
        }
        else if(ang2 > 0)    {
            return ang2 - ang1;
        }
    }
}

complex<double> rotate_it(complex<double> a, double theta) {///rotate vector by theta
    return a * polar(1.0, theta);
}

void draw_fpop_ins_1()  {
    string ss =
    "Furthest pair of points between n points can be found if we take convex hull of those points and find diameter of it. "
    "diameter of a convex polygon means furthest distance between any two vertex of that polygon. "
    "\n\nA linear O(n) solution exists with rotating calipers technique. "
    "We put the polygon between the two jaws of a caliper, and tighten. "
    "We then rotate the calipers a full circle around the polygon, while keeping the caliper tight at all times. "
    "The answer is then the maximum distance between the two jaws at any point during the procedure. "
    "That two jaws are two parallel supports called antipodal and lines between antipodal are only to be considered to find diameter"
    "\n\n\n                         plot atleast 3 non-linear points on white space"
    "\n                                                      Then"
    "\n                                        Click Enter To Continue"
    ;
    print_doc(0.01, 0.9, ss, 75, 0);

    ///cout << "hi2\n";
}

void draw_fpop_sim_1()  {
    string ss =
    "\n\nblue and red lines are two antipodal"
    "\ndark violet line is current distance which is between two antipodal and is considered for diameter"
    "\nGreen line is big diameter measured up to date"
    "\n\n\n                                        Click Enter To Continue"
    ;
    print_doc(0.01, 0.9, ss, 70, 0);

    ///cout << "hi2\n";
}

void draw_fpop()    {
    draw_right_region();
    glColor3f(0.0, 0.0, 0.0);
    ///drawing points
    glBegin(GL_POINTS);
        FOR(i, 0, (int)vv.size() - 1)   {
            glVertex3f(vv[i].F, vv[i].S, 0);
        }
    glEnd();
    ///drawing the polygon
    glColor3f(0.0, 0.0, 0.0);

    FOR(i, 0, (int)hull.size() - 1)   {
        if(i == hull.size() - 1)    {
            glBegin(GL_LINES);
                glVertex3f(hull[i].F, hull[i].S, 0);
                glVertex3f(hull[0].F, hull[0].S, 0);
            glEnd();
            continue;
        }
        glBegin(GL_LINES);
            glVertex3f(hull[i].F, hull[i].S, 0);
            glVertex3f(hull[i + 1].F, hull[i + 1].S, 0);
        glEnd();
    }



    if(ent_pressed == 0)    {
        glColor3f(0.0, 0.0, 0.0);
        draw_fpop_ins_1();
    }
    else if(ent_pressed == 1)    {
        glColor3f(col(204), col(0), col(204));///dark violet
        if(fpop_ans_pnt1 != -1)    {///drawing current diameter between antinodal points
            glBegin(GL_LINES);
                glVertex3f(hull[fpop_pos1].F, hull[fpop_pos1].S, 0);
                glVertex3f(hull[fpop_pos2].F, hull[fpop_pos2].S, 0);
            glEnd();
        }

        if(dis({hull[fpop_pos1].F, hull[fpop_pos1].S}, {hull[fpop_pos2].F, hull[fpop_pos2].S}) > fpop_ans)    {
            fpop_ans = dis({hull[fpop_pos1].F, hull[fpop_pos1].S}, {hull[fpop_pos2].F, hull[fpop_pos2].S});
            fpop_ans_pnt1 = fpop_pos1;
            fpop_ans_pnt2 = fpop_pos2;
        }

        ///drawing 1st support of antinodal lines
        glColor3f(0.0, 0.0, 1.0);
        glBegin(GL_LINES);///transform support back to hull[i]
            glVertex3f(0 + hull[fpop_pos1].F, 0 + hull[fpop_pos1].S, 0);
            glVertex3f(fpop_support1.F + hull[fpop_pos1].F, fpop_support1.S + hull[fpop_pos1].S, 0);
        glEnd();
        ///drawing 2nd support of antinodal lines
        glColor3f(1.0, 0.0, 0.0);
        glBegin(GL_LINES);///transform support back to hull[i]
            glVertex3f(0 + hull[fpop_pos2].F, 0 + hull[fpop_pos2].S, 0);
            glVertex3f(fpop_support2.F + hull[fpop_pos2].F, fpop_support2.S + hull[fpop_pos2].S, 0);
        glEnd();

        ///fpop_ans = max(fpop_ans, dis({hull[fpop_pos1].F, hull[fpop_pos1].S}, {hull[fpop_pos2].F, hull[fpop_pos2].S}));


        int nxpos1 = (fpop_pos1 + 1) % hull.size();
        PDD trans1 = {hull[nxpos1].F - hull[fpop_pos1].F, hull[nxpos1].S - hull[fpop_pos1].S};///transform hull[pos1] to (0, 0)
        double angle1 = get_angle({fpop_support1.F, fpop_support1.S}, {trans1.F, trans1.S});

        /****************************************************************/
        int nxpos2 = (fpop_pos2 + 1) % hull.size();
        PDD trans2 = {hull[nxpos2].F - hull[fpop_pos2].F, hull[nxpos2].S - hull[fpop_pos2].S};///transform hull[pos2] to (0, 0)
        double angle2 = get_angle({fpop_support2.F, fpop_support2.S}, {trans2.F, trans2.S});
        /**
        if(1)  {///eita bad dibo
            complex<double> c1 = {fpop_support1.F, fpop_support1.S};
            complex<double> c2 = {trans1.F, trans1.S};
            cout << "angles1 : " << arg(c1) * 180 / pi << " " << arg(c2) * 180 / pi << "\n";
        }
        if(1)  {///eita bad dibo
            complex<double> c1 = {fpop_support2.F, fpop_support2.S};
            complex<double> c2 = {trans2.F, trans2.S};
            cout << "angles2 : " << arg(c1) * 180 / pi << " " << arg(c2) * 180 / pi << "\n";
        }
        */

        ///double mn_angle;

        double mn_angle = min(abs(angle1), abs(angle2));
        mn_angle *= -1;///as it should turn cw


        if(abs(angle1) < abs(angle2))    {
            fpop_pos1 = nxpos1;
            ///mn_angle = angle1;
        }
        else    {
            fpop_pos2 = nxpos2;
            ///mn_angle = angle2;
        }


        ///cout << "angle : "  << angle1 * 180 / pi << " " << angle2 * 180 / pi << " : " << mn_angle * 180 / pi << "\n";
        ///cout << "pos " << fpop_pos1 << " " << fpop_pos2 << "\n\n";
        /****************************************************************/

        complex<double> yo1 = rotate_it({fpop_support1.F, fpop_support1.S}, mn_angle);
        ///complex<double> yo1 = rotate_it({support1.F, support1.S}, angle1);
        fpop_support1 = {yo1.real(), yo1.imag()};



        /****************************************************************/

        complex<double> yo2 = rotate_it({fpop_support2.F, fpop_support2.S}, mn_angle);
        ///complex<double> yo2 = rotate_it({support2.F, support2.S}, angle2);
        fpop_support2 = {yo2.real(), yo2.imag()};

        fpop_tot_rot += (-1 * mn_angle);
        cout << "total rotation : " << (fpop_tot_rot * 180 / pi) << " " << (mn_angle * 180 / pi) << "\n";
        if(fpop_tot_rot > 2.0 * pi)    {///pi is enough ... but used 2 * pi anyway
            ent_pressed = 2;
        }

        draw_right_region();
        glColor3f(0.0, 0.0, 0.0);
        draw_fpop_sim_1();
    }
    else if(ent_pressed == 2)   {
        cout << fpop_ans << " " << fpop_brute_ans << "\n";
    }

    glColor3f(0.0, 1.0, 0.0);
    if(fpop_ans_pnt1 != -1)    {///drawing max diameter till date
        glBegin(GL_LINES);
            glVertex3f(hull[fpop_ans_pnt1].F, hull[fpop_ans_pnt1].S, 0);
            glVertex3f(hull[fpop_ans_pnt2].F, hull[fpop_ans_pnt2].S, 0);
        glEnd();
    }

    /**
    ///drawing support in (0, 0)
    glBegin(GL_LINES);
        glVertex3f(0, 0, 0);
        glVertex3f(support1.F, support1.S, 0);
    glEnd();
    */

    ///drawing (0, 0) point
    /**
    glColor3f(0.0, 1.0, 0.0);
    glBegin(GL_POINTS);
        glVertex3f(0, 0, 0);
    glEnd();
    */

    return;
}
///furthest pair of points end


///minimum area-enclosing rectangle of points start

PDD maer_mn_x = {1000000000, -1};///point whose value in x axis is minimum
int maer_mn_xpos = -1;///position of point whose value in x axis is minimum
PDD maer_mx_x = {-1000000000, -1};///point whose value in x axis is max
int maer_mx_xpos = -1;///position of point whose value in x axis is max

PDD maer_mn_y = {-1, 1000000000};///
int maer_mn_ypos = -1;
PDD maer_mx_y = {-1, -1000000000};
int maer_mx_ypos = -1;

PDD maer_support1 = {0, 1};///this is 90 degree up, support for leftmost mn_x
PDD maer_support2 = {0, -1};///this is 90 degree down, support for rightmost mx_x
PDD maer_support3 = {0, 1};///this is 180 degree right, support for upmost mx_y
PDD maer_support4 = {0, -1};///this is 180 degree left, support for upmost mn_y
int maer_pos1 = -1;
int maer_pos2 = -1;
int maer_pos3 = -1;
int maer_pos4 = -1;
void init_maer() {
    maer_mn_x = {1000000000, -1};///point whose value in x axis is minimum
    maer_mn_xpos = -1;///position of point whose value in x axis is minimum
    maer_mx_x = {-1000000000, -1};///point whose value in x axis is max
    maer_mx_xpos = -1;///position of point whose value in x axis is max

    maer_mn_y = {-1, 1000000000};///
    maer_mn_ypos = -1;
    maer_mx_y = {-1, -1000000000};
    maer_mx_ypos = -1;

    maer_support1 = {0, 1};///this is 90 degree up, support for leftmost mn_x
    maer_support2 = {0, -1};///this is 90 degree down, support for rightmost mx_x
    maer_support3 = {1, 0};///this is 180 degree right, support for upmost mx_y
    maer_support4 = {-1, 0};///this is 180 degree left, support for upmost mn_y

    maer_pos1 = -1;maer_pos2 = -1;maer_pos3 = -1;maer_pos4 = -1;

    /**
    maer_brute_ans = 0;fpop_ans = 0;
    maer_tot_rot = 0;///total rotation
    maer_ans_pnt1 = -1;fpop_ans_pnt1 = -1;
    */
}

void draw_maer_ins_1()  {
    string ss =
    "minimum area enclosing rectangle  between n points can be found if we take convex hull of those points and apply rotating calipers to it"
    "\n\nWe place two calipers around convex polygon perpendicular to each other and they form a rectangle. "
    "For both calipers, we put the polygon between the two jaws of a caliper, and tighten. "
    "We then rotate the calipers a full circle around the polygon, while keeping the calipers tight at all times. "
    "The answer is then the minimum area of a rectangle formed by two calipers at any point during the procedure. "
    "\n\n\n                         plot atleast 4 non-linear points on white space"
    "\n                                                      Then"
    "\n                                        Click Enter To Continue"
    ;
    print_doc(0.01, 0.9, ss, 72, 0);

    ///cout << "hi2\n";
}

void draw_maer_sim_1()  {
    string ss =
    "\n\nblue lines are parallel to each other"
    "\nred lines are parallel to each other"
    "\nblue and red lines are perpendicular to each other"
    "\ntogether they form a rectangle to be considered for minimum area enclosing rectangle"
    ""
    "\n\n\n                                        Click Enter To Continue"
    ;
    print_doc(0.01, 0.9, ss, 70, 0);

    ///cout << "hi2\n";
}

void draw_maer()    {
    draw_right_region();
    glColor3f(0.0, 0.0, 0.0);
    ///drawing points
    glBegin(GL_POINTS);
        FOR(i, 0, (int)vv.size() - 1)   {
            glVertex3f(vv[i].F, vv[i].S, 0);
        }
    glEnd();
    ///drawing the polygon
    glColor3f(0.0, 0.0, 0.0);

    FOR(i, 0, (int)hull.size() - 1)   {
        if(i == hull.size() - 1)    {
            glBegin(GL_LINES);
                glVertex3f(hull[i].F, hull[i].S, 0);
                glVertex3f(hull[0].F, hull[0].S, 0);
            glEnd();
            continue;
        }
        glBegin(GL_LINES);
            glVertex3f(hull[i].F, hull[i].S, 0);
            glVertex3f(hull[i + 1].F, hull[i + 1].S, 0);
        glEnd();
    }



    if(ent_pressed == 0)    {
        glColor3f(0.0, 0.0, 0.0);
        draw_maer_ins_1();
    }
    else if(ent_pressed == 1)    {
        ///glColor3f(col(204), col(0), col(204));///dark violet



        complex<double> vec_a;
        complex<double> vec_b;
        complex<double> vec_c, vec_d;///parametric equation
        double t1 = -11;///scaling factor 1
        double t2 = 11;///scaling factor 2
        ///drawing 1st support of antinodal lines

        ///transform support back to hull[i]
        vec_a = {0 + hull[maer_pos1].F, 0 + hull[maer_pos1].S};
        vec_b = {maer_support1.F + hull[maer_pos1].F, maer_support1.S + hull[maer_pos1].S};
        vec_c = vec_a + (t1 * (vec_b - vec_a));///parametric equation
        vec_d = vec_a + (t2 * (vec_b - vec_a));///parametric equation
        glColor3f(0.0, 0.0, 1.0);
        glBegin(GL_LINES);
            glVertex3f(vec_c.real(), vec_c.imag(), 0);
            glVertex3f(vec_d.real(), vec_d.imag(), 0);
        glEnd();
        ///drawing 2nd support of antinodal lines
        vec_a = {0 + hull[maer_pos2].F, 0 + hull[maer_pos2].S};
        vec_b = {maer_support2.F + hull[maer_pos2].F, maer_support2.S + hull[maer_pos2].S};
        vec_c = vec_a + (t1 * (vec_b - vec_a));///parametric equation
        vec_d = vec_a + (t2 * (vec_b - vec_a));///parametric equation
        glColor3f(0.0, 0.0, 1.0);
        glBegin(GL_LINES);
            glVertex3f(vec_c.real(), vec_c.imag(), 0);
            glVertex3f(vec_d.real(), vec_d.imag(), 0);
        glEnd();
        ///drawing 3rd support of antinodal lines
        vec_a = {0 + hull[maer_pos3].F, 0 + hull[maer_pos3].S};
        vec_b = {maer_support3.F + hull[maer_pos3].F, maer_support3.S + hull[maer_pos3].S};
        vec_c = vec_a + (t1 * (vec_b - vec_a));///parametric equation
        vec_d = vec_a + (t2 * (vec_b - vec_a));///parametric equation
        glColor3f(1.0, 0.0, 0.0);
        glBegin(GL_LINES);
            glVertex3f(vec_c.real(), vec_c.imag(), 0);
            glVertex3f(vec_d.real(), vec_d.imag(), 0);
        glEnd();
        ///drawing 4th support of antinodal lines
        vec_a = {0 + hull[maer_pos4].F, 0 + hull[maer_pos4].S};
        vec_b = {maer_support4.F + hull[maer_pos4].F, maer_support4.S + hull[maer_pos4].S};
        vec_c = vec_a + (t1 * (vec_b - vec_a));///parametric equation
        vec_d = vec_a + (t2 * (vec_b - vec_a));///parametric equation
        glColor3f(1.0, 0.0, 0.0);
        glBegin(GL_LINES);
            glVertex3f(vec_c.real(), vec_c.imag(), 0);
            glVertex3f(vec_d.real(), vec_d.imag(), 0);
        glEnd();

        ///fpop_ans = max(fpop_ans, dis({hull[fpop_pos1].F, hull[fpop_pos1].S}, {hull[fpop_pos2].F, hull[fpop_pos2].S}));

        ///measuring angle1
        int nxpos1 = (maer_pos1 + 1) % hull.size();
        PDD trans1 = {hull[nxpos1].F - hull[maer_pos1].F, hull[nxpos1].S - hull[maer_pos1].S};///transform hull[pos1] to (0, 0)
        double angle1 = get_angle({maer_support1.F, maer_support1.S}, {trans1.F, trans1.S});

        ///measuring angle2
        int nxpos2 = (maer_pos2 + 1) % hull.size();
        PDD trans2 = {hull[nxpos2].F - hull[maer_pos2].F, hull[nxpos2].S - hull[maer_pos2].S};///transform hull[pos2] to (0, 0)
        double angle2 = get_angle({maer_support2.F, maer_support2.S}, {trans2.F, trans2.S});

        ///measuring angle3
        int nxpos3 = (maer_pos3 + 1) % hull.size();
        PDD trans3 = {hull[nxpos3].F - hull[maer_pos3].F, hull[nxpos3].S - hull[maer_pos3].S};///transform hull[pos2] to (0, 0)
        double angle3 = get_angle({maer_support3.F, maer_support3.S}, {trans3.F, trans3.S});

        ///measuring angle4
        int nxpos4 = (maer_pos4 + 1) % hull.size();
        PDD trans4 = {hull[nxpos4].F - hull[maer_pos4].F, hull[nxpos4].S - hull[maer_pos4].S};///transform hull[pos2] to (0, 0)
        double angle4 = get_angle({maer_support4.F, maer_support4.S}, {trans4.F, trans4.S});


        double mn_angle;

        vector< pair<double, int> > ang_pos;///contains angle and position of it's support
        ang_pos.PB({abs(angle1), 1});
        ang_pos.PB({abs(angle2), 2});
        ang_pos.PB({abs(angle3), 3});
        ang_pos.PB({abs(angle4), 4});

        sort(ang_pos.begin(), ang_pos.end());
        mn_angle = -1 * ang_pos[0].F;

        /**
        double mn_angle = min(abs(angle1), abs(angle2));
        mn_angle *= -1;///as it should turn cw
        */
        if(ang_pos[0].S == 1)    {
            maer_pos1 = nxpos1;
        }
        else if(ang_pos[0].S == 2)    {
            maer_pos2 = nxpos2;
        }
        else if(ang_pos[0].S == 3)    {
            maer_pos3 = nxpos3;
        }
        else if(ang_pos[0].S == 4)    {
            maer_pos4 = nxpos4;
        }

        /**
        if(abs(angle1) < abs(angle2))    {
            maer_pos1 = nxpos1;
            ///mn_angle = angle1;
        }
        else    {
            maer_pos2 = nxpos2;
            ///mn_angle = angle2;
        }
        */


        ///cout << "angle : "  << angle1 * 180 / pi << " " << angle2 * 180 / pi << " : " << mn_angle * 180 / pi << "\n";
        ///cout << "pos " << fpop_pos1 << " " << fpop_pos2 << "\n\n";

        ///rotate 1st support with minimum angle
        complex<double> yo1 = rotate_it({maer_support1.F, maer_support1.S}, mn_angle);
        ///complex<double> yo1 = rotate_it({support1.F, support1.S}, angle1);
        maer_support1 = {yo1.real(), yo1.imag()};

        ///rotate 2nd support with minimum angle
        complex<double> yo2 = rotate_it({maer_support2.F, maer_support2.S}, mn_angle);
        ///complex<double> yo2 = rotate_it({support2.F, support2.S}, angle2);
        maer_support2 = {yo2.real(), yo2.imag()};

        ///rotate 3rd support with minimum angle
        complex<double> yo3 = rotate_it({maer_support3.F, maer_support3.S}, mn_angle);
        ///complex<double> yo2 = rotate_it({support2.F, support2.S}, angle2);
        maer_support3 = {yo3.real(), yo3.imag()};

        ///rotate 4th support with minimum angle
        complex<double> yo4 = rotate_it({maer_support4.F, maer_support4.S}, mn_angle);
        ///complex<double> yo2 = rotate_it({support2.F, support2.S}, angle2);
        maer_support4 = {yo4.real(), yo4.imag()};

        draw_right_region();
        glColor3f(0.0, 0.0, 0.0);
        draw_maer_sim_1();
    }

    return;
}
///minimum area-enclosing rectangle of points end

void drawScene() {
    ///cout << "hi\n";
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity(); //Reset the drawing perspective
	glMatrixMode(GL_MODELVIEW);
	glPointSize(5);

    if(page_no == 1)    {
        draw_1st_page();
    }
    else if(page_no == 2)   {
        draw_2nd_page();
    }
    else if(page_no == 3)   {
        draw_3rd_page();
    }
    else if(page_no == 4)   {
        if(p3_menu_no == vec_rep_idx)    {///vector representation
            draw_vec_rep();
        }
        else if(p3_menu_no == scale_rotate_vec_idx)   {///scaling and rotating a vector
            draw_scale_rotate_vec();
        }
        else if (p3_menu_no == par_rep_idx)  {///parametric representation of a line
            draw_par_rep();
        }
        else if (p3_menu_no == line_int_idx)  {///intersection point between parametric lines
            draw_line_int();
        }
        else if(p3_menu_no == conv_idx)    {///convex hull
            drawConvexHull();
        }
        else if(p3_menu_no == cpop_idx)    {///cpop
            drawClosestPair();
        }
        else if(p3_menu_no == fpop_idx)    {///fpop
            draw_fpop();
        }
        else if(p3_menu_no == maer_idx)    {///minimum area enclosing rectangle
            draw_maer();
        }
    }
    else if(page_no == 5)   {///about page
        draw_about_page();
    }
    else if(page_no == 6)   {///reference page
        draw_ref_page();
    }
	glutSwapBuffers();
}

void keyboard(unsigned char key, int x, int y) {
//find key codes: https://www.cambiaresearch.com/articles/15/javascript-char-codes-key-codes

    if(key == 13)    {///means enter pressed
        if(page_no == 1)    {
            page_no = 2;
        }
        else if(page_no == 2)    {
            if(p2_menu_no == 1)    {///menu
                page_no = 3;
            }
            else if(p2_menu_no == 2)    {///about
                page_no = 5;
                drawScene();
            }
            else if(p2_menu_no == 3)    {///reference
                page_no = 6;
                drawScene();
            }
            else if(p2_menu_no == 4)   {///exit
                exit(0);
            }
        }
        else if(page_no == 3)    {///geometry menu
            ent_pressed = 0;
            page_no = 4;
            if(p3_menu_no == vec_rep_idx)    {///vector representation
                ///take input please
                vv.clear();
                take_inp = 1;
                ///initialize anything related to vector representation here
            }
            else if (p3_menu_no == scale_rotate_vec_idx)  {///scaling and rotating a vector
                init_scale_rotate_vec(5, pi / 2, 50);
            }
            else if (p3_menu_no == par_rep_idx)  {///parametric representation of a line
                ///take input please
                vv.clear();
                take_inp = 1;
            }
            else if (p3_menu_no == line_int_idx)  {///intersection point between parametric lines
                ///take input please
                vv.clear();
                take_inp = 1;
            }
            else if(p3_menu_no == conv_idx)    {///convex hull
                ///take input please
                vv.clear();
                take_inp = 1;
                ///initialize anything related to convex hull here
                init_convex_hull();
            }
            else if(p3_menu_no == cpop_idx)   {///closest pair of points

                ///take input please
                vv.clear();
                take_inp = 1;
                ///initialize anything related to cpop here
                cpop_ara_y.clear();
                while(!cpop_valid.empty()){cpop_valid.pop();}
                cpop_dd = 1000000000;
                cpop_pos = -1;
            }
            else if(p3_menu_no == fpop_idx)   {///furthest pair of points
                ///take input please
                vv.clear();
                take_inp = 1;
                ///initialize anything related to convex hull here
                init_convex_hull();
                ///initialize anything related to fpop here
                init_fpop();

            }
            else if(p3_menu_no == maer_idx)   {///minimum area enclosing rectangle
                ///take input please
                vv.clear();
                take_inp = 1;
                ///initialize anything related to convex hull here
                init_convex_hull();
                ///initialize anything related to fpop here
                init_maer();
            }
            else if(p3_menu_no == back_idx)   {///go back
                page_no = 2;
            }
        }
        else if(page_no == 4)    {
            if(p3_menu_no == vec_rep_idx)    {///vector representation
                if(vv.size() == 0)    {
                    return;
                }
                if(ent_pressed == 0)    {
                    while(vv.size() > 1) {
                        vv.pop_back();
                    }
                    ///stop taking input please
                    take_inp = 0;
                    ent_pressed = 1;
                }
                else if(ent_pressed == 1)    {
                    page_no = 3;
                }
            }
            else if (p3_menu_no == scale_rotate_vec_idx)  {///scaling and rotating a vector
                if(ent_pressed == 0)    {
                    ent_pressed = 1;
                }
                else if(ent_pressed == 1)  {
                    spd = mn_spd;
                }
                else if(ent_pressed == 2)  {
                    page_no = 3;
                    spd = fst_spd;
                }
            }
            else if (p3_menu_no == par_rep_idx)  {///parametric representation of a line
                if(vv.size() < 2)    {
                    return;
                }
                if(ent_pressed == 0)    {
                    while(vv.size() > 2) {
                        vv.pop_back();
                    }
                    ///stop taking input please
                    take_inp = 0;
                    ent_pressed = 1;
                }
                else if(ent_pressed == 1)  {
                    page_no = 3;
                }
            }
            else if (p3_menu_no == line_int_idx)  {///intersection point between parametric lines
                if(ent_pressed == 0)    {
                    if(vv.size() < 2)    {
                        return;
                    }
                    while(vv.size() > 2) {
                        vv.pop_back();
                    }
                    ent_pressed = 1;
                }
                else if(ent_pressed == 1)  {
                    if(vv.size() < 4)    {
                        return;
                    }
                    while(vv.size() > 4) {
                        vv.pop_back();
                    }
                    ///stop taking input please
                    take_inp = 0;
                    ent_pressed = 2;
                }
                else if(ent_pressed == 2)  {
                    page_no = 3;
                }
            }
            else if(p3_menu_no == conv_idx)    {///convex hull
                if(vv.size() < 3)    {
                    return;
                }
                if(ent_pressed == 0)    {
                    ///stop taking input please
                    take_inp = 0;
                    spd = sim_spd;
                    ent_pressed = 1;
                }
                else if(ent_pressed == 1)   {

                    run_convex_hull();
                    ent_pressed = 2;
                    conv_pos2 = conv_ara2.size();

                }
                else if(ent_pressed == 2)   {///upper hull simulation
                    ///if someone preses enter simulation should be over
                    spd = mn_spd;
                }
                else if(ent_pressed == 3)   {///upper hull simulation done
                    ent_pressed = 4;
                }
                else if(ent_pressed == 4)   {///lower hull simulation
                    spd = mn_spd;
                }
                else if(ent_pressed == 5)   {///lower hull simulation done
                    page_no = 3;
                    spd = fst_spd;
                }
            }
            else if(p3_menu_no == cpop_idx)   {///closest pair of points
                if(ent_pressed == 0)    {
                    if(vv.size() < 2)    {
                        return;
                    }
                    ///stop taking input please
                    take_inp = 0;
                    spd = sim_spd;
                    ///do initializing after taking input
                    cpop_ara = vv;
                    ///sorted by x from small to large then sorted similarly for y
                    sort(cpop_ara.begin(), cpop_ara.end());
                    cpop_tmp_dd = 1000000000;
                    ///brute_force
                    FOR(i, 0, (int)cpop_ara.size() - 1)   {
                        FOR(j, i + 1, (int)cpop_ara.size() - 1)   {
                            cpop_tmp_dd = min(cpop_tmp_dd, dis({cpop_ara[i].F, cpop_ara[i].S}, {cpop_ara[j].F, cpop_ara[j].S}));
                        }
                    }
                    ent_pressed = 1;
                }
                else if(ent_pressed == 1)   {///simulation of cpop
                    spd = mn_spd;
                }
                else if(ent_pressed == 2)   {///simulation done
                    page_no = 3;
                    spd = fst_spd;
                }
            }
            else if(p3_menu_no == fpop_idx)   {///furthest pair of points
                if(vv.size() < 3)    {
                    return;
                }
                if(ent_pressed == 0)    {
                    ///stop taking input please
                    take_inp = 0;
                    run_convex_hull();
                    ///do initialize after taking input
                    FOR(i, 0, (int)hull.size() - 1)   {
                        glVertex3f(hull[i].F, hull[i].S, 0);
                        if(hull[i].F < fpop_mn_x.F)    {
                            fpop_mn_x = hull[i];
                            fpop_mn_xpos = i;
                        }
                        if(hull[i].F > fpop_mx_x.F)    {
                            fpop_mx_x = hull[i];
                            fpop_mx_xpos = i;
                        }
                    }
                    fpop_pos1 = fpop_mn_xpos;
                    fpop_pos2 = fpop_mx_xpos;

                    fpop_brute_ans = 0;

                    FOR(i, 0, (int)hull.size() - 1)   {
                        FOR(j, i + 1, (int)hull.size() - 1)   {
                            fpop_brute_ans = max(fpop_brute_ans, dis({hull[i].F, hull[i].S}, {hull[j].F, hull[j].S}));
                        }
                    }
                    cout << "max : " << fpop_brute_ans << "\n";

                    spd = sim_spd;
                    ent_pressed = 1;
                }
                else if(ent_pressed == 1)    {
                    spd = mn_spd;
                }
                else if(ent_pressed == 2)    {
                    page_no = 3;
                    spd = fst_spd;
                }
            }
            else if(p3_menu_no == maer_idx)   {///furthest pair of points
                if(vv.size() < 4)    {
                    return;
                }
                if(ent_pressed == 0)    {
                    ///stop taking input please
                    take_inp = 0;
                    run_convex_hull();

                    ///do initialize after taking input
                    FOR(i, 0, (int)hull.size() - 1)   {
                        ///glVertex3f(hull[i].F, hull[i].S, 0);
                        if(hull[i].F < maer_mn_x.F)    {
                            maer_mn_x = hull[i];
                            maer_mn_xpos = i;
                        }
                        if(hull[i].F > maer_mx_x.F)    {
                            maer_mx_x = hull[i];
                            maer_mx_xpos = i;
                        }
                    }
                    maer_pos1 = maer_mn_xpos;
                    maer_pos2 = maer_mx_xpos;

                    FOR(i, 0, (int)hull.size() - 1)   {
                        ///glVertex3f(hull[i].F, hull[i].S, 0);
                        if(hull[i].S < maer_mn_y.S)    {
                            maer_mn_y = hull[i];
                            maer_mn_ypos = i;
                        }
                        if(hull[i].S > maer_mx_y.S)    {
                            maer_mx_y = hull[i];
                            maer_mx_ypos = i;
                        }
                    }
                    maer_pos3 = maer_mx_ypos;
                    maer_pos4 = maer_mn_ypos;
                    spd = sim_spd;
                    ent_pressed = 1;
                }
                else if(ent_pressed == 1)    {
                    page_no = 3;
                    spd = fst_spd;
                }
            }
        }
        else if(page_no == 5 || page_no == 6)   {///about page, reference page
            page_no = 2;
        }

    }
    return;
}
void specialKeys(int key, int x, int y) {
    if(key == GLUT_KEY_UP)    {///up arrow
        ///cout << "hi1\n";
        if(page_no == 2)    {
            p2_menu_no = max(p2_menu_min, p2_menu_no - 1);
        }
        if(page_no == 3)    {
            p3_menu_no = max(p3_menu_min, p3_menu_no - 1);
        }
    }
    if(key == GLUT_KEY_DOWN)    {///down arrow
        if(page_no == 2)    {
            ///cout << "hi2\n";
            p2_menu_no = min(p2_menu_max, p2_menu_no + 1);
        }
        if(page_no == 3)    {
            p3_menu_no = min(p3_menu_max, p3_menu_no + 1);
        }
    }
}


int main(int argc, char** argv) {
    /***/
    /**
    vv.PB({-0.554688, 0.075});
    vv.PB({-0.295312, 0.419444});
    vv.PB({-0.282813, -0.102778});
    */
	//Initialize GLUT
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(800, 400);

	//Create the window
	glutCreateWindow("Geometry Simulator");

	glutFullScreen();
	initializeTexture();
	//Set handler functions
	glutDisplayFunc(drawScene);

	glutSpecialFunc(specialKeys); //Special Key Handler
	glutKeyboardFunc(keyboard);   //Basic keyboard key handler
	glutMouseFunc(mouse);         //Mouse Handler

	glutTimerFunc(spd, update, 0); //Add a timer
	/// glutIdleFunc(idle_func);

	glutMainLoop();
	return 0;
}
