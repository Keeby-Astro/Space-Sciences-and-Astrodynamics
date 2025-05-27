#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

/*
    This program presents the graphical solution of the motion of three
    bodies in the plane for data provided in the input definitions below.

    MATLAB's ode45 Runge-Kutta solver is used.

    G                           - gravitational constant (km^3/kg/s^2)
    t0, tf                      - initial and final times (s)
    m1, m2, m3                  - masses of the three bodies (kg)
    m                           - total mass (kg)
    X1,Y1; X2,Y2; X3,Y3         - coordinates of the three masses (km)
    VX1,VY1; VX2,VY2; VX3,VY3   - velocity components of the three
                                    masses (km/s)
    XG, YG                      - coordinates of the center of mass (km)
    y0                          - column vector of the initial conditions
    t                           - column vector of times at which the solution
                                    was computed
    y                           - matrix, the columns of which contain the
                                    position and velocity components evaluated at
                                    the times t(:):
                                    y(:,1) , y(:, 2) = X1(:), Y1(:)
                                    y(:,3) , y(:, 4) = X2(:), Y2(:)
                                    y(:,5) , y(:, 6) = X3(:), Y3(:)
                                    y(:,7) , y(:, 8) = VX1(:), VY1(:)
                                    y(:,9) , y(:,10) = VX2(:), VY2(:)
                                    y(:,11), y(:,12) = VX3(:), VY3(:)
    
    User h-functions required: none
    User subfunctions required: rates, plotit
*/
static const double G  = 6.67259e-20;

//...Input data:
static const double m1 = 1.e29;
static const double m2 = 1.e29;
static const double m3 = 1.e29;
static const double m  = m1 + m2 + m3;

static const double t0 = 0.0;     // initial time (s)
static const double tf = 67000.0; // final time (s)

static const double X1_init  = 0.0;      
static const double Y1_init  = 0.0;       
static const double X2_init  = 300000.0;  
static const double Y2_init  = 0.0;       
static const double X3_init  = 2.0*X2_init; 
static const double Y3_init  = 0.0;       

static const double VX1_init = 0.0;   
static const double VY1_init = 0.0;   
static const double VX2_init = 250.0; 
static const double VY2_init = 250.0; 
static const double VX3_init = 0.0;   
static const double VY3_init = 0.0;
//...End input data   

std::vector<double> rates(double t, const std::vector<double> &y);
void plotit(const std::vector<double> &tspan, 
            const std::vector<std::vector<double>> &Y);

//...Pass the initial conditions and time interval to ode45, which
// calculates the position and velocity of each particle at discrete
// times t, returning the solution in the column vector y. ode45 uses
// the subfunction 'rates' below to evaluate the accelerations at each
// integration time step.
void rk4_solve(std::vector<double> &tspan, std::vector<std::vector<double>> &Y, double dt)
{
    size_t N = static_cast<size_t>((tf - t0) / dt) + 1;
    tspan.resize(N);
    Y.resize(N, std::vector<double>(12, 0.0));

    // Initialize time
    tspan[0] = t0;

    // Initialize state vector
    Y[0][0]  = X1_init;
    Y[0][1]  = Y1_init;
    Y[0][2]  = X2_init;
    Y[0][3]  = Y2_init;
    Y[0][4]  = X3_init;
    Y[0][5]  = Y3_init;
    Y[0][6]  = VX1_init;
    Y[0][7]  = VY1_init;
    Y[0][8]  = VX2_init;
    Y[0][9]  = VY2_init;
    Y[0][10] = VX3_init;
    Y[0][11] = VY3_init;

    // Loop for each step
    for (size_t i = 0; i < N - 1; ++i)
    {
        double current_t = tspan[i];
        const std::vector<double> &current_y = Y[i];

        std::vector<double> k1 = rates(current_t, current_y);

        std::vector<double> yMid1(12, 0.0);
        for (int j = 0; j < 12; j++) {
            yMid1[j] = current_y[j] + 0.5*dt*k1[j];
        }
        std::vector<double> k2 = rates(current_t + 0.5*dt, yMid1);

        std::vector<double> yMid2(12, 0.0);
        for (int j = 0; j < 12; j++) {
            yMid2[j] = current_y[j] + 0.5*dt*k2[j];
        }
        std::vector<double> k3 = rates(current_t + 0.5*dt, yMid2);

        std::vector<double> yEnd(12, 0.0);
        for (int j = 0; j < 12; j++) {
            yEnd[j] = current_y[j] + dt*k3[j];
        }
        std::vector<double> k4 = rates(current_t + dt, yEnd);

        std::vector<double> next_y(12, 0.0);
        for (int j = 0; j < 12; j++) {
            next_y[j] = current_y[j] + (dt/6.0)*(k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j]);
        }

        Y[i+1] = next_y;
        tspan[i+1] = current_t + dt;
    }
}

int main()
{
    double dt = 10.0; 

    std::vector<double> t;
    std::vector<std::vector<double>> Y;

    // Integrate the system
    rk4_solve(t, Y, dt);

    plotit(t, Y);

    return 0;
}

std::vector<double> rates(double /*t*/, const std::vector<double> &y)
{
    double X1  = y[0];
    double Y1  = y[1];
    double X2  = y[2];
    double Y2  = y[3];
    double X3  = y[4];
    double Y3  = y[5];
    double VX1 = y[6];
    double VY1 = y[7];
    double VX2 = y[8];
    double VY2 = y[9];
    double VX3 = y[10];
    double VY3 = y[11];

    //...Equations C.8:
    double R12 = std::pow(std::hypot(X2 - X1, Y2 - Y1), 3);
    double R13 = std::pow(std::hypot(X3 - X1, Y3 - Y1), 3);
    double R23 = std::pow(std::hypot(X3 - X2, Y3 - Y2), 3);

    // Equations C.7:
    double AX1 = G*m2*(X2 - X1)/R12 + G*m3*(X3 - X1)/R13;
    double AY1 = G*m2*(Y2 - Y1)/R12 + G*m3*(Y3 - Y1)/R13;

    double AX2 = G*m1*(X1 - X2)/R12 + G*m3*(X3 - X2)/R23;
    double AY2 = G*m1*(Y1 - Y2)/R12 + G*m3*(Y3 - Y2)/R23;

    double AX3 = G*m1*(X1 - X3)/R13 + G*m2*(X2 - X3)/R23;
    double AY3 = G*m1*(Y1 - Y3)/R13 + G*m2*(Y2 - Y3)/R23;

    std::vector<double> dydt(12, 0.0);
    dydt[0]  = VX1;
    dydt[1]  = VY1;
    dydt[2]  = VX2;
    dydt[3]  = VY2;
    dydt[4]  = VX3;
    dydt[5]  = VY3;
    dydt[6]  = AX1;
    dydt[7]  = AY1;
    dydt[8]  = AX2;
    dydt[9]  = AY2;
    dydt[10] = AX3;
    dydt[11] = AY3;

    return dydt;
}

void plotit(const std::vector<double> &tspan, 
            const std::vector<std::vector<double>> &Y)
{
    size_t lastIndex = tspan.size() - 1;
    double X1 = Y[lastIndex][0];
    double Y1 = Y[lastIndex][1];
    double X2 = Y[lastIndex][2];
    double Y2 = Y[lastIndex][3];
    double X3 = Y[lastIndex][4];
    double Y3 = Y[lastIndex][5];

    double VX1 = Y[lastIndex][6];
    double VY1 = Y[lastIndex][7];
    double VX2 = Y[lastIndex][8];
    double VY2 = Y[lastIndex][9];
    double VX3 = Y[lastIndex][10];
    double VY3 = Y[lastIndex][11];

    // Locate the center of mass at the final time step:
    double XG = (m1*X1 + m2*X2 + m3*X3) / m;
    double YG = (m1*Y1 + m2*Y2 + m3*Y3) / m;

    // Print results
    std::cout << std::fixed << std::setprecision(3);

    std::cout << "\n--------------------------------------------------\n";
    std::cout << "Final time: " << tspan[lastIndex] << " s\n";

    std::cout << "\n-- Final positions (inertial frame) --\n";
    std::cout << "Body 1: (X1, Y1) = (" << X1 << ", " << Y1 << ") km\n";
    std::cout << "Body 2: (X2, Y2) = (" << X2 << ", " << Y2 << ") km\n";
    std::cout << "Body 3: (X3, Y3) = (" << X3 << ", " << Y3 << ") km\n";

    std::cout << "\n-- Final velocities (inertial frame) --\n";
    std::cout << "Body 1: (VX1, VY1) = (" << VX1 << ", " << VY1 << ") km/s\n";
    std::cout << "Body 2: (VX2, VY2) = (" << VX2 << ", " << VY2 << ") km/s\n";
    std::cout << "Body 3: (VX3, VY3) = (" << VX3 << ", " << VY3 << ") km/s\n";

    std::cout << "\n-- Center of mass --\n";
    std::cout << "XG = " << XG << " km\n";
    std::cout << "YG = " << YG << " km\n";

    // Coordinates of each particle relative to the center of mass at the final time:
    double X1G = X1 - XG; 
    double Y1G = Y1 - YG;
    double X2G = X2 - XG; 
    double Y2G = Y2 - YG;
    double X3G = X3 - XG; 
    double Y3G = Y3 - YG;

    std::cout << "\n-- Positions relative to center of mass --\n";
    std::cout << "Body 1: (X1G, Y1G) = (" << X1G << ", " << Y1G << ") km\n";
    std::cout << "Body 2: (X2G, Y2G) = (" << X2G << ", " << Y2G << ") km\n";
    std::cout << "Body 3: (X3G, Y3G) = (" << X3G << ", " << Y3G << ") km\n";

    std::cout << "--------------------------------------------------\n\n";
}
