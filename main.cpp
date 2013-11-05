#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#define pi atan(1.0)*4

using namespace std;

/*Programmet:
 *Klassen Planet inneholder info om planeter
 *Klassen Solver inneholder runge-kutta algoritmen
 *main() oppretter array av planetobjekter og sender denne til
 *konstruktøren i Solverklassen.
 *Solverklassen har hovdeløkken i konstruktøren, hvor den kaller
 *rk4() som igjen kaller f() for k-verdiene.
 *De deriverte blir beregnet i f().
 */


class Planet {
    double x0, y0, vx0, vy0; //Initial coordinates and velocities
    double mass;//kg
    double distance; //Distance from origo(sun) in astrominical units
    char* name;//Name of planet

public:
        Planet(double mass, double distance, char* name, double velocity){
        this->mass = mass;
        this->name = name;
        this->x0 = distance;
        this->y0 = 0.0;
        this->vx0 = 0.0;
        this->vy0 = velocity;//kmt_to_au(velocity);//average orbital velocity in au/yr
    }
    double get_x0() {
        return x0;
    }
    double get_y0() {
        return y0;
    }
    double get_vx0() {
        return vx0;
    }
    double get_vy0() {
        return vy0;
    }
    double get_mass() {
        return mass;
    }

    /*Parameter: velocity in km/h
     *Return: velocity in au/y*
     */
    /*double kmt_to_au(double v){
        double au = 149.6*pow(10,6);
        double y = 3600*24*365;
        return (double) v*y/au;
    }*/
};

/*Denne klassen inneholder algoritmen/løsningen av runge kutta metoden.
 *Konstruktøren tar som input: n antall planeter, peker til array av planetobjekter.
 *
 *Hovedølkken:
 *Oppretter en array av lengde 4*antall planter
 *Initialiserer denne: q = {x0,y0,vx0,vy0,x1,y1,vx1,vy1,...}
 *Lager en while-løkke som løper så lenge siste tidssteg ikke er nådd.
 *Kaller rk4-metoden.
 *Øker tidssteget og gjentar prosedyren.
 *
 *RK5 ALGORITMEN (Generelt):
 *
 *k1=hf(t_i,y_i)                  Økningen basert på begynnelsen av intervallet, ved dy/dt (Euler's metode)
 *k2=hf(t_i+h/2,y_i+k1/2)         Økningen basert på midtpunktet av intervallet, ved dy/dt+1/2hk1
 *k3=hf(t_i+h/2,y_i+k2/2)         Økningen basert på midtpunktet av intervallet, ved dy/dt+1/2hk2
 *k4=hf(t_i+h,y_i+k3)             Økningen basert på enden av intervallet, ved dy/dt+hk3
 *y_i+1 = y_i+1/6(k1+2k2+2k3+k4)
 *
 *f
 */
class Solver {
    int n;
    double h, t_i, t_f, x_0, y_0,E_0;
    double *q;
    ofstream file;
public:
    Solver(int n,Planet **p) {
        this->n = n;
        t_i = 0;             //t=0
        t_f = 40;//20        //t=julian year
        h = (t_f-t_i)/1000;  //steplength
        q = new double[4*n];
        init_q(p,n);         //Putting the initial values into vector q
        while(t_i<t_f){
        rk4(p);              //Runge-Kutta
        t_i +=h;
        }
        file.close();
       }
public:
    int init_q(Planet **p,int n) {
        int j = 0;
        for(int i = 0; i<4*n; i=i+4) {
            q[i] = p[j]->get_x0();
            q[i+1] = p[j]->get_y0();
            q[i+2] = p[j]->get_vx0();
            q[i+3] = p[j]->get_vy0();
            j++;
        }
        return 0;
    }
    /*ALGORITMEN rk4:
     *Oppretter peker til 5 arrayer k1,k2,k3,k4,temp
     *setter k1 = returarrayen til metoden f
     *setter k2 = returarrayen til metoden f, men jobber på temp[] som inneholder q+h*k1/2
     *setter k3 = returarrayen til metoden f, men jobber på temp[] som inneholder q+h*k2/2
     *setter k4 = returarrayen til metoden f, men jobber på temp[] som inneholder q+hk3
     *
     *løper gjennom hele q[] og setter q[i] = q[i]+h/6*(k1[i]+2k2[i]+2k3[i]+k4[i]
     */
public:
    int rk4(Planet **p) {
        double *k1, *k2, *k3, *k4,*temp;
        temp = new double[4*n];

        k1 = f(q,p);
        for(int i = 0; i<4*n; i++){
            temp[i] = q[i] + k1[i]*h/ 2.0;
        }
        k2 = f(temp,p);
        for(int i = 0; i<4*n; i++){
            temp[i] = q[i] + k2[i]*h/ 2.0;
        }
        k3 = f(temp,p);
        for(int i = 0; i<4*n; i++) {
            temp[i] = q[i] + k3[i]*h;
        }
        k4 = f(temp,p);

        for(int i = 0; i<4*n; i++) {
            q[i] =q[i]+h/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
        }

        for(int i = 0; i<4*n; i+=4) {
            cout<<q[i+0]<<" "<<q[i+1]<<" "<<q[i+2]<<" "<<q[i+3]<<endl;
        }
        //Skriver data til fil
        atf(q,n);
        return 0;
    }
    /*Oppretter plass i minnet for k-arrayene
     *setter k0 = dx/dt, k1 = dy/dt.
     *Hvis i = 0, har vi med solen å gjøre. Da settes disse til 0,
     *siden vi antar solen står i ro.
     *Så vil vi sette k2 = dvx/dt, k3 = dvy/dt
     *Løper også gjennom planetene og summerer akselerasjonen
     *ved å beregne kreftene fra hver planet. Men hopper over den planeten vi jobber på. (j!=i)
     *Så oppdaterer vi i-løkken, og begynner på de 4 neste veridene.
     */
    double* f(double* q,Planet **p){
        double* k = new double[4*n];
        double r = 0;

        for (int i = 0; i < 4*n; i+=4){
            k[i+0] = q[i+2];
            k[i+1] = q[i+3];
            if(i == 0) {
                k[i+2] = 0.0;
                k[i+3] = 0.0;
            }
            else {
                k[i+2] = 0.0;
                k[i+3] = 0.0;
                for(int j = 0; j<4*n; j+=4){
                    if(i!= j){
                        // Kalkulerer kraften F fra j på i
                        double x =q[i]-q[j];
                        double y = q[i+1]-q[j+1];
                        //r = sqrt(q[i+0]*q[i+0]+q[i+1]*q[i+1]);
                        r = sqrt(x*x + y*y);
                        k[i+2] -= 4*pi*pi*p[j/4]->get_mass()/pow(r,3)*x;
                        k[i+3] -= 4*pi*pi*p[j/4]->get_mass()/pow(r,3)*y;
                    }
                }
            }
        }
        return k;
    }

    /*For debugging*/
    int print_array(double *temp, int n){
        cout<<"Array: "<<endl;
        for(int i = 0; i<4*n; i++){
            cout<<temp[i]<<endl;
        }
        return 0;
    }

   /*Array to file
    *Skriver arrayen til fil. Merk: Banen må oppdateres for kjøring på annen maskin.
    */
    int atf(double *temp,int n) {
        if(!file.is_open())
            //file.open("/uio/hume/student-u32/andehus/Project3/data.dat",ios::trunc);
            file.open("/home/anders/Git/glowing-meme/data.dat",ios::trunc);
        for(int i = 4; i<4*n; i+=4) {
            file<<temp[i]<<" "<<temp[i+1]<<" "<<temp[i+2]<<" "<<temp[i+3]<<" ";
        }
        file<<endl;
        return 0;
    }

};

/*Lager objekter av klassen planet. n er antall planeter dvs. lengden på Planet arrayen.
 *Massen er i solmasser og hastigheten i AU³/yr²
 *Konstruktøren tar: Planet(double mass, double distance, char* name, double velocity_y0)
 *For å kjøre for flere planeter: Kommenter vekk planeten, evt. oppdater index,
 *og oppdater n til antall planeter.
 */
int main()
{
    int n = 2;
    Planet **p = new Planet*[n];
    p[0] = new Planet(1.0,0,"sun",0);
    p[1] = new Planet(3.0*pow(10,-6),1,"earth",2*pi);//sqrt(2)*1.999*pi forsvinner
    //p[2] = new Planet(1000.0*pow(10,-3), 5.20, "jupiter",0.5*pi);//5.20*11.86*2*pi
    //p[2] = new Planet(9.5*pow(10,-4), 5.20, "jupiter",5.20*1.09*2*pi);
    //p[3] = new Planet(1.66*pow(10,-7),0.39, "mercury",0.39*0.32*2*pi);
    //p[4] = new Planet(2.66*pow(10,-6), 0.72, "venus",0.72*1.59*2*pi); //solar masses velocity= 0.75 AU*1.59Julian year*2*pi
    //p[5] = new Planet(3.3*pow(10,-7),1.52, "mars",1.52*1.88*2*pi);
    //p[6] = new Planet(2.85*pow(10,-4), 9.54,"saturn",9.54*29.46*2*pi);
    //p[7] = new Planet(4.35*pow(10,-5), 19.19,"uranus",19.19*1.01*2*pi);
    //p[8] = new Planet(5.0*pow(10,-5),30.06,"neptune",30.06*1.0*2*pi);
    //p[9] = new Planet(6.5*pow(10,-9),39.53, "pluto",39.53*1.0*2*pi);
    Solver* s = new Solver(n,p);
    return 0;
}

