#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>
#include <string>
#include <sstream>

#include "cntr.hpp"
#include "read_inputfile.hpp"


using namespace std;

#define CFUNC cntr::function<double>
#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>
#define CPLX complex<double>

#define FERMION -1
#define BOSON 1
#define USE_DYSON_MAT 1


namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}


void set_particle_hole_symmetry(int &tstp1,cntr::herm_matrix<double> &A,cntr::herm_matrix<double> &B);
double check_particle_hole_symmetry(int &tstp1,double err, cntr::herm_matrix<double> &A,cntr::herm_matrix<double> &B);





////////////////////////////////////////////////
// Dyson equation for phonon Green's function //
////////////////////////////////////////////////
// <<Matsubara,timestep>>-------------------------------------
void step_D(int tstp,GREEN &D0, GREEN &D, GREEN &Pi, GREEN &D0_Pi, GREEN &Pi_D0,
double beta, double h, int SolverOrder){
//set D0_Pi=-D0*Pi
cntr::convolution_timestep(tstp,D0_Pi,D0,Pi,integration::I<double>(SolverOrder),beta,h);
D0_Pi.smul(tstp,-1.0);

cntr::convolution_timestep(tstp,Pi_D0,Pi,D0, integration::I<double>(SolverOrder),beta,h);
Pi_D0.smul(tstp,-1.0);


//solve [1-D0*Pi]*D=[1+D0_Pi]*D=D0
if(tstp==-1){
    //cntr::vie2_mat(D,D0_Pi,Pi_D0,D0,beta,SolverOrder);
     cntr::vie2_mat(D,D0_Pi,Pi_D0,D0,beta);

}else {
    cntr::vie2_timestep(tstp,D,D0_Pi,Pi_D0,D0,integration::I<double>(SolverOrder),beta,h);
}

}
// <<bootstrap>>----------------------------------------------
void start_D(GREEN &D0, GREEN &D, GREEN &Pi, GREEN &D0_Pi, GREEN &Pi_D0,
double beta, double h, int SolverOrder){
//set D0_Pi=-D0*Pi
for(int n=0; n<=SolverOrder; n++){
  cntr::convolution_timestep(n,D0_Pi,D0,Pi,integration::I<double>(SolverOrder),beta,h);
  D0_Pi.smul(n,-1.0);

  cntr::convolution_timestep(n,Pi_D0,Pi,D0,integration::I<double>(SolverOrder),beta,h);
  Pi_D0.smul(n,-1.0);
}

//solve [1-D0*Pi]*D=[1+D0_Pi]*D=D0
cntr::vie2_start(D,D0_Pi,Pi_D0,D0,integration::I<double>(SolverOrder),beta,h);
}


// << Matsubara and Timestep >> ---------------------------------------
void Sigma_Mig(int tstp, GREEN &G, GREEN &Sigma, GREEN &D0, GREEN &D,
	      GREEN &Pi, GREEN &D0_Pi, GREEN &Pi_D0, CFUNC &g_el_ph, double beta,
	      double h, int SolverOrder){
assert(G.size1()==1);
int Norb=G.size1();//Nambu index=1
int Ntau=G.ntau();

GREEN_TSTP gGg(tstp,Ntau,Norb,FERMION);
G.get_timestep(tstp,gGg);//copy time step from G
gGg.right_multiply(g_el_ph,1.0);
gGg.left_multiply(g_el_ph,1.0);

//Get Pi(t,t') = -i2g^2 G(t,t')G(t',t)

GREEN_TSTP Pol(tstp,Ntau,1,BOSON);//=i gGg(t,t') G(t',t)

Pi.set_timestep_zero(tstp);

Bubble1(tstp,Pol,0,0,gGg,0,0,G,0,0);
Pi.incr_timestep(tstp,Pol,-2.0);

//Solve the dyson equation for D = D_0 + D_0*Pi*D for the time step
step_D(tstp, D0, D, Pi, D0_Pi, Pi_D0, beta, h, SolverOrder);

//Get Sig(t,t')=ig^2 D(t,t') G(t,t') 

Bubble2(tstp,Sigma,0,0,D,0,0,gGg,0,0);
}

// << Bootstrapping >> ----------------------------------------------
void Sigma_Mig(GREEN &G, GREEN &Sigma, GREEN &D0, GREEN &D, GREEN &Pi,
	      GREEN &D0_Pi, GREEN &Pi_D0, CFUNC &g_el_ph, double beta, double h,
	      int SolverOrder){
// this is the starting version of the routine
assert(G.size1()==1);
int Norb=G.size1();//Nambu index=1
int Ntau=G.ntau();

GREEN gGg(SolverOrder,Ntau,Norb,FERMION);
GREEN Pol(SolverOrder,Ntau,1,BOSON); //=i gGg(t,t') G(t',t)

for(int n=0; n<=SolverOrder; n++){
    gGg.set_timestep(n,G);//copy time step from G
    gGg.right_multiply(n,g_el_ph);
    gGg.left_multiply(n,g_el_ph);

    //Get Pi(t,t') = -i2g^2 G(t,t')G(t',t)
    Pi.set_timestep_zero(n);

    Bubble1(n,Pol,0,0,gGg,0,0,G,0,0);
    Pi.incr_timestep(n,Pol,-2.0);
}

//Solve the dyson equiation for D=D_0+D_0*Pi*D for first time steps
start_D(D0, D, Pi, D0_Pi, Pi_D0, beta, h, SolverOrder);

//Get Sig(t,t')=ig^2 D(t,t')G(t,t')
for(int n=0; n<=SolverOrder; n++) Bubble2(n,Sigma,0,0,D,0,0,gGg,0,0);
}


template <typename T>
std::string ToString(T val)
{
    std::stringstream stream;
    stream << val;
    return stream.str();
}

double sin_square(double t,double Amp, double T,double shift){
//T is the period of the sin square function
//double shift=1000;
//double h=0.02;
double gamma=M_PI/T;
if (t>=shift and t<= (T+shift)) return Amp*sin(gamma*(t-shift))*sin(gamma*(t-shift));
else return 0.0;
}



template <class myclass> void convert_and_print(myclass& G, string name){
ostringstream os;
string filename;
os.str("");
os.clear();
os << name;
filename= os.str();
const char *filename_c = filename.c_str();
cout<<"name = "<<filename_c<<endl;
G.print_to_file(filename_c);
}



double gaussian_function(double amp, double mu, double sigma,double t){
return amp*exp(-(t-mu)*(t-mu)/(sigma*sigma));
}



double electric_field(double t,double Omega, double Amp_f, double delay, double time_duration){
        //double delay(25.0),time_duration(5.0);
        return sin(Omega*t)*gaussian_function(Amp_f,delay,time_duration,t);
}





void updateX_equilibrium(double omega0,double g_Xn,double Omega,double mixing, double& X_err,complex <double> & XA, double Amp_f, complex<double> nA, int iter){

if (iter==0) cout<<"omega0 "<<omega0<<" g_Xn "<<g_Xn<<endl;

cout << "XA before equilibrium update : " << XA<<" iter "<<iter<<endl;
//cout << "XB before equilibrium update : " << latt.X_(1,1)<<" iter "<<iter<<endl;
//cdmatrix Xold=latt.X_;
complex <double> Xold=XA;
//XA.real() = 1./omega0*(sqrt(2)*g_Xn*(1.0-2.*nA.real())+electric_field(0.0,Omega,Amp_f));
XA.real() = 1./omega0*sqrt(2.)*g_Xn*(1.0-2.*nA.real());
//latt.X_(1,1)=-latt.X_(0,0); //I impose by hand that X_B = -X_A
//latt.X_=mixing*Xold+(1.0-mixing)*latt.X_;
XA.real()=mixing*Xold.real()+(1.0-mixing)*XA.real();
//X_err=abs(latt.X_(0,0)-Xold(0,0))+abs(latt.X_(1,1)-Xold(1,1));
X_err=abs(XA.real()-Xold.real());
//cout<<" electric field f @equilibrium = "<<electric_field(0.0,Omega,Amp_f)<<endl;
//cout << "XA after equilibrium update : " << latt.X_(0,0)<<" iter "<<iter<<endl;
//cout << "XB after equilibrium update : " << latt.X_(1,1)<<" iter "<<iter<<endl;
}




cdouble fx(double omega0,cdouble pA){
	return omega0*pA.real();
}

cdouble fp(double t,cdouble xA,cdouble pA, cdouble nA,double omega0,double g_ph,double Omega,double Amp_f, double gamma_boson, double shift_f, double T_f){
	return - omega0*xA.real() - sqrt(2.)*g_ph*(2.*nA.real()-1.0) -gamma_boson*pA.real() - sqrt(2.)* electric_field(t, Omega,Amp_f,shift_f,T_f);
}

void rk4(cntr::function <double> nA, cdouble x,cdouble p, cdouble &xx, cdouble &pp,double omega0,double g_ph, double Omega, double Amp_f,int tstp,double h,int iter,double gamma_boson, double shift_f, double T_f){

	cdouble k1x(0.,0.), k2x(0.,0.), k3x(0.,0.), k4x(0.,0.), k1p(0.,0.), k2p(0.,0.), k3p(0.,0.), k4p(0.,0.), n0(0.,0.), n1(0.,0.), n2(0., 0.), n3(0.,0.), n31(0.0);

	k1x = fx(omega0,p);
	n0=nA[tstp-1].real();
	k1p = fp(h*tstp,x,p,n0,omega0,g_ph,Omega,Amp_f,gamma_boson,shift_f,T_f);

	k2x = fx(omega0,p+h*k1p*0.5);
	//n1=poly(tstp,h,nA,h*tstp+h*0.5);
	n1=(nA[tstp-1].real()+nA[tstp].real())*0.5;
	k2p = fp(h*tstp+h*0.5,x+h*k1x*0.5,p+h*k1p*0.5,n1,omega0,g_ph,Omega,Amp_f,gamma_boson, shift_f,T_f);

	k3x = fx(omega0,p+h*k2p*0.5);
	n2=n1;
	k3p = fp(h*tstp+h*0.5,x+h*k2x*0.5,p+h*k2p*0.5,n2,omega0,g_ph,Omega,Amp_f,gamma_boson,shift_f,T_f);

	k4x = fx(omega0,p+h*k3p);
	//n3=poly(tstp,h,nA,h*tstp+h);
	n31=nA[tstp].real();
	k4p = fp(h*tstp+h,x+h*k3x,p+h*k3p, n31,omega0,g_ph,Omega,Amp_f,gamma_boson, shift_f,T_f);

	xx = x + h*(k1x+2.0*k2x+2.0*k3x+k4x)/6.0;
	pp = p + h*(k1p+2.0*k2p+2.0*k3p+k4p)/6.0;

	cout<<"tstp = "<<tstp<<" iter = "<<iter<<endl;
	cout<<"n at tstp "<<tstp-1<<" = "<<n0<<endl;//<<"; n1_poly = "<<n1<<"; n poly at tstp "<<tstp<<" = "<<n3<<"; n3 from GA at tstp"<<n31<<endl;
	cout<<"X at tstp "<<tstp-1<<" = "<<x<<endl;//"; xx at tstp "<<xx<<endl;
	cout<<"P at tstp "<<tstp-1<<" = "<<p<<endl;//"; pp at tstp "<<pp<<endl;
	cout<<endl;
}


class ohmic_dos {  // Ohmic DOS, BANDWIDTH 2V
  public:
    double hi_;
    double lo_;
        ohmic_dos(){ lo_=0.001;hi_=10.0;} //questi sono gli estremi di integrazione. Consideriamo solo frequenze positive dato che DOS per bosoni e` definita solo per freq positive (questo e`anche il motivo per cui x sempre >0 quindi non devo considerare il modulo nella DOS). Poi lo spettro e` la DOS pero` estesa anche a freq negative in maniera che sia una funzione DISPARI
    double operator()(double x){
	//double omega_ph=2.0;
	//double alpha=0.5;
	double omega_ph=0.2;
	double alpha=2.*0.25;//alpha should be 1.0 (this is the value in the Boltzmann) code, but I set it to 2 since it turns out that a factor 0.5 is introduced somewhere in the green_equilibrium_boson routine. alhpa MUST BE 2 X THE VALUE IN BOLTZMANN CODE
      double arg=alpha*x/(omega_ph*omega_ph)*exp(-x/omega_ph);
      return ( arg < 0 ? 0.0 : arg);

    }
};

class fermionic_bump {  
  public:
    double hi_;
    double lo_;
        fermionic_bump(){ lo_=-5.0;hi_=5.0;}
    double operator()(double x){
	double shift=2.5, Amp=1.0/M_PI, T=4.0, gamma=M_PI/T;
      double arg, omega0(shift);
	if (x<0) {omega0=-shift;}
	if((x>=omega0-T/2.) and (x<=omega0+T/2.) )
             {arg=Amp*cos(gamma*(x-omega0))*cos(gamma*(x-omega0));}
	else arg=0.0;
      return ( arg < 0 ? 0.0 : arg);
    }
};


class dos_flat_fermionic_bath {  // Fermionic DOS, const DOS
  public:
    double hi_;
    double lo_;
        dos_flat_fermionic_bath(){lo_=-2;hi_=2;} //integro tra -2 e 2
    double operator()(double x){
        double W=8;//the BW is from -2 to 2 (Bethe Lattice), ma metto da -4 a 4 per avere flat BW tra -2 e 2
        double Gamma=0.01;//0.001;
      double arg=x;
      if( abs(x) <= W/2) arg=Gamma; //Gamma is the DOS(w) = 0.005
        else arg=0 ;
      return ( arg < 0 ? 0.0 : arg);
    }
};


class bethedos {
  public:
    /** \brief <b> Higher edge of the density of states </b> */
    double hi_;
    /** \brief <b> Lower edge of the density of states </b> */
    double lo_;
    /** \brief <b> Hopping integral and 4V corresponds to the bandwidth </b> */
    double V_;
    bethedos() {
        V_ = 1;
        lo_ = -2;
        hi_ = 2;
    }
    double operator()(double x) {
        double arg = 4.0 * V_ * V_ - x * x;
        double num = V_ * V_ * 3.14159265358979323846 * 2;
        return (arg < 0 ? 0.0 : sqrt(arg) / num);
    }
};


void green_equilibrium_mat_bethe_2(cntr::herm_matrix<double> &G,double beta){
   bethedos dos;
   green_equilibrium_mat(G,dos,beta);
}

//*****************************************************************************************************
void prepare_band1(double &h, double &b, cntr::herm_matrix<double> &A){

   int nt=A.nt();
   double mu;
   CPLX data,img;
   img.real()=0.0;
   img.imag()=1.0;
   mu=b;
   double pi = acos(-1);
   for (int i=0;i<=nt;i++){ 
     for (int j=0;j<=i;j++){
           data=-1.0*img*exp(img*mu*(i*h-j*h))*exp((-1.0*(i*h-j*h)*(i*h-j*h))/4.0*pi);
           A.set_ret(i,j,data);
       }
    }

    for (int i=0;i<=nt;i++){
     for (int j=0;j<=i;j++){
           data=img*exp(img*mu*(i*h-j*h))*exp((-1.0*(i*h-j*h)*(i*h-j*h))/4.0*pi);
           A.set_les(j,i,data);
       }
    }

}

//******************************************************************************************
void prepare_band2(double &h, double &b, cntr::herm_matrix<double> &A){

   int nt=A.nt();
   double mu;
   CPLX data,img;
   img.real()=0.0;
   img.imag()=1.0;
   mu=b;
   double pi = acos(-1);
   for (int i=0;i<=nt;i++){ 
     for (int j=0;j<=i;j++){
           data=-1.0*img*exp(img*mu*(i*h-j*h))*exp((-1.0*(i*h-j*h)*(i*h-j*h))/4.0*pi);
           A.set_ret(i,j,data);
       }
    }

    for (int i=0;i<=nt;i++){  
     for (int j=0;j<=i;j++){
           data.real()=0.0;
           data.imag()=0.0;
           A.set_les(j,i,data);
       }
    }

}

//*******************************************************************************************************

double check_particle_hole_symmetry(int &tstp1,double err, cntr::herm_matrix<double> &A,cntr::herm_matrix<double> &B){

             int ntau=A.ntau();
             int tstp=tstp1;
             cntr::herm_matrix<double> gtmp(tstp,ntau,1,-1);
             set_particle_hole_symmetry(tstp,gtmp,B);
             err=cntr::distance_norm2(tstp,gtmp,A);
             if(err==0){
               cout << "P-H symmetry satisfied " << err << endl;
             }
             else{
               cout << "P-H symmetry is not satisfied: error " << err << endl;
             }

             return err;

}

/******************************************************************************/

void set_particle_hole_symmetry(int &tstp1,cntr::herm_matrix<double> &A,cntr::herm_matrix<double> &B){
    // A is the to particle hole symmetric version of B
    int ntau=A.ntau();
    int tstp=tstp1;
    // checks ..
    if(tstp==-1){
        for(int i=0;i<=ntau;i++){ // Matsubara: G_B(tau) = G_A(\beta-tau)
            CPLX data;
            B.get_mat(-i+ntau,data);
            A.set_mat(i,data);
        }
    }
    else {
       for (int i=0;i<=tstp;i++){  // Retarted: G_B(t,t') = -cong((G_A(t,t'))) 
           CPLX data;
           B.get_ret(tstp,i,data);
           A.set_ret(tstp,i,-conj(data));
       }

       for (int i=0;i<=ntau;i++){  // tv: G_B(t,\tau) = -cong((G_A(t,\beta-tau))) 
           CPLX data;
           B.get_tv(tstp,-i+ntau,data);
           A.set_tv(tstp,i,-conj(data));
       }
       for (int i=0;i<=tstp;i++){  // less: G_B(t,t') = -((G_A(t,t')))_Retarted + cong((G_A(t,t)))_less  
           CPLX data1,data2;
//           B.get_les(i,tstp,data1);
//           B.get_ret(tstp,i,data2);
           B.get_gtr(tstp,i,data1);
           A.set_les(i,tstp,-data1);  //less: G_B(t,t') = conj(G_A(t,t'))_greater = - G_A(t',t)_greater
//           A.set_les(i,tstp,-data1+conj(data2));
       }



    }

}

/**************************************************************************/

void sigma(int tstp,cntr::herm_matrix<double> &S,cntr::herm_matrix<double> &G,cntr::function<double> &U){
        // PI(t,t')=ii*G(t,t')*G(t',t)
        cntr::herm_matrix_timestep<double> Pi(tstp,G.ntau(),1,1);
        cntr::Bubble1(tstp,Pi,G,G);
        // Sigma(t,t')=ii*Pi(t,t')*G(t,t')
        cntr::Bubble2(tstp,S,G,Pi);
        // Sigma(t,t') -> U(t)*Sigma(t,t')
        S.left_multiply(tstp,U);
        // Sigma(t,t') -> Sigma(t,t')*U(t')
        S.right_multiply(tstp,U);
        // Sigma(t,t') -> -Sigma(t,t')
        S.smul(tstp,-1.0);
//        cout << U[tstp] << endl;
}

/****************************************************************************/
int main(int argc,char *argv[]){
	cntr::herm_matrix<double> G_A,Delta_A,G_temp_A,F_A,Sigma_A,Sigma_fermionic_flat_bath,Sigma_fermionic_bump,Sigma_bath; // For G
        cntr::herm_matrix<double> G_tilda_A,G_tilda_temp_A; // For g_tilda
	cntr::herm_matrix<double> G_A_equi,Delta_A_equi,Delta_A_old,F_A_equi,Sigma_A_equi,G_tilda_A_equi;



	GREEN D,Pi,D0,D0_Pi,Pi_D0,dtemp,Sigma_GD1; //For the phonon propagator
        cntr::herm_matrix<double> g0,gF; // for solving Hartree propagators
	cntr::function<double> f,U,V_func,nA_vector,el_field_vector,f_equi; // single-time function  on  contour
	CFUNC g_elph_t,J_hop_t, zero_t;
	int nt,ntau,size,sig,tstp,kt,n_iter,previous_present,Uc1,Uc2,n_fourier,jump_in_time,param_number;
	double beta,h,err,err_ph,X_err=10.0,mixing,errmax,mu_A,D1_double,D0_double,U0,U1,A_A,e0,b;
        bool matsubara_converged=false;
        cdouble n_A,data;
        char fname1[500],fname2[500];
        int nsweep=6;
	//double bath_beta;
	//double g=sqrt(2.5*0.2);
	double T_Vt,Amp_Vt,shift_Vt,Amp_f,shift_f,T_f,omega0,g_Xn,g_ph_bath,Omega;
	complex <double> ii(0.0,1.0);
	double gamma_boson=0.001; //damping of the oscillation due to the interaction with the bath of bosons (same value as in the Boltzmann code)
	double hopping=1.0,delta_hopping=0.0;
	int Sigma_GD0_int=0, Sigma_GD1_int=0;

	         try{
        ////////////////////////////////////////////////////////////////////////////////////    
        // (II) READ INPUT
          {
                if(argc<2) throw("COMMAND LINE ARGUMENT MISSING");
                find_param(argv[1],"__nt=",nt);
                find_param(argv[1],"__ntau=",ntau);
                find_param(argv[1],"__beta=",beta);
                //find_param(argv[1],"__bath_beta=",bath_beta);
                find_param(argv[1],"__h=",h); cout<<"h= "<<h<<endl;
                find_param(argv[1],"__size=",size);
                find_param(argv[1],"__sig=",sig);
                find_param(argv[1],"__kt=",kt);cout<<"kt= "<<kt<<endl;
                find_param(argv[1],"__D1=",D1_double);
                find_param(argv[1],"__D0=",D0_double);
                find_param(argv[1],"__U0=",U0);
                find_param(argv[1],"__n_iter=",n_iter);
                find_param(argv[1],"__errmax=",errmax);
                find_param(argv[1],"__previous_present=",previous_present);
                find_param(argv[1],"__U1=",U1);
                find_param(argv[1],"__Uc1=",Uc1);
                find_param(argv[1],"__Uc2=",Uc2);
                find_param(argv[1],"__delta_hopping=",delta_hopping);
                //find_param(argv[1],"__t=",t);
                find_param(argv[1],"__b_gauss=",b);
                find_param(argv[1],"__gamma_boson=",gamma_boson);
                //find_param(argv[1],"__hopping=",hopping);
                find_param(argv[1],"__param=",param_number);
                find_param(argv[1],"__n_fourier=",n_fourier);
                find_param(argv[1],"__jump_in_time=",jump_in_time);
		find_param(argv[1],"__T_Vt=",T_Vt);cout<<"T_Vt = "<<T_Vt<<endl; //T_Vt is in hopping times
  		find_param(argv[1],"__Amp_Vt=",Amp_Vt);
  		find_param(argv[1],"__shift_Vt=",shift_Vt);cout<<"shift_Vt = "<<shift_Vt<<endl; //shift_Vt is in hopping times
		find_param(argv[1], "__g_Xn=", g_Xn);cout<<"g_Xn = "<<g_Xn<<endl;
		find_param(argv[1], "__g_ph_bath=", g_ph_bath);cout<<"g_ph_bath = "<<g_ph_bath<<endl;
		find_param(argv[1], "__omega0=", omega0); cout<<"omega0 = "<<omega0<<endl;
		find_param(argv[1], "__Omega=", Omega); cout<<"Omega = "<<Omega<<endl;
		find_param(argv[1], "__Amp_f=", Amp_f); cout<<"Amp_f = "<<Amp_f<<endl;
		find_param(argv[1], "__shift_f=", shift_f);cout<<"shift_f = "<<shift_f<<endl; //shift_f is in hopping times
		find_param(argv[1], "__T_f=", T_f);cout<<"T_f = "<<T_f<<endl; //T_f is in hopping times
		//find_param(argv[1], "__mixing=", mixing);cout<<"mixing = "<<mixing<<endl;
		find_param(argv[1], "__Sigma_GD0=", Sigma_GD0_int); cout<<"Sigma_GD0 = "<<Sigma_GD0_int<<endl;
		find_param(argv[1], "__Sigma_GD1=", Sigma_GD1_int); cout<<"Sigma_GD1 = "<<Sigma_GD1_int<<endl;
		//find_param(argv[1], "__fermionic_bump=", fermionic_bump_int); cout<<"fermionic_bump = "<<fermionic_bump_int<<endl;
		find_param(argv[1], "__mixing=", mixing);cout<<"mixing = "<<mixing<<endl;

          }

        tstp=-1;
        A_A= 0.0;
	cout<<"nt = "<<nt<<endl;
	cout<<"Amp_Vt "<<Amp_Vt<<endl;
	cout<<"T_Vt "<<T_Vt<<endl;
	cout<<"shift_Vt "<<shift_Vt<<endl; 
	cout<<"Amp_f "<<Amp_f<<endl;
	cout<<"T_f "<<T_f<<endl;
	cout<<"shift_f "<<shift_f<<endl; 
	cout<<"hopping "<<hopping<<endl;
	cout<<"gamma_boson "<<gamma_boson<<endl; 
	cout<<"delta hopping "<< delta_hopping<<endl;
	cout<<"U0 "<<U0<<endl;
	cout<<"U1 "<<U1<<endl;
	cout<<"Uc1 "<<Uc1<<endl;
	cout<<"Uc2 "<<Uc2<<endl;


	string param_num=ToString(param_number);
	string direct_out="/nfs/tfkp00/picanoan/Desktop/C++/codes/realtime/IPT_code/param_"+param_num+"/";
       string corpus = "_beta_"+patch::to_string(beta)+"_h_"+patch::to_string(h) +"_T_Vt_"+patch::to_string(T_Vt)+"_Amp_Vt_"+patch::to_string(Amp_Vt) +"_shift_Vt_"+patch::to_string(shift_Vt)+"_nt_"+patch::to_string(nt);



                /******************************************************/
                
                // Step 0: Intialization

                /******************************************************/
             {

                G_A=cntr::herm_matrix<double>(nt,ntau,size,sig);
                Delta_A=cntr::herm_matrix<double>(nt,ntau,size,sig);
                G_tilda_A=cntr::herm_matrix<double>(nt,ntau,size,sig);
                //G_tilda_B=cntr::herm_matrix<double>(nt,ntau,size,sig);
                G_temp_A=cntr::herm_matrix<double>(nt,ntau,size,sig);
                Sigma_A=cntr::herm_matrix<double>(nt,ntau,size,sig);
                //Sigma_B=cntr::herm_matrix<double>(nt,ntau,size,sig);
                Sigma_bath=cntr::herm_matrix<double>(nt,ntau,size,sig);
                Sigma_GD1=cntr::herm_matrix<double>(nt,ntau,size,sig);
                //BAND1=cntr::herm_matrix<double>(nt,ntau,size,sig);
                //BAND2=cntr::herm_matrix<double>(nt,ntau,size,sig);
                //G_new=cntr::herm_matrix<double>(-1,ntau,size,sig); 
                F_A=cntr::herm_matrix<double>(nt,ntau,size,sig);
                Sigma_fermionic_flat_bath=cntr::herm_matrix<double>(nt,ntau,size,sig); //flat fermionic bath
		Sigma_fermionic_bump=cntr::herm_matrix<double>(nt,ntau,size,sig); //bump fermionic bath

		//equilibrium GF
                G_A_equi=cntr::herm_matrix<double>(nt,ntau,size,sig);
                Delta_A_equi=cntr::herm_matrix<double>(nt,ntau,size,sig);
                Delta_A_old=cntr::herm_matrix<double>(nt,ntau,size,sig);
                G_tilda_A_equi=cntr::herm_matrix<double>(nt,ntau,size,sig);
                Sigma_A_equi=cntr::herm_matrix<double>(nt,ntau,size,sig);
                F_A_equi=cntr::herm_matrix<double>(nt,ntau,size,sig);

	


		//phonon
		D = GREEN(nt,ntau,size,BOSON);
		Pi = GREEN(nt,ntau,size,BOSON);
		D0 = GREEN(nt,ntau,size,BOSON);
		D0_Pi = GREEN(nt,ntau,size,BOSON);
		Pi_D0 = GREEN(nt,ntau,size,BOSON);



                g0=cntr::herm_matrix<double>(-1,ntau,size,sig);
                gF=cntr::herm_matrix<double>(-1,ntau,size,sig);
                //chi=cntr::herm_matrix<double>(nt,ntau,1,1);
                //CHI=cntr::herm_matrix<double>(nt,ntau,1,1);
                f=cntr::function<double>(nt,size);
                f_equi=cntr::function<double>(nt,size);
                //f1=cntr::function<double>(nt,size);
                U=cntr::function<double>(nt,size);
		V_func=cntr::function<double>(nt,size);
                //A_A1=cntr::function<double>(nt,size);
                //hyb1=cntr::function<double>(nt,size);
                //XA=cntr::function<double>(nt,size);
                //PA=cntr::function<double>(nt,size);
		nA_vector=cntr::function<double>(nt,size);
		//nB_vector=cntr::function<double>(nt,size);
		el_field_vector=cntr::function<double>(nt,size);
                g_elph_t = CFUNC(nt,size);
                J_hop_t = CFUNC(nt,size);
                zero_t = CFUNC(nt,size);


		
                //this part is used for the quench. I select a time window tstp \to [Uc1,Uc2] where U is U1; outside from this window U is U2
                for(int i=-1;i<=nt;i++){
                    if((i>=Uc1 and i<=Uc2 ))  U[i]=U1;
                      else U[i]=U0;
		 }

        // 1) I include here the flat fermionic bath
	Sigma_fermionic_flat_bath.clear();
        dos_flat_fermionic_bath dos_fermion; //fermionic bath DOS
        green_equilibrium(Sigma_fermionic_flat_bath,dos_fermion,beta,h); 
	//string Sigma_fermionic_flat_bath_str=direct_out+"Sigma_fermionic_flat_bath.out";
        //convert_and_print(Sigma_fermionic_flat_bath,Sigma_fermionic_flat_bath_str);
        //end of the fermonic bath
 

	// 2) For the coupling with the fermionic bump bath I use the CFUNC V_func	
 	for (int i=-1;i<=nt;i++){
	if(i==-1) V_func[i]=0.0;
	else	  V_func[i]=sin_square(i*h,Amp_Vt,T_Vt,shift_Vt);   //shift_Vt is in hopping times, like in the Boltzmann code
	}

        fermionic_bump dos_fermion_bump; //fermionic bath DOS
        green_equilibrium(Sigma_fermionic_bump,dos_fermion_bump, -beta,h); //with this I have D_fermion(t,t')that I can always use from now on. -beta is to have a fermi(-|beta|omega)
	//Sigma_fermionic_bump(t,t')=V(t)G_fermonic_bump(t,t')[V(t')]*	
	for(int tstp=-1;tstp<=nt;tstp++){
	Sigma_fermionic_bump.left_multiply(tstp,V_func);
	Sigma_fermionic_bump.right_multiply(tstp,V_func);
	}





	// For the excitations through the electric field
 	for (int i=-1;i<=nt;i++){
	if (i==-1) el_field_vector[i]=0.0;
	else el_field_vector[i]=electric_field(i*h,Omega,Amp_f,shift_f,T_f);
	}

	//Initialization of some excitations
	std::vector<double> dEl_Ph_g(nt+2,0.0);
	std::vector<double> dHopping(nt+2,0.0);
	double sigma_hopping=sqrt(0.16);
	double mu_hopping= 5.*sigma_hopping;
 	for (int i=0;i<=nt;i++) dHopping[i]=delta_hopping*exp(-(i*h-mu_hopping)*(i*h-mu_hopping)/(2.*sigma_hopping*sigma_hopping));
	
        cdmatrix Iu = cdmatrix::Identity(1,1);
        for(int it=-1 ; it<=nt ; it++) J_hop_t.set_value(it,hopping*(1.+dHopping[it+1])*Iu);
        for(int it=-1 ; it<=nt ; it++) 	g_elph_t.set_value(it,(g_Xn+dEl_Ph_g[it+1])*Iu);
        for(int it=-1 ; it<=nt ; it++) 	zero_t.set_value(it,0.0*Iu);


         
        // 2)  Here I include here the bosonic bath, D0 that is used to dissipate energy (I call it D0). It gives a contribution to the Sigma equal to g*g*G(t,t')*D0_ph(t,t')
	// I have two options: one is the one I have implemented
        //ohmic_dos dos; //ohmic bath DOS
        //green_equilibrium_ohmic_bath(D0,dos,beta,h); //with this I have D(t,t')that I can always use from now on; I can find this function in "./realtime/non_eq_codes/cntr/cntr_equilibrium.hpp" 
	 // the other one is the one from the library
	 // initialize D0
	 D0.clear();
         cntr::green_single_pole_XX(D0,omega0,beta,h);  // I initialize D0(t,t')/2
         for(int it=-1 ; it<=nt ; it++) D0.smul(it,2.0);

	//string D0_str = direct_out+"D0.out";
        //convert_and_print(D0,D0_str);
        //end of the bosonic bath
       
        //mu_A = -D1_double; // U/2 shift have taken in function f. So it is U/2 + D
	 mu_A=0.;
             
               cout << "Intialization done for functions" << endl;

           }  // Initialization done


              /************************************************************/

              //Step:1 Solve the equilibrium problem for tstp=-1.


             /***********************************************************/
             tstp=-1;

                if(previous_present == 0){  // this is the non-restarting procedure 

		     //cntr::green_equilibrium_bethe(G_A,beta,h);
                     green_equilibrium_mat_bethe_2(G_A,beta);
                     Delta_A.set_timestep(tstp,G_A);
		     Delta_A.right_multiply(tstp,J_hop_t);
		     Delta_A.left_multiply(tstp,J_hop_t);

                     G_temp_A.set_timestep(tstp,G_A);

                     n_A=0.50;
		     //G_A.print_to_file("G_test.out");


                     for (int i=0;i<=n_iter;i++){
		       f[tstp]=0.; //we are at half-filling
		        if (U[tstp] != 0.0) {
                        // SOLVE (id/dt + mu - Delta - S_Hartree)^{-1} = G_tilda that I need to calculate the IPT Sigma_A

                        #if USE_DYSON_MAT==1
                        // (i) DYSON MAT:
                        cntr::dyson_mat(G_tilda_A,Delta_A,mu_A,f,beta);
			//cntr::dyson_mat(G_A,F_A,mu,S_Hartree,beta);
                        #else
                        // (ii) ALTERNATIVE: solving Dyson by
                        // g0=(id/dt + mu - S_Hartree)^{-1}
                        // gF=g0*Delta
                        // (1-g0*Delta)*G_tilda=g0
                        e0=f[tstp].real();
                        cntr::green_from_eps(g0,mu_A,&e0,beta,h);
                        cntr::convolution_matsubara(gF,g0,Delta_A,integration::I<double>(kt),beta);
                        gF.smul(tstp,-1.0);
                        cntr::vie2_mat_sweep(G_tilda_A,gF,gF,g0,beta,integration::I<double>(kt),nsweep);
                         #endif // use dyson mat

			// I start constructing the self energy
                         sigma(tstp,Sigma_A,G_tilda_A,U);  //Sigma_A contains the IPT self-energy
		        }   //enf of if over U>0
		       else {
                       Sigma_A.right_multiply(tstp,zero_t);
		       }
                      //Here I add the phonon and fermion bath
		      Sigma_bath.right_multiply(tstp,zero_t); // I put Sigma_bath at the given tstp to zero
                      if (i==0) Sigma_bath.incr_timestep(tstp,Sigma_fermionic_flat_bath);

		      if (Sigma_GD0_int ==1 ){
	              cntr::Bubble2(tstp,Sigma_bath,G_A,D0); // Sigma_bath contains the self-energy due to the interaction with the equilibrium phonon bath
	              Sigma_bath.smul(tstp,g_ph_bath*g_ph_bath); //I do not need to multiply by ii since this is already in Bubble2
		      }
		      Sigma_bath.incr_timestep(tstp,Sigma_fermionic_bump); //I add to Sigma_bath the term due to the fermionic bath at negative temperature
	              // Now I calculate Sigma_GD1 due to the interaction between the electrons and the phonons 
		      if (Sigma_GD1_int ==1){
		      Sigma_Mig(tstp, G_A, Sigma_GD1, D0, D, Pi, D0_Pi, Pi_D0, g_elph_t, beta, h, kt);
                      Sigma_bath.incr_timestep(tstp,Sigma_GD1); // I add to Sigma_bath to the Sigma_GD1
		      }
                      Sigma_A.incr_timestep(tstp,Sigma_bath); // I add Sigma_bath to the IPT Sigma_A

                      F_A.set_timestep(tstp,Sigma_A); // F=Sigma on timestep tstp
                      F_A.incr_timestep(tstp,Delta_A); // F += Delta on timestep tstp


		      //------------I am now ready to solve the Dyson eqn for the full G ---------------//
                     // SOLVE (id/dt + mu - Delta -Sigma- S_Hartree)^{-1} = G
                     #if USE_DYSON_MAT==1
                      // (i) DYSON MAT:
                        cntr::dyson_mat(G_A,F_A,mu_A,f,beta);
                     #else
                        // (ii) ALTERNATIVE: solving Dyson by
                        // g0=(id/dt + mu - S_Hartree)^{-1}
                        // gF=g0*F
                        // (1-g0*F)*G=g0
                        e0=f[tstp].real();
                        cntr::green_from_eps(g0,mu_A,&e0,beta,h);
                        cntr::convolution_matsubara(gF,g0,F_A,integration::I<double>(kt),beta);
                        gF.smul(tstp,-1.0);
                        cntr::vie2_mat_sweep(G_A,gF,gF,g0,beta,integration::I<double>(kt),nsweep);
                     #endif // use dyson mat
                     n_A=G_A.density_matrix(tstp);
 		     nA_vector[tstp]=n_A;
                     Delta_A.set_timestep(tstp,G_A);
		     Delta_A.right_multiply(tstp,J_hop_t);
                     Delta_A.left_multiply(tstp,J_hop_t);

		     //cout << "iteration : " << i << endl;
                     cout<<"n_A "<<2.*n_A.real()<<endl;     

                     if(err<errmax and i>2){
                           matsubara_converged=true;
                           break;
                     }

		     //implement the mixing	
		     if(i>=1){
			G_A.smul(tstp,1.0-mixing);
			G_A.incr_timestep(tstp,G_temp_A,mixing);
	             }

		     //calculate the error
		     err=cntr::distance_norm2(tstp,G_A,G_temp_A);
                     cout << "iteration : " << i << " |G_A - G_A_temp|= " << err << endl;

		     //G_A becomes G_temp_A for the next iteration
                     G_temp_A.set_timestep(tstp,G_A);
                     
                     } //end of matsubara iteration
                
                if(!matsubara_converged){
                        std::cerr << "cannot converge matsubara; end here " << endl;
                        // should end here ....
			return 1;
                }

               
		cout << "matsubara done at tstp = " <<tstp<< endl;
                
                }
 
               else {  //previous_present != 0 branch (I guess this is a restarting routine)

                sprintf(fname1,"%s/G_A.out",argv[2]);
                sprintf(fname2,"%s/f.out",argv[2]);
                
                G_A.read_from_file(fname1);  // In principle u should read Delta and F_A 
                f.read_from_file(fname2);
                G_temp_A.set_timestep(tstp,G_A);
                Delta_A.set_timestep(tstp,G_A);
		Delta_A.right_multiply(tstp,J_hop_t);
                Delta_A.left_multiply(tstp,J_hop_t);


                cout << "reading and allocation done succesesfully from folder" << endl;
                for (int i=0;i<=n_iter;i++){

	                if (U[tstp] != 0.0)
                	{

                  #if USE_DYSON_MAT==1
                        // (i) DYSON MAT:
                         cntr::dyson_mat(G_tilda_A,Delta_A,mu_A,f,beta);
//                        cntr::dyson_mat(G_A,F_A,mu,S_Hartree,beta);
                  #else
                        // (ii) ALTERNATIVE: solving Dyson by
                        // g=(id/dt + mu - S_Hartree)^{-1}
                        // gF=g0*Delta
                        // (1-g0*Delta)*G_tilda=g0
                        e0=f[tstp].real();//((G_tilda_A.density_matrix(tstp).real()-0.5)*U[tstp].real());// f[tstp].real();
                        cntr::green_from_eps(g0,mu_A,&e0,beta,h);
                        cntr::convolution_matsubara(gF,g0,Delta_A,integration::I<double>(kt),beta);
                        gF.smul(tstp,-1.0);
                        cntr::vie2_mat_sweep(G_tilda_A,gF,gF,g0,beta,integration::I<double>(kt),nsweep);
                  #endif
                  n_A=G_A.density_matrix(tstp);
                  sigma(tstp,Sigma_A,G_tilda_A,U);
		  }  // end if cycle over U
		       else {
                       Sigma_A.right_multiply(tstp,zero_t);
		       }
		        Sigma_bath.right_multiply(tstp,zero_t); // I put Sigma_bath at the given tstp to zero
	                Sigma_bath.incr_timestep(tstp,Sigma_fermionic_flat_bath);		
			if (Sigma_GD0_int ==1){
			cntr::Bubble2(tstp,Sigma_bath,G_A,D0);
                        Sigma_bath.smul(tstp,g_ph_bath*g_ph_bath); //I do not need to multiply by ii since this is already in Bubble2
			 }
		         Sigma_bath.incr_timestep(tstp,Sigma_fermionic_bump);
			 
			 if (Sigma_GD1_int ==1){		 
                         Sigma_Mig(tstp, G_A, Sigma_GD1, D0, D, Pi, D0_Pi, Pi_D0, g_elph_t, beta, h, kt);
                         Sigma_bath.incr_timestep(tstp,Sigma_GD1); // I add to Sigma_bath to the Sigma_GD1
			 }

                         Sigma_A.incr_timestep(tstp,Sigma_bath);
			
                  F_A.set_timestep(tstp,Sigma_A);
                  F_A.incr_timestep(tstp,Delta_A);

		   #if USE_DYSON_MAT==1
                      // (i) DYSON MAT:
                        cntr::dyson_mat(G_A,F_A,mu_A,f,beta);
                   #else
                        // (ii) ALTERNATIVE: solving Dyson by
                        // g=(id/dt + mu - S_Hartree)^{-1}
                        // gF=g0*Delta
                        // (1-g0*Delta)*G_tilda=g0
                        e0=f[tstp].real();
                        cntr::green_from_eps(g0,mu_A,&e0,beta,h);
                        cntr::convolution_matsubara(gF,g0,F_A,integration::I<double>(kt),beta);
                        gF.smul(tstp,-1.0);
                        cntr::vie2_mat_sweep(G_A,gF,gF,g0,beta,integration::I<double>(kt),nsweep);
                   #endif // use dyson mat
                   Delta_A.set_timestep(tstp,G_A);
		   Delta_A.right_multiply(tstp,J_hop_t);
                   Delta_A.left_multiply(tstp,J_hop_t);

		     if(i>=1){
			G_A.smul(tstp,1.0-mixing);
			G_A.incr_timestep(tstp,G_temp_A,mixing);
	             }


                  err=cntr::distance_norm2(tstp,G_A,G_temp_A);
                  cout << "iteration : " << i << " |G_A - G_A_temp|= " << err << endl;
                  G_temp_A.set_timestep(tstp,G_A);


		  if(err<errmax && i>2){
                       matsubara_converged=true;
                       break;
                  }
		  f[tstp]=0.;
                         }

                  if(!matsubara_converged){
                        std::cerr << "cannot cnverge matsubara; end here " << endl;
                        // should end here ....
			return 1;
                  }

		 cout << "Matsubara done" << endl;
               
               }   // loop over for Matsubara
               
               /************************************************/
              
               // Step: 2 Solve the time dependent problem from nt=0...5 time steps.

               /************************************************/ 
                for(int i=0;i<=n_iter;i++){

                          cntr::dyson_start(G_A,mu_A,f,F_A,integration::I<double>(kt),beta,h);

                          for(tstp=0;tstp<=kt;tstp++) {
            	           Delta_A.set_timestep(tstp,G_A);
                           Delta_A.right_multiply(tstp,J_hop_t);
                           Delta_A.left_multiply(tstp,J_hop_t);
                           nA_vector[tstp]=G_A.density_matrix(tstp).real();
			   f[tstp]=0.;
  			  }
         

                        err=0.0;
                        for(tstp=0;tstp<=kt;tstp++) err += cntr::distance_norm2(tstp,G_A,G_temp_A);
                        cout << "START iteration : " << i << " |G_A - G_A_temp|= " << err << endl;
                        if(err<errmax and i>2){
                                matsubara_converged=true;
                                break;
                        }
                  
                        cntr::dyson_start(G_tilda_A,mu_A,f,Delta_A,integration::I<double>(kt),beta,h);           
                    
                        for(tstp=0;tstp<=kt;tstp++){

		           if (U[tstp] != 0.) {sigma(tstp,Sigma_A,G_tilda_A,U);}
			   else { Sigma_A.right_multiply(tstp,zero_t);}

                           //************ Here I could add the Sigma due to the baths****// 
			   Sigma_bath.right_multiply(tstp,zero_t);
 	                   //Sigma_bath.incr_timestep(tstp,Sigma_fermionic_flat_bath);
			   if (Sigma_GD0_int ==1 ){                              
                           cntr::Bubble2(tstp,Sigma_bath,G_A,D0);
                           Sigma_bath.smul(tstp,g_ph_bath*g_ph_bath); //I do not need to multiply by ii since this is already in Bubble2             
			    }
			   if (Sigma_GD1_int == 1){
			   Sigma_Mig(G_A, Sigma_GD1, D0, D, Pi, D0_Pi, Pi_D0, g_elph_t, beta, h, kt);
                           Sigma_bath.incr_timestep(tstp,Sigma_GD1); // I add to Sigma_bath to the Sigma_GD1
			   }
                           Sigma_bath.incr_timestep(tstp,Sigma_fermionic_bump);
                           Sigma_A.incr_timestep(tstp,Sigma_bath);

                           F_A.set_timestep(tstp,Sigma_A); 
                           F_A.incr_timestep(tstp,Delta_A);

                        }  //end of for cycle over tstp on bootstrap

                        for(tstp=0;tstp<=kt;tstp++){
                           G_temp_A.set_timestep(tstp,G_A);
                        }       
                    } //end of cycle over DMFT-iter on bootstrap

		   for (tstp=0;tstp<=kt;tstp++){
		     cout<<"------- tstp = "<<tstp<<" n_A "<<2.*nA_vector[tstp].real()<<"-------------"<<endl;     
		     }

                /************************************************************************/
 
                // Step 3: Solve the time dependent problem from kt+1 to nt time steps. 

               /*************************************************************************/

		//new part

		//I start the time propagation from the same inital condition for each DMFT iteration
		//I put the equilibrium results into the equilibrium GF
		for (int tstp=-1;tstp<=kt;tstp++){
			G_A_equi.set_timestep(tstp,G_A);
			Sigma_A_equi.set_timestep(tstp,Sigma_A);
			Delta_A_equi.set_timestep(tstp,Delta_A);
			//Delta_A_old.set_timestep(tstp,Delta_A);
			F_A_equi.set_timestep(tstp,F_A);
			f_equi[tstp]=f[tstp];
			} 			

		//I propagate at equilibrium in order to get Delta_A_equilibrium for each timestep
                for(tstp=kt+1;tstp<=nt;tstp++){
                        // extrapoaltion of Delta by one timestep:
                        cntr::extrapolate_timestep(tstp-1,Delta_A_equi,integration::I<double>(kt));
                        cntr::extrapolate_timestep(tstp-1,Sigma_A,integration::I<double>(kt));                        
                        // iteration on the timestep, fixed iteration number, no error check!
                        f[tstp]=f[tstp-1]; //guess value
                        cout << " ... at timestep " << tstp <<endl;

			for(int i=0;i<=n_iter;i++){

				cout<<"iter "<<i<<endl;
                                G_temp_A.set_timestep(tstp,G_A);
                                F_A.set_timestep(tstp,Sigma_A); // F=Sigma on timestep tstp
                                F_A.incr_timestep(tstp,Delta_A_equi);
                                cntr::dyson_timestep(tstp,G_A,mu_A,f,F_A,integration::I<double>(kt),beta,h);
                                n_A=G_A.density_matrix(tstp);
			        nA_vector[tstp]=n_A;
				f[tstp]=0.;
                                Delta_A_equi.set_timestep(tstp,G_A);
				Delta_A_equi.right_multiply(tstp,J_hop_t);
                                Delta_A_equi.left_multiply(tstp,J_hop_t);
				if (U[tstp] != 0. ){
                                cntr::dyson_timestep(tstp,G_tilda_A,mu_A,f,Delta_A_equi,integration::I<double>(kt),beta,h);
                                sigma(tstp,Sigma_A,G_tilda_A,U);
				}
				else { Sigma_A.right_multiply(tstp,zero_t);} 

				err = cntr::distance_norm2(tstp,G_A,G_temp_A);
				cout << "iteration : " << i << " |G_A - G_A_temp|= " << err << endl;
				if(err<errmax and i>2 ){
				matsubara_converged=true;
				break;
				}
				// G_temp_A.set_timestep(tstp,G_A);
			
                        } //end of cycle over DMFT iterations
			
		     cout<<"------- equilibrium tstp = "<<tstp<<" n_A "<<2.*nA_vector[tstp].real()<<"-------------"<<endl;  
                 }   //end of for cycle over all the tstp
		//end of propagation at equilibrium--- now I have Delta_A_equi for each t,t'

		//print the total energy for iter 0
                string energy_name=direct_out+"temp_total_energy"+corpus+"_iter_"+patch::to_string(0)+".out";
                FILE *fout2=fopen(energy_name.c_str(),"w");
                cdouble x_A,y_A;
                cdouble occu_t_A;
                kt=5;
                for(tstp=-1;tstp<=nt;tstp++){
                occu_t_A=G_A.density_matrix(tstp);
                cntr::convolution_density_matrix(tstp,&x_A,Sigma_A,G_A,integration::I<double>(kt),beta,h);
                cntr::convolution_density_matrix(tstp,&y_A,Delta_A_equi,G_A,integration::I<double>(kt),beta,h);
                fprintf(fout2,"%i ",tstp);
                fprintf(fout2, "%.13g ",occu_t_A.real());
		fprintf(fout2,"%.13g ",x_A.real());
		fprintf(fout2,"%.13g ",y_A.real());
		//fprintf(fout2,"%.13g ",U[tstp].real()*occu_t_A.real()*occu_t_A.real());
		//fprintf(fout2,"%.13g ",U[tstp].real());
                //fprintf(fout2, "%.13g ",XA[tstp].real());
                //fprintf(fout2, "%.13g ",PA[tstp].real());
                fprintf(fout2, "%.13g ",V_func[tstp].real()); //coupling with the fermionic bump bath
                //fprintf(fout2, "%.13g ",el_field_vector[tstp].real()); //electric field coupled to the system
                fprintf(fout2,"\n");
                }
                fclose(fout2);



		//I start the DMFT iterations
		for(int i=1;i<=n_iter;i++){

			 cout<<endl;
			 cout<<"iter "<<i<<endl;

			//For each GF, I start the time prop from the equilibrium values
			for (int tstp=-1;tstp<=kt;tstp++){
				G_A.set_timestep(tstp,G_A_equi);
				Sigma_A.set_timestep(tstp,Sigma_A_equi);
				Delta_A.set_timestep(tstp,Delta_A_equi);
				//Delta_A_old.set_timestep(tstp,Delta_A_equi);
				F_A.set_timestep(tstp,F_A_equi);
				f[tstp]=f_equi[tstp];
				} 		

			for (int tstp=kt+1;tstp<=nt;tstp++){
				//if (i==1) Delta_A.set_timestep(tstp,Delta_A_equi);
				//else Delta_A.set_timestep(tstp,Delta_A_old);
				Delta_A.set_timestep(tstp,G_A);
				Delta_A.right_multiply(tstp,J_hop_t);
				Delta_A.left_multiply(tstp,J_hop_t);
				} 		




			for(int tstp=kt+1;tstp<=nt;tstp++){
				// extrapoaltion of Delta by one timestep:
				//cntr::extrapolate_timestep(tstp-1,Delta_A,integration::I<double>(kt));
				// extrapoaltion of Sigma by one timestep:
				cntr::extrapolate_timestep(tstp-1,Sigma_A,integration::I<double>(kt));
				
				// iteration on the timestep, fixed iteration number, no error check!
				//f[tstp]=f[tstp-1]; //guess value
				f[tstp]=0.;
				//cout << " ... at timestep " << tstp << endl;

				G_temp_A.set_timestep(tstp,G_A);
				F_A.set_timestep(tstp,Sigma_A); // F=Sigma on timestep tstp
				F_A.incr_timestep(tstp,Delta_A);
				cntr::dyson_timestep(tstp,G_A,mu_A,f,F_A,integration::I<double>(kt),beta,h);
				n_A=G_A.density_matrix(tstp);
				nA_vector[tstp]=n_A;
				f[tstp]=0.;
				//Delta_A.set_timestep(tstp,G_A);
				//Delta_A.right_multiply(tstp,J_hop_t);
				//Delta_A.left_multiply(tstp,J_hop_t);
				if (U[tstp] != 0. ){
				cntr::dyson_timestep(tstp,G_tilda_A,mu_A,f,Delta_A,integration::I<double>(kt),beta,h);
				sigma(tstp,Sigma_A,G_tilda_A,U);
				}
				else { Sigma_A.right_multiply(tstp,zero_t);} 

			    //************ Here I add the Sigma due to the baths****// 
			   Sigma_bath.right_multiply(tstp,zero_t);
			   //Sigma_bath.incr_timestep(tstp,Sigma_fermionic_flat_bath);
			   if (Sigma_GD0_int ==1){
			   cntr::Bubble2(tstp,Sigma_bath,G_A,D0);
			   Sigma_bath.smul(tstp,g_ph_bath*g_ph_bath);
			   }
			   if (Sigma_GD1_int ==1){
			   Sigma_Mig(tstp, G_A, Sigma_GD1, D0, D, Pi, D0_Pi, Pi_D0, g_elph_t, beta, h, kt);
			   Sigma_bath.incr_timestep(tstp,Sigma_GD1); // I add to Sigma_bath to the Sigma_GD1
			   }
			   Sigma_bath.incr_timestep(tstp,Sigma_fermionic_bump);
			   Sigma_A.incr_timestep(tstp,Sigma_bath);

			   //*************End of part where I add Sigma due to the baths*********//
			   cout<<"------- tstp = "<<tstp<<" iter = "<<i<<" n_A "<<2.*nA_vector[tstp].real()<<"-----------"<<endl;      

			//update Delta...
			//Delta_A.set_timestep(tstp,G_A);
			//Delta_A.right_multiply(tstp,J_hop_t);
			//Delta_A.left_multiply(tstp,J_hop_t);
			//... and put the value in Delta_old
			//Delta_A_old.set_timestep(tstp,Delta_A);

			 } //end of for cycle over all the tstp



			//when I reach convergence, the G at iteration n is the same of the G at n-1, for each timestep
	                err=0.0;
                        for(tstp=kt+1;tstp<=nt;tstp++) err += cntr::distance_norm2(tstp,G_A,G_temp_A);
                        cout << "|G_A - G_A_temp|= " << err << endl;
                        if(err<errmax and i>2){
                                matsubara_converged=true;
                                break;
                        }
        

		//print the total energy for each iteration
                string energy_name=direct_out+"temp_total_energy"+corpus+"_iter_"+patch::to_string(i)+".out";
                FILE *fout2=fopen(energy_name.c_str(),"w");
                cdouble x_A,y_A;
                cdouble occu_t_A;
                kt=5;
                for(tstp=-1;tstp<=nt;tstp++){
                occu_t_A=G_A.density_matrix(tstp);
                cntr::convolution_density_matrix(tstp,&x_A,Sigma_A,G_A,integration::I<double>(kt),beta,h);
                cntr::convolution_density_matrix(tstp,&y_A,Delta_A,G_A,integration::I<double>(kt),beta,h);
                fprintf(fout2,"%i ",tstp);
                fprintf(fout2, "%.13g ",occu_t_A.real());
		fprintf(fout2,"%.13g ",x_A.real());
		fprintf(fout2,"%.13g ",y_A.real());
		//fprintf(fout2,"%.13g ",U[tstp].real()*occu_t_A.real()*occu_t_A.real());
		//fprintf(fout2,"%.13g ",U[tstp].real());
                //fprintf(fout2, "%.13g ",XA[tstp].real());
                //fprintf(fout2, "%.13g ",PA[tstp].real());
                fprintf(fout2, "%.13g ",V_func[tstp].real()); //coupling with the fermionic bump bath
                //fprintf(fout2, "%.13g ",el_field_vector[tstp].real()); //electric field coupled to the system
                fprintf(fout2,"\n");
                }
                fclose(fout2);


		}//end of loop over DMFT-iter
		//end new part



		//old code
		/*
                for(tstp=kt+1;tstp<=nt;tstp++){
                        // extrapoaltion of Delta by one timestep:
                        cntr::extrapolate_timestep(tstp-1,Delta_A,integration::I<double>(kt));
                        cntr::extrapolate_timestep(tstp-1,Sigma_A,integration::I<double>(kt));
                        
                        // iteration on the timestep, fixed iteration number, no error check!
                        f[tstp]=f[tstp-1]; //guess value
                        cout << " ... at timestep " << tstp << endl;

			for(int i=0;i<=n_iter;i++){
                                G_temp_A.set_timestep(tstp,G_A);
                                F_A.set_timestep(tstp,Sigma_A); // F=Sigma on timestep tstp
                                F_A.incr_timestep(tstp,Delta_A);
                                cntr::dyson_timestep(tstp,G_A,mu_A,f,F_A,integration::I<double>(kt),beta,h);
                                n_A=G_A.density_matrix(tstp);
			        nA_vector[tstp]=n_A;
				f[tstp]=0.;
                                Delta_A.set_timestep(tstp,G_A);
				Delta_A.right_multiply(tstp,J_hop_t);
                                Delta_A.left_multiply(tstp,J_hop_t);
				if (U[tstp] != 0. ){
                                cntr::dyson_timestep(tstp,G_tilda_A,mu_A,f,Delta_A,integration::I<double>(kt),beta,h);
                                sigma(tstp,Sigma_A,G_tilda_A,U);
				}
				else { Sigma_A.right_multiply(tstp,zero_t);} 

                 	Sigma_bath.right_multiply(tstp,zero_t);
                           //Sigma_bath.incr_timestep(tstp,Sigma_fermionic_flat_bath);
			    if (Sigma_GD0_int ==1) {
                           cntr::Bubble2(tstp,Sigma_bath,G_A,D0);
                           Sigma_bath.smul(tstp,g_ph_bath*g_ph_bath);
				}

			    if (Sigma_GD1_int ==1) {
                           Sigma_Mig(tstp, G_A, Sigma_GD1, D0, D, Pi, D0_Pi, Pi_D0, g_elph_t, beta, h, kt);
                           Sigma_bath.incr_timestep(tstp,Sigma_GD1); // I add to Sigma_bath to the Sigma_GD1
			   }
			Sigma_bath.incr_timestep(tstp,Sigma_fermionic_bump);

			Sigma_A.incr_timestep(tstp,Sigma_bath);
			err = cntr::distance_norm2(tstp,G_A,G_temp_A);
			cout << "iteration : " << i << " |G_A - G_A_temp|= " << err << endl;
				if(err<errmax and i>2 ){
				matsubara_converged=true;
				break;
				}
				// G_temp_A.set_timestep(tstp,G_A);
			
                        } //end of cycle over DMFT iterations

		     cout<<"------- tstp = "<<tstp<<" --------"<<endl;
                     cout<<"n_A "<<2.*nA_vector[tstp].real()<<endl;      
                 }   //end of for cycle over all the tstp
		*/

		   //string G_A_str=direct_out+"G_A"+corpus+".out";
		   string G_A_str=direct_out+"G_A.out";
		   string Sigma_A_str=direct_out+"Sigma_A.out";
		   string f_str=direct_out+"f.out";
		   string U_str=direct_out+"U.out";
		   string Delta_A_str=direct_out+"Delta_A.out";
		   string Sigma_fermionic_bump_str=direct_out+"Sigma_fermionic_bump.out";
		   string D0_str=direct_out+"D0.out";
		   string Sigma_bath_str=direct_out+"Sigma_bath.out";

                   convert_and_print(G_A,G_A_str);
                   //convert_and_print(D0,D0_str);
                   //convert_and_print(Sigma_bath,Sigma_bath_str);
                   //convert_and_print(Sigma_fermionic_bump,Sigma_fermionic_bump_str);

		   cout<<"g = "<<g_Xn<<endl;

//                   G_B.print_to_file("G_B.out");
//                   G_tilda_A.print_to_file("G_tilda_A.out");
                   convert_and_print(Sigma_A,Sigma_A_str);
                   //convert_and_print(f,f_str);
                   //convert_and_print(U,U_str);
                   //convert_and_print(Delta_A,Delta_A_str);
		   //convert_and_print(Sigma_fermionic_bump,Sigma_fermionic_bump_str);
//             } // stop the tstp>5 calculations.

             } // try
        catch(char *message){
            cerr << "exception\n**** " << message << " ****" << endl;
                cerr << "input_file [ --test ]\n" << endl;
        }
        catch(...){
            cerr << "unspecified exception " << endl;
                cerr << "\n input_file [ --test ]\n" << endl;
        }




             return 0;

              }


              

