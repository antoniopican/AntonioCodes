//Author: Antonio Picano

#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <vector>
#include <cmath>
#include <sstream>
#include <omp.h>
#include "./find_param.h"
#include "./ness_decls.hpp"
#include <random>
#include <chrono>
#include <sstream>

using namespace std;
using namespace ness;



double gaussian_function(double amp, double mu, double sigma,double t){
return amp*exp(-(t-mu)*(t-mu)/(sigma*sigma));
}

double electric_field(double t,double Omega,double Amp_f,double time_duration,double delay){
	return gaussian_function(Amp_f,delay,time_duration,t);

}  

double sin_square(double t,double Amp, double T,double shift);


void print_to_file(vector <int>& time,vector <double>& occ_numberA,vector <double>& half_kin_energyA,vector <double> & pot_energyA,vector <double>& occ_numberB,vector <double>& half_kin_energyB,vector <double> & pot_energyB,vector <double>& f,vector <double>& XA , vector <double> & PA, vector <double>& XB, vector <double> & PB,std::string filename,int interval,int it, int & first_time, vector <double> gamma, vector <double> delta_omega, vector <double>& dW_A, vector <double> & dW_B, vector <double> & XXA, vector <double> & XPA,  vector <double> & PPA,  vector <double> & XXB, vector <double> & XPB,  vector <double> & PPB,  vector <double> & XXXA,  vector <double> & XXXB){
	int start,stop;
	std::ofstream outputFile;
	outputFile.open(filename.c_str(),std::ios_base::app); //open in append mode
	cout<<"**************************"<<endl;
	cout<<" I am inside print_to_file; it = "<<it<<"; interval =  "<<interval<<"; first_time = "<<first_time<<endl;
	cout<<"**************************"<<endl;
	if(first_time==1){
	 outputFile << "# time" << "\t" << "nA" << "\t"<< "half_KinA" <<"\t"<<"PotA" << "\t" << "nB" << "\t"<< "half_KinB" <<"\t"<<"PotB"<<"\t"<<"el_field"<<"\t"<<"XA"<<"\t"<<"PA"<<"\t"<<"XB"<<"\t"<<"PB"<<"\t"<<"gamma"<<"\t"<<"domega"<<"\t"<<"dW_A"<<"\t"<<"dW_B"<<"\t"<<"XX_A"<<"\t"<<"XP_A"<<"\t"<<"PP_A"<<"\t"<<"XX_B"<<"\t"<<"XP_B"<<"\t"<<"PP_B"<<"\t"<<"XXX_A"<<"\t"<<"XXX_B"<<std::endl;
	start=0;
	stop=interval;
	}
	else {start=it-interval;stop=it;}
	for (int i = start; i<stop; i++)
	{
		outputFile <<std::setprecision(5) << time[i] << "\t" << 2.*occ_numberA[i] << "\t"<<half_kin_energyA[i] <<"\t"<<pot_energyA[i]<<"\t" << 2.*occ_numberB[i] << "\t"<< half_kin_energyB[i] <<"\t"<<pot_energyB[i]<<"\t"<<f[i]<<"\t"<<XA[i]<<"\t"<<PA[i]<<"\t"<<XB[i]<<"\t"<<PB[i]<<"\t"<<gamma[i]<<"\t"<<delta_omega[i]<<"\t"<<dW_A[i]<<"\t"<<dW_B[i]<< "\t"<<XXA[i]<<"\t"<<XPA[i]<<"\t"<<PPA[i]<<"\t"<<XXB[i]<<"\t"<<XPB[i]<<"\t"<<PPB[i]<<"\t"<<XXXA[i]<<"\t"<<XXXB[i]<<std::endl;
	}	
	outputFile.close();
	first_time=0;
}
vector<double> get_gamma_and_delta_omega(GF & Polar, double g_ph,double omega0){
	vector <double> result(2,0.0);
	cplx deriv_Polariz_Retard_in_zero(0.0,0.0);
	deriv_Polariz_Retard_in_zero=1./(12.*Polar.dgrid_)*(-Polar.Retarded[2](0,0)+8.*Polar.Retarded[1](0,0)-8.*Polar.Retarded[Polar.ngrid_-1](0,0)+Polar.Retarded[Polar.ngrid_-2](0,0));
	result[0]=2.*g_ph*g_ph*omega0*deriv_Polariz_Retard_in_zero.imag();
	//result[1]=-2.*g_ph*g_ph*Polar.Retarded[0](0,0).real();	
	result[1]=0.;	
	return result;
}


void print_to_file_single_traj(vector <int>& time,int traj, vector<vector<double>>& occ_numberA, vector <vector <double>>& half_kin_energyA, vector<vector <double>> & pot_energyA, vector<vector <double>>& occ_numberB, vector<vector <double>>& half_kin_energyB, vector <vector<double>> & pot_energyB, vector <double>& f,vector <vector <double>> & XA , vector<vector <double>> & PA, vector <vector <double>> & XB, vector <vector <double>> & PB, std::string filename,int interval,int it, int & first_time, vector <vector <double>> & gamma, vector <vector <double>> & delta_omega, vector <vector <double>> & dW_A, vector <vector <double>> & dW_B) {
	int start,stop;
	std::ofstream outputFile;
	outputFile.open(filename.c_str(),std::ios_base::app); //open in append mode
	cout<<"**************************"<<endl;
	cout<<" I am inside print_to_file; it = "<<it<<"; interval =  "<<interval<<"; first_time = "<<first_time<<endl;
	cout<<"**************************"<<endl;
	if(first_time==1){
	 outputFile << "# time" << "\t" << "nA" << "\t"<< "half_KinA" <<"\t"<<"PotA" << "\t" << "nB" << "\t"<< "half_KinB" <<"\t"<<"PotB"<<"\t"<<"el_field"<<"\t"<<"XA"<<"\t"<<"PA"<<"\t"<<"XB"<<"\t"<<"PB"<<"\t"<<"gamma"<<"\t"<<"domega"<<"\t"<<"dW_A"<<"\t"<<"dW_B"<<std::endl;
	start=0;
	stop=interval;
	}
	else {start=it-interval+1;stop=it;}
	for (int i = start; i<=stop; i++)
	{
		outputFile <<std::setprecision(5) << time[i] << "\t" << 2.*occ_numberA[traj][i] << "\t"<<half_kin_energyA[traj][i] <<"\t"<<pot_energyA[traj][i]<<"\t" << 2.*occ_numberB[traj][i] << "\t"<< half_kin_energyB[traj][i] <<"\t"<<pot_energyB[traj][i]<<"\t"<<f[i]<<"\t"<<XA[traj][i]<<"\t"<<PA[traj][i]<<"\t"<<XB[traj][i]<<"\t"<<PB[traj][i]<<"\t"<<gamma[traj][i]<<"\t"<<delta_omega[traj][i]<<"\t"<<dW_A[traj][i]<<"\t"<<dW_B[traj][i]<<std::endl;
	}	
	outputFile.close();
	first_time=0;
}

void print_to_file_single_traj_restart(vector <int>& time,int traj, vector<vector<double>>& occ_numberA, vector <vector <double>>& half_kin_energyA, vector<vector <double>> & pot_energyA, vector<vector <double>>& occ_numberB, vector<vector <double>>& half_kin_energyB, vector <vector<double>> & pot_energyB, vector <double>& f,vector <vector <double>> & XA , vector<vector <double>> & PA, vector <vector <double>> & XB, vector <vector <double>> & PB, std::string filename,int interval,int it, vector <vector <double>> & gamma, vector <vector <double>> & delta_omega, vector <vector <double>> & dW_A, vector <vector <double>> & dW_B, int time_restart) {
	//int start,stop;
	std::ofstream outputFile;
	outputFile.open(filename.c_str(),std::ios_base::app); //open in append mode
	cout<<"**************************"<<endl;
	cout<<" I am inside print_to_file; it = "<<it<<"; interval =  "<<interval<<endl;
	cout<<"**************************"<<endl;

	for (int i = it-interval+1; i <= it; i++)
	{
		outputFile <<std::setprecision(5) << time[i-time_restart-1] << "\t" << 2.*occ_numberA[traj][i-time_restart-1] << "\t"<<half_kin_energyA[traj][i-time_restart-1] <<"\t"<<pot_energyA[traj][i-time_restart-1]<<"\t" << 2.*occ_numberB[traj][i-time_restart-1] << "\t"<< half_kin_energyB[traj][i-time_restart-1] <<"\t"<<pot_energyB[traj][i-time_restart-1]<<"\t"<<f[i-time_restart-1]<<"\t"<<XA[traj][i-time_restart-1]<<"\t"<<PA[traj][i-time_restart-1]<<"\t"<<XB[traj][i-time_restart-1]<<"\t"<<PB[traj][i-time_restart-1]<<"\t"<<gamma[traj][i-time_restart-1]<<"\t"<<delta_omega[traj][i-time_restart-1]<<"\t"<<dW_A[traj][i-time_restart-1]<<"\t"<<dW_B[traj][i-time_restart-1]<<std::endl;
	}	
	outputFile.close();

}



void mixing_G(GF & Gloc,GF & G_old,double mixing){	
 for(int w=0;w<Gloc.ngrid_;w++){
 Gloc.Retarded[w].noalias()=Gloc.Retarded[w]*(1-mixing)+G_old.Retarded[w]*mixing;
 Gloc.Lesser[w].noalias()=Gloc.Lesser[w]*(1-mixing)+G_old.Lesser[w]*mixing;
}
}

template<class Matrix>
void DensityMatrix(const GF &G, Matrix &M){
    cplx ii(0.0,1.0);
    M.resize(G.size1_,G.size2_);
    M.setZero();
	for(long i=0; i<G.ngrid_; i++) M+=G.Lesser[i];
	M=M*G.dgrid_/(2.0*M_PI*ii);
}


//to be adapted to the 2x2 case
template<class Matrix>
void KinEnergy(GF G, GF Delta, Matrix &M){
    cplx minus_ii(0.0,-1.0);
    M.resize(G.size1_,G.size2_);
    M.setZero();
    GF Delta_G=G;
    Delta_G.clear();
    Delta=G;// in bethe lattice Delta=G; if I do not impose this condition, Delta and G differ a bit due to the fact that G was updated at the end of the last dmft iter and then the kinetic energy is not perfectly real
    Delta_G.left_multiply(Delta,1.0);
    //cdmatrix integral_Gret(1,1);
    //integral_Gret.setZero();
    for(long w=0; w<G.ngrid_; w++) {
        Delta_G.Lesser[w]= Delta.Retarded[w]*G.Lesser[w] + Delta.Lesser[w] * G.Retarded[w].adjoint();
	Delta_G.Retarded[w] =Delta.Retarded[w] * G.Retarded[w]; //even if not needed
	M+=Delta_G.Lesser[w];
	}
    M=M*G.dgrid_*minus_ii/(2.0*M_PI);
}

//to be adapted to the 2x2 case
template<class Matrix>
void PotEnergy(GF G, GF Sigma, Matrix &M){
    cplx minus_ii(0.0,-1.0);
    M.resize(G.size1_,G.size2_);
    M.setZero();
    GF Sigma_G=G;
    Sigma_G.clear();
    Sigma_G.left_multiply(Sigma,1.0);
    for(long w=0; w<G.ngrid_; w++) {
        Sigma_G.Lesser[w]= Sigma.Retarded[w]*G.Lesser[w] + Sigma.Lesser[w] * G.Retarded[w].adjoint();
	Sigma_G.Retarded[w] =Sigma.Retarded[w] * G.Retarded[w]; //even if not needed
	M+=Sigma_G.Lesser[w];
	}
    M=M*G.dgrid_*minus_ii/(2.0*M_PI);
}



// for debug:
void outputGF_label(GF &G,std::string name1,int idx,std::string name2){
    ostringstream os;
    string filename;
    os.str("");
    os.clear();
    os << name1 << idx << name2;
    filename= os.str();
    outputGF_new(G,filename);
}

void set_matrixelement(GF &Gout,int i1,int i2,GF &Gin, int j1,int j2){
    for(int i=0;i<Gin.ngrid_;i++){
        Gout.Lesser[i](i1,i2)=Gin.Lesser[i](j1,j2);
        Gout.Retarded[i](i1,i2)=Gin.Retarded[i](j1,j2);
    }
}


void G_into_Delta(GF & Gloc, GF & Delta){
    GF GA(Gloc.dgrid_,Gloc.ngrid_,1);
    GF GB(Gloc.dgrid_,Gloc.ngrid_,1);
    set_matrixelement(GA,0,0,Gloc,0,0);
    set_matrixelement(GB,0,0,Gloc,1,1);

    GF DeltaA(Gloc.dgrid_,Gloc.ngrid_,1);
    GF DeltaB(Gloc.dgrid_,Gloc.ngrid_,1);
 
    for(long w = 0; w < Gloc.ngrid_; w++){
	DeltaA.Retarded[w]=GB.Retarded[w];
	DeltaA.Lesser[w]=GB.Lesser[w];
	DeltaB.Retarded[w]=GA.Retarded[w];
	DeltaB.Lesser[w]=GA.Lesser[w];
    	}

    Delta.clear();
    set_matrixelement(Delta,0,0,DeltaA,0,0);
    set_matrixelement(Delta,1,1,DeltaB,0,0);
    }



int TestPositivity(GF &G){
cout<<"Sono in Test Positivity"<<endl;
	int el_size = G.el_size_;
	for(long i=0; i<G.ngrid_; i++){

		if(G.Retarded[i](0,0).imag() > 0.0 || -G.Lesser[i](0,0).imag()/(2.0*G.Retarded[i](0,0).imag())>1.0+1E-12){
                        printf("error: %.15g %.15g\n", -G.Lesser[i](0,0).imag()/(2.0*G.Retarded[i](0,0).imag()), G.Retarded[i](0,0).imag());
                        return i;
                }
	
	}
	return -1;
}



void updateSE(double U, GF &G, GF &locSE, fftw_complex *input_GL, fftw_complex *input_GG,
                fftw_complex *output_GL, fftw_complex *output_GG, fftw_plan &plan_GL, fftw_plan &plan_GG)
{

        cout<<"I enter into updateSE!!!"<<endl;
        cdmatrix zero2cd=cdmatrix::Zero(2,2);
        long N=G.ngrid_;
        long Nft=3*N;
        int Nhalf = N / 2;
        cplx imag_one(0,1);
        int test=-1;

        if((test=TestPositivity(G))>=0){
                std::cout << "G < 0 at f=" << test << std::endl;
        }
         cout<<"I enter into updateSE!!!"<<endl;
        std::vector< cdmatrix> GL(Nft,zero2cd);
        std::vector< cdmatrix > GG(Nft,zero2cd);
        cdmatrix tempr, templ;
         cout<<"I enter into updateSE!!!"<<endl;
        for(int row=0; row<2; row++){
                for(int col=0; col<2; col++){
                        int ind = row * 2 + col;
                        for(long i=0; i<Nhalf; i++){
                                input_GL[i][0]=0.0;
                                input_GG[i][0]=0.0;
                                input_GL[i][1] = G.Lesser[i](row,col).imag();
                                input_GG[i][1] = G.Lesser[i](row,col).imag() + 2.0 * G.Retarded[i](row,col).imag();
                        }
                        for(long i=Nhalf; i<2*N+Nhalf; i++){
                                input_GL[i][0]=0.0;
                                input_GL[i][1]=0.0;
                                input_GG[i][0]=0.0;
                                input_GG[i][1]=0.0; //G>=G<+2i Im(GR)
                        }

                        for(long i=2*N + Nhalf; i<Nft; i++){
                                input_GL[i][0]=0.0;
                                input_GG[i][0]=0.0;

                                input_GL[i][1]=G.Lesser[i - 2*N](row,col).imag();
                                input_GG[i][1]=(G.Lesser[i - 2*N](row,col).imag()+2.0*G.Retarded[i - 2*N](row,col).imag()); //G>=G<+2i Im(GR)
                        }

                        fftw_execute(plan_GL);
                        fftw_execute(plan_GG);
                        for(long i=0; i<Nft; i++){
                                GL[i](row,col)=output_GL[i][0]+imag_one*output_GL[i][1];
                                GG[i](row,col)=output_GG[i][0]+imag_one*output_GG[i][1];
                        }
                }
        }

 cout<<"I enter into updateSE!!!"<<endl;

        double Usq=U*U;
        double w=Usq/(double) (Nft) * G.dgrid_ * G.dgrid_/ (4.0 * M_PI * M_PI); //Normierung
        double sign = -1.0;
        double cross_sign = 1.0;

        //diagonal terms
        //(0,0)
        for(long i=0; i<Nft; i++){
                //GL=U*GL*GL*GG
                //input_GG for GL because of reversed FT
                input_GG[i][0]=w* ( GL[i](0,0)*GL[i](1,1)*GG[i](1,1)).real();
               input_GG[i][1]=w* ( GL[i](0,0)*GL[i](1,1)*GG[i](1,1)).imag();
                input_GL[i][0]=w* ( GG[i](0,0)*GG[i](1,1)*GL[i](1,1)).real();
                input_GL[i][1]=w* ( GG[i](0,0)*GG[i](1,1)*GL[i](1,1)).imag();

                input_GG[i][0]+=sign * w* ( GL[i](0,1)*GG[i](1,1)*GL[i](1,0)).real();
                input_GG[i][1]+=sign * w* ( GL[i](0,1)*GG[i](1,1)*GL[i](1,0)).imag();
                input_GL[i][0]+=sign * w* ( GG[i](0,1)*GL[i](1,1)*GG[i](1,0)).real();
                input_GL[i][1]+=sign * w* ( GG[i](0,1)*GL[i](1,1)*GG[i](1,0)).imag();

        }

        //print2fftw_cvectors(Nft, input_GL,input_GG, "tlocSE.out");
        fftw_execute(plan_GL);
        fftw_execute(plan_GG);

        for(long i=0; i<N/2; i++){
                //*(locSE.p_les(N/2 + i, 0, 0)) =  (output_GG[i][0]+imag_one*output_GG[i][1]);
                *(locSE.p_les(i, 0, 0)) =  (output_GG[i][0]+imag_one*output_GG[i][1]);
        }
        for(long i=2*N+N/2; i<N*3; i++){
                //*(locSE.p_les(i - 2*N - N/2, 0, 0)) =  (output_GG[i][0]+imag_one*output_GG[i][1]);
                *(locSE.p_les(i - 2*N, 0, 0)) =  (output_GG[i][0]+imag_one*output_GG[i][1]);
        }
        for(long i=0; i<Nft; i++){
                double theta = 1.0;
                if (i > N+N/2) theta = 0.0;
                if (i == 0)
                {
                        theta = 0.5;
                        input_GL[i][0] = (input_GL[i][0] - input_GG[i][0]) * theta;
                        input_GL[i][1] = (input_GL[i][1] - input_GG[i][1]) * theta;
                }
                else
                {
                        input_GL[i][0] = (input_GL[i][0] - input_GG[Nft - i][0]) * theta;
                        input_GL[i][1] = (input_GL[i][1] - input_GG[Nft - i][1]) * theta;
                }

        }
        fftw_execute(plan_GL);
        for(long i=0; i<N/2; i++){
                *(locSE.p_ret(i, 0, 0)) = output_GL[i][0] + imag_one * output_GL[i][1];
                *(locSE.p_ret(i + N/2, 0, 0)) = output_GL[i + 2*N + N/2][0] + imag_one * output_GL[i + 2*N + N/2][1];

        }
 cout<<"I enter into updateSE!!!"<<endl;
        //(1,1)
        for(long i=0; i<Nft; i++){
                //GL=U*GL*GL*GG
                //input_GG for GL because of reversed FT
                input_GG[i][0]=w*( GL[i](1,1)*GL[i](0,0)*GG[i](0,0)).real();
                input_GG[i][1]=w*( GL[i](1,1)*GL[i](0,0)*GG[i](0,0)).imag();
                input_GL[i][0]=w*( GG[i](1,1)*GG[i](0,0)*GL[i](0,0)).real();
                input_GL[i][1]=w*( GG[i](1,1)*GG[i](0,0)*GL[i](0,0)).imag();

                input_GG[i][0]+=sign * w*( GL[i](1,0)*GG[i](0,0)*GL[i](0,1)).real();
                input_GG[i][1]+=sign * w*( GL[i](1,0)*GG[i](0,0)*GL[i](0,1)).imag();
                input_GL[i][0]+=sign * w*( GG[i](1,0)*GL[i](0,0)*GG[i](0,1)).real();
                input_GL[i][1]+=sign * w*( GG[i](1,0)*GL[i](0,0)*GG[i](0,1)).imag();
        }

        fftw_execute(plan_GL);
       fftw_execute(plan_GG);
 cout<<"I enter into updateSE!!!"<<endl;
        for(long i=0; i<N/2; i++){
                //*(locSE.p_les(N/2+i,1,1))=(output_GG[i][0]+imag_one*output_GG[i][1]);
                *(locSE.p_les(i,1,1))=(output_GG[i][0]+imag_one*output_GG[i][1]);

        }
        for(long i=2*N+N/2; i<N*3; i++){
                //*(locSE.p_les(i-2*N-N/2,1,1))=(output_GG[i][0]+imag_one*output_GG[i][1]);
                *(locSE.p_les(i-2*N,1,1))=(output_GG[i][0]+imag_one*output_GG[i][1]);
        }
        for(long i=0; i<Nft; i++){
                double theta = 1.0;
                //if (i < N + N/2) theta = 0.0;
                if (i > N+N/2) theta = 0.0;
                if (i == 0)
                {
                        theta = 0.5;
                        input_GL[i][0] = (input_GL[i][0] - input_GG[i][0]) * theta;
                        input_GL[i][1] = (input_GL[i][1] - input_GG[i][1]) * theta;
                }
                else
                {
                        input_GL[i][0] = (input_GL[i][0] - input_GG[Nft - i][0]) * theta;
                        input_GL[i][1] = (input_GL[i][1] - input_GG[Nft - i][1]) * theta;
                }
        }
        fftw_execute(plan_GL);
        for(long i=0; i<N/2; i++){
                //*(locSE.p_ret(N/2+i,1,1))=output_GL[i][0] + imag_one * output_GL[i][1]; // GR=0.5*i(GG-GL)
                *(locSE.p_ret(i,1,1))=output_GL[i][0] + imag_one * output_GL[i][1]; // GR=0.5*i(GG-GL)
                //*(locSE.p_ret(i,1,1))=output_GL[i + 2*N + N/2][0] + imag_one * output_GL[i + 2*N + N/2][1]; // GR=0.5*i(GG-GL)
                *(locSE.p_ret(i+N/2,1,1))=output_GL[i + 2*N + N/2][0] + imag_one * output_GL[i + 2*N + N/2][1]; // GR=0.5*i(GG-GL)
        }
 cout<<"I enter into updateSE!!!"<<endl;
        //(0,1)
        for(long i=0; i<Nft; i++){
                //GL=U*GL*GL*GG
                //input_GG for GL because of reversed FT
                input_GG[i][0]=w* ( GL[i](0,1)*GL[i](1,0)*GG[i](0,1)).real() * cross_sign;
                input_GG[i][1]=w* ( GL[i](0,1)*GL[i](1,0)*GG[i](0,1)).imag() * cross_sign;
                input_GL[i][0]=w* ( GG[i](0,1)*GG[i](1,0)*GL[i](0,1)).real() * cross_sign;
                input_GL[i][1]=w* ( GG[i](0,1)*GG[i](1,0)*GL[i](0,1)).imag() * cross_sign;

                input_GG[i][0]+=sign * w* ( GL[i](0,0)*GG[i](0,1)*GL[i](1,1)).real() * cross_sign;
                input_GG[i][1]+=sign * w* ( GL[i](0,0)*GG[i](0,1)*GL[i](1,1)).imag() * cross_sign;
                input_GL[i][0]+=sign * w* ( GG[i](0,0)*GL[i](0,1)*GG[i](1,1)).real() * cross_sign;
                input_GL[i][1]+=sign * w* ( GG[i](0,0)*GL[i](0,1)*GG[i](1,1)).imag() * cross_sign;
        }

        fftw_execute(plan_GL);
        fftw_execute(plan_GG);

        for(long i=0; i<N/2; i++){
                //*(locSE.p_les(N/2+i,0,1))=(output_GG[i][0]+imag_one*output_GG[i][1]);
                *(locSE.p_les(i,0,1))=(output_GG[i][0]+imag_one*output_GG[i][1]);
        }
        for(long i=2*N+N/2; i<N*3; i++){
                //*(locSE.p_les(i-2*N-N/2,0,1))=(output_GG[i][0]+imag_one*output_GG[i][1]);
                *(locSE.p_les(i-2*N,0,1))=(output_GG[i][0]+imag_one*output_GG[i][1]);
        }
        for(long i=0; i<Nft; i++){
                double theta = 1.0;
                //if (i < N + N/2) theta = 0.0;
                if (i > N+N/2) theta = 0.0;
                if (i == 0)
               {
                        theta = 0.5;
                        input_GL[i][0] = (input_GL[i][0] - input_GG[i][0]) * theta;
                        input_GL[i][1] = (input_GL[i][1] - input_GG[i][1]) * theta;
                }
                else
                {
                        input_GL[i][0] = (input_GL[i][0] - input_GG[Nft - i][0]) * theta;
                        input_GL[i][1] = (input_GL[i][1] - input_GG[Nft - i][1]) * theta;
                }

        }
        fftw_execute(plan_GL);
        for(long i=0; i<N/2; i++){
                //*(locSE.p_ret(N/2+i,0,1))=output_GL[i][0] + imag_one * output_GL[i][1]; // GR=0.5*i(GG-GL)
                *(locSE.p_ret(i,0,1))=output_GL[i][0] + imag_one * output_GL[i][1]; // GR=0.5*i(GG-GL)
                //*(locSE.p_ret(i,0,1))=output_GL[i + 2*N + N/2][0] + imag_one * output_GL[i + 2*N + N/2][1]; // GR=0.5*i(GG-GL)
                *(locSE.p_ret(i+N/2,0,1))=output_GL[i + 2*N][0] + imag_one * output_GL[i + 2*N][1]; // GR=0.5*i(GG-GL)
        }
        //(1,0)
        for(long i=0; i<Nft; i++){
                //GL=U*GL*GL*GG
                //input_GG for GL because of reversed FT
                input_GG[i][0]=w*( GL[i](1,0)*GL[i](0,1)*GG[i](1,0)).real() * cross_sign;
                input_GG[i][1]=w*( GL[i](1,0)*GL[i](0,1)*GG[i](1,0)).imag() * cross_sign;
                input_GL[i][0]=w*( GG[i](1,0)*GG[i](0,1)*GL[i](1,0)).real() * cross_sign;
                input_GL[i][1]=w*( GG[i](1,0)*GG[i](0,1)*GL[i](1,0)).imag() * cross_sign;

                input_GG[i][0]+=sign * w*( GL[i](1,1)*GG[i](1,0)*GL[i](0,0)).real() * cross_sign;
                input_GG[i][1]+=sign * w*( GL[i](1,1)*GG[i](1,0)*GL[i](0,0)).imag() * cross_sign;
                input_GL[i][0]+=sign * w*( GG[i](1,1)*GL[i](1,0)*GG[i](0,0)).real() * cross_sign;
                input_GL[i][1]+=sign * w*( GG[i](1,1)*GL[i](1,0)*GG[i](0,0)).imag() * cross_sign;
        }

        fftw_execute(plan_GL);
        fftw_execute(plan_GG);

        for(long i=0; i<N/2; i++){
                //*(locSE.p_les(N/2+i,1,0))=(output_GG[i][0]+imag_one*output_GG[i][1]);
                *(locSE.p_les(i,1,0))=(output_GG[i][0]+imag_one*output_GG[i][1]);
        }
        for(long i=2*N+N/2; i<N*3; i++){
                //*(locSE.p_les(i - 2*N - N/2,1,0))=(output_GG[i][0]+imag_one*output_GG[i][1]);
                *(locSE.p_les(i - 2*N,1,0))=(output_GG[i][0]+imag_one*output_GG[i][1]);
        }
        for(long i=0; i<Nft; i++){
                double theta = 1.0;
                //if (i < N + N/2) theta = 0.0;
                if (i > N+N/2) theta = 0.0;
                if (i == 0)
                {
                        theta = 0.5;
                        input_GL[i][0] = (input_GL[i][0] - input_GG[i][0]) * theta;
                        input_GL[i][1] = (input_GL[i][1] - input_GG[i][1]) * theta;
                }
                else
                {
                        input_GL[i][0] = (input_GL[i][0] - input_GG[Nft - i][0]) * theta;
                        input_GL[i][1] = (input_GL[i][1] - input_GG[Nft - i][1]) * theta;
                }
        }
        fftw_execute(plan_GL);
        for(long i=0; i<N/2; i++){
                //*(locSE.p_ret(N/2+i,1,0))=output_GL[i][0] + imag_one * output_GL[i][1]; // GR=0.5*i(GG-GL)
                *(locSE.p_ret(i,1,0))=output_GL[i][0] + imag_one * output_GL[i][1]; // GR=0.5*i(GG-GL)
                //*(locSE.p_ret(i,1,0))=output_GL[i + 2*N + N/2][0] + imag_one * output_GL[i + 2*N + N/2][1]; // GR=0.5*i(GG-GL)
                *(locSE.p_ret(i+N/2,1,0))=output_GL[i + 2*N][0] + imag_one * output_GL[i + 2*N][1]; // GR=0.5*i(GG-GL)
       }
        if((test=TestPositivity(locSE))>=0){
                std::cout << "locSE < 0 at f=" << test << std::endl;
        }

}


void ForceImag(GF &G){
	cdmatrix tmp;
    for(long i=0; i< G.ngrid_; i++){
        tmp=G.Lesser[i];
        G.Lesser[i]=0.5*(tmp-tmp.adjoint());
	}
}


double fermi(double w,double beta){
    double arg=w*beta;
    if(arg*arg>10000.0){
        return (w>0 ? 0.0 : 1.0);
    }else{
        return 1.0/(exp(arg)+1.0);
    }
}

double fermi_negative_beta(double w,double beta){
    double arg=w*beta;
    if(arg*arg>10000.0){
        return (w>0 ? 1.0 : 0.0);
    }else{
        return 1.0/(exp(-arg)+1.0);
    }
}


double fermi_modified(double w,double beta, double fcutoff){
    double arg=w*beta;
    double T1 = fcutoff/5.0;      //this is the period of the function sin^2 that I put on top of the equilibrium Fermi-dirac distribution. Negative frequencies 
    double T2 = fcutoff/5.0;           //positive frequencies
    double A = 1.0;            //this is the amplitude of the  Fermi-dirac distribution
    double gamma1 = M_PI * (1.0/T1); 
    double gamma2 = M_PI * (1.0/T2);
    double shift = fcutoff/8.0;    //this is the shift I have to apply to sin^2
    double shift_prime = fcutoff/2.0; 
    double shift_second = fcutoff/2.0; 
    double Amp1 = 0.5;  //0.5
    double Amp_prime = 0.0;
    double Amp_second = 0.0;
    double Amp2 = Amp1; //0.5
    if(arg*arg>10000.0){
        return (w>0 ? 0.0 : 1.0);
    }else{
      
           if(w >= shift && w <= (shift + T2 ) ){
               return (A * 1.0/(exp(arg)+1.0) +  Amp2 * sin(gamma2 * (w-shift)) *  sin(gamma2 * (w-shift)) +   Amp_prime * sin(gamma2 * (w-shift_prime)) *  sin(gamma2 * (w-shift_prime)) +  Amp_second * sin(gamma2 * (w-shift_second)) *  sin(gamma2 * (w-shift_second))     ) ;
               }
	   
           else if(w <= -shift && w >= (-shift - T1 ) ){
              return (A * 1.0/(exp(arg)+1.0) -  Amp1 * sin(gamma1 * (w+shift)) *  sin(gamma1 * (w+shift)) -   Amp_prime * sin(gamma1 * (w+shift_prime)) *  sin(gamma1 * (w+shift_prime)) -  Amp_second * sin(gamma1 * (w+shift_second)) *  sin(gamma1 * (w+shift_second))     ) ;
               }
	   
          else return A * 1.0/(exp(arg)+1.0);
    }
}

void force_equi(GF &G,double beta){
     for(long w = 0; w < G.ngrid_; w++){
        double omega=G.grid_[w];
	//cout<<"omega = "<<omega<<endl;
        double fw=fermi(omega,beta);
        G.Lesser[w]=fw*(G.Retarded[w].adjoint()-G.Retarded[w]);
        //cout << G.Lesser[w] << endl;
    }
}

void setBath_diag(double Jcoupl, double chemPot, double beta, GF &Bath){
	cplx imag_one(0,1);
	// asser bath is square
    cdmatrix Retarded(Bath.size1_,Bath.size1_);
    Retarded.setZero();
	for(int j=0;j<Bath.size1_;j++) Retarded(j,j)=-imag_one*Jcoupl;
	for(long w=0; w<Bath.ngrid_; w++) Bath.Retarded[w]=Retarded;
    force_equi(Bath,beta);
}

void force_distribution_bump(GF &G,double beta,double shift,double A){
     for(long w = 0; w < G.ngrid_; w++){
        double omega=G.grid_[w];
        double fw=fermi(omega,beta);
        fw+=exp(-(omega-shift)*(omega-shift))*A-exp(-(omega+shift)*(omega+shift))*A;
        G.Lesser[w]=fw*(G.Retarded[w].adjoint()-G.Retarded[w]);
        //cout << G.Lesser[w] << endl;
      }
}



void set_phonon_bath(double omega_ph,double beta, GF &Phonon, int traj){
cplx ii(0.0,1.0);
//I define the DOS of the phonon bath A(w),bose(w), with which I calculate Phonon Les and Great
vector <double> A(Phonon.ngrid_,0.0);
vector <double> bose(Phonon.ngrid_,0.0);
cdmatrix Lesser(Phonon.size1_,Phonon.size1_), Retarded(Phonon.size1_,Phonon.size1_);
Lesser.setZero();
Retarded.setZero();
double omega = 0.0, DOS_pos = 0.0, DOS_neg = 0.0, DOS_zero = 0.0;
for (long w=0;w<Phonon.ngrid_;w++){
//cout<<"w = "<<w<<endl;
omega=Phonon.grid_[w];
double coeff=0.25;//if I want that \int from 0 to \inf |DOS(w)| dw=1, I should impose coeff=1.0; 

	if (omega>0){	
			//cout<<"omega = "<<omega<<endl;
			A[w]= coeff*omega/(omega_ph*omega_ph)*exp(-omega/omega_ph);
			DOS_pos+=A[w];
			bose[w]=1.0/(exp(omega*beta)-1);
			for(int j=0;j<Phonon.size1_;j++) Lesser(j,j)=-ii*2.*M_PI*A[w]*bose[w];  //the minus sign is because we are dealing with bosons
			for(int j=0;j<Phonon.size1_;j++) Retarded(j,j)=-ii*M_PI*A[w];

			Phonon.Lesser[w]=Lesser;
			Phonon.Retarded[w]=Retarded;
					}

	if(omega<0){
			//cout<<"omega = "<<omega<<endl;
			A[w]= coeff*omega/(omega_ph*omega_ph)*exp(-(-omega)/omega_ph);
			DOS_neg+=A[w];
			bose[w]=1.0/(exp((omega)*beta)-1);
			for(int j=0;j<Phonon.size1_;j++) Lesser(j,j)=-ii*2.*M_PI*A[w]*bose[w];
			for(int j=0;j<Phonon.size1_;j++) Retarded(j,j)=-ii*M_PI*A[w];

			Phonon.Lesser[w]=Lesser;
			Phonon.Retarded[w]=Retarded;
					}

        if(omega==0){ //I make lim w->0 of A(w) and b(w)
                        A[w]= coeff*omega/(omega_ph*omega_ph);
			DOS_zero+=A[w];
                        //bose[w]=1.0/(omega*beta);
                        for(int j=0;j<Phonon.size1_;j++) Lesser(j,j)=-ii*2.*M_PI*coeff/(omega_ph*omega_ph*beta);  //this is -ii*2.*M_PI*A[w]*b[w] in the limit of A[w],b[w] for w->0
                        for(int j=0;j<Phonon.size1_;j++) Retarded(j,j)=-ii*M_PI*A[w];
                        Phonon.Lesser[w]=Lesser;
                        Phonon.Retarded[w]=Retarded;
                        //cout<<"Lesser= "<<Lesser.imag()<<endl;
                        }


}

//outputGF_new(Phonon,"temp_Phonon.out");

DOS_pos*=Phonon.dgrid_;
DOS_neg*=Phonon.dgrid_;
DOS_zero*=Phonon.dgrid_;	
if (traj==0){
cout<<endl;
cout<<"***** Phonon bath ****************"<<endl;
cout<<"DOS_neg = "<<DOS_neg<<endl;
cout<<"DOS_pos = "<<DOS_pos<<endl;
cout<<"DOS_zero = "<<DOS_zero<<endl;
cout<<"DOS_pos - DOS_neg = "<<DOS_pos - DOS_neg <<endl;
cout<<"DOS_pos_00 + DOS_neg_00 + DOS_zero_00 = "<<DOS_pos + DOS_neg + DOS_zero <<endl;
cout<<"2 x ( 1 - (DOS_pos + DOS_neg + DOS_zero) ) = "<< 2.0 * (1.0 - (DOS_pos + DOS_neg + DOS_zero) ) <<endl;
cout<<endl;
}



}


void set_fermion_bath_bump(GF &Bath,double V,double beta, double shift){

	cplx ii(0.0,1.0);
	//double shift=1.5, T=3.0, Amp=1.0/T;
	double  T=3.0, Amp=1.0/T;
	double DOS(0.0),f(1.0),omega0(shift), gamma(M_PI/T);

	cplx ret;
	cdmatrix Retarded(Bath.size1_,Bath.size1_);
	Retarded.setZero();

	cplx les;
	cdmatrix Lesser(Bath.size1_,Bath.size1_);
	Lesser.setZero();

	double DOS_neg(0.0),DOS_pos(0.0),DOS_zero(0.0),f_neg(0.0),f_pos(0.0),f_zero(0.0),occ_pos(0.0),occ_neg(0.0),occ_zero(0.0);
	//DOS_neg.setZero(),DOS_pos.setZero(),DOS_zero.setZero(),f_neg.setZero(),f_pos.setZero(),f_zero.setZero();


	for (long w=0;w<Bath.ngrid_;w++){
	double omega=Bath.grid_[w];
	f=fermi_negative_beta(omega,beta); //fermi-dirac at negative temperature

	if (omega<0) {omega0=-shift;}
	if( (omega>=omega0-T/2) and (omega<=omega0+T/2) ) DOS=Amp*cos(gamma*(omega-omega0))*cos(gamma*(omega-omega0));
		
	ret=-ii*M_PI*DOS;
	les=ii*2.*M_PI*DOS*f;


	for (int j=0;j<Bath.size1_;j++){
	Retarded(j,j)=ret;
	Lesser(j,j)=les;
	}

	Bath.Retarded[w]=Retarded;
	Bath.Lesser[w]=Lesser;
	//make the numerical check
	if (omega<0){ 
	DOS_neg+=DOS*Bath.dgrid_;
	f_neg+=f*Bath.dgrid_;
	occ_neg+=f*DOS*Bath.dgrid_;
	}
	if (omega>0){
	DOS_pos+=DOS*Bath.dgrid_;
	f_pos+=f*Bath.dgrid_;
	occ_pos+=f*DOS*Bath.dgrid_;
	}
	if (omega==0){
	DOS_zero+=DOS*Bath.dgrid_;
	f_zero+=f*Bath.dgrid_;
	occ_zero+=f*DOS*Bath.dgrid_;
	}

	}
	Bath.smul(V*V);


}


void setBath_Bethe_diag(double Jcoupl, double wcut, double beta, GF &Bath){
	// semielitic DOS in range [-2*wcut,2*wcut]
    // asser bath is square
    cplx ret;
    cdmatrix Retarded(Bath.size1_,Bath.size1_);
    Retarded.setZero();
	for(long w=0; w<Bath.ngrid_; w++){
        double omega=Bath.grid_[w];
        double om1=omega/wcut;
        double arg=om1*om1-4.0;
        if(arg>0){
            if(om1>0){
                ret=cplx((om1-sqrt(arg))/(2.0*wcut),0.0);
            }else{
                ret=cplx((om1+sqrt(arg))/(2.0*wcut),0.0);
            }
        }else{
                ret=cplx(om1/(2.0*wcut),-sqrt(-arg)/(2.0*wcut));
        }
        for(int j=0;j<Bath.size1_;j++) Retarded(j,j)=ret;
        Bath.Retarded[w]=Retarded;
    }
    force_equi(Bath,beta);
    Bath.smul(Jcoupl*wcut);
}

//time-dependent lattice functions
class lattice_func{
public:
    int ntheta_;
    double dtheta_;
    double grid_spacing_;
	long numOfFreq_;
    int size_; 
    vector<GF> Gk_;
    GF Sigma_,Gloc_,Bath_,G1_,G2_,Delta_,Floc_,D0_Phonon_,Sigma_G_D0,Sigma_fermion_bump_,Gweiss_,Polariz_;
    vector<GF> Fk_; // this is wasteful -- I use uly the Lesser part of Fk.
    vector<double> wk_;
    vector<cdmatrix> hk_;
    cdmatrix Hartree_,hloc_,X_,P_;
    ///////////////////////////////////////////////////////////////
    lattice_func(){
        ntheta_=0;
        dtheta_=0.0;
        grid_spacing_=0;
        numOfFreq_=0;
        size_=2;//1 ; we are now in symmetry broken phase
    }
    lattice_func(int ntheta,int numOfFreq,double grid_spacing){
        //assert(ntheta>1 and (ntheta%2)==0);
	assert(ntheta==0);  //in order to not waste memory, since the code is not k-depe
    	double norm=0.0;
        ntheta_=ntheta;
        dtheta_=M_PI/(2.0*ntheta);//Not anymore M_PI/(ntheta): symmetry broken phase
        grid_spacing_=grid_spacing;
        numOfFreq_=numOfFreq;
        size_=2; //1 we are now in symmetry broken phase
        Gk_.resize(ntheta_+1);
        Fk_.resize(ntheta_+1);
        wk_.resize(ntheta_+1);
        hk_.resize(ntheta_+1);
        for(int i=0;i<=ntheta_;i++){
            cdmatrix tmp(2,2); //symmetry broken phase
            double theta=i*dtheta_;
            double eps=-2.0*cos(theta);
            double simpson=(i==0 || i==ntheta_ ? 1.0 : (1+(i%2))*2.0);
            wk_[i]=sin(theta)*sin(theta)*dtheta_*(simpson/3.0);
            norm+=wk_[i];
            Gk_[i]=GF(grid_spacing_, numOfFreq_,size_);
            Fk_[i]=GF(grid_spacing_, numOfFreq_,size_);
            tmp.setZero();
	    tmp(0,1)=eps; 
	    tmp(1,0)=eps; 
            hk_[i]=tmp;
        }
        // renormalize weights to sum_k wk =1
        for(int i=0;i<=ntheta_;i++) wk_[i]=wk_[i]/norm;
        Sigma_=GF(grid_spacing_, numOfFreq_,size_);
        Gloc_=GF(grid_spacing_, numOfFreq_,size_);
        Floc_=GF(grid_spacing_, numOfFreq_,size_);
        G1_=GF(grid_spacing_, numOfFreq_,size_);
        G2_=GF(grid_spacing_, numOfFreq_,size_);
        Delta_=GF(grid_spacing_, numOfFreq_,size_);
        Bath_=GF(grid_spacing_, numOfFreq_,size_);
        Sigma_fermion_bump_=GF(grid_spacing_, numOfFreq_,size_);
	D0_Phonon_=GF(grid_spacing_, numOfFreq_,size_);
	Sigma_G_D0=GF(grid_spacing_, numOfFreq_,size_);
        Gweiss_=GF(grid_spacing_, numOfFreq_,size_);
        Polariz_=GF(grid_spacing_, numOfFreq_,size_);
        //Hartree_=cdmatrix(1,1);
        //Hartree_.setZero(); /// zero for now ...
        Hartree_=cdmatrix(size_,size_);
        Hartree_.setZero();
        hloc_=cdmatrix(size_,size_);
        hloc_.setZero();
        X_=cdmatrix(size_,size_);
        X_.setZero();
        P_=cdmatrix(size_,size_);
        P_.setZero();
	//we want to split hk in hloc +dhk
    }
    ~lattice_func(){
        // nothing to do
    }
    lattice_func(const lattice_func &ll){
        // copy content
        ntheta_=ll.ntheta_;
        dtheta_=ll.dtheta_;
        grid_spacing_=ll.grid_spacing_;
        numOfFreq_=ll.numOfFreq_;
        size_=ll.size_;
           Gk_.resize(ntheta_+1);
        Fk_.resize(ntheta_+1);
        wk_.resize(ntheta_+1);
        hk_.resize(ntheta_+1);
        for(int i=0;i<=ntheta_;i++){
            Gk_[i]=ll.Gk_[i];
            Fk_[i]=ll.Fk_[i];
            wk_[i]=ll.wk_[i];
            hk_[i]=ll.hk_[i];
        }
        Sigma_=ll.Sigma_;
        Gloc_=ll.Gloc_;
        Gweiss_=ll.Gweiss_;
        Floc_=ll.Floc_;
        Bath_=ll.Bath_;
        Sigma_fermion_bump_=ll.Sigma_fermion_bump_;
	D0_Phonon_=ll.D0_Phonon_;
	Sigma_G_D0=ll.Sigma_G_D0;
	Polariz_=ll.Polariz_;
        Delta_=ll.Delta_;
        G1_=ll.G1_;
        G2_=ll.G2_;
        Hartree_=ll.Hartree_;
	hloc_=ll.hloc_;
	hk_=ll.hk_;
	X_=ll.X_;
	P_=ll.P_;
    }
    lattice_func& operator=(const lattice_func &ll){
        ntheta_=ll.ntheta_;
        dtheta_=ll.dtheta_;
        grid_spacing_=ll.grid_spacing_;
        numOfFreq_=ll.numOfFreq_;
        size_=ll.size_;
           Gk_.resize(ntheta_+1);
        Fk_.resize(ntheta_+1);
        wk_.resize(ntheta_+1);
        hk_.resize(ntheta_+1);
        for(int i=0;i<=ntheta_;i++){
            Gk_[i]=ll.Gk_[i];
            Fk_[i]=ll.Fk_[i];
            wk_[i]=ll.wk_[i];
            hk_[i]=ll.hk_[i];
        }
        Sigma_=ll.Sigma_;
        Gloc_=ll.Gloc_;
        Gweiss_=ll.Gweiss_;
        Floc_=ll.Floc_;
        Bath_=ll.Bath_;
	Sigma_fermion_bump_=ll.Sigma_fermion_bump_;
	D0_Phonon_=ll.D0_Phonon_;
	Sigma_G_D0=ll.Sigma_G_D0;
	Polariz_=ll.Polariz_;
        Delta_=ll.Delta_;
        G1_=ll.G1_;
        G2_=ll.G2_;
        Hartree_=ll.Hartree_;
        hloc_=ll.hloc_;
        hk_=ll.hk_;
	X_=ll.X_;
	P_=ll.P_;
        return *this;
    }
    void get_gloc(void){
        Gloc_.clear();
        for(int i=0;i<=ntheta_;i++) Gloc_.incr(Gk_[i],wk_[i]);
    }

    void get_gloc_nonkdepe(void){
 	cdmatrix one(size_,size_);
	one.setZero();
        for (int j=0;j<size_;j++) one(j,j)=1.0;
        Gloc_.clear();
	for(long w = 0; w < Gloc_.ngrid_; w++){
        Gloc_.Retarded[w].noalias() =  (Gloc_.grid_[w]*one-Hartree_-hloc_-Sigma_.Retarded[w]-Delta_.Retarded[w]).inverse();
	Gloc_.Lesser[w].noalias() = Gloc_.Retarded[w]*(Sigma_.Lesser[w]+Delta_.Lesser[w])*Gloc_.Retarded[w].adjoint();
	}
	}

    void get_gloc_nonkdepe_ret(void){

  	cdmatrix one(size_,size_);
	one.setZero();
	cdmatrix DOS_pos(size_,size_), DOS_neg(size_,size_), DOS_zero(size_,size_);
	DOS_pos.setZero(),DOS_neg.setZero(),DOS_zero.setZero();

        for (int j=0;j<size_;j++) one(j,j)=1.0;
        Gloc_.clear();
	for(long w = 0; w < Gloc_.ngrid_; w++){
        Gloc_.Retarded[w].noalias() =  ( Gloc_.grid_[w]*one - Hartree_ - hloc_ - Sigma_.Retarded[w] - Delta_.Retarded[w] ).inverse();
	double omega=Gloc_.grid_[w];
	if (omega<0){ 
	DOS_neg+=-1./M_PI*Gloc_.Retarded[w].imag()*Gloc_.dgrid_;
	}
	if (omega>0){
	DOS_pos+=-1./M_PI*Gloc_.Retarded[w].imag()*Gloc_.dgrid_;

	}
	if (omega==0){
	DOS_zero+=-1./M_PI*Gloc_.Retarded[w].imag()*Gloc_.dgrid_;
	}

	}
    }


 void get_floc(void){
        Floc_.clear();
        for(int i=0;i<=ntheta_;i++) Floc_.incr(Fk_[i],wk_[i]);
    }


    void get_gloc_moments(void){
        G1_.clear();
        G2_.clear();
        GF gtmp;
        for(int i=0;i<=ntheta_;i++){
            gtmp=Gk_[i];
            gtmp.left_multiply(hk_[i],1.0);
            G1_.incr(gtmp,wk_[i]);
            gtmp.right_multiply(hk_[i],1.0);
            G2_.incr(gtmp,wk_[i]);
        }
        // make G1 and G2 diagonal
        GF Gzero(Gloc_.dgrid_,Gloc_.ngrid_,1);
        Gzero.clear();
        set_matrixelement(G1_,0,1,Gzero,0,0);
        set_matrixelement(G1_,1,0,Gzero,0,0);
        set_matrixelement(G2_,0,1,Gzero,0,0);
        set_matrixelement(G2_,1,0,Gzero,0,0);

    }
    void get_delta(void){
        // (1+G1) * Delta = G2, assuming that sum_k h_k=0
        //cdmatrix one(2,1);
        //cdmatrix tmp(2,2);
        //one.setZero();
        //one(0,0)=1.0;
        //one(1,1)=1.0;       
	cdmatrix one(size_,size_);
	one.setZero();
        for (int j=0;j<size_;j++) one(j,j)=1.0;
	cdmatrix tmp(size_,size_);

        for(long w = 0; w < Gloc_.ngrid_; w++){
            tmp=(one+G1_.Retarded[w]).inverse();
            // retarded
            Delta_.Retarded[w].noalias() = tmp*G2_.Retarded[w];
            // lesser:  DL + G1R DL + G1L DA = G2L
            Delta_.Lesser[w].noalias() = tmp*(G2_.Lesser[w] - G1_.Lesser[w]*(Delta_.Retarded[w].adjoint()));
        }
	Delta_.incr(Bath_);
      // now Gloc_=(w - Delta - Sigma)^{-1}, local terms added later.
	//Delta_.Lesser[w].noalias() = Gloc_.Retarded[w].inverse() * Gloc_.Lesser[w] * (Gloc_.Retarded[w].adjoint()).inverse()  - Sigma_.Lesser[w];

    	}
    

    void G_into_Delta(void){

    GF GA(Gloc_.dgrid_,Gloc_.ngrid_,1);
    GF GB(Gloc_.dgrid_,Gloc_.ngrid_,1);
    set_matrixelement(GA,0,0,Gloc_,0,0);
    set_matrixelement(GB,0,0,Gloc_,1,1);

    GF DeltaA(Gloc_.dgrid_,Gloc_.ngrid_,1);
    GF DeltaB(Gloc_.dgrid_,Gloc_.ngrid_,1);
 
    for(long w = 0; w < Gloc_.ngrid_; w++){
	DeltaA.Retarded[w]=GB.Retarded[w];
	DeltaA.Lesser[w]=GB.Lesser[w];
	DeltaB.Retarded[w]=GA.Retarded[w];
	DeltaB.Lesser[w]=GA.Lesser[w];
    	}

    Delta_.clear();
    set_matrixelement(Delta_,0,0,DeltaA,0,0);
    set_matrixelement(Delta_,1,1,DeltaB,0,0);

    }
//---------------case 1): the starting distribution function is the equilibrium one-------------------------------------------------------------------------------//
  
    void init_distribution_equi(int i,double beta){
        //cdmatrix one(1,1);
        //one(0,0)=1.0;
	cdmatrix one(size_,size_);
	one.setZero();
        for (int j=0;j<size_;j++) one(j,j)=1.0;
        for(long w = 0; w < Gk_[i].ngrid_; w++){
           double fw=fermi(Gk_[i].grid_[w],beta);
           Fk_[i].Lesser[w]=one*fw;
        }
    }
    void init_distribution_equi(double beta){
        for(int i=0;i<=ntheta_;i++) init_distribution_equi(i,beta);
        }

    void init_distribution_equi_nonkdepe(double beta){
        //cdmatrix one(1,1);
	//one.setZero();
        //one(0,0)=1.0;
	cdmatrix one(size_,size_);
	one.setZero();
        for (int j=0;j<size_;j++) one(j,j)=1.0;
	Floc_.clear();
        for(long w = 0; w < Gloc_.ngrid_; w++){
           double fw=fermi(Gloc_.grid_[w],beta);
           Floc_.Lesser[w]=one*fw;
        }
    }

      void set_Floc_from_Gloc(void){
	cdmatrix one(size_,size_);
	one.setZero();
        for (int j=0;j<size_;j++) one(j,j)=1.0;
	Floc_.clear();
	for(int j=0;j<Gloc_.size1_;j++) {
        for(long w = 0; w < Gloc_.ngrid_; w++){
           //double fw=fermi(Gloc_.grid_[w],beta);
           //Floc_.Lesser[w]=one*fw;
	   Floc_.Lesser[w](j,j)=-0.5*Gloc_.Lesser[w](j,j).imag()/Gloc_.Retarded[w](j,j).imag();
            }
	}


     }

//--------------- end case 1)-------------------------------------------------------------------------------//


//---------------case 2): the starting distribution function is a nonequilibrium one and is put by hand by myself-------------------------------------------------------------------------------//

    void init_distribution_nonequi(int i,double beta,double fcutoff){
        //cdmatrix one(1,1);
        //one(0,0)=1.0;
	cdmatrix one(size_,size_);
	one.setZero();
        for (int j=0;j<size_;j++) one(j,j)=1.0;
        for(long w = 0; w < Gk_[i].ngrid_; w++){
	  double fw=fermi_modified(Gk_[i].grid_[w],beta,fcutoff);
           Fk_[i].Lesser[w]=one*fw;
        }
    }
    void init_distribution_nonequi(double beta, double fcutoff){
          for(int i=0;i<=ntheta_;i++) init_distribution_nonequi(i,beta,fcutoff);
        }

//--------------- end case 2)-------------------------------------------------------------------------------//


  
  
//---------------case 3) : the input f are all the nk at a given time step comping from the realtime code-------------------------------------------------------------------------------//
      void init_distribution_nonequi_eps(double beta,double fcutoff,int time_to_boltzmann){
                   //cdmatrix one(1,1);
        	   //one(0,0)=1.0;
		   cdmatrix one(size_,size_);
		   one.setZero();
        	   for (int j=0;j<size_;j++) one(j,j)=1.0;

                   double fw;
                   ostringstream os;
		   string filename;
		   os.str("");
		   os.clear();
		   //os << "nk_" << time_to_boltzmann<<".txt";
		   os<<"Init_F.txt";
		   filename= os.str();
		   const char *filename_c = filename.c_str();
		   fstream myfile(filename_c, std::ios_base::in);
      		   // FILE *fout1=fopen(filename_c,"w");		
		   //fstream myfile("Init_F.txt", std::ios_base::in);       //this is in case I have as input file Init_F.txt
	
		   for(int i=0;i<=ntheta_;i++){ 
		     myfile>>fw;

			for(long w = 0; w < Gk_[i].ngrid_; w++){
			  // myfile>>fw;
			   Fk_[i].Lesser[w]=one*fw;
			}
			//fclose(fout1);
    		   }
}

	void load_init_distribution_local(int tstp){
	//cdmatrix one(1,1);
	//one(0,0)=1.0;	
	cdmatrix one(size_,size_);
	one.setZero();
	for (int j=0;j<size_;j++) one(j,j)=1.0;
	double fw;
	ostringstream os;
	string filename;
	os.str("");
	os.clear();
	os<<"Init_F"<<tstp<<".txt";
	filename= os.str();
	const char *filename_c = filename.c_str();
	fstream myfile(filename_c, std::ios_base::in);
	for(long w = 0; w < Floc_.ngrid_; w++){
	myfile>>fw;
	Floc_.Lesser[w]=one*fw;
	}
	}


//-------------------------------------end of case 3) -------------------------------------------------------------------------------------//
 
    void set_distribution_function(void){
            for(long w = 0; w < Floc_.ngrid_; w++){
                Gloc_.Lesser[w].noalias()=Floc_.Lesser[w]*(Gloc_.Retarded[w].adjoint())-Gloc_.Retarded[w]*Floc_.Lesser[w];
            }
		cdmatrix R(2,2);
		R.setZero();
		DensityMatrix(Gloc_,R);
    }

    void set_distribution_function_kdepe(void){
        for(int i=0;i<=ntheta_;i++){
            for(long w = 0; w < Gk_[i].ngrid_; w++){
                Gk_[i].Lesser[w]=Fk_[i].Lesser[w]*(Gk_[i].Retarded[w].adjoint())-Gk_[i].Retarded[w]*Fk_[i].Lesser[w];
            }
        }
    }

    void set_Hartree(double U){
 	Eigen::Matrix2cd R(2,2);
        Hartree_.setZero();
        DensityMatrix(Gloc_,R);
        Hartree_(0,0)=(R(0,0).real()-0.5)*U;//also mu is included
        Hartree_(1,1)=(R(1,1).real()-0.5)*U;
    }


   
   void set_hloc_diag(double ea,double eb,double g){
	hloc_(0,0)=ea+sqrt(2)*g*X_(0,0);
	hloc_(1,1)=eb+sqrt(2)*g*X_(1,1);
    }

    void set_hk_diag(double ea,double eb){
        for(int i=0;i<=ntheta_;i++){
	    //do nothing
            //hk_[i](0,0)=ea;
            //hk_[i](1,1)=eb;
        }
    }


    void set_sigma(GF &Sa){
        Sigma_.clear();
        set_matrixelement(Sigma_,0,0,Sa,0,0);
    }
    void output_all(std::string filename){
        ostringstream os;
        string Gfile;
        for(int i=0;i<=ntheta_;i++){
            os.str("");
            os.clear();
            os << filename << "_k" << i << ".out";
            Gfile= os.str();
            outputGF_new(Gk_[i], Gfile);
        }
        os.str("");
        os.clear();
        os << filename << "temp_Gloc.out";
        Gfile= os.str();
        outputGF_new(Gloc_, Gfile);
        os.str("");
        os.clear();
        os << filename << "temp_Sigma.out";     
        Gfile= os.str();
        outputGF_new(Sigma_,Gfile);      //remember that with Sigma_ we mean Sigma_interaction, whilst with Bath we mean Sigma_Bath
    }
}; //end of the lattice_func class


void avg_G_into_Delta(vector <lattice_func> &latt,GF & Delta_avg,vector <GF>& GA,vector <GF>& GB, int num_of_traj){	

    GF DeltaA_atomic(Delta_avg.dgrid_,Delta_avg.ngrid_,1);
    GF DeltaB_atomic(Delta_avg.dgrid_,Delta_avg.ngrid_,1);

    DeltaA_atomic.clear();
    DeltaB_atomic.clear();
   
   
    //prendo gli elemento 00 e 11 da tutte le Gloc e le metto in un vettore di 1D GA e GB. La dim del vettore e` pari al numero di traiettorie
    #pragma omp parallel for
    for (int i=0;i<num_of_traj;i++)
    {
    	    set_matrixelement(GA[i],0,0,latt[i].Gloc_,0,0);
    	    set_matrixelement(GB[i],0,0,latt[i].Gloc_,1,1);
    }


    for (int i=0;i<num_of_traj;i++)
    {
	//faccio la media delle GA e GB e le copio in DeltaA e DeltaB
        for(long w = 0; w < Delta_avg.ngrid_; w++){
	DeltaA_atomic.Retarded[w]+=GB[i].Retarded[w]/num_of_traj;
	DeltaA_atomic.Lesser[w]+=GB[i].Lesser[w]/num_of_traj;
	DeltaB_atomic.Retarded[w]+=GA[i].Retarded[w]/num_of_traj;
	DeltaB_atomic.Lesser[w]+=GA[i].Lesser[w]/num_of_traj;
    	}
	
     }

    

   //ogni Delta 2x2 (per ogni traiettoria) e` data dal valor medio delle GA e GB 
    Delta_avg.clear();
    set_matrixelement(Delta_avg,0,0,DeltaA_atomic,0,0);
    set_matrixelement(Delta_avg,1,1,DeltaB_atomic,0,0);

    //At the end I copy Delta_avg into all the latt[i].Delta_
    #pragma omp parallel for
    for (int i=0;i<num_of_traj;i++) {
    latt[i].Delta_ = Delta_avg;
    }	
}

void updateSE1(double U, GF &Gloc, GF &Sloc, fft_solver &solver)
{
    GF GA(Gloc.dgrid_,Gloc.ngrid_,1);
    GF GB(Gloc.dgrid_,Gloc.ngrid_,1);
    GF locSEB(Gloc.dgrid_,Gloc.ngrid_,1);
    GF locSEA(Gloc.dgrid_,Gloc.ngrid_,1);
    set_matrixelement(GA,0,0,Gloc,0,0);
    set_matrixelement(GB,0,0,Gloc,1,1);
    GF tGA(GA.dgrid_, GA.ngrid_ * 3 / 2 + 1, GA.size1_, time_gf);
    GF tlocSEA(GA.dgrid_, GA.ngrid_ * 3 / 2 + 1, GA.size1_, time_gf);
    GF tGB(GB.dgrid_, GB.ngrid_ * 3 / 2 + 1, GB.size1_, time_gf);
    GF tlocSEB(GB.dgrid_, GB.ngrid_ * 3 / 2 + 1, GB.size1_, time_gf);
    solver.to_time(tGA, GA, 1);
    solver.to_time(tGB, GB, 1);
    for(long t = 0; t < tGA.ngrid_; t++)
        {
                tlocSEA.Lesser[t](0,0) = tGA.Lesser[t](0,0) * (-tGA.Greater(t).adjoint()(0,0)) * tGA.Lesser[t](0,0);
                tlocSEA.Retarded[t](0,0) = tGA.Greater(t)(0,0) * (-tGA.Lesser[t].adjoint()(0,0)) * tGA.Greater(t)(0,0)-tlocSEA.Lesser[t](0,0);
                tlocSEB.Lesser[t](0,0) = tGB.Lesser[t](0,0) * (-tGB.Greater(t).adjoint()(0,0)) * tGB.Lesser[t](0,0);
                tlocSEB.Retarded[t](0,0) = tGB.Greater(t)(0,0) * (-tGB.Lesser[t].adjoint()(0,0)) * tGB.Greater(t)(0,0)-tlocSEB.Lesser[t](0,0);
        }
    tlocSEA.smul(U * U);
        tlocSEA.reset_grid(tGA.dgrid_);
    tlocSEB.smul(U * U);
        tlocSEB.reset_grid(tGB.dgrid_);
    tlocSEA.Retarded[0] *= 0.5;
        tlocSEA.Retarded[tlocSEA.ngrid_ - 1] *= 0.0;
    tlocSEB.Retarded[0] *= 0.5;
        tlocSEB.Retarded[tlocSEB.ngrid_ - 1] *= 0.0;
    solver.to_freq(locSEA, tlocSEA);
        solver.to_freq(locSEB, tlocSEB);
    ForceImag(locSEB);
    	ForceImag(locSEA);
    Sloc.clear();
    set_matrixelement(Sloc,0,0,locSEA,0,0);
    set_matrixelement(Sloc,1,1,locSEB,0,0);
   
}

void add_Sigma_phonon_to_Sigma2_ipt(GF &Sloc, GF& Sigma_Phonon){
    Sloc.incr(Sigma_Phonon);
}

void add_Sigma_fermion_bump_to_Sigma2_ipt(GF &Sloc, GF& Sigma_fermion){
    Sloc.incr(Sigma_fermion);
}


void set_Sigma_Phonon_G_D(double g_el_ph, GF & Gloc, GF & D0,GF & Sigma_ph,fft_solver &solver){
     cplx ii(0.0,1.0);

     GF GA(Gloc.dgrid_,Gloc.ngrid_,1);
     GF GB(Gloc.dgrid_,Gloc.ngrid_,1);
     GF locSEB(Gloc.dgrid_,Gloc.ngrid_,1);
     GF locSEA(Gloc.dgrid_,Gloc.ngrid_,1);
     GF D0B(Gloc.dgrid_,Gloc.ngrid_,1);
     GF D0A(Gloc.dgrid_,Gloc.ngrid_,1);
     set_matrixelement(GA,0,0,Gloc,0,0);
     set_matrixelement(GB,0,0,Gloc,1,1);
     set_matrixelement(D0A,0,0,D0,0,0);
     set_matrixelement(D0B,0,0,D0,1,1);

     GF tGA(GA.dgrid_, GA.ngrid_ * 3 / 2 + 1, GA.size1_, time_gf);
     GF tlocSEA(GA.dgrid_, GA.ngrid_ * 3 / 2 + 1, GA.size1_, time_gf);
     GF tD0A(GA.dgrid_, GA.ngrid_ * 3 / 2 + 1, GA.size1_, time_gf);
     GF tGB(GB.dgrid_, GB.ngrid_ * 3 / 2 + 1, GB.size1_, time_gf);
     GF tlocSEB(GB.dgrid_, GB.ngrid_ * 3 / 2 + 1, GB.size1_, time_gf);
     GF tD0B(GB.dgrid_, GB.ngrid_ * 3 / 2 + 1, GB.size1_, time_gf);

     solver.to_time(tGA, GA, 1);
     solver.to_time(tGB, GB, 1);
     solver.to_time(tD0A, D0A, 1);
     solver.to_time(tD0B, D0B, 1);


     for(long t = 0; t < tGA.ngrid_; t++)
     {
        tlocSEA.Lesser[t](0,0) =  tGA.Lesser[t](0,0) * tD0A.Lesser[t](0,0);
        tlocSEA.Retarded[t](0,0) = tGA.Greater(t)(0,0) * tD0A.Greater(t)(0,0) - tlocSEA.Lesser[t](0,0); 
        tlocSEB.Lesser[t](0,0) = tGB.Lesser[t](0,0) * tD0B.Lesser[t](0,0);
        tlocSEB.Retarded[t](0,0) =tGB.Greater(t)(0,0) * tD0B.Greater(t)(0,0) - tlocSEB.Lesser[t](0,0); 
     }


     tlocSEA.smul(g_el_ph * g_el_ph);
     tlocSEA.smul(ii);
     tlocSEA.reset_grid(tGA.dgrid_);
     tlocSEB.smul(g_el_ph * g_el_ph);
     tlocSEB.smul(ii);
     tlocSEB.reset_grid(tGB.dgrid_);
     tlocSEA.Retarded[0] *= 0.5;
     tlocSEA.Retarded[tlocSEA.ngrid_ - 1] *= 0.0;
     tlocSEB.Retarded[0] *= 0.5;
     tlocSEB.Retarded[tlocSEB.ngrid_ - 1] *= 0.0;


     solver.to_freq(locSEA, tlocSEA);
     solver.to_freq(locSEB, tlocSEB);


     //ForceImag(locSEA);
     //ForceImag(locSEB);  //I do not know if I need to impose imag Sigma also for Sigma_phonon
     Sigma_ph.clear();


     set_matrixelement(Sigma_ph,0,0,locSEA,0,0);
     set_matrixelement(Sigma_ph,1,1,locSEB,0,0);

}


void get_gweiss_from_Delta(GF & Gloc, GF & Gweiss, GF & Sigma,GF & Delta, cdmatrix & Hartree,cdmatrix & hloc){
    // solve Gweiss = (w -Delta)^{-1}.
    cdmatrix one(2,2);
    one.setZero();
    one(0,0)=1.0;
    one(1,1)=1.0;
    for(long w = 0; w < Gweiss.ngrid_; w++){
        Gweiss.Retarded[w].noalias() = (Gweiss.grid_[w]*one - Hartree- hloc- Delta.Retarded[w]).inverse();
        //Gweiss.Retarded[w].noalias() = (Gloc.Retarded[w].inverse() + Sigma.Retarded[w]).inverse();
        Gweiss.Lesser[w].noalias() = Gweiss.Retarded[w] * Delta.Lesser[w]* Gweiss.Retarded[w].adjoint();
    }
}
class sigma_solver{
public:
    fft_solver solver_;
    sigma_solver(long numOfFreq):
    solver_(numOfFreq, FFTW_ESTIMATE)
    {
        // nothing more to do
    }
    void sigma2_bold(lattice_func &latt,double U){
        latt.get_gloc();
        updateSE1(U,latt.Gloc_,latt.Sigma_,solver_);
	add_Sigma_phonon_to_Sigma2_ipt(latt.Sigma_,latt.Sigma_G_D0);
	add_Sigma_fermion_bump_to_Sigma2_ipt(latt.Sigma_,latt.Sigma_fermion_bump_);
    }
    void sigma2_ipt(vector <lattice_func> &latt, double U,double g_ph_bath, double V,int it,int iter, int num_of_traj, int traj_num){
        GF Gweiss=latt[traj_num].Gloc_;

	G_into_Delta(latt[traj_num].Gloc_,latt[traj_num].Delta_);  //at equilibrium 
	latt[traj_num].Sigma_.clear();

        if (U>0){ 
           get_gweiss_from_Delta(latt[traj_num].Gloc_,Gweiss,latt[traj_num].Sigma_,latt[traj_num].Delta_,latt[traj_num].Hartree_,latt[traj_num].hloc_);
	   updateSE1(U,Gweiss,latt[traj_num].Sigma_,solver_);
	}
	if (g_ph_bath !=0.){
	   set_Sigma_Phonon_G_D(g_ph_bath,latt[traj_num].Gloc_,latt[traj_num].D0_Phonon_,latt[traj_num].Sigma_G_D0,solver_);
	   add_Sigma_phonon_to_Sigma2_ipt(latt[traj_num].Sigma_,latt[traj_num].Sigma_G_D0);
	}

	if (V>0) add_Sigma_fermion_bump_to_Sigma2_ipt(latt[traj_num].Sigma_,latt[traj_num].Sigma_fermion_bump_);
	if (it==0 and iter<=1) latt[traj_num].Sigma_.incr(latt[traj_num].Bath_);

    }
};

void get_polarizability(GF& Gloc,GF& Polariz ,fft_solver & solver){
cplx ii(0.0,1.0);
GF GA(Gloc.dgrid_,Gloc.ngrid_,1);
GF GB(Gloc.dgrid_,Gloc.ngrid_,1);
GF PolB(Polariz.dgrid_,Polariz.ngrid_,1);
GF PolA(Polariz.dgrid_,Polariz.ngrid_,1);
set_matrixelement(GA,0,0,Gloc,0,0);
set_matrixelement(GB,0,0,Gloc,1,1);

GF tGA(GA.dgrid_, GA.ngrid_ * 3 / 2 + 1, GA.size1_, time_gf);
GF tPolA(PolA.dgrid_, PolA.ngrid_ * 3 / 2 + 1, PolA.size1_, time_gf);
GF tGB(GB.dgrid_, GB.ngrid_ * 3 / 2 + 1, GB.size1_, time_gf);
GF tPolB(PolB.dgrid_, PolB.ngrid_ * 3 / 2 + 1, PolB.size1_, time_gf);
solver.to_time(tGA, GA, 1);
solver.to_time(tGB, GB, 1);

    for(long t = 0; t < tGA.ngrid_; t++)
        {
		tPolA.Retarded[t](0,0) = tGA.Retarded[t](0,0) * (-tGA.Greater(t).adjoint()(0,0)) + tGA.Greater(t)(0,0) * tGA.Retarded[t].adjoint()(0,0); //Pi_R(t)=G_R(t)*G_>(-t)+G_>(t)*G_A(-t)
                tPolB.Retarded[t](0,0) = tGB.Retarded[t](0,0) * (-tGB.Greater(t).adjoint()(0,0)) + tGB.Greater(t)(0,0) * tGB.Retarded[t].adjoint()(0,0); //Pi_R(t)=G_R(t)*G_>(-t)+G_>(t)*G_A(-t)

                //N.B.: in the following I write in the Lesser component of the polarizability what de facto is the Keldysh component
                tPolA.Lesser[t](0,0) =  tGA.Lesser[t](0,0) * (-tGA.Greater(t).adjoint()(0,0))+tGA.Greater(t)(0,0)*(-tGA.Lesser[t].adjoint()(0,0));    ; //Pi_K(t)=G_<(t)*G_>(-t)+G_>(t)*G_<(-t) 
                tPolB.Lesser[t](0,0) =  tGB.Lesser[t](0,0) * (-tGB.Greater(t).adjoint()(0,0)) +  tGB.Greater(t)(0,0)*(-tGB.Lesser[t].adjoint()(0,0));    ; //Pi_K(t)=G_<(t)*G_>(-t)+G_>(t)*G_<(-t) 


        }

    tPolA.smul(ii); //ii and not ii*0.5 since there is the sum over \sigma
    tPolA.reset_grid(tGA.dgrid_);
        tPolB.smul(ii);
        tPolB.reset_grid(tGB.dgrid_);
    tPolA.Retarded[0] *= 0.5;
    tPolA.Retarded[tPolA.ngrid_ - 1] *= 0.0;
        tPolB.Retarded[0] *= 0.5;
        tPolB.Retarded[tPolB.ngrid_ - 1] *= 0.0;
    solver.to_freq(PolA, tPolA);
    solver.to_freq(PolB, tPolB);
    //ForceImag(locSEA);
    //ForceImag(locSEB);  //I do not know if I need to impose imag Sigma also for Sigma_phonon
    Polariz.clear();
    set_matrixelement(Polariz,0,0,PolA,0,0);
    set_matrixelement(Polariz,1,1,PolB,0,0);
}

void set_Delta_retarded(GF & Delta,GF &Floc, sigma_solver &imp){

		cdmatrix ii(1,1); ii.setZero();ii(0,0)=cplx(0.0,1.0);

		cdmatrix DOS(1,1),f(1,1),DOS_pos(1,1),DOS_neg(1,1),DOS_zero(1,1),occ_pos(1,1),occ_neg(1,1),occ_zero(1,1);
		DOS_pos.setZero(),DOS_neg.setZero(),DOS_zero.setZero(),occ_pos.setZero(),occ_neg.setZero(),occ_zero.setZero();



	    double omega=0.0;
	    Delta.clear();
	    GF DA(Delta.dgrid_,Delta.ngrid_,1);
	    GF DB(Delta.dgrid_,Delta.ngrid_,1);
	    for(long w = 0; w < Delta.ngrid_; w++){
			omega=Delta.grid_[w];
			if (omega>=-2.0 and omega<=2.0){
			   DA.Retarded[w]=M_PI*(-1.)   /   (2.*M_PI)*sqrt(4.0-omega*omega)  *    ii;
			   DB.Retarded[w]=M_PI*(-1.)/(2.*M_PI)*sqrt(4.0-omega*omega)*ii;
			  }

			DOS=-1./M_PI*DA.Retarded[w].imag();
			
			if (omega<0){ 
			DOS_neg+=DOS.real()*Delta.dgrid_;
			}

			
			if (omega>0){
			DOS_pos+=DOS.real()*Delta.dgrid_;
			}
			if (omega==0){
			DOS_zero+=DOS.real()*Delta.dgrid_;
			}

	    }
		cout<<endl;
		cout<<"DA before"<<endl;
		cout<<"DOS_neg = "<<DOS_neg<<endl;
		cout<<"DOS_pos = "<<DOS_pos<<endl;
		cout<<"DOS_zero = "<<DOS_zero<<endl;
		cout<<"DOS_pos - DOS_neg = "<<DOS_pos -DOS_neg <<endl;
		cout<<"DOS_pos + DOS_neg + DOS_zero = "<<DOS_pos + DOS_neg + DOS_zero <<endl;
		cout<<endl;

	    Delta.clear();
	    set_matrixelement(Delta,0,0,DA,0,0);
	    set_matrixelement(Delta,1,1,DB,0,0);

}




double dmft_iteration(int iter,vector <lattice_func> &latt,sigma_solver & imp,int ipt,double U,int ret,double mixing,int tstp,string flag,double g,double V,int num_of_traj,int traj_num,double eA,double eB, double g_ph_bath){
    double err;
    GF G_old=latt[traj_num].Gloc_;
    latt[traj_num].set_Hartree(U);
    latt[traj_num].set_hloc_diag(eA,eB,g);
    imp.sigma2_ipt(latt,U,g_ph_bath,V,tstp,iter,num_of_traj,traj_num);
    if(ret==1) {	
	latt[traj_num].get_gloc_nonkdepe_ret();
        latt[traj_num].set_distribution_function();
    }else{
       	latt[traj_num].get_gloc_nonkdepe();
    }
    mixing_G(latt[traj_num].Gloc_,G_old,mixing);
    err = GF2norm(latt[traj_num].Gloc_, G_old);
    
    return err;
}

//to be adapted to the 2x2 case
void get_ISigma_F(lattice_func &latt_new, lattice_func &latt){
    // assuming that retarded dyson has been solved
    // I = Sigma< + Sigma_R*F-F*Sigma_A
    // N.B. inside Sigma I put: the IPT Sigma and the Sigma from the three baths (that eas already added in the DMFT iteration in the function updateSE1 ), and the Delta since the code is not anymore kdependent

	latt_new.Floc_.clear();
        for(int w=0;w<latt.Floc_.ngrid_;w++){
           cdmatrix Rtmp=latt.Sigma_.Retarded[w]+latt.Delta_.Retarded[w];
            cdmatrix Ltmp=latt.Sigma_.Lesser[w]+latt.Delta_.Lesser[w];
	     latt_new.Floc_.Lesser[w]=Ltmp+Rtmp*latt.Floc_.Lesser[w]-latt.Floc_.Lesser[w]*(Rtmp.adjoint());	    
        }
	
}

void updateX_equilibrium_no_noise(double omega0,double g_Xn,double Omega,lattice_func &latt, int iter,double mixing, double& X_err, double U,double Amp_f,double T_f,double shift_f){
cdmatrix one(2,2);
one.setZero();
one(0,0)=1.0;
one(1,1)=1.0;
cdmatrix R(2,2);
R.setZero();
DensityMatrix(latt.Gloc_,R); 
cdmatrix el_field(2,2);
el_field.setZero();
el_field(0,0) = electric_field(0.0,Omega,Amp_f,T_f,shift_f);
el_field(1,1) = - electric_field(0.0,Omega,Amp_f,T_f,shift_f);
cdmatrix Xold=latt.X_;
latt.X_(0,0) = sqrt(2.)*g_Xn*(1.0-2.*R(0,0).real())/omega0;
//latt.X_(0,0) = sqrt(2.)*g_Xn*(1.0-2.*R(0,0).real())/(omega0-2.*g_Xn*g_Xn*latt.Polariz_.Retarded[0](0,0).real());
latt.X_(1,1) = -latt.X_(0,0);
latt.X_=mixing*Xold+(1.0-mixing)*latt.X_;
X_err=abs(latt.X_(0,0)-Xold(0,0))+abs(latt.X_(1,1)-Xold(1,1));
}


void dFdt(lattice_func &latt_new_atomic,vector <lattice_func> &latt,sigma_solver &imp,double U,double errmax,int Niter,int ipt,std::string flag,double mixing, int ret,int tstp, double g_Xn,double V,int num_traj,int traj_num,double eA,double eB, double g_ph_bath){
    double err=10.0;
    int i=0;

   //cout<<"U = "<<U<<endl;
    while(i<Niter and err>errmax){
      err=dmft_iteration(i,latt,imp,ipt,U,ret,mixing,tstp,flag,g_Xn,V,num_traj,traj_num,eA,eB,g_ph_bath);
      if (traj_num==0) cout << "time " << flag << " -- iter " << i << " err " << err << endl;
        i++;
    }
    //after convergence, I calculate the scattering integral 
    get_ISigma_F(latt_new_atomic,latt[traj_num]);
    
}

void dFdt_new(lattice_func & latt, GF & Delta_avg, sigma_solver &imp, double U,int ipt,std::string flag,double mixing, int ret, double g_Xn,double V,int num_traj,int traj_num,double eA,double eB, vector <double>& err, vector <int>& iter, double g_ph_bath){

	GF G_old=latt.Gloc_;
 
	latt.set_Hartree(U);
	latt.set_hloc_diag(eA,eB,g_Xn);
  
		GF Gweiss=latt.Gloc_;
		//avg_G_into_Delta(ltmp,Delta_avg,num_of_traj);
		latt.Sigma_.clear();

		if (U>0){ 
		   get_gweiss_from_Delta(latt.Gloc_,Gweiss,latt.Sigma_,Delta_avg,latt.Hartree_,latt.hloc_);
		   updateSE1(U,Gweiss,latt.Sigma_,imp.solver_);
		}
		if (g_Xn >0){
		   set_Sigma_Phonon_G_D(g_ph_bath,latt.Gloc_,latt.D0_Phonon_,latt.Sigma_G_D0,imp.solver_);
		   add_Sigma_phonon_to_Sigma2_ipt(latt.Sigma_,latt.Sigma_G_D0);
		}

		if (V>0) add_Sigma_fermion_bump_to_Sigma2_ipt(latt.Sigma_,latt.Sigma_fermion_bump_);


	//latt.set_Hartree(U);
	//latt.set_hloc_diag(eA,eB,g_Xn);
	if (ret==1){
	latt.get_gloc_nonkdepe_ret();
	latt.set_distribution_function();
	}
	else{
	latt.get_gloc_nonkdepe();
	}
	mixing_G(latt.Gloc_,G_old,mixing);
	err[traj_num] = GF2norm(latt.Gloc_, G_old);
        if (traj_num==0) cout <<"traj "<<traj_num<<" time " << flag << " -- iter " << iter[traj_num] << " err " << err[traj_num] << endl;
	iter[traj_num]++;
   
}



void set_les(lattice_func &lout,lattice_func &lin){
        for(int w=0;w<lin.Bath_.ngrid_;w++){
            lout.Gloc_.Lesser[w]=lin.Gloc_.Lesser[w];
        }
}

void incr_les(lattice_func &lout,lattice_func &lin,double a){
        for(int w=0;w<lin.Bath_.ngrid_;w++){
            lout.Gloc_.Lesser[w]+=a*lin.Gloc_.Lesser[w];
        }
}

void smul_les(lattice_func &lout,cplx a){
        for(int w=0;w<lout.Bath_.ngrid_;w++){
            lout.Gloc_.Lesser[w]*=a;
        }
}

void set_les_F(lattice_func &lout,lattice_func &lin){
        for(int w=0;w<lin.Bath_.ngrid_;w++){
            lout.Floc_.Lesser[w]=lin.Floc_.Lesser[w];
        }
}

void incr_les_F(lattice_func &lout,lattice_func &lin,double a){
        for(int w=0;w<lin.Bath_.ngrid_;w++){
            lout.Floc_.Lesser[w]+=a*lin.Floc_.Lesser[w];
        }
}

void smul_les_F(lattice_func &lout,cplx a){
        for(int w=0;w<lout.Floc_.ngrid_;w++){
            lout.Floc_.Lesser[w]*=a;
        }
}

void set_les_k(lattice_func &lout,lattice_func &lin){
    for(int i=0;i<=lin.ntheta_;i++){
        for(int w=0;w<lin.Bath_.ngrid_;w++){
            lout.Gk_[i].Lesser[w]=lin.Gk_[i].Lesser[w];
        }
    }
}
void incr_les_k(lattice_func &lout,lattice_func &lin,double a){
    for(int i=0;i<=lin.ntheta_;i++){
        for(int w=0;w<lin.Bath_.ngrid_;w++){
            lout.Gk_[i].Lesser[w]+=a*lin.Gk_[i].Lesser[w];
        }
    }
}
void smul_les_k(lattice_func &lout,cplx a){
    for(int i=0;i<=lout.ntheta_;i++){
        for(int w=0;w<lout.Bath_.ngrid_;w++){
            lout.Gk_[i].Lesser[w]*=a;
        }
    }
}
void set_les_F_k(lattice_func &lout,lattice_func &lin){
    for(int i=0;i<=lin.ntheta_;i++){
        for(int w=0;w<lin.Bath_.ngrid_;w++){
            lout.Fk_[i].Lesser[w]=lin.Fk_[i].Lesser[w];
        }
    }
}
void incr_les_F_k(lattice_func &lout,lattice_func &lin,double a){
    for(int i=0;i<=lin.ntheta_;i++){
        for(int w=0;w<lin.Bath_.ngrid_;w++){
            lout.Fk_[i].Lesser[w]+=a*lin.Fk_[i].Lesser[w];
        }
    }
}
void smul_les_F_k(lattice_func &lout,cplx a){
    for(int i=0;i<=lout.ntheta_;i++){
        for(int w=0;w<lout.Bath_.ngrid_;w++){
            lout.Fk_[i].Lesser[w]*=a;
        }
    }
}


void set_X(lattice_func &lout, lattice_func &lin){
     lout.X_=lin.X_;
}

void incr_X(lattice_func &lout,lattice_func &lin,double a){        
            lout.X_+=a*lin.X_;            
}

void smul_X(lattice_func &lout,cplx a){
            lout.X_*=a;
}

void set_P(lattice_func &lout, lattice_func &lin){
	lout.P_=lin.P_;
}
void incr_P(lattice_func &lout,lattice_func &lin,double a){        
            lout.P_+=a*lin.P_;            
}

void incr_P_noise(lattice_func &lout,double incr, double dW_A, double dW_B){
        //lout.P_+=a*lin.P_;            
        lout.P_(0,0) += incr*dW_A;
        lout.P_(1,1) += incr*dW_B;
	//cout<<"dW_A = "<<dW_A<<" dW_B = "<<dW_B<<endl;
}


void smul_P(lattice_func &lout,cplx a){
            lout.P_*=a;
}

void fx(lattice_func &lout, lattice_func& lin,double omega0);

void set_X_fx(lattice_func &lout, lattice_func &lin, double omega0){
	fx(lout,lin,omega0);
}


void fx(lattice_func &lout, lattice_func& lin,double omega0){
	 lout.X_=omega0*lin.P_;
	//lout.X_=lin.P_;
}

void fp(lattice_func &lout,lattice_func &lin,double el_field,double omega0,double g_ph,double gamma_phonon_bath);

void set_P_fp(lattice_func &lout, lattice_func &lin, double el_field,double omega0,double g_ph,double gamma_phonon_bath,sigma_solver & imp){
	get_polarizability(lin.Gloc_,lin.Polariz_,imp.solver_);
	fp(lout,lin,el_field,omega0,g_ph,gamma_phonon_bath);
}


void fp(lattice_func &lout,lattice_func &lin, double el_field,double omega0,double g_ph, double gamma_phonon_bath){
	cplx ii(0.0,1.0);
	cdmatrix one(2,2);
	one.setZero();
	one(0,0)=1.0;
	one(1,1)=1.0;
	cdmatrix R(2,2);
	R.setZero();
	cdmatrix electric_field(2,2);
	electric_field.setZero();
	electric_field(0,0)= el_field;
	electric_field(1,1)= - el_field;
	DensityMatrix(lin.Gloc_,R);
	cplx deriv_Polariz_Retard_in_zero_00(0.0,0.0),  deriv_Polariz_Retard_in_zero_11(0.0,0.0);

	deriv_Polariz_Retard_in_zero_00=1./(12.*lin.Polariz_.dgrid_)*(-lin.Polariz_.Retarded[2](0,0)+8.*lin.Polariz_.Retarded[1](0,0)-8.*lin.Polariz_.Retarded[lin.Polariz_.ngrid_-1](0,0)+lin.Polariz_.Retarded[lin.Polariz_.ngrid_-2](0,0));
	deriv_Polariz_Retard_in_zero_11=1./(12.*lin.Polariz_.dgrid_)*(-lin.Polariz_.Retarded[2](1,1)+8.*lin.Polariz_.Retarded[1](1,1)-8.*lin.Polariz_.Retarded[lin.Polariz_.ngrid_-1](1,1)+lin.Polariz_.Retarded[lin.Polariz_.ngrid_-2](1,1));

	lout.P_(0,0) = -omega0*lin.X_(0,0)-(2.*g_ph*g_ph*omega0*deriv_Polariz_Retard_in_zero_00.imag()+gamma_phonon_bath)*lin.P_(0,0) - sqrt(2)*g_ph*(2.*R(0,0).real()-one(0,0)) -sqrt(2)*el_field;
	lout.P_(1,1) = -omega0*lin.X_(1,1)-(2.*g_ph*g_ph*omega0*deriv_Polariz_Retard_in_zero_11.imag()+gamma_phonon_bath)*lin.P_(1,1) - sqrt(2)*g_ph*(2.*R(1,1).real()-one(0,0)) +sqrt(2)*el_field;

}




double sin_ramp(double t,double tramp,double x_i,double x_f){
    if(t<=0) return x_i;
    if(t>=tramp) return x_f;
    double s=sin(t/tramp*M_PI*0.5);
    return x_i+(x_f-x_i)*s*s;
}

double sin_square(double t,double Amp, double T,double shift){
//T is the period of the sin square function
double gamma=M_PI/T;
if (t>=shift and t<=(T+shift)) return Amp*sin(gamma*(t-shift))*sin(gamma*(t-shift));
else return 0.0;
}

void display_density(int tstp, int i, double GF_err,int ret, lattice_func &latt,double el_field){
cdmatrix R(2,2);
R.setZero();
DensityMatrix(latt.Gloc_,R);
//cout << "R:\n" << R << endl;
cout<<"inside display density with thread "<<omp_get_thread_num()<<endl;
cout << "tstp = " << tstp << endl;
cout << "iter = " << i << endl;
cout<< " GF err = " << GF_err<<endl;
cout<<"el_field =  "<<el_field<<endl;
cout << " nA = " << 2.*R(0,0);
cout << " nB = " << 2.*R(1,1); 
cout << " n_tot = " << 2.*(R(0,0)+R(1,1));
cout << " 2-n_tot = " << 2.0 - 2.*(R(0,0)+R(1,1));
cout << " n_diff = " << 2.*(R(0,0)-R(1,1));
cout << endl;
cout << " XA = " << latt.X_(0,0);
cout << " XB = " << latt.X_(1,1); 
cout << " XA + XB  = " << latt.X_(0,0)+latt.X_(1,1); 
cout << " XA - XB  = " << latt.X_(0,0)-latt.X_(1,1); 
cout << endl;
cout << " PA = " << latt.P_(0,0);
cout << " PB = " << latt.P_(1,1); 
cout<<endl;
cout <<"--------------------------------------------------"<< endl;
}


void RK_step(vector <lattice_func> & ltmp,vector <lattice_func> & latt,  vector <lattice_func> & lk1, vector <lattice_func> & lk2,double h,double mult_h, vector <double> & dW_A, vector <double> & dW_B,int NiterR, GF& Delta_avg,vector <GF> & GA, vector <GF> & GB,vector<sigma_solver*> imp_new_vector, double U, double errmax,int ipt,std::string flag,double mixing, int ret, double g_Xn,double V,int num_traj, double eA,double eB, vector <double>& error_vector, vector <int>& iter_vector,double omega0,complex <double> mult,double  f_electric_field, double beta, double gamma_phonon_bath,double g_ph_bath, double omega_fermion){	
#pragma omp parallel for
for (int i=0;i<num_traj;i++){
//I calculate k2 (h/2):
ltmp[i]=latt[i];
incr_les_F(ltmp[i],lk1[i],h*mult_h); // F_ltmp=F_t0+h/2*I_lk1; where I_lk1=I[F_t0] 
set_fermion_bath_bump(ltmp[i].Sigma_fermion_bump_,V,beta,omega_fermion);
incr_X(ltmp[i],lk1[i],h*mult_h);
incr_P(ltmp[i],lk1[i],h*mult_h);
incr_P_noise(ltmp[i],mult_h*(1./sqrt(omega0)),dW_A[i],dW_B[i]);
}
// Now that I have updated the functions, I calculate Delta as the averge of Gloc
int min_iter=0;
double max_error=10.0;
fill(error_vector.begin(), error_vector.end(), max_error);
fill(iter_vector.begin(), iter_vector.end(),min_iter);
avg_G_into_Delta(ltmp,Delta_avg,GA,GB,num_traj);
while(min_iter<NiterR and max_error>errmax){
//I calculate Delta_avg
//avg_G_into_Delta(ltmp,Delta_avg,GA,GB,num_traj);
#pragma omp parallel for
for (int i=0;i<num_traj;i++) {
dFdt_new(ltmp[i],Delta_avg,*imp_new_vector[i],U,ipt,flag,mixing,ret,g_Xn,V,num_traj,i,eA,eB,error_vector,iter_vector,g_ph_bath); // do the dmft iter and then calculate I_lk2[F_ltmp]
}
min_iter = *min_element(iter_vector.begin(), iter_vector.end());
max_error= *max_element(error_vector.begin(), error_vector.end());
}
cout<<"min iter "<<min_iter<<endl;
cout<<"max error "<<max_error<<endl;
//after convergence, I calculate the scattering integral 
#pragma omp parallel for
for (int i=0;i<num_traj;i++){	
get_ISigma_F(lk2[i],ltmp[i]);
smul_les_F(lk2[i],mult);
set_X_fx(lk2[i],ltmp[i],omega0);
set_P_fp(lk2[i],ltmp[i],f_electric_field,omega0,g_Xn,gamma_phonon_bath,*imp_new_vector[i]);
}
}

void outputXP(cdmatrix &X,cdmatrix &P, std::string filename){
	std::ofstream outputFile;
	outputFile.open(filename.c_str());
	outputFile << "# X00" << "\t" << "X11" << "\t"<< "P00" << "\t"<< "P11" << std::endl;
	outputFile <<std::setprecision(30) << X(0,0).real() << "\t" << X(1,1).real() << "\t"<< P(0,0).real() << "\t"<< P(1,1).real()<< std::endl;
	outputFile.close();
}


void readGF(GF &G, std::string filename){
	std::ifstream readFile;
	readFile.open(filename.c_str());
	double freq=0.0;
	double grid_sp=G.ngrid_;
	long N=G.ngrid_;
	cplx imag_one(0.0, 1.0);
	std::string line;
	cdmatrix Rr(2,2), Lr(2,2), Ri(2,2), Li(2,2);

	//for (long i = 0; i < N; i++)
	int i = 0;
	while(std::getline(readFile, line))
	{
		if (i >= N) 
		{
			break;
			std::cout << "wrong file!" << std::endl;
		}
		if (line[0] == '#') continue;
		std::istringstream iss(line);
		iss >> freq >> Rr(0,0) >> Ri(0,0) >> Li(0,0) >> Lr(0,0) >> 
		Rr(0,1) >> Ri(0,1) >> Li(0,1) >> Lr(0,1) >> 
		Rr(1,0) >> Ri(1,0) >> Li(1,0) >> Lr(1,0) >> 
		Rr(1,1) >> Ri(1,1) >> Li(1,1) >> Lr(1,1);
		//freq+=grid_sp;
		//G.set_Retarded(i, Rr + imag_one * Ri);
		//G.set_Lesser(i, Lr + imag_one * Li);
		G.Retarded[i] = Rr + imag_one * Ri;
		G.Lesser[i] = Lr + imag_one * Li;
		//if(i == 0) std::cout << G.Retarded[i](0,0) << '\t' << Rr(0,0) << '\t' << Ri(0,0) << '\t' << imag_one << std::endl;
		i++;
	}
	std::cout << i << " lines read from data file" <<filename<<std::endl;
	
	readFile.close();
}

void readXP(cdmatrix &X, cdmatrix &P, std::string filename){
	std::ifstream readFile;
	readFile.open(filename.c_str());
	std::string line;

	int i = 0;
	while(std::getline(readFile, line))
	{
		if (line[0] == '#') continue;
		std::istringstream iss(line);
		iss >> X(0,0) >> X(1,1) >> P(0,0) >> P(1,1);
		i++;
	}
	std::cout << i << " lines read from data file" <<filename<<std::endl;
	
	readFile.close();
}



int main(int argc, char** argv){

  omp_set_num_threads(omp_get_max_threads()); 
  cout<<"The number of threads is "<<omp_get_max_threads()<<endl;
  double mult_traj;
  double start_time(0.),end_time(0.),diff_time(0.);
  double start_time_print(0.),end_time_print(0.);
  double start_time_k1(0.),end_time_k1(0.);
  double start_time_k2(0.),end_time_k2(0.);
  double start_time_k3(0.),end_time_k3(0.);
  double start_time_k4(0.),end_time_k4(0.);
  double start_time_k5(0.),end_time_k5(0.);

  double start_delta_k2(0.), end_delta_k2(0.);
  double start_delta_k3(0.), end_delta_k3(0.);
  double start_delta_k4(0.), end_delta_k4(0.);
  double start_delta_k5(0.), end_delta_k5(0.);

  int restart, time_restart_init, time_restart_end;
  int mult_print=2; //I do not print just at the end of the tstp iterations but when it % ((tstp_end-tstp_start)/mult_print)==0
  int latt_output_all;
  double grid_spacing=0.0005; //these two values will be modified after the 
  long numOfFreq=40000;  //reading of param_boltz.in
  int it =0;
   int ret;
  int interval=100; //intervallo per stampare energia e GF,1000
  int interval_energy=500; //100
  int  first_time=1; 
  int tmp_int;
  double freq_cutoff=10.0; //this value will be modified according to the param file
  double beta, beta_ph;
  double chemPot=0.000;
  double mixing;
  long numOfIt=100;
  double U,U_i,U_f,tramp;
  int ntheta,time_to_boltzmann;
  int nt; 
  double Jbath,h;
  double lerr = 0.0; //error at equilibrium, given by the param file
  double err=1e-6; //error during time prop
  int Niter = 10,NiterR;
  //int i=0;
  int ipt;
  int outputf,param;
  double Amp_Vt,T_Vt, shift_Vt=0; //excitation through the coupling with the fermionic bath. T_Vt is already multiplied by h
  double Amp_f,T_f,shift_f=0; //coupling with the electric field; T_f is already multiplied by h
  char fname[100];
  double omega0,g_Xn,g_ph_bath; //frequency of the phonon mode, el-ph coupling
  int num_G=1;
  int val;
  double omega_fermion(0.0); //central frequency of the fermionic bath 
  double eAseed,XA=0.0,XB=0.0,PA=0.0,PB=0.0,Omega;
  double f_electric_field=0.0;
  double gamma_phonon_bath(0.0);	


  cdmatrix one(2,2);
  one.setZero();
  one(0,0)=1.0;
  one(1,1)=1.0;

  // At t=0 X and P are 0, but it is present a static field eA on A and -eA on B. This creates an asimmetry between nA and nB. Once I get these nA and nB - I am still at t=0 -, I update X(t=0).

    //double mean(0.0),std_A(0.0),std_B(0.0),dW_A(0.0),dW_B(0.0);
    default_random_engine generator;
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());

    ostringstream os;
    string filename;
	{
		if (argv[1] == NULL) 
		{
			std::cout << "no input file!" << std::endl;
			return 1;
		}
	find_param(argv[1], "__numOfFreq",tmp_int);
        numOfFreq=tmp_int;cout<<"Num of Freq "<<numOfFreq<<endl;
        find_param(argv[1], "__freq_cutoff", freq_cutoff);cout<<"freq cutoff "<<freq_cutoff<<endl;

        find_param(argv[1], "__ntheta", ntheta);
	//ntheta=2;
        find_param(argv[1], "__U_i", U_i);
        find_param(argv[1], "__U_f", U_f);
        U=U_i;
        find_param(argv[1], "__ret", ret);
        find_param(argv[1], "__tramp", tramp);	
 	find_param(argv[1], "__Jbath", Jbath);
	find_param(argv[1], "__beta", beta);
	find_param(argv[1], "__beta_ph", beta_ph);
        find_param(argv[1], "__Niter", Niter);
	find_param(argv[1], "__NiterR", NiterR);
	find_param(argv[1], "__mixing", mixing);
	find_param(argv[1], "__err", lerr);
        find_param(argv[1], "__latt_output_all", latt_output_all);
        grid_spacing=freq_cutoff*2.0/numOfFreq;
	//Jbath=2.0*grid_spacing;
	cout<<"Grid_spacing = "<<grid_spacing<<endl;
        find_param(argv[1], "__h", h);
        find_param(argv[1], "__nt", nt);
	find_param(argv[1],"__time_to_boltzmann=",time_to_boltzmann);
        find_param(argv[1], "__outputf", outputf);
        find_param(argv[1], "__ipt", ipt);
	find_param(argv[1],"__Amp_Vt",Amp_Vt);
	find_param(argv[1],"__Amp_f",Amp_f);
	find_param(argv[1],"__T_f",T_f);
	find_param(argv[1],"__shift_f",shift_f);
	find_param(argv[1],"__T_Vt",T_Vt); 
	find_param(argv[1],"__shift_Vt",shift_Vt);
	find_param(argv[1],"__param",param);
        find_param(argv[1], "__eAseed", eAseed);
        find_param(argv[1], "__g_Xn", g_Xn);
        find_param(argv[1], "__g_ph_bath", g_ph_bath);
	find_param(argv[1], "__omega0", omega0);
        find_param(argv[1], "__Omega", Omega);
        find_param(argv[1], "__omega_fermion", omega_fermion);
        find_param(argv[1], "__num_G", num_G);
        find_param(argv[1], "__gamma_phonon_bath", gamma_phonon_bath);
	find_param(argv[1], "__mult_traj",mult_traj);
        find_param(argv[1], "__restart",restart);
        find_param(argv[1], "__time_restart_init",time_restart_init);
        find_param(argv[1], "__time_restart_end",time_restart_end);

        //find_param(argv[1], "__epsBseed", epsBseed);


	    if (argc >= 3)
	    {
    		std::istringstream iss( argv[2] );
//		int val;
		if (iss >> val)
		{
		    // Conversion successful
		    cout<<"traject number "<<val<<endl;

		}
	    }

	}
	cout<<"number of G = "<<num_G<<endl;
	cout<<"eAseed = "<<eAseed<<endl;
	cout<<"g_Xn = "<<g_Xn<<endl;
        cout<<"g_ph_bath = "<<g_ph_bath<<endl;
	cout<<"beta = "<<beta<<endl;
	cout<<"omega0 = "<<omega0<<endl;
	cout<<"Omega = "<<Omega<<endl;
	cout<<"omega_fermion = "<<omega_fermion<<endl;
	cout<<"shift_Vt = "<<shift_Vt<<endl; //in terms of hopping time, and not time_steps
	cout<<"Amp_Vt = "<<Amp_Vt<<endl;
	cout<<"T_Vt = "<<T_Vt<<endl; //in terms of hopping time, and not time_steps
	cout<<"shift_f = "<<shift_f<<endl;
	cout<<"Amp_f = "<<Amp_f<<endl;
	cout<<"T_f = "<<T_f<<endl;
	cout<<"gamma_phonon_bath = "<<gamma_phonon_bath<<endl;
	cout<<"ntheta = "<<ntheta<<endl;
        cout<<"restart = "<<restart<<endl;
        cout<<"time_restart_init = "<<time_restart_init<<endl;
        cout<<"time_restart_end = "<<time_restart_end<<endl;



	string direct_out="/home/vault/mptf/mptf007h/boltzmann_CDW/TEST_ANTONIO/param_"+std::to_string(param)+"/";
	string direct_restart="/home/vault/mptf/mptf007h/boltzmann_CDW/TEST_ANTONIO/restart/";

	//string direct_out="/nfs/tfkp00/picanoan/Desktop/C++/codes/boltzmann/standalone_ness_mio/CDW/param_"+std::to_string(param)+"/";
	//string direct_restart="/nfs/tfkp00/picanoan/Desktop/C++/codes/boltzmann/standalone_ness_mio/CDW/restart/";

	cout<<"The output directory of the Boltz code is "<<direct_out<<endl;
	cout<<"The restart directory of the Boltz code is "<<direct_restart<<endl;

	ostringstream name_block;
        name_block.str("");
        name_block.clear();

	name_block<<"_U_"<<U_i<<"_g_el_ph_"<<g_Xn<<"_beta_"<<beta<<"_omega0_"<<omega0<<"_Omega_"<<Omega<<"_h_"<<h<<"_T_Vt_"<<T_Vt<<"_Amp_Vt_"<<Amp_Vt<<"_shift_Vt_"<<shift_Vt<<"_T_f_"<<T_f<<"_Amp_f_"<<Amp_f<<"_shift_f_"<<shift_f<<"_gamma_nuovo_"<<gamma_phonon_bath<<"_g_ph_bath_"<<g_ph_bath<<"_adim";


	cout<<"mult traj "<<mult_traj<<endl;
	int num_traj= (int) (mult_traj*omp_get_max_threads());
        cout<<"num traj "<<num_traj<<endl;
	vector <int>  first_time_vector(num_traj,1); 

	double mean(0.0);
	
        vector <double> std_A(num_traj,0.0), std_B(num_traj,0.0), dW_A(num_traj,0.0),dW_B(num_traj,0.);
        //vector <double> gamma_and_domega(2,0.0);
   	vector <int> time(nt+1,0);
	vector <double> f_vector(nt+1,0.);

	vector<vector<double>> occ_numberA, half_kin_energyA, pot_energyA, occ_numberB, half_kin_energyB, pot_energyB, XA_vector, PA_vector, XB_vector, PB_vector,gamma_vector, delta_omega_vector, dW_A_vector, dW_B_vector; //XXA_vector, XXXA_vector,  XXXB_vector, XPA_vector, PPA_vector, XXB_vector,PPB_vector,XPB_vector;

    lattice_func latt_atomic(ntheta,numOfFreq,grid_spacing); 
    vector <lattice_func> latt,ltmp,lk1,lk2,lk3,lk4;

    sigma_solver imp(numOfFreq);
    vector<sigma_solver*> imp_new_vector;
    sigma_solver* p = NULL;

    int size=2;	
    GF G_avg(grid_spacing, numOfFreq,size), Delta_avg(grid_spacing, numOfFreq,size), Pol_avg(grid_spacing, numOfFreq,size);

    for (int i=0;i<num_traj;i++)
    {
	    latt.push_back(latt_atomic);
	    ltmp.push_back(latt_atomic);
	    lk1.push_back(latt_atomic);
	    lk2.push_back(latt_atomic);
	    lk3.push_back(latt_atomic);
	    lk4.push_back(latt_atomic);
	    p = new sigma_solver(numOfFreq);
            imp_new_vector.push_back(p);

    }

   
     //I need a container for all the latt[i].Gloc_ 
     GF GA_atomic(Delta_avg.dgrid_,Delta_avg.ngrid_,1);
     GF GB_atomic(Delta_avg.dgrid_,Delta_avg.ngrid_,1);
     GA_atomic.clear();
     GB_atomic.clear();

     vector <GF> GA, GB;
     for (int i=0;i<num_traj;i++){
     GA.push_back(GA_atomic);
     GB.push_back(GB_atomic);
     }

    double V;
    complex <double> mult=cplx(0.0,-1.0); 

	vector <double> error_vector(num_traj,10.);
	double max_error(10.);
	vector <int> iter_vector(num_traj,0);
	int min_iter(0);	



    if (restart==0){
    cout<<"restart "<<restart<<endl;


    U=U_i;
    V=0.;
    f_electric_field=0.;

	occ_numberA.resize(num_traj, std::vector<double>(nt+1, 0.));
        half_kin_energyA.resize(num_traj, std::vector<double>(nt+1, 0.));
        pot_energyA.resize(num_traj, std::vector<double>(nt+1, 0.));
        occ_numberB.resize(num_traj, std::vector<double>(nt+1, 0.));
        half_kin_energyB.resize(num_traj, std::vector<double>(nt+1, 0.));
        pot_energyB.resize(num_traj, std::vector<double>(nt+1, 0.));
        XA_vector.resize(num_traj, std::vector<double>(nt+1, 0.));
        PA_vector.resize(num_traj, std::vector<double>(nt+1, 0.));
        XB_vector.resize(num_traj, std::vector<double>(nt+1, 0.));
        PB_vector.resize(num_traj, std::vector<double>(nt+1, 0.));
        gamma_vector.resize(num_traj, std::vector<double>(nt+1, 0.));
        delta_omega_vector.resize(num_traj, std::vector<double>(nt+1, 0.));
        dW_A_vector.resize(num_traj, std::vector<double>(nt+1, 0.));
        dW_B_vector.resize(num_traj, std::vector<double>(nt+1, 0.));

        time.resize(nt+1, 0);
        f_vector.resize(nt+1, 0.);




 
    for (int i=0;i<num_traj;i++){
    setBath_diag(Jbath,chemPot,beta,latt[i].Bath_);
    set_phonon_bath(omega0,beta,latt[i].D0_Phonon_,i);
    set_fermion_bath_bump(latt[i].Sigma_fermion_bump_,V,beta,omega_fermion);
    }
    
    for (int i=0;i<num_traj;i++)
    {
    double eA=eAseed;
    double eB=-eAseed;

    latt[i].init_distribution_equi_nonkdepe(beta);
    latt[i].Gloc_.clear();
    latt[i].Polariz_.clear();
    latt[i].Sigma_.clear();
    latt[i].Sigma_.incr(latt[i].Bath_);       
    latt[i].Sigma_.incr(latt[i].Sigma_fermion_bump_);
    latt[i].Delta_.clear();
   //latt[i].set_Hartree(U);
    latt[i].X_(0,0)=0.;
    latt[i].X_(1,1)=0.;
    latt[i].P_(0,0)=0.;
    latt[i].P_(1,1)=0.;
 
    latt[i].set_hloc_diag(eA,eB,g_Xn);


    if (ret==1) {
	    latt[i].get_gloc_nonkdepe_ret(); 
            latt[i].set_distribution_function();
    }
    else latt[i].get_gloc_nonkdepe();
    }

        //equilibrium iteration
	start_time=omp_get_wtime();
	{
	double GF_err(10.), X_err(10.); //variables declared inside the parallel region are private
	int iter(0);	
	double eA=eAseed;
        double eB=-eAseed;
	int i=0;
	while( (GF_err > lerr or X_err>lerr) and  iter<Niter){
        if (iter>=1) {eA=0.0; eB=0.0;}
	GF_err=dmft_iteration(iter,latt,imp,ipt,U,ret,mixing,it,"R1_",g_Xn,V,num_traj,i,eA,eB,g_ph_bath);
	get_polarizability(latt[i].Gloc_,latt[i].Polariz_,imp.solver_);
	updateX_equilibrium_no_noise(omega0,g_Xn,Omega,latt[i],iter,mixing,X_err,U,Amp_f,T_f,shift_f);	
	display_density(it,iter,GF_err,ret,latt[i],f_electric_field);
	iter++;
    	}
	cout<<"------------------- end of EQUILIBRIUM iter for traj "<<i<<"; X "<<latt[i].X_(0,0)<<endl;
	//latt[i].Sigma_.incr(latt[i].Bath_,-1.0);  //I remove the fermionic bath       
	}
        end_time=omp_get_wtime();
        cout<<"equiilibrium ; execution time: "<<(end_time-start_time)*1000<<endl;


	
	cout<<endl;
	cout<<"------------------- end of EQUILIBRIUM iteration -----------------"<<endl;
	cout<<endl;	

	// I copy the functions I have calculated for traj 0 at equilibrium into the others
	for (int i=1;i<num_traj;i++) latt[i]=latt[0]; 


        os.str("");
        os.clear();
        os <<direct_out<<"temp_G_equilibrium"<<name_block.str()<<"_num_traj_"<<num_traj<<".out";
        filename= os.str();
        outputGF_new(latt[omp_get_thread_num()].Gloc_,filename);
	cout<<"------  I am printing "<<filename<<"--------------"<<endl;

        os.str("");
        os.clear();
        os <<direct_out<<"temp_Polariz_equilibrium"<<name_block.str()<<"_num_traj_"<<num_traj<<".out";
        filename= os.str();
	outputGF_new(latt[omp_get_thread_num()].Polariz_,filename);
	cout<<"------  I am printing "<<filename<<"--------------"<<endl;

        os.str("");
        os.clear();
        os <<direct_out<<"temp_Delta_equilibrium"<<name_block.str()<<"_num_traj_"<<num_traj<<".out";
        filename= os.str();
	outputGF_new(latt[omp_get_thread_num()].Delta_,filename);
	cout<<"------  I am printing "<<filename<<"--------------"<<endl;



	double eA=0.0; //variables declared outside parallel loop are shared but this is fine since they will not be modified by any thread
	double eB=0.0;
       
	//return 0;
	
	cout << "it = " << it << " t0 = "  << it*h << " U = " << U <<" nt = "<<nt<<"  h = "<<h<<endl;

	time[it]=it;
	f_vector[it]=electric_field(0.0,Omega,Amp_f,T_f,shift_f);

	
	//In case I want to save all the trajectories and not just the mean values
	#pragma omp parallel for
	for (int i=0;i<num_traj;i++){
	cdmatrix R(2,2);      R.setZero();
	cdmatrix R_Ekin(2,2); R_Ekin.setZero();
	cdmatrix R_Epot(2,2); R_Epot.setZero();
	cdmatrix Xold(2,2);   Xold.setZero();
	cdmatrix Pold(2,2);   Pold.setZero();
	//Calculate Ekin and Epot, to be adapted to the 2x2 case
	KinEnergy(latt[i].Gloc_,latt[i].Delta_,R_Ekin);
	PotEnergy(latt[i].Gloc_,latt[i].Sigma_,R_Epot);
	DensityMatrix(latt[i].Gloc_,R);
	occ_numberA[i][it]=R(0,0).real();	 
	half_kin_energyA[i][it]=R_Ekin(0,0).real(); 
	pot_energyA[i][it]=R_Epot(0,0).real();//+U*(R(0,0).real()-0.5); 
	occ_numberB[i][it]=R(1,1).real();	 
	half_kin_energyB[i][it]=R_Ekin(1,1).real(); 
	pot_energyB[i][it]=R_Epot(1,1).real();//+U*(R(0,0).real()-0.5); 
	XA_vector[i][it]=latt[i].X_(0,0).real();
	PA_vector[i][it]=latt[i].P_(0,0).real();
	XB_vector[i][it]=latt[i].X_(1,1).real();
	PB_vector[i][it]=latt[i].P_(1,1).real();
	vector <double> gamma_and_domega(2,0.0);
	gamma_and_domega=get_gamma_and_delta_omega(latt[i].Polariz_,g_Xn,omega0);
	gamma_vector[i][it]=gamma_and_domega[0];
	delta_omega_vector[i][it]=gamma_and_domega[1];
        dW_A_vector[i][it]=dW_A[i];
        dW_B_vector[i][it]=dW_B[i];
	}
	
     for(it=1;it<=nt;it++){

        double t0=(it-1)*h;
	start_time=omp_get_wtime();

	start_time_k1=omp_get_wtime();
	#pragma omp parallel for
	for (int i=0;i<num_traj;i++){

	get_polarizability(latt[i].Gloc_,latt[i].Polariz_,imp_new_vector[i]->solver_);
	//std_A[i]=noise_coeff*sqrt(h)*g_Xn*sqrt(abs(latt[i].Polariz_.Lesser[0](0,0).imag()));
        //std_B[i]=noise_coeff*sqrt(h)*g_Xn*sqrt(abs(latt[i].Polariz_.Lesser[0](1,1).imag())); 
	std_A[i]=sqrt(h)*sqrt(omega0*g_Xn*g_Xn*abs(latt[i].Polariz_.Lesser[0](0,0).imag())+2.*gamma_phonon_bath/beta);
	std_B[i]=sqrt(h)*sqrt(omega0*g_Xn*g_Xn*abs(latt[i].Polariz_.Lesser[0](1,1).imag())+2.*gamma_phonon_bath/beta);
	std::normal_distribution <double> gaussian_distribution_A(mean,std_A[i]);
        std::normal_distribution <double> gaussian_distribution_B(mean,std_B[i]);
        dW_A[i]=gaussian_distribution_A(generator);
        dW_B[i]=gaussian_distribution_B(generator);
        //dW_A[i]=0.5;
        //dW_B[i]=0.5;


	//dF/dt=I[F]	
        //I calculate k1: 
	ltmp[i]=latt[i];
	f_electric_field=electric_field(t0,Omega,Amp_f,T_f,shift_f);
        get_ISigma_F(lk1[i],ltmp[i]); // I put in lk1 the scattering integral at t0
	smul_les_F(lk1[i],mult);
	set_X_fx(lk1[i],ltmp[i],omega0);
	set_P_fp(lk1[i],ltmp[i],f_electric_field,omega0,g_Xn,gamma_phonon_bath,*imp_new_vector[i]);
	}
	end_time_k1=omp_get_wtime();

        //I calculate k2 (h/2):   	
	start_time_k2=omp_get_wtime();
	U=sin_ramp(t0+0.5*h,tramp,U_i,U_f);
	f_electric_field=electric_field(t0+0.5*h,Omega,Amp_f,T_f,shift_f);
	V=sin_square(t0+0.5*h,Amp_Vt,T_Vt,shift_Vt); // t -> t +h/2
	RK_step(ltmp,latt,lk1,lk2,h,0.5,dW_A,dW_B,NiterR,Delta_avg,GA,GB,imp_new_vector,U,lerr,ipt,"R2_",mixing,ret,g_Xn,V,num_traj,eA,eB,error_vector,iter_vector,omega0,mult,f_electric_field,beta,gamma_phonon_bath,g_ph_bath,omega_fermion);
	end_time_k2=omp_get_wtime();

        //I calculate k3 (h/2):   	
	start_time_k3=omp_get_wtime();
	U=sin_ramp(t0+0.5*h,tramp,U_i,U_f);
	f_electric_field=electric_field(t0+0.5*h,Omega,Amp_f,T_f,shift_f);
	V=sin_square(t0+0.5*h,Amp_Vt,T_Vt,shift_Vt); // t -> t +h/2
	RK_step(ltmp,latt,lk2,lk3,h,0.5,dW_A,dW_B,NiterR,Delta_avg,GA,GB,imp_new_vector,U,lerr,ipt,"R3_",mixing,ret,g_Xn,V,num_traj,eA,eB,error_vector,iter_vector,omega0,mult,f_electric_field,beta,gamma_phonon_bath,g_ph_bath,omega_fermion);
	end_time_k3=omp_get_wtime();

        //I calculate k4 (h):   	
	start_time_k4=omp_get_wtime();
	U=sin_ramp(t0+1.0*h,tramp,U_i,U_f);
	f_electric_field=electric_field(t0+1.0*h,Omega,Amp_f,T_f,shift_f);
	V=sin_square(t0+1.0*h,Amp_Vt,T_Vt,shift_Vt); // t -> t +h/2
	RK_step(ltmp,latt,lk3,lk4,h,1.0,dW_A,dW_B,NiterR,Delta_avg,GA,GB,imp_new_vector,U,lerr,ipt,"R4_",mixing,ret,g_Xn,V,num_traj,eA,eB,error_vector,iter_vector,omega0,mult,f_electric_field,beta,gamma_phonon_bath,g_ph_bath,omega_fermion);
	end_time_k4=omp_get_wtime();

	//now I consider the next tstp
	start_time_k5=omp_get_wtime();
	V=sin_square(t0+1.0*h,Amp_Vt,T_Vt,shift_Vt); // t -> t +h
        U=sin_ramp(t0+1.0*h,tramp,U_i,U_f);
	f_electric_field=electric_field(t0+1.0*h,Omega,Amp_f,T_f,shift_f);
	#pragma omp parallel for
	for (int i=0;i<num_traj;i++){	
	set_fermion_bath_bump(latt[i].Sigma_fermion_bump_,V,beta,omega_fermion);
	//update F
        set_les_F(ltmp[i],lk1[i]);
        incr_les_F(ltmp[i],lk2[i],2.0);
        incr_les_F(ltmp[i],lk3[i],2.0);
        incr_les_F(ltmp[i],lk4[i],1.0);
        incr_les_F(latt[i],ltmp[i],h/6.0);//update F_latt = F_latt+ h/6 * (I_lk1+2* I_lk2 + 2 * I_lk3+ I_lk4)
  	//F(t0+h)= F(t0)+ h/6 * (I_lk1+2* I_lk2 + 2 * I_lk3+ I_lk4  )

	//update X
        set_X(ltmp[i],lk1[i]);
        incr_X(ltmp[i],lk2[i],2.0);
        incr_X(ltmp[i],lk3[i],2.0);
        incr_X(ltmp[i],lk4[i],1.0);
	incr_X(latt[i],ltmp[i],h/6.0);

	//update P
        set_P(ltmp[i],lk1[i]);
        incr_P(ltmp[i],lk2[i],2.0);
        incr_P(ltmp[i],lk3[i],2.0);
        incr_P(ltmp[i],lk4[i],1.0);
        incr_P(latt[i],ltmp[i],h/6.0);
	incr_P_noise(latt[i],1.0*(1./sqrt(omega0)),dW_A[i],dW_B[i]);
	}


       	// Now that I have updated the functions, I make the DMFT iter at tstp+h
	min_iter=0;
	max_error=10.0;
        fill(error_vector.begin(), error_vector.end(), max_error);
        fill(iter_vector.begin(), iter_vector.end(),min_iter);
	avg_G_into_Delta(latt,Delta_avg,GA,GB,num_traj);
	while(min_iter<NiterR and max_error>lerr){
	//I calculate Delta_avg
	start_delta_k5=omp_get_wtime();
	//avg_G_into_Delta(latt,Delta_avg,GA,GB,num_traj);
	end_delta_k5=omp_get_wtime();
	#pragma omp parallel for
	for (int i=0;i<num_traj;i++) {
	dFdt_new(latt[i],Delta_avg,*imp_new_vector[i],U,ipt,"R5_",mixing,ret,g_Xn,V,num_traj,i,eA,eB,error_vector,iter_vector,g_ph_bath); // do the dmft iter and then calculate I_lk2[F_ltmp]
	}
	min_iter = *min_element(iter_vector.begin(), iter_vector.end());
	max_error= *max_element(error_vector.begin(), error_vector.end());
	}
	cout<<"min iter "<<min_iter<<endl;
	cout<<"max error "<<max_error<<endl;
	end_time_k5=omp_get_wtime();


	time[it]=it;
	f_vector[it]=f_electric_field;
	cout << "it = " << it << " t0 = "  << it*h << " U = " << U <<" V "<<V<<" nt = "<<nt<<"  h = "<<h<<endl;
        display_density(it,0,max_error,ret,latt[0],f_electric_field);

	start_time_print=omp_get_wtime();
	#pragma omp parallel for
	for (int i=0;i<num_traj;i++) {
	get_polarizability(latt[i].Gloc_,latt[i].Polariz_,imp_new_vector[i]->solver_); //I update the polarizability
	cdmatrix R(2,2);      R.setZero();
	cdmatrix R_Ekin(2,2); R_Ekin.setZero();
	cdmatrix R_Epot(2,2); R_Epot.setZero();
	KinEnergy(latt[i].Gloc_,latt[i].Delta_,R_Ekin);
	PotEnergy(latt[i].Gloc_,latt[i].Sigma_,R_Epot);
	DensityMatrix(latt[i].Gloc_,R);
	occ_numberA[i][it]=R(0,0).real();	 
	half_kin_energyA[i][it]=R_Ekin(0,0).real(); 
	pot_energyA[i][it]=R_Epot(0,0).real();
	occ_numberB[i][it]=R(1,1).real();	 
	half_kin_energyB[i][it]=R_Ekin(1,1).real(); 
	pot_energyB[i][it]=R_Epot(1,1).real();
	XA_vector[i][it]=latt[i].X_(0,0).real();
	PA_vector[i][it]=latt[i].P_(0,0).real();
	XB_vector[i][it]=latt[i].X_(1,1).real();
	PB_vector[i][it]=latt[i].P_(1,1).real();    
	vector <double> gamma_and_domega(2,0.0);
	gamma_and_domega=get_gamma_and_delta_omega(latt[i].Polariz_,g_Xn,omega0);
	gamma_vector[i][it]=gamma_and_domega[0];
	delta_omega_vector[i][it]=gamma_and_domega[1];
        dW_A_vector[i][it]=dW_A[i];
        dW_B_vector[i][it]=dW_B[i];	
      	} 
	end_time_print=omp_get_wtime();

        if(it%interval_energy==0 ){   

	#pragma omp parallel for
	for (int i=0;i<num_traj;i++){
	ostringstream os;
    	string filename;
	os.str("");
        os.clear();
        os<<direct_out<<"temp_occupation_plus_energy"<<name_block.str()<<"_traj_num_"<<i<<".out";
        filename= os.str();

	print_to_file_single_traj(time,i,occ_numberA,half_kin_energyA,pot_energyA,occ_numberB,half_kin_energyB,pot_energyB,f_vector,XA_vector,PA_vector,XB_vector,PB_vector,filename,interval_energy,it,first_time_vector[i],gamma_vector,delta_omega_vector,dW_A_vector,dW_B_vector);//,XXA_vector,XPA_vector,PPA_vector,XXB_vector,XPB_vector,PPB_vector,XXXA_vector,XXXB_vector);
	}
	}

	if (it%interval==0){
	//Then I calculate the average of Gloc, Delta, Pol, over the different trajectories at that given time it
	G_avg.clear();
	//Delta_avg.clear();
	Pol_avg.clear();

      	for (int i=0;i<num_traj;i++){
	 G_avg.incr(latt[i].Gloc_,1./num_traj);
	 //Delta_avg.incr(latt[i].Delta_,1./num_traj);
	 Pol_avg.incr(latt[i].Polariz_,1./num_traj);
	}
        os.str("");
        os.clear();
        os <<direct_out<<"temp_G_t" << it <<name_block.str()<<"_num_traj_"<<num_traj<<".out";
        filename= os.str();
        outputGF_new(G_avg,filename);
	cout<<"------  I am printing "<<filename<<"--------------"<<endl;

	
        os.str("");
        os.clear();
	os <<direct_out<<"temp_Delta_t" << it <<name_block.str()<<"_num_traj_"<<num_traj<<".out";
        filename= os.str();
        //outputGF_new(Delta_avg,filename);
	cout<<"------  I am printing "<<filename<<"--------------"<<endl;


        os.str("");
        os.clear();
        os <<direct_out<<"temp_Polariz_t" << it <<name_block.str()<<"_num_traj_"<<num_traj<<".out";
        filename= os.str();
        outputGF_new(Pol_avg,filename);
	cout<<"------  I am printing "<<filename<<"--------------"<<endl;

        }  //end of the printing loop
	

	end_time=omp_get_wtime();
        cout<<"tstp "<<it<<"; execution time: "<<(end_time-start_time)*1000<<endl;
	cout<<"tstp "<<it<<"; time to print in percentage : "<<(end_time_print-start_time_print)/(end_time-start_time)*100<<endl;
	cout<<"tstp "<<it<<"; time to k1 in percentage : "<<(end_time_k1-start_time_k1)/(end_time-start_time)*100<<endl;
	cout<<"tstp "<<it<<"; time to k2 in percentage : "<<(end_time_k2-start_time_k2)/(end_time-start_time)*100<<endl;
	cout<<"tstp "<<it<<"; time to k3 in percentage : "<<(end_time_k3-start_time_k3)/(end_time-start_time)*100<<endl;
	cout<<"tstp "<<it<<"; time to k4 in percentage : "<<(end_time_k4-start_time_k4)/(end_time-start_time)*100<<endl;
	cout<<"tstp "<<it<<"; time to k5 in percentage : "<<(end_time_k5-start_time_k5)/(end_time-start_time)*100<<endl;
	cout<<endl;

	if (   it % ( nt / mult_print  )==0){
	for (int i=0;i<num_traj;i++){
        os.str("");
        os.clear();
        os <<direct_restart<<"temp_G_t" << it <<name_block.str()<<"_traj_num_"<<i<<".out";
        filename= os.str();
        outputGF(latt[i].Gloc_,filename);
	cout<<"------  I am printing "<<filename<<"--------------"<<endl;


        os.str("");
        os.clear();
        os <<direct_restart<<"temp_XP_t"<< it <<name_block.str()<<"_traj_num_"<<i<<".out";
        filename= os.str();
        outputXP(latt[i].X_,latt[i].P_,filename);
	cout<<"------  I am printing "<<filename<<"--------------"<<endl<<endl;

	}
	}

        }  //end of cycle over tstp	

	//At the end of the restart==0 condition, I need to print all the GF I may need to restart
	for (int i=0;i<num_traj;i++){
        os.str("");
        os.clear();
        os <<direct_restart<<"temp_G_t" << nt <<name_block.str()<<"_traj_num_"<<i<<".out";
        filename= os.str();
        outputGF(latt[i].Gloc_,filename);
	cout<<"------  I am printing "<<filename<<"--------------"<<endl;

        os.str("");
        os.clear();
        os <<direct_restart<<"temp_XP_t" << nt <<name_block.str()<<"_traj_num_"<<i<<".out";
        filename= os.str();
        outputXP(latt[i].X_,latt[i].P_,filename);
	cout<<"------  I am printing "<<filename<<"--------------"<<endl<<endl;

	}


	}//end restart==0 condition

	else{//restart==1 condition

	cout<<"restart "<<restart<<endl;

        double t0=(time_restart_init)*h; //if time_restart_init coincides with nt, this is the last tstp calculated in restart==0
	U=sin_ramp(t0,tramp,U_i,U_f);
	f_electric_field=electric_field(t0,Omega,Amp_f,T_f,shift_f);
	V=sin_square(t0,Amp_Vt,T_Vt,shift_Vt); 
	cout<<"U : "<<U<<", f : "<<f_electric_field<<", V : "<<V<<", g_Xn : "<<g_Xn<<endl;
	double eA=0.;
	double eB=0.;	

	occ_numberA.resize(num_traj, std::vector<double>( time_restart_end-time_restart_init +1, 0.));
        half_kin_energyA.resize(num_traj, std::vector<double>( time_restart_end-time_restart_init  +1, 0.));
        pot_energyA.resize(num_traj, std::vector<double>( time_restart_end-time_restart_init  +1, 0.));
        occ_numberB.resize(num_traj, std::vector<double>( time_restart_end-time_restart_init +1, 0.));
        half_kin_energyB.resize(num_traj, std::vector<double>(  time_restart_end-time_restart_init  +1, 0.));
        pot_energyB.resize(num_traj, std::vector<double>( time_restart_end-time_restart_init  +1, 0.));
        XA_vector.resize(num_traj, std::vector<double>(   time_restart_end-time_restart_init  +1, 0.));
        PA_vector.resize(num_traj, std::vector<double>(   time_restart_end-time_restart_init  +1, 0.));
        XB_vector.resize(num_traj, std::vector<double>(   time_restart_end-time_restart_init  +1, 0.));
        PB_vector.resize(num_traj, std::vector<double>(   time_restart_end-time_restart_init  +1, 0.));
        gamma_vector.resize(num_traj, std::vector<double>(  time_restart_end-time_restart_init+1, 0.));
        delta_omega_vector.resize(num_traj, std::vector<double>(  time_restart_end-time_restart_init +1, 0.));
        dW_A_vector.resize(num_traj, std::vector<double>(  time_restart_end-time_restart_init  +1, 0.));
        dW_B_vector.resize(num_traj, std::vector<double>(  time_restart_end-time_restart_init  +1, 0.));
        time.resize( time_restart_end-time_restart_init  +1, 0);
        f_vector.resize(time_restart_end-time_restart_init  +1, 0);



	//As first thing in the restarting procedure, I upload all the Gloc, and I calculate the Delta_avg
	for (int i=0;i<num_traj;i++){
        os.str("");
        os.clear();
        os <<direct_restart<<"temp_G_t" <<time_restart_init <<name_block.str()<<"_traj_num_"<<i<<".out";
        filename= os.str();
        readGF(latt[i].Gloc_,filename);
	cout<<"------  I am reading "<<filename<<"--------------"<<endl;

	latt[i].set_Floc_from_Gloc();

        os.str("");
        os.clear();
        os <<direct_restart<<"temp_XP_t" << time_restart_init <<name_block.str()<<"_traj_num_"<<i<<".out";
        filename= os.str();
        readXP(latt[i].X_,latt[i].P_,filename);
	cout<<"------  I am reading "<<filename<<"--------------"<<endl<<endl;
	}

	//Amp_f=0;//aggiunto

	//Then I calculate Delta
	avg_G_into_Delta(latt,Delta_avg,GA,GB,num_traj);
	display_density(time_restart_init,0,max_error,ret,latt[0],f_electric_field);

	//latt[0].X_(0,0)*= -1.;
	//latt[0].P_(0,0)*= -1.;

	#pragma omp parallel for
	for (int i=0;i<num_traj;i++){

   	 setBath_diag(Jbath,chemPot,beta,latt[i].Bath_); //even if I do not need in the restart section the flat fermionic bath
    	 set_phonon_bath(omega0,beta,latt[i].D0_Phonon_,i);
         set_fermion_bath_bump(latt[i].Sigma_fermion_bump_,V,beta,omega_fermion);
    	 latt[i].set_Hartree(U);
         latt[i].set_hloc_diag(eA,eB,g_Xn);	
	 GF Gweiss=latt[i].Gloc_;
	 latt[i].Sigma_.clear();

	if (U>0){ 
	   get_gweiss_from_Delta(latt[i].Gloc_,Gweiss,latt[i].Sigma_,Delta_avg,latt[i].Hartree_,latt[i].hloc_);
	   updateSE1(U,Gweiss,latt[i].Sigma_,imp_new_vector[i]->solver_);
	}

	if (g_Xn>0){
	   set_Sigma_Phonon_G_D(g_ph_bath,latt[i].Gloc_,latt[i].D0_Phonon_,latt[i].Sigma_G_D0,imp_new_vector[i]->solver_);
	   add_Sigma_phonon_to_Sigma2_ipt(latt[i].Sigma_,latt[i].Sigma_G_D0);
	}

	if (V>0) add_Sigma_fermion_bump_to_Sigma2_ipt(latt[i].Sigma_,latt[i].Sigma_fermion_bump_);

	if (ret==1){
	latt[i].get_gloc_nonkdepe_ret();
	latt[i].set_distribution_function();
	}
	else{
	latt[i].get_gloc_nonkdepe();
	}

       	}

	cout<<"restart end - restart init = "<< time_restart_end-time_restart_init<<endl;
	cout<<"*****************"<<endl;
	//display_density(time_restart_init,0,max_error,ret,latt[0],f_electric_field);

        for(int it=time_restart_init+1;it<=time_restart_end;it++){

        double t0=(it-1)*h;
	cout<<"it "<<it<<" t0 "<<t0<<endl;
	start_time=omp_get_wtime();
	start_time_k1=omp_get_wtime();
	cout<<"t0 = "<<t0<<endl;
	
	#pragma omp parallel for
	for (int i=0;i<num_traj;i++){
	get_polarizability(latt[i].Gloc_,latt[i].Polariz_,imp_new_vector[i]->solver_);
	//std_A[i]=noise_coeff*sqrt(h)*g_Xn*sqrt(abs(latt[i].Polariz_.Lesser[0](0,0).imag()));
        //std_B[i]=noise_coeff*sqrt(h)*g_Xn*sqrt(abs(latt[i].Polariz_.Lesser[0](1,1).imag())); 
	std_A[i]=sqrt(h)*sqrt(omega0*g_Xn*g_Xn*abs(latt[i].Polariz_.Lesser[0](0,0).imag())+2.*gamma_phonon_bath/beta);
	std_B[i]=sqrt(h)*sqrt(omega0*g_Xn*g_Xn*abs(latt[i].Polariz_.Lesser[0](1,1).imag())+2.*gamma_phonon_bath/beta);
	std::normal_distribution <double> gaussian_distribution_A(mean,std_A[i]);
        std::normal_distribution <double> gaussian_distribution_B(mean,std_B[i]);
        dW_A[i]=gaussian_distribution_A(generator);
        dW_B[i]=gaussian_distribution_B(generator);
        //dW_A[i]=0.5;
        //dW_B[i]=0.5;


	//dF/dt=I[F]	
        //I calculate k1: 
	ltmp[i]=latt[i];
	f_electric_field=electric_field(t0,Omega,Amp_f,T_f,shift_f);
        get_ISigma_F(lk1[i],ltmp[i]); // I put in lk1 the scattering integral at t0
	smul_les_F(lk1[i],mult);
	set_X_fx(lk1[i],ltmp[i],omega0);
	set_P_fp(lk1[i],ltmp[i],f_electric_field,omega0,g_Xn,gamma_phonon_bath,*imp_new_vector[i]);
	}
	end_time_k1=omp_get_wtime();
	//display_density(it,0,max_error,ret,lk1[0],f_electric_field);


        //I calculate k2 (h/2):   	
	start_time_k2=omp_get_wtime();
	U=sin_ramp(t0+0.5*h,tramp,U_i,U_f);
	f_electric_field=electric_field(t0+0.5*h,Omega,Amp_f,T_f,shift_f);
	V=sin_square(t0+0.5*h,Amp_Vt,T_Vt,shift_Vt); // t -> t +h/2
	RK_step(ltmp,latt,lk1,lk2,h,0.5,dW_A,dW_B,NiterR,Delta_avg,GA,GB,imp_new_vector,U,lerr,ipt,"R2_",mixing,ret,g_Xn,V,num_traj,eA,eB,error_vector,iter_vector,omega0,mult,f_electric_field,beta,gamma_phonon_bath,g_ph_bath,omega_fermion);
	end_time_k2=omp_get_wtime();
	//display_density(it,0,max_error,ret,lk2[0],f_electric_field);


        //I calculate k3 (h/2):   	
	start_time_k3=omp_get_wtime();
	U=sin_ramp(t0+0.5*h,tramp,U_i,U_f);
	f_electric_field=electric_field(t0+0.5*h,Omega,Amp_f,T_f,shift_f);
	V=sin_square(t0+0.5*h,Amp_Vt,T_Vt,shift_Vt); // t -> t +h/2
	RK_step(ltmp,latt,lk2,lk3,h,0.5,dW_A,dW_B,NiterR,Delta_avg,GA,GB,imp_new_vector,U,lerr,ipt,"R3_",mixing,ret,g_Xn,V,num_traj,eA,eB,error_vector,iter_vector,omega0,mult,f_electric_field,beta,gamma_phonon_bath,g_ph_bath,omega_fermion);
	end_time_k3=omp_get_wtime();
	//display_density(it,0,max_error,ret,lk3[0],f_electric_field);


        //I calculate k4 (h):   	
	start_time_k4=omp_get_wtime();
	U=sin_ramp(t0+1.0*h,tramp,U_i,U_f);
	f_electric_field=electric_field(t0+1.0*h,Omega,Amp_f,T_f,shift_f);
	V=sin_square(t0+1.0*h,Amp_Vt,T_Vt,shift_Vt); // t -> t +h/2
	RK_step(ltmp,latt,lk3,lk4,h,1.0,dW_A,dW_B,NiterR,Delta_avg,GA,GB,imp_new_vector,U,lerr,ipt,"R4_",mixing,ret,g_Xn,V,num_traj,eA,eB,error_vector,iter_vector,omega0,mult,f_electric_field,beta,gamma_phonon_bath,g_ph_bath,omega_fermion);
	end_time_k4=omp_get_wtime();
	//display_density(it,0,max_error,ret,lk4[0],f_electric_field);


	//now I consider the next tstp
	start_time_k5=omp_get_wtime();
	V=sin_square(t0+1.0*h,Amp_Vt,T_Vt,shift_Vt); // t -> t +h
        U=sin_ramp(t0+1.0*h,tramp,U_i,U_f);
	f_electric_field=electric_field(t0+1.0*h,Omega,Amp_f,T_f,shift_f);
	#pragma omp parallel for
	for (int i=0;i<num_traj;i++){	
	set_fermion_bath_bump(latt[i].Sigma_fermion_bump_,V,beta,omega_fermion);
	//update F
        set_les_F(ltmp[i],lk1[i]);
        incr_les_F(ltmp[i],lk2[i],2.0);
        incr_les_F(ltmp[i],lk3[i],2.0);
        incr_les_F(ltmp[i],lk4[i],1.0);
        incr_les_F(latt[i],ltmp[i],h/6.0);//update F_latt = F_latt+ h/6 * (I_lk1+2* I_lk2 + 2 * I_lk3+ I_lk4)
  	//F(t0+h)= F(t0)+ h/6 * (I_lk1+2* I_lk2 + 2 * I_lk3+ I_lk4  )

	//update X
        set_X(ltmp[i],lk1[i]);
        incr_X(ltmp[i],lk2[i],2.0);
        incr_X(ltmp[i],lk3[i],2.0);
        incr_X(ltmp[i],lk4[i],1.0);
	incr_X(latt[i],ltmp[i],h/6.0);

	//update P
        set_P(ltmp[i],lk1[i]);
        incr_P(ltmp[i],lk2[i],2.0);
        incr_P(ltmp[i],lk3[i],2.0);
        incr_P(ltmp[i],lk4[i],1.0);
        incr_P(latt[i],ltmp[i],h/6.0);
	incr_P_noise(latt[i],1.0*(1./sqrt(omega0)),dW_A[i],dW_B[i]);
	}


       	// Now that I have updated the functions, I make the DMFT iter at tstp+h
	min_iter=0;
	max_error=10.0;
        fill(error_vector.begin(), error_vector.end(), max_error);
        fill(iter_vector.begin(), iter_vector.end(),min_iter);
	avg_G_into_Delta(latt,Delta_avg,GA,GB,num_traj);
	while(min_iter<NiterR and max_error>lerr){
	//I calculate Delta_avg
	start_delta_k5=omp_get_wtime();
	//avg_G_into_Delta(latt,Delta_avg,GA,GB,num_traj);
	end_delta_k5=omp_get_wtime();
	#pragma omp parallel for
	for (int i=0;i<num_traj;i++) {
	dFdt_new(latt[i],Delta_avg,*imp_new_vector[i],U,ipt,"R5_",mixing,ret,g_Xn,V,num_traj,i,eA,eB,error_vector,iter_vector,g_ph_bath); // do the dmft iter and then calculate I_lk2[F_ltmp]
	}
	min_iter = *min_element(iter_vector.begin(), iter_vector.end());
	max_error= *max_element(error_vector.begin(), error_vector.end());
	}
	cout<<"min iter "<<min_iter<<endl;
	cout<<"max error "<<max_error<<endl;
	end_time_k5=omp_get_wtime();

	time[it-time_restart_init-1]=it;
	f_vector[it-time_restart_init-1]=f_electric_field;
	cout << "it = " << it << " t0 = "  << t0 << " U = " << U <<" V "<<V<<" nt = "<<nt<<"  h = "<<h<<endl;
        display_density(it,0,max_error,ret,latt[0],f_electric_field);


	start_time_print=omp_get_wtime();
	#pragma omp parallel for
	for (int i=0;i<num_traj;i++) {
	get_polarizability(latt[i].Gloc_,latt[i].Polariz_,imp_new_vector[i]->solver_); //I update the polarizability
	cdmatrix R(2,2);      R.setZero();
	cdmatrix R_Ekin(2,2); R_Ekin.setZero();
	cdmatrix R_Epot(2,2); R_Epot.setZero();
	KinEnergy(latt[i].Gloc_,latt[i].Delta_,R_Ekin);
	PotEnergy(latt[i].Gloc_,latt[i].Sigma_,R_Epot);
	DensityMatrix(latt[i].Gloc_,R);
	occ_numberA[i][it-time_restart_init-1]=R(0,0).real();	 
	half_kin_energyA[i][it-time_restart_init-1]=R_Ekin(0,0).real(); 
	pot_energyA[i][it-time_restart_init-1]=R_Epot(0,0).real();
	occ_numberB[i][it-time_restart_init-1]=R(1,1).real();	 
	half_kin_energyB[i][it-time_restart_init-1]=R_Ekin(1,1).real(); 
	pot_energyB[i][it-time_restart_init-1]=R_Epot(1,1).real();
	XA_vector[i][it-time_restart_init-1]=latt[i].X_(0,0).real();
	PA_vector[i][it-time_restart_init-1]=latt[i].P_(0,0).real();
	XB_vector[i][it-time_restart_init-1]=latt[i].X_(1,1).real();
	PB_vector[i][it-time_restart_init-1]=latt[i].P_(1,1).real();    
	vector <double> gamma_and_domega(2,0.0);
	gamma_and_domega=get_gamma_and_delta_omega(latt[i].Polariz_,g_Xn,omega0);
	gamma_vector[i][it-time_restart_init-1]=gamma_and_domega[0];
	delta_omega_vector[i][it-time_restart_init-1]=gamma_and_domega[1];
        dW_A_vector[i][it-time_restart_init-1]=dW_A[i];
        dW_B_vector[i][it-time_restart_init-1]=dW_B[i];	
      	} 
	end_time_print=omp_get_wtime();

        if(it%interval_energy==0 ){   

	#pragma omp parallel for
	for (int i=0;i<num_traj;i++){
	ostringstream os;
    	string filename;
	os.str("");
        os.clear();
        os<<direct_out<<"temp_occupation_plus_energy"<<name_block.str()<<"_traj_num_"<<i<<".out";
        filename= os.str();

	print_to_file_single_traj_restart(time,i,occ_numberA,half_kin_energyA,pot_energyA,occ_numberB,half_kin_energyB,pot_energyB,f_vector,XA_vector,PA_vector,XB_vector,PB_vector,filename,interval_energy,it,gamma_vector,delta_omega_vector,dW_A_vector,dW_B_vector,time_restart_init);//,XXA_vector,XPA_vector,PPA_vector,XXB_vector,XPB_vector,PPB_vector,XXXA_vector,XXXB_vector);
	}
	}

       	if(it%interval==0 ){   

	//Then I calculate the average of Gloc, Delta, Pol, over the different trajectories at that given time it
	G_avg.clear();
	//Delta_avg.clear();
	Pol_avg.clear();

      	for (int i=0;i<num_traj;i++){
	 G_avg.incr(latt[i].Gloc_,1./num_traj);
	 //Delta_avg.incr(latt[i].Delta_,1./num_traj);
	 Pol_avg.incr(latt[i].Polariz_,1./num_traj);
	}
        os.str("");
        os.clear();
        os <<direct_out<<"temp_G_t" << it <<name_block.str()<<"_num_traj_"<<num_traj<<".out";
        filename= os.str();
        outputGF_new(G_avg,filename);
	cout<<"------  I am printing "<<filename<<"--------------"<<endl;

	
        os.str("");
        os.clear();
	//os <<direct_out<<"temp_Delta_t" << it <<name_block.str()<<"_num_traj_"<<num_traj<<".out";
        filename= os.str();
        outputGF_new(Delta_avg,filename);
	cout<<"------  I am printing "<<filename<<"--------------"<<endl;


        os.str("");
        os.clear();
        os <<direct_out<<"temp_Polariz_t" << it <<name_block.str()<<"_num_traj_"<<num_traj<<".out";
        filename= os.str();
        outputGF_new(Pol_avg,filename);
	cout<<"------  I am printing "<<filename<<"--------------"<<endl;

        }  //end of the printing loop
	

	end_time=omp_get_wtime();
        cout<<"tstp "<<it<<"; execution time: "<<(end_time-start_time)*1000<<endl;
	cout<<"tstp "<<it<<"; time to print in percentage : "<<(end_time_print-start_time_print)/(end_time-start_time)*100<<endl;
	cout<<"tstp "<<it<<"; time to k1 in percentage : "<<(end_time_k1-start_time_k1)/(end_time-start_time)*100<<endl;
	cout<<"tstp "<<it<<"; time to k2 in percentage : "<<(end_time_k2-start_time_k2)/(end_time-start_time)*100<<endl;
	cout<<"tstp "<<it<<"; time to k3 in percentage : "<<(end_time_k3-start_time_k3)/(end_time-start_time)*100<<endl;
	cout<<"tstp "<<it<<"; time to k4 in percentage : "<<(end_time_k4-start_time_k4)/(end_time-start_time)*100<<endl;
	cout<<"tstp "<<it<<"; time to k5 in percentage : "<<(end_time_k5-start_time_k5)/(end_time-start_time)*100<<endl;
	cout<<endl;

	if (   (it-time_restart_init) % ( (time_restart_end-time_restart_init) / mult_print  )==0){
	for (int i=0;i<num_traj;i++){
        os.str("");
        os.clear();
        os <<direct_restart<<"temp_G_t" << it <<name_block.str()<<"_traj_num_"<<i<<".out";
        filename= os.str();
        outputGF(latt[i].Gloc_,filename);
	cout<<"------  I am printing "<<filename<<"--------------"<<endl;


        os.str("");
        os.clear();
        os <<direct_restart<<"temp_XP_t"<< it <<name_block.str()<<"_traj_num_"<<i<<".out";
        filename= os.str();
        outputXP(latt[i].X_,latt[i].P_,filename);
	cout<<"------  I am printing "<<filename<<"--------------"<<endl<<endl;

	}
	}

	}//end iteration over tstp for restart==1 branch

	
	//At the end of the itertions over tstp for the restart==1 condition, I need to print all the GF I may need to restart another time
	for (int i=0;i<num_traj;i++){
        os.str("");
        os.clear();
        os <<direct_restart<<"temp_G_t" << time_restart_end <<name_block.str()<<"_traj_num_"<<i<<".out";
        filename= os.str();
        outputGF(latt[i].Gloc_,filename);
	cout<<"------  I am printing "<<filename<<"--------------"<<endl;


        os.str("");
        os.clear();
        os <<direct_restart<<"temp_XP_t"<<time_restart_end<<name_block.str()<<"_traj_num_"<<i<<".out";
        filename= os.str();
        outputXP(latt[i].X_,latt[i].P_,filename);
	cout<<"------  I am printing "<<filename<<"--------------"<<endl<<endl;

	}
	

	}//end restart branch

  	
       //delate the content of imp_new_veector
      for (vector<sigma_solver*>::iterator pObj = imp_new_vector.begin();
        pObj != imp_new_vector.end(); ++pObj) {
        delete *pObj; // Note that this is deleting what pObj points to,
                    // which is a pointer
          }


     imp_new_vector.clear();
 	
     fftw_cleanup();

		return 0;


}
/*--------Plotting------*/

void outputGF(GF &G, std::string filename){
	std::ofstream outputFile;
	outputFile.open(filename.c_str());
	outputFile << "# frequency" << "\t" << "Spectrum" << "\t"<< "Occupation" << "\t"<< "Retarded Real" << "\t"<< "Retarded imag"<< "\t"<< "Lesser imag"<< "\t"<< "Lesser real"<< std::endl;
	long N=G.ngrid_;

	for (long w = 0; w < N; w++)
	{
		int i = w * G.el_size_;
		outputFile <<std::setprecision(30) << G.grid_[w] << "\t" << G.g_ret[i].real() << "\t"<< G.g_ret[i].imag() << "\t"<< G.g_les[i].imag() << "\t"<< G.g_les[i].real() << "\t"
			   << G.g_ret[i+1].real() << "\t"<< G.g_ret[i+1].imag() << "\t"<< G.g_les[i+1].imag() << "\t"<< G.g_les[i+1].real() << "\t" 
			   << G.g_ret[i+2].real() << "\t"<< G.g_ret[i+2].imag() << "\t"<< G.g_les[i+2].imag() << "\t"<< G.g_les[i+2].real() << "\t"
			   << G.g_ret[i+3].real() << "\t"<< G.g_ret[i+3].imag() << "\t"<< G.g_les[i+3].imag() << "\t"<< G.g_les[i+3].real() << std::endl;
	}
	
	outputFile.close();
}


 void print2vectors(std::vector<double> &V1, std::vector<double> &V2, std::string filename){
	std::ofstream outputFile;
	outputFile.open(filename.c_str());
	for(long i=0; i< V1.size(); i++){
		outputFile << V1[i] << "\t" << V2[i] << std::endl;
	}
	outputFile.close();
}

void print2fftw_cvectors(int n, fftw_complex *V1, fftw_complex *V2, std::string filename){
	std::ofstream outputFile;
	outputFile.open(filename.c_str());
	for(long i=0; i< n; i++){
		outputFile << V1[i][0] << '\t' << V1[i][1] << "\t" << V2[i][0] << '\t' << V2[i][1] << std::endl;
	}
	outputFile.close();
}

 void print2doubles(double &V1, double &V2, std::string filename){
	std::ofstream outputFile;
	outputFile.open(filename.c_str(), std::ofstream::out | std::ofstream::app);
	
	outputFile << V1 << "\t" << V2 << std::endl;
	
	outputFile.close();
}
 void print2doubles(double &V1, double &&V2, std::string filename){
	std::ofstream outputFile;
	outputFile.open(filename.c_str(), std::ofstream::out | std::ofstream::app);
	
	outputFile << V1 << "\t" << V2 << std::endl;
	
	outputFile.close();
}
void printSpectrum(std::vector<double> &freq, std::vector<double> &pot, std::vector<std::vector<double> > &spec, std::string filename){
	std::ofstream outputFile;
	outputFile.open(filename.c_str());
	outputFile << "#" << "var" << "\t";
	for(long p=0; p< pot.size(); p++){
		outputFile << pot[p] << "\t";// << V2[i] << endl;
	}
	outputFile << std::endl;
	
	for(long f=0; f<freq.size(); f++){
		outputFile << freq[f] << "\t";
		for(long p=0; p<pot.size(); p++){
			outputFile << spec[p][f] << "\t";
		}
		outputFile << std::endl;
	}
	outputFile.close();
}			

void Spectrum(GF &G, std::vector<cdmatrix> &Spectrum){
	Spectrum.clear();

	for(long i=0; i<G.ngrid_; i++){
		cdmatrix temp;
		//G.get_Retarded(i, temp);
		temp = G.Retarded[i];
		//Spectrum.push_back(-temp.imag()/M_PI);
		
	}
}
