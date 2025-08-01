#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <vector>
#include <cmath>
#include <sstream>
#include "omp.h"
#include "./find_param.h"
#include "./ness_decls.hpp"

using namespace std;
using namespace ness;

double sin_square(double t,double Amp, double T,double shift);
void print_to_file(vector <int>& time,vector <double>& occ_number,vector <double>& kin_energy, vector <double> & pot_energy,std::string filename,int interval,int it,int & first_time, vector <vector <double>> & matrix_G );
//void print_to_file_equilibrium(vector <int>& time,vector <double>& occ_number,vector <double>& kin_energy, vector <double> & pot_energy,std::string filename,int interval,int it,int first_time,int numOfFreq);
void print_to_file_matrix_G(vector <int>& time,vector <vector <double>> & matrix_G,std::string filename,int interval,int it, int & first_time);
void print_to_file_matrix_der_F(vector <int>& time,vector <vector <double>> & matrix_G,std::string filename,int interval,int it, int & first_time);



double get_derivative_in_omega_zero(double delta_omega, double f_omega_2, double f_omega_1, double f_omega_minus_1, double f_omega_minus_2){
        double deriv=0.;
        deriv=1./(12.*delta_omega)*(-f_omega_2+8.*f_omega_1-8.*f_omega_minus_1+f_omega_minus_2);
        //deriv_Polariz_Retard_in_zero=1./(12.*G.dgrid_)*(-G.Retarded[2](i,i)+8.*G.Retarded[1](i,i)-8.*G.Retarded[G.ngrid_-1](i,i)+G.Retarded[G.ngrid_-2](i,i));
return deriv;
}

void fill_matrix(vector <vector <double>> & matrix, GF & G,int tstp){
//in the matrix I put Im_G_ret_omega0, Im_G_les_omega0, F_G_omega0, derivative_omega0_Im_Gret, derivative_omega0_Im_Gles, derivative_omega0_F_G

matrix[0][tstp]= G.Retarded[0](0,0).imag(); //Im_G_ret at omega=0
matrix[1][tstp]= G.Lesser[0](0,0).imag(); //Im_G_les at omega=0
matrix[2][tstp]=-0.5*G.Lesser[0](0,0).imag()/G.Retarded[0](0,0).imag(); //F_G at omega=0



matrix[3][tstp]=get_derivative_in_omega_zero(G.dgrid_, G.Retarded[2](0,0).imag(), G.Retarded[1](0,0).imag(), G.Retarded[G.ngrid_-1](0,0).imag(), G.Retarded[G.ngrid_-2](0,0).imag()); //Im_G_ret

matrix[4][tstp]=get_derivative_in_omega_zero(G.dgrid_, G.Lesser[2](0,0).imag(), G.Lesser[1](0,0).imag(), G.Lesser[G.ngrid_-1](0,0).imag(), G.Lesser[G.ngrid_-2](0,0).imag()); //Im_G_lesser

matrix[5][tstp]=get_derivative_in_omega_zero(G.dgrid_,-0.5*G.Lesser[2](0,0).imag()/G.Retarded[2](0,0).imag(),-0.5*G.Lesser[1](0,0).imag()/G.Retarded[1](0,0).imag(),-0.5*G.Lesser[G.ngrid_-1](0,0).imag()/G.Retarded[G.ngrid_-1](0,0).imag(),-0.5*G.Lesser[G.ngrid_-2](0,0).imag()/ G.Retarded[G.ngrid_-2](0,0).imag()); //Im_G_lesser

}


void fill_matrix_der_F(vector <vector <double>> & matrix, GF & G,int tstp, const vector <double> & freq_array){

        int ind=0;
        for(int i=0; i<G.ngrid_/2;i++){
                double w=G.grid_[i];
                //if (w==0 or (w>0 and w<=10 and round(w)==w)){
                auto it = std::find(freq_array.begin(), freq_array.end(), w);
                if (it != freq_array.end()) {
                        int ind_pp=i+2;
                        int ind_p=i+1;
                        int ind_m=i-1;
                        int ind_mm=i-2;
                        if (w==0){
                                ind_m=G.ngrid_-1;
                                ind_mm=G.ngrid_-2;
                                }

                        matrix[ind][tstp]=get_derivative_in_omega_zero(G.dgrid_,-0.5*G.Lesser[ind_pp](0,0).imag()/G.Retarded[ind_pp](0,0).imag(),-0.5*G.Lesser[ind_p](0,0).imag()/G.Retarded[ind_p](0,0).imag(),-0.5*G.Lesser[ind_m](0,0).imag()/G.Retarded[ind_m](0,0).imag(),-0.5*G.Lesser[ind_mm](0,0).imag()/ G.Retarded[ind_mm](0,0).imag()); //der_F
                        //cout<<"tstp "<<tstp<<"; w = "<<w<<"; der_F = "<<matrix[ind][tstp]<<endl;
                        ind++;
                }

        }
}





void mixing_G(GF & Gloc,GF & G_old,double mixing){	
 for(int w=0;w<Gloc.ngrid_;w++){
 Gloc.Retarded[w].noalias()=Gloc.Retarded[w]*(1-mixing)+G_old.Retarded[w]*mixing;
 Gloc.Lesser[w].noalias()=Gloc.Lesser[w]*(1-mixing)+G_old.Lesser[w]*mixing;
}
}
void output(string Gname,GF & G,int tstp,int iter=0, string flag="R5"){
char fname[100];
sprintf (fname,"temp_%s_t_%d_iter_%d_%s.out",Gname.c_str(),tstp,iter,flag.c_str());
outputGF_new(G,fname);
}

template<class Matrix>
void DensityMatrix(const GF &G, Matrix &M){
    cplx ii(0.0,1.0);
    M.resize(G.size1_,G.size2_);
    M.setZero();
	for(long i=0; i<G.ngrid_; i++) {
	M+=G.Lesser[i];
	}
	M=M*G.dgrid_/(2.0*M_PI*ii);

}

template<class Matrix>
void KinEnergy(GF G, GF Delta, Matrix &M,int tstp){
    cplx minus_ii(0.0,-1.0);
    M.resize(G.size1_,G.size2_);
    M.setZero();
    GF Delta_G=G;
    Delta_G.clear();
    Delta=G;// in bethe lattice Delta=G; if I do not impose this condition, Delta and G differ a bit due to the fact that G was updated at the end of the last dmft iter and then the kinetic energy is not perfectly real
    Delta_G.left_multiply(Delta,1.0);
    cdmatrix integral_Gret(1,1);
    integral_Gret.setZero();
    for(long w=0; w<G.ngrid_; w++) {
        Delta_G.Lesser[w]= Delta.Retarded[w]*G.Lesser[w] + Delta.Lesser[w] * G.Retarded[w].adjoint();
	Delta_G.Retarded[w] =Delta.Retarded[w] * G.Retarded[w]; //even if not needed
	M+=Delta_G.Lesser[w];
        integral_Gret+=G.Retarded[w];
	}
    M=M*G.dgrid_*minus_ii/(2.0*M_PI);
    integral_Gret=integral_Gret*G.dgrid_;

     cout<<"Integral Gret = "<<integral_Gret.real()<<" "<<integral_Gret.imag()<<endl;
    //output("Delta_G",Delta_G,tstp);
    //output("G",G,tstp);
    //output("D",Delta,tstp);
}

template<class Matrix>
void PotEnergy(GF G, GF Sigma, Matrix &M,int tstp){
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
    //output("Sigma_G",Sigma_G,tstp);
    //output("G",G,tstp);
    //output("S",Sigma,tstp);
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
    // Gout(i1,i2) = Gin(j1,j2)
    for(int i=0;i<Gin.ngrid_;i++){
        Gout.Lesser[i](i1,i2)=Gin.Lesser[i](j1,j2);
        Gout.Retarded[i](i1,i2)=Gin.Retarded[i](j1,j2);
    }
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
			//if(row == 0 && col == 1){print2fftw_cvectors(Nft, output_GL,output_GG, "tG.out");}
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
		//if(i < 100) std::cout << i << GL[i](0,0) * 0.0005 / (2 * M_PI) << '\t' << GL[i](1,1) * 0.0005 / (2 * M_PI) << '\t' << GG[i](1,1) * 0.0005 / (2 * M_PI) <<std::endl;
		
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
	//std::cout << input_GL[0][0] + imag_one * input_GL[0][1] << std::endl;
	//std::cout << input_GL[0] << std::endl;
	fftw_execute(plan_GL);
	for(long i=0; i<N/2; i++){
		//*(locSE.p_ret(N/2 + i, 0, 0)) = output_GL[i][0] + imag_one * output_GL[i][1]; 
		//*(locSE.p_ret(i, 0, 0)) = output_GL[i + 2*N + N/2][0] + imag_one * output_GL[i + 2*N + N/2][1];
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
	//calculateRetarded_SelfE(locSE);
	//ForceImag(locSE);
	if((test=TestPositivity(locSE))>=0){
		std::cout << "locSE < 0 at f=" << test << std::endl;
	}
	
}







void ForceImag(GF &G){
	cdmatrix tmp;
    for(long i=0; i< G.ngrid_; i++){
        tmp=G.Lesser[i];
        G.Lesser[i]=0.5*(tmp-tmp.adjoint());
		//G.g_les[i * el_size] = G.g_les[i * el_size].imag() * cplx(0,1);
		//G.g_les[i * el_size + 3] = G.g_les[i * el_size + 3].imag() * cplx(0,1);
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




/*
double fermi_modified(double w,double beta, double fcutoff){
    double arg=w*beta;
    double T = fcutoff/2;      //this is the period of the function sin^2 that I put on top of the equilibrium Fermi-dirac distribution 
    double A = 0.2;            //this is the amplitude of the function sin^2 that I put on top of the equilibrium Fermi-dirac distribution
    double gamma = M_PI * (1/T);
    double shift = fcutoff/5;    //this is the shift I have to apply to sin^2
    double Amp = 0.6;
    if(arg*arg>10000.0){
        return (w>0 ? 0.0 : 1.0);
    }else{
        return Amp * 1.0/(exp(-arg)+1.0);
    }
}

*/


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


void set_Sigma_diss(double gamma, GF &G, GF &Bath){

    //cdmatrix one(Bath.size1_,Bath.size1_);
    //one.setIdentity();
    cplx one(1.0,0.);

    // I initialize to zero 
    cdmatrix Retarded(Bath.size1_,Bath.size1_);
    Retarded.setZero();
    cdmatrix Lesser(Bath.size1_,Bath.size1_);
    Lesser.setZero();

    //I calculate the value of G_ret and G_les at t=0
	cdmatrix M_ret(Bath.size1_,Bath.size1_);
    M_ret.setZero();
	cdmatrix M_les(Bath.size1_,Bath.size1_);
    M_les.setZero();

    
	for(long i=0; i<G.ngrid_; i++) {
    	M_les += G.Lesser[i];
    	M_ret += G.Retarded[i];
	}

	M_les=M_les*G.dgrid_/(2.0*M_PI);
	M_ret=M_ret*G.dgrid_/(2.0*M_PI);

    
    Retarded = one * gamma * M_ret;
    Lesser = one * gamma * M_les;

    //Set all the values of Retarded(w) and Lesser(w) equal to the value just calculated
	for(long w=0; w<Bath.ngrid_; w++){
        Bath.Retarded[w] = Retarded;
        Bath.Lesser[w] = Lesser;
    }

    //cout<<"Ret = "<<Retarded<<", Les = "<<Lesser<<endl; 

}





void set_phonon_bath(double omega_ph,double beta, GF &Phonon){
cplx ii(0.0,1.0);
//I define the DOS of the phonon bath A(w),bose(w), with which I calculate Phonon Les and Great
vector <double> A(Phonon.ngrid_,0.0);
vector <double> bose(Phonon.ngrid_,0.0);
cdmatrix Lesser(Phonon.size1_,Phonon.size1_), Retarded(Phonon.size1_,Phonon.size1_);
Lesser.setZero();
Retarded.setZero();
double omega = 0.0;
for (long w=0;w<Phonon.ngrid_;w++){
//cout<<"w = "<<w<<endl;
omega=Phonon.grid_[w];
double coeff=0.25;//if I want that \int from 0 to \inf |DOS(w)| dw=1, I should impse coeff=1.0 
//use omega > and <0
	//if (w>0 and w<= (Phonon.ngrid_/2-1)){ // se sono per omega strettamente positive	
	if (omega>0){	
			//cout<<"omega = "<<omega<<endl;
			A[w]= coeff*omega/(omega_ph*omega_ph)*exp(-omega/omega_ph);
			bose[w]=1.0/(exp(omega*beta)-1);
			for(int j=0;j<Phonon.size1_;j++) Lesser(j,j)=-ii*2.*M_PI*A[w]*bose[w];  //the minus sign is because we are dealing with bosons
			for(int j=0;j<Phonon.size1_;j++) Retarded(j,j)=-ii*M_PI*A[w];

			Phonon.Lesser[w]=Lesser;
			Phonon.Retarded[w]=Retarded;
					}

	//if (w> (Phonon.ngrid_/2-1)){ // se sono per omega strettamente negative
			if(omega<0){
			//cout<<"omega = "<<omega<<endl;
			A[w]= coeff*omega/(omega_ph*omega_ph)*exp(-(-omega)/omega_ph);
			bose[w]=1.0/(exp((omega)*beta)-1);
			for(int j=0;j<Phonon.size1_;j++) Lesser(j,j)=-ii*2.*M_PI*A[w]*bose[w];
			for(int j=0;j<Phonon.size1_;j++) Retarded(j,j)=-ii*M_PI*A[w];

			Phonon.Lesser[w]=Lesser;
			Phonon.Retarded[w]=Retarded;
					}

        if(omega==0){ //I make lim w->0 of A(w) and b(w)
                        A[w]= coeff*omega/(omega_ph*omega_ph);
                        //bose[w]=1.0/(omega*beta);
                        for(int j=0;j<Phonon.size1_;j++) Lesser(j,j)=-ii*2.*M_PI*coeff/(omega_ph*omega_ph*beta);  //this is -ii*2.*M_PI*A[w]*b[w] in the limit of A[w],b[w] for w->0
                        for(int j=0;j<Phonon.size1_;j++) Retarded(j,j)=-ii*M_PI*A[w];
                        Phonon.Lesser[w]=Lesser;
                        Phonon.Retarded[w]=Retarded;
                        cout<<"Lesser= "<<Lesser.imag()<<endl;
                        }

}
//outputGF_new(Phonon,"temp_Phonon.out");
}



void set_fermion_bath_bump(GF &Bath,double V){
cplx ii(0.0,1.0);
double shift=2.5, Amp=1.0/M_PI, T=4.0;
double DOS(0.0),f(1.0),omega0(shift), gamma(M_PI/T);

cplx ret;
cdmatrix Retarded(Bath.size1_,Bath.size1_);
Retarded.setZero();

cplx les;
cdmatrix Lesser(Bath.size1_,Bath.size1_);
Lesser.setZero();

	for (long w=0;w<Bath.ngrid_;w++){
	double omega=Bath.grid_[w];
	if (omega<0) {omega0=-shift; f=0.0; }
	//if(omega>=0){
		if( (omega>=omega0-T/2) and (omega<=omega0+T/2) ){
		DOS=Amp*cos(gamma*(omega-omega0))*cos(gamma*(omega-omega0));
		}
		//else if (abs(omega)<=T/2) DOS=Amp*cos(gamma*omega)*cos(gamma*omega);
		else DOS=0.0;

	//}

	ret=-ii*M_PI*DOS;
	les=ii*2.*M_PI*DOS*f;
	
		for (int j=0;j<Bath.size1_;j++){
		Retarded(j,j)=ret;
		Lesser(j,j)=les;
		}

	Bath.Retarded[w]=Retarded;
	Bath.Lesser[w]=Lesser;
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
    int size_; //=2
    double omega_ph_;
    //double lambda_ = 2.5; //0.0; //lambda should be 0.5 or 1 or 2.5
    //double g_el_ph_ =sqrt(lambda_*omega_ph_);   
    vector<GF> Gk_;
    GF Sigma_,Gloc_,Bath_,G1_,G2_,Delta_,Floc_,Phonon_,Sigma_Phonon_,Bath_fermion_bump_,Gweiss_,Sigma_plus_Delta_,Sigma_diss_,Sigma_2U_;
    vector<GF> Fk_; // this is wasteful -- I use uly the Lesser part og Fk
    vector<double> wk_;
    vector<cdmatrix> hk_;
    cdmatrix Hartree_;
    ///////////////////////////////////////////////////////////////
    lattice_func(){
        ntheta_=0;
        dtheta_=0.0;
        grid_spacing_=0;
        numOfFreq_=0;
	omega_ph_=0.;
        size_=1;
        // sizes of vectors are zero
    }
    lattice_func(int ntheta,int numOfFreq,double grid_spacing){
        assert(ntheta>1 && (ntheta%2)==0);
        double norm=0.0;
        ntheta_=ntheta;
        dtheta_=M_PI/(ntheta);
        grid_spacing_=grid_spacing;
        numOfFreq_=numOfFreq;
	omega_ph_=0.2;
        size_=1;
        Gk_.resize(ntheta_+1);
        Fk_.resize(ntheta_+1);
        wk_.resize(ntheta_+1);
        hk_.resize(ntheta_+1);
        for(int i=0;i<=ntheta_;i++){
            cdmatrix tmp(1,1);
            double theta=i*dtheta_;

	    //modified by AP
	    /*
	    double simpson;
	    if (i==0 or i== ntheta) simpson=0.0;
	    else if (i==1 or i==ntheta-1) simpson=55.0/24.0;
	    else if (i==2 or i==ntheta-2) simpson=-1.0/6.0;
	    else if (i==3 or i==ntheta-3) simpson=11.0/8.0;
	    else simpson=1.0;
	    double eps=8.0/M_PI*(theta*0.5 - 0.125* sin(4.0*theta )) -2.0;    //1) case one: max of derivative is 8/M_PI
	    wk_[i]=dtheta_*simpson* (  sin(2*theta)*sin(2*theta) *   sqrt ( 4.0 -  eps*eps) );  // integral(f) = h * Somma_su_tutte_le_i di  (simpson[i]*f[i])
	    */
	    
            double eps=-2.0*cos(theta);
            double simpson=(i==0 || i==ntheta_ ? 1.0 : (1+(i%2))*2.0);
            wk_[i]=sin(theta)*sin(theta)*dtheta_*(simpson/3.0);
	    


            norm+=wk_[i];
            Gk_[i]=GF(grid_spacing_, numOfFreq_,size_);
            Fk_[i]=GF(grid_spacing_, numOfFreq_,size_);
            tmp.setZero();
            tmp(0,0)=eps;
            hk_[i]=tmp;
        }
        // renormalize weights to sum_k wk =1
        for(int i=0;i<=ntheta_;i++) wk_[i]=wk_[i]/norm;
        Sigma_=GF(grid_spacing_, numOfFreq_,size_);
        Sigma_diss_=GF(grid_spacing_, numOfFreq_,size_);
        Sigma_2U_=GF(grid_spacing_, numOfFreq_,size_);
        Sigma_plus_Delta_=GF(grid_spacing_, numOfFreq_,size_);
        Gloc_=GF(grid_spacing_, numOfFreq_,size_);
        Floc_=GF(grid_spacing_, numOfFreq_,size_);
        G1_=GF(grid_spacing_, numOfFreq_,size_);
        G2_=GF(grid_spacing_, numOfFreq_,size_);
        Delta_=GF(grid_spacing_, numOfFreq_,size_);
        Bath_=GF(grid_spacing_, numOfFreq_,size_);
        Bath_fermion_bump_=GF(grid_spacing_, numOfFreq_,size_);
	Phonon_=GF(grid_spacing_, numOfFreq_,size_);
	Sigma_Phonon_=GF(grid_spacing_, numOfFreq_,size_);
        Gweiss_=GF(grid_spacing_, numOfFreq_,size_);
        Hartree_=cdmatrix(1,1);
        Hartree_.setZero(); /// zero for now ...
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
	omega_ph_=ll.omega_ph_;
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
        Sigma_diss_=ll.Sigma_diss_;
        Sigma_2U_=ll.Sigma_2U_;
        Sigma_plus_Delta_=ll.Sigma_plus_Delta_;
        Gloc_=ll.Gloc_;
        Gweiss_=ll.Gweiss_;
        Floc_=ll.Floc_;
        Bath_=ll.Bath_;
        Bath_fermion_bump_=ll.Bath_fermion_bump_;
	Phonon_=ll.Phonon_;
	Sigma_Phonon_=ll.Sigma_Phonon_;
        Delta_=ll.Delta_;
        G1_=ll.G1_;
        G2_=ll.G2_;
        Hartree_=ll.Hartree_;
    }
    lattice_func& operator=(const lattice_func &ll){
        ntheta_=ll.ntheta_;
        dtheta_=ll.dtheta_;
        grid_spacing_=ll.grid_spacing_;
        numOfFreq_=ll.numOfFreq_;
        size_=ll.size_;
	omega_ph_=ll.omega_ph_;
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
        Sigma_diss_=ll.Sigma_diss_;
        Sigma_2U_=ll.Sigma_2U_;
        Sigma_plus_Delta_=ll.Sigma_plus_Delta_;
        Gloc_=ll.Gloc_;
        Gweiss_=ll.Gweiss_;
        Floc_=ll.Floc_;
        Bath_=ll.Bath_;
	Bath_fermion_bump_=ll.Bath_fermion_bump_;
	Phonon_=ll.Phonon_;
	Sigma_Phonon_=ll.Sigma_Phonon_;
        Delta_=ll.Delta_;
        G1_=ll.G1_;
        G2_=ll.G2_;
        Hartree_=ll.Hartree_;
        return *this;
    }
    void get_gloc(void){
        Gloc_.clear();
        for(int i=0;i<=ntheta_;i++) Gloc_.incr(Gk_[i],wk_[i]);
    }

    void get_gloc_nonkdepe(void){
 	cdmatrix one(1,1);
        one(0,0)=1.0;
        Gloc_.clear();
	for(long w = 0; w < Gloc_.ngrid_; w++){
        Gloc_.Retarded[w].noalias() =  (Gloc_.grid_[w]*one -Hartree_- Sigma_.Retarded[w]-Delta_.Retarded[w]).inverse();
        //Gloc_.Lesser[w].noalias()=Floc_.Lesser[w]*Gloc_.Retarded[w].adjoint()-Gloc_.Retarded[w]*Floc_.Lesser[w];
	Gloc_.Lesser[w].noalias() = Gloc_.Retarded[w] *(Sigma_.Lesser[w]+Delta_.Lesser[w])* Gloc_.Retarded[w].adjoint();
	}
	}
    void get_gloc_nonkdepe_ret(void){
 	cdmatrix one(1,1);
        one(0,0)=1.0;
        Gloc_.clear();
	for(long w = 0; w < Gloc_.ngrid_; w++){
        Gloc_.Retarded[w].noalias() =  (Gloc_.grid_[w]*one -Hartree_- Sigma_.Retarded[w]-Delta_.Retarded[w]).inverse();
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
    }
    
    void get_delta(void){
	Delta_=Gloc_;
        //cdmatrix R(1,1);
        //R.setZero();
        //DensityMatrix(Delta_,R);
        //cout<<"nAup from Delta = "<<R(0,0).real()<<endl;
	}
    /*
    void get_delta(void){
        // (1+G1) * Delta = G2, assuming that sum_k h_k=0
        cdmatrix one(1,1);
        cdmatrix tmp(1,1);
        one(0,0)=1.0;
        //for(long w = 0; w < Gloc_.ngrid_; w++){
	
            //tmp=(one+G1_.Retarded[w]).inverse();
            // retarded
            //Delta_.Retarded[w].noalias() = tmp*G2_.Retarded[w];
            // lesser:  DL + G1R DL + G1L DA = G2L
            //Delta_.Lesser[w].noalias() = tmp*(G2_.Lesser[w] - G1_.Lesser[w]*(Delta_.Retarded[w].adjoint()));
        //}
      // now Gloc_=(w - Delta - Sigma)^{-1}, local terms added later.
	
	 //Delta_.Lesser[w].noalias() = Gloc_.Retarded[w].inverse() * Gloc_.Lesser[w] * (Gloc_.Retarded[w].adjoint()).inverse()  - Sigma_.Lesser[w];
	Delta_.Retarded[w]=Gloc_.Retarded[w];
	Delta_.Lesser[w]=Gloc_.Lesser[w];
    	}
    }
    */
	
    void lattice_dyson_ret(void){
        // solve *retarded* lattice dyson equation on each K-point
        // Gk = (idt + mu - Hk - Sigma)^{-1}
        // where Hk and Sigma has been set.
        // does not touch lesser component
        cdmatrix one(1,1);
        one(0,0)=1.0;
	#pragma omp parallel for
        for(int i=0;i<=ntheta_;i++){
		#pragma omp parallel for
            for(long w = 0; w < Gk_[i].ngrid_; w++){
                Gk_[i].Retarded[w].noalias() = (Gk_[i].grid_[w]*one - hk_[i] - Hartree_ - Sigma_.Retarded[w]-Bath_.Retarded[w]-Sigma_Phonon_.Retarded[w]-Bath_fermion_bump_.Retarded[w]).inverse();              //I took out Bath from the formula
               //   Gk_[i].Retarded[w].noalias() = (Gk_[i].grid_[w]*one - hk_[i] - Hartree_ - Sigma_.Retarded[w]).inverse(); 
                //cout << "i " << i << " w : " << w << " " <<Gk_[i].Retarded[w] << endl;
            }
        }
    }
    void lattice_dyson_equi(double beta){
        lattice_dyson_ret();
        for(int i=0;i<=ntheta_;i++) force_equi(Gk_[i],beta);
    }
    void lattice_dyson(void){
        // solve lattice dyson equation on each K-point
        // Gk = (idt + mu - Hk - Sigma)^{-1}
        // where Hk and Sigma has been set.
        cdmatrix one(1,1);
        one(0,0)=1.0;
	#pragma omp parallel for	
        for(int i=0;i<=ntheta_;i++){
		#pragma omp parallel for
            for(long w = 0; w < Gk_[i].ngrid_; w++){
                Gk_[i].Retarded[w].noalias() = (Gk_[i].grid_[w]*one - hk_[i] - Hartree_ - Sigma_.Retarded[w]-Bath_.Retarded[w]-Sigma_Phonon_.Retarded[w]-Bath_fermion_bump_.Retarded[w]).inverse();
                Gk_[i].Lesser[w].noalias() = Gk_[i].Retarded[w] * (Sigma_.Lesser[w] + Bath_.Lesser[w]+Sigma_Phonon_.Lesser[w]+Bath_fermion_bump_.Lesser[w])* Gk_[i].Retarded[w].adjoint();
            }
        }
    }


  
    void getG_ret(int h){  
        cdmatrix one(1,1);
        one(0,0)=1.0;
	cdmatrix two(1,1);
        one(0,0)=2.0;
	cdmatrix h_matrix(1,1);
        one(0,0)=h;
	cplx imag_one(0,1);
        for(int i=0;i<=ntheta_;i++){
            for(long w = 0; w < Gk_[i].ngrid_; w++){
               // Gk_[i].Retarded[w].noalias() = (Gk_[i].grid_[w]*one - hk_[i] - Hartree_ - Sigma_.Retarded[w]-Bath_.Retarded[w]).inverse();              //I took out Bath from the formula
	      Gk_[i].Retarded[w].noalias() = Gk_[i].Retarded[w] * (one - (two * h_matrix * imag_one) * (-Gk_[i].grid_[w]*one + hk_[i] - Hartree_ - Sigma_.Retarded[w])) - (two * h_matrix * imag_one) ; 
                //cout << "i " << i << " w : " << w << " " <<Gk_[i].Retarded[w] << endl;

	       Gk_[i].Lesser[w]=Fk_[i].Lesser[w]*(Gk_[i].Retarded[w].adjoint())-Gk_[i].Retarded[w]*Fk_[i].Lesser[w];
            }
        }
  }

//---------------case 1): the starting distribution function is the equilibrium one-------------------------------------------------------------------------------//
  
    void init_distribution_equi(int i,double beta){
        cdmatrix one(1,1);
        one(0,0)=1.0;
        for(long w = 0; w < Gk_[i].ngrid_; w++){
           double fw=fermi(Gk_[i].grid_[w],beta);
           Fk_[i].Lesser[w]=one*fw;
        }
    }
    void init_distribution_equi(double beta){
        for(int i=0;i<=ntheta_;i++) init_distribution_equi(i,beta);
        }

    void init_distribution_equi_nonkdepe(double beta){
        cdmatrix one(1,1);
	one.setZero();
        one(0,0)=1.0;
	Floc_.clear();
        for(long w = 0; w < Gloc_.ngrid_; w++){
           double fw=fermi(Gloc_.grid_[w],beta);
           Floc_.Lesser[w]=one*fw;
        }
    }

//--------------- end case 1)-------------------------------------------------------------------------------//


//---------------case 2): the starting distribution function is a nonequilibrium one and is put by hand by myself-------------------------------------------------------------------------------//
  
  void init_distribution_nonequi(int i,double beta,double fcutoff){
        cdmatrix one(1,1);
        one(0,0)=1.0;
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
        cdmatrix one(1,1);
        one(0,0)=1.0;	
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
cdmatrix one(1,1);
one(0,0)=1.0;	
double fw;
ostringstream os;
string filename;
os.str("");
os.clear();
os<<"Init_F"<<tstp<<".txt";
filename= os.str();
const char *filename_c = filename.c_str();
fstream myfile(filename_c, std::ios_base::in);
// FILE *fout1=fopen(filename_c,"w");		
//fstream myfile("Init_F.txt", std::ios_base::in);       //this is in case I have as input file Init_F.txt
//for(int i=0;i<=ntheta_;i++){ 
//myfile>>fw;
for(long w = 0; w < Floc_.ngrid_; w++){
myfile>>fw;
Floc_.Lesser[w]=one*fw;
}
//fclose(fout1);
//}
}


//-------------------------------------end of case 3) -------------------------------------------------------------------------------------//

  

  /*  
//---------------case 3) in which the input file is Init_F.txt-------------------------------------------------------------------------------//

//I add the next two functions 


    void init_distribution_nonequi(double beta,double fcutoff){
        for(int i=0;i<=ntheta_;i++) init_distribution_nonequi(i,beta,fcutoff);
        }

   void init_distribution_nonequi(int i,double beta, double fcutoff){
        cdmatrix one(1,1);
        one(0,0)=1.0;
	double fw;
	fstream myfile("Init_F.txt", std::ios_base::in);       //this is in case I have as input file Init_F.txt
	
        for(long w = 0; w < Gk_[i].ngrid_; w++){
	  //double fw=fermi_modified(Gk_[i].grid_[w],beta,fcutoff);
	   myfile>>fw;
           Fk_[i].Lesser[w]=one*fw;
        }

    }

//--------------------------------------------------end of case 3)---------------------------------------------------------------------------//

//----------------------------------------------------------------------------------------------//  
*/
  
    void set_distribution_function(void){
            for(long w = 0; w < Floc_.ngrid_; w++){
                Gloc_.Lesser[w].noalias()=Floc_.Lesser[w]*(Gloc_.Retarded[w].adjoint()-Gloc_.Retarded[w]);
            }
    }
    void get_hartree(double U){
        Hartree_(0,0)=0.0;
	cdmatrix R(2,2);
        DensityMatrix(Gloc_,R); //here I resize R to 1x1 matrix
	//Hartree_(0,0)=U*(R(0,0).real()-0.5);
	Hartree_(0,0)=0.0;
	//What I found when I got the code
        //Hartree_(0,0)=(R(1,1).real()-0.5)*U;
        //Hartree_(1,1)=(R(0,0).real()-0.5)*U;

    }
    void set_hk_diag(double ea,double eb){
        // donothing
        /*for(int i=0;i<=ntheta_;i++){
            hk_[i](0,0)=ea;
            hk_[i](1,1)=eb;
        }*/
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


void updateSE1(double U, const GF &GA, GF &locSEA, fft_solver &solver, GF & Bath, GF & Sigma_Phonon, GF & Bath_fermion_bump, GF & Sigma_diss, GF & Sigma2U)
{
	GF tGA(GA.dgrid_, GA.ngrid_ * 3 / 2 + 1, GA.size1_, time_gf), tlocSEA(GA.dgrid_, GA.ngrid_ * 3 / 2 + 1, GA.size1_, time_gf);
	solver.to_time(tGA, GA, 1);
    
    for(long t = 0; t < tGA.ngrid_; t++)
	{
		tlocSEA.Lesser[t](0,0) = tGA.Lesser[t](0,0) * (-tGA.Greater(t).adjoint()(0,0)) * tGA.Lesser[t](0,0);
		tlocSEA.Retarded[t](0,0) = tGA.Greater(t)(0,0) * (-tGA.Lesser[t].adjoint()(0,0)) * tGA.Greater(t)(0,0)-tlocSEA.Lesser[t](0,0);
	}
	tlocSEA.smul(U * U);
	tlocSEA.reset_grid(tGA.dgrid_);
	tlocSEA.Retarded[0] *= 0.5;
	tlocSEA.Retarded[tlocSEA.ngrid_ - 1] *= 0.0;
	solver.to_freq(locSEA, tlocSEA);
	ForceImag(locSEA);
    Sigma2U = locSEA;
    locSEA.incr(Bath);       
    locSEA.incr(Sigma_Phonon);
	locSEA.incr(Bath_fermion_bump);       
	locSEA.incr(Sigma_diss);       
	ForceImag(locSEA);
}

void get_Sigma_ph(double g_ph_bath,const GF & GA,const GF & Phonon,GF & Sigma_ph,fft_solver &solver){
cplx ii(0.0,1.0);
//GF tGA(GA.dgrid_, GA.ngrid_ * 3 / 2 + 1, GA.size1_, time_gf);
//GF tPhonon(Phonon.dgrid_, Phonon.ngrid_ * 3 / 2 + 1, Phonon.size1_, time_gf);
//GF tSigma_ph(Phonon.dgrid_, Phonon.ngrid_ * 3 / 2 + 1, Phonon.size1_, time_gf);
GF tGA(GA.dgrid_, GA.ngrid_ * 3 / 2 + 1, GA.size1_, time_gf);
GF tPhonon(GA.dgrid_, GA.ngrid_ * 3 / 2 + 1, GA.size1_, time_gf);
GF tSigma_ph(GA.dgrid_,GA.ngrid_ * 3 / 2 + 1, GA.size1_, time_gf);
solver.to_time(tGA, GA, 1);  //I get time G lesser and retarded
solver.to_time(tPhonon,Phonon,1);
for (long t = 0 ; t< tGA.ngrid_ ; t++ ){
tSigma_ph.Lesser[t](0,0) = tGA.Lesser[t](0,0) * tPhonon.Lesser[t](0,0);
tSigma_ph.Retarded[t](0,0)= tGA.Greater(t)(0,0) * tPhonon.Greater(t)(0,0) - tSigma_ph.Lesser[t](0,0); 
}
tSigma_ph.smul(g_ph_bath*g_ph_bath);
tSigma_ph.reset_grid(tGA.dgrid_); 
tSigma_ph.smul(ii); //why this factor?
tSigma_ph.Retarded[0]*=0.5;
tSigma_ph.Retarded[tSigma_ph.ngrid_-1]*=0.0;
solver.to_freq(Sigma_ph,tSigma_ph);
ForceImag(Sigma_ph); 
}


void get_gweiss_from_Delta(GF & Gloc, GF & Gweiss, GF & Sigma,GF & Delta, cdmatrix & Hartree){
    // solve Gweiss = (w -Delta)^{-1}. Local terms added TODO
    cdmatrix one(1,1);
    one(0,0)=1.0;
    for(long w = 0; w < Gweiss.ngrid_; w++){
        Gweiss.Retarded[w].noalias() = (Gweiss.grid_[w]*one - Hartree- Delta.Retarded[w]).inverse();
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
    void sigma2_bold(lattice_func &latt,double U,double gamma){
        latt.get_gloc();
        set_Sigma_diss(gamma,latt.Gloc_,latt.Sigma_diss_);
        updateSE1(U,latt.Gloc_,latt.Sigma_,solver_,latt.Bath_,latt.Sigma_Phonon_,latt.Bath_fermion_bump_,latt.Sigma_diss_,latt.Sigma_2U_);
    }
    void sigma2_ipt(lattice_func &latt,double U, double g_ph_bath,double gamma){
        
        //N.B. this is sigma2_ipt when I want to calculate Sigma with Gweiss 
	    //This was for the k-dependent code
        //GF Gweiss=latt.Gloc_;
        //latt.get_gloc_moments();
        //latt.get_delta();
        //get_gweiss_from_Delta(Gweiss,latt.Delta_);
        //updateSE1(U,Gweiss,latt.Sigma_,solver_);

        GF Gweiss=latt.Gloc_;
        latt.get_delta();
        get_gweiss_from_Delta(latt.Gloc_,Gweiss,latt.Sigma_,latt.Delta_,latt.Hartree_);
	    get_Sigma_ph(g_ph_bath,latt.Gloc_,latt.Phonon_,latt.Sigma_Phonon_,solver_);	
        set_Sigma_diss(gamma,latt.Gloc_,latt.Sigma_diss_);
        updateSE1(U,Gweiss,latt.Sigma_,solver_,latt.Bath_,latt.Sigma_Phonon_,latt.Bath_fermion_bump_,latt.Sigma_diss_,latt.Sigma_2U_);
/*  
        //N.B. this is sigma2_ipt when I want to calculate Sigma using the full Green function 

        updateSE1(U,latt.Gloc_,latt.Sigma_,solver_);

*/
    }
};

double dmft_iteration(int iter,lattice_func &latt,sigma_solver &imp,int ipt,double U,int ret,double mixing,int tstp,string flag, double g_ph_bath, double gamma ){
    double err; 
    GF G_old=latt.Gloc_;
   
    if(ipt){
        imp.sigma2_ipt(latt,U,g_ph_bath,gamma);
    }else{
        imp.sigma2_bold(latt,U,gamma);
    }
    if(ret==1) {
	latt.get_hartree(U);	
	latt.get_gloc_nonkdepe_ret();
        latt.set_distribution_function();
    }else{
	latt.get_hartree(U);
       	latt.get_gloc_nonkdepe();


    }
    mixing_G(latt.Gloc_,G_old,mixing);
    
    //output("GA",latt.Gloc_,tstp,iter,flag);
    //output("Sigma",latt.Sigma_,tstp,iter,flag);
    //output("Delta",latt.Delta_,tstp,iter,flag);
    //output("Gold",G_old,tstp,iter,flag);
    //output("F",latt.Floc_,tstp,iter,flag);
    
    err = GF2norm(latt.Gloc_, G_old);
    return err;
}

void get_ISigma_F(lattice_func &latt_new,lattice_func &latt){
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

void dFdt(lattice_func &latt_new,lattice_func &latt,sigma_solver &imp,double U,double errmax,int Niter,int ipt,std::string flag,double mixing, int ret,int tstp, double g_ph_bath, double gamma){
    double err=10.0;
    int i=0;
   
   cout<<"U = "<<U<<endl;
    while(i<Niter && err>errmax){
      err=dmft_iteration(i,latt,imp,ipt,U,ret,mixing,tstp,flag,g_ph_bath,gamma);
        cout << "time " << flag << " -- iter " << i << " err " << err << endl;
        i++;
    }
    get_ISigma_F(latt_new,latt);
   // latt_new.output_all("test");                    // I do not want to print anything at the moment
}
void set_les(lattice_func &lout,lattice_func &lin){
    //for(int i=0;i<=lin.ntheta_;i++){
        for(int w=0;w<lin.Bath_.ngrid_;w++){
            lout.Gloc_.Lesser[w]=lin.Gloc_.Lesser[w];
        }
    //}
}
void incr_les(lattice_func &lout,lattice_func &lin,double a){
    //for(int i=0;i<=lin.ntheta_;i++){
        for(int w=0;w<lin.Bath_.ngrid_;w++){
            lout.Gloc_.Lesser[w]+=a*lin.Gloc_.Lesser[w];
        }
    //}
}
void smul_les(lattice_func &lout,cplx a){
    //for(int i=0;i<=lout.ntheta_;i++){
        for(int w=0;w<lout.Bath_.ngrid_;w++){
            lout.Gloc_.Lesser[w]*=a;
        }
    //}
}
void set_les_F(lattice_func &lout,lattice_func &lin){
    //for(int i=0;i<=lin.ntheta_;i++){
        for(int w=0;w<lin.Bath_.ngrid_;w++){
            lout.Floc_.Lesser[w]=lin.Floc_.Lesser[w];
        }
    //}
}
void incr_les_F(lattice_func &lout,lattice_func &lin,double a){
    //for(int i=0;i<=lin.ntheta_;i++){
        for(int w=0;w<lin.Bath_.ngrid_;w++){
            lout.Floc_.Lesser[w]+=a*lin.Floc_.Lesser[w];
        }
    //}
}
void smul_les_F(lattice_func &lout,cplx a){
    //for(int i=0;i<=lout.ntheta_;i++){
        for(int w=0;w<lout.Floc_.ngrid_;w++){
            lout.Floc_.Lesser[w]*=a;
        }
    //}
}
double sin_ramp(double t,double tramp,double x_i,double x_f){
  
   // cout<<" t0 = "<<t<< " tramp = "<<tramp<<" Ui = "<<x_i<<"U_f"<<x_f<<endl;   
 
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





int main(int argc, char** argv){
    int latt_output_all;
	double grid_spacing=0.0005; //these two values will be modified after the 
	long numOfFreq=40000;  //reading of param_boltz.in
  int it =0, ret=1,interval=1,first_time=1,first_time_matrix_G=1, first_time_matrix_Sigma=1, first_time_matrix_Delta=1,first_time_matrix_Sigma_plus_Delta=1,first_time_matrix_der_F=1 ; 
  int tmp_int;
  double err_F;
  double freq_cutoff=10.0;
  double g_ph_bath=0.0;
  double beta,beta_ph;
  double chemPot=0.000;
  double mixing;
  long numOfIt=100;
  double U=-2.5,U_i,U_f,tramp,gamma_init=0.0, gamma_final=0.0, gamma_final_dmft=0.0;
  int ntheta,time_to_boltzmann;
  int nt; 
  double Jbath,h;
  double lerr = 0.0;
  int Niter = 10,NiterR;
  int i=0;
  int ipt;
   int outputf,param;
	double Amp_Vt,T_Vt;int shift_Vt=0;
	char fname[100];
	cdmatrix R(1,1);
	cdmatrix R_Ekin(1,1);
	cdmatrix R_Epot(1,1);


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
	ntheta=2;
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
	cout<<"Grid_spacing = "<<grid_spacing<<endl;
        find_param(argv[1], "__h", h);
        find_param(argv[1], "__nt", nt);
	find_param(argv[1],"__time_to_boltzmann=",time_to_boltzmann);
        find_param(argv[1], "__outputf", outputf);
        find_param(argv[1], "__ipt", ipt);
	find_param(argv[1],"__Amp_Vt",Amp_Vt);
	find_param(argv[1],"__T_Vt",T_Vt); 
	find_param(argv[1],"__shift_Vt",shift_Vt);
	find_param(argv[1],"__g_ph_bath",g_ph_bath);
	find_param(argv[1],"__interval_plot",interval); interval=nt;
	find_param(argv[1],"__gamma_init",gamma_init);
	find_param(argv[1],"__gamma_final",gamma_final);
	find_param(argv[1],"__param",param);
 	}
	
	cout<<"Amp_Vt = "<<Amp_Vt<<endl;
	cout<<"T_Vt = "<<T_Vt<<endl;
	cout<<"shift_Vt = "<<shift_Vt<<endl;
	cout<<"g_ph_bath = "<<g_ph_bath<<endl;
	cout<<"interval = "<<interval<<endl;
	cout<<"gamma_init = "<<gamma_init<<endl;
	cout<<"gamma_final = "<<gamma_final<<endl;


    string direct_out="/nfs/tfkp00/picanoan/Desktop/C++/codes/boltzmann/standalone_ness_mio/diss/param_"+std::to_string(param)+"/";
    cout<<"The output directory of the Boltz code is "<<direct_out<<endl;
    vector <int> time(nt+1,0);
    vector <double> occ_number(nt+1,0.0), kin_energy(nt+1,0.0), pot_energy(nt+1,0.0) ;

    int row_mat=6;
    vector<vector<double>> matrix_G(row_mat , vector<double> (nt+1, 0.));
    vector<vector<double>> matrix_Sigma(row_mat , vector<double> (nt+1, 0.));
    vector<vector<double>> matrix_Delta(row_mat , vector<double> (nt+1, 0.));
    vector<vector<double>> matrix_Sigma_plus_Delta(row_mat , vector<double> (nt+1, 0.));
    //This is when I want to plot multiple frequencies
    std::vector<double> freq_to_plot;
    //for (double i = 0.0; i <= 20.0; i += 0.5) freq_to_plot.push_back(i);
    //This is when I want to plot only w=0
    for (double i = 0.0; i <= 20.0; i += 100) freq_to_plot.push_back(i);
    vector<vector<double>> matrix_der_F(freq_to_plot.size() , vector<double> (nt+1, 0.));

    ostringstream block_name;
    block_name.str("");
    block_name.clear();
    block_name<<"_U_"<<U_i<<"_beta_"<<beta<<"_h_"<<h<<"_T_Vt_"<<T_Vt<<"_Amp_Vt_"<<Amp_Vt<<"_tstp_max_"<<nt<<"_gamma_init_"<<gamma_init<<"_gamma_final_"<<gamma_final;
    cout<<"block_name = "<<block_name.str()<<endl;


    ostringstream block_name_reduced;
    block_name_reduced.str("");
    block_name_reduced.clear();
    block_name_reduced<<"_U_"<<U_i<<"_beta_"<<beta<<"_gamma_init_"<<gamma_init;
    cout<<"block_name = "<<block_name_reduced.str()<<endl;



    //T_Vt*=h;
    //shift_Vt*=h;
    double V=0.0;
 
    lattice_func latt(2,numOfFreq,grid_spacing);
    // for RK:
    lattice_func ltmp(2,numOfFreq,grid_spacing);
    lattice_func lk1(2,numOfFreq,grid_spacing);
    lattice_func lk2(2,numOfFreq,grid_spacing);
    lattice_func lk3(2,numOfFreq,grid_spacing);
    lattice_func lk4(2,numOfFreq,grid_spacing);

    setBath_diag(Jbath,chemPot,beta_ph,latt.Bath_);
    set_phonon_bath(latt.omega_ph_,beta,latt.Phonon_);
    set_fermion_bath_bump(latt.Bath_fermion_bump_,V);
    sigma_solver imp(numOfFreq);
    i=0;
    double GF_err = 10.0;
    // start with some local greens function, sigma=0
     latt.init_distribution_equi_nonkdepe(beta);
    //latt.init_distribution_nonequi(beta,freq_cutoff);
    // latt.init_distribution_nonequi_eps(beta,freq_cutoff,time_to_boltzmann);

    complex <double> mult=cplx(0.0,-1.0);
    latt.Gloc_.clear();
    latt.Sigma_.clear();
    latt.Sigma_.incr(latt.Bath_);       
    latt.Sigma_.incr(latt.Bath_fermion_bump_);
    latt.Delta_.clear();
    latt.get_hartree(U);

	 if (ret==1) {latt.get_gloc_nonkdepe_ret(); 
                      latt.set_distribution_function();
      			}
        else latt.get_gloc_nonkdepe();

        DensityMatrix(latt.Gloc_,R);
        cout << " first iter "<<",  err= " << GF_err << " dens_tot= " << R(0,0) <<endl;


        setBath_diag(0.0,chemPot,beta_ph,latt.Bath_); //set the bath to zero after the very first iteration

	while(GF_err > lerr && i < Niter){
	cout << "------------ Iteration: " << i << " -------------" << endl;     
	GF_err=dmft_iteration(i,latt,imp,ipt,U,ret,mixing,it,"R1_",g_ph_bath,gamma_init);       
        DensityMatrix(latt.Gloc_,R);
        cout << "iter= " << i << " err= " << GF_err << " dens_tot= " << R(0,0) <<endl;
	cout<< "ret = "<<ret<<endl;
        i++;
        cout <<"--------------------------------------------------"<< endl;
    	}

	mixing=0.;

      cout << "it = " << it << " t0 = "  << it*h << " U = " << U <<" nt = "<<nt<<"  h = "<<h<<" V fermionic bump bath "<<V<<endl<<endl;
	 //outputGF_new(latt.Sigma_, "selfA.out");
   	 //outputGF_new(latt.Gloc_, "temp_GA.out");
    	 //outputGF_new(imp.Gweiss_A_, "GweissA.out");
         //outputGF_new(latt.Sigma_Phonon_, "Sigma_Phonon.out");
	 //outputGF_new(latt.Phonon_, "Phonon_Bath.out");

	//Calculate Ekin
	KinEnergy(latt.Gloc_,latt.Delta_,R_Ekin,0);
	cout << "EkinA = " << R_Ekin(0,0).real() << " " << R_Ekin(0,0).imag()<<endl;
	PotEnergy(latt.Gloc_,latt.Sigma_2U_,R_Epot,0);
	cout << "Epot = " << R_Epot(0,0).real() << " " << R_Epot(0,0).imag()<<endl;
	cout<<endl;
	
        os.str("");
        os.clear();
        os <<direct_out<<"temp_G_t_0"<<block_name.str()<<".out";
        filename= os.str();
        outputGF_new(latt.Gloc_,filename);

        os.str("");
        os.clear();
        os <<direct_out<<"temp_G_t_0"<<block_name_reduced.str()<<".out";
        filename= os.str();
        outputGF_new(latt.Gloc_,filename);





        os.str("");
        os.clear();
        os <<direct_out<<"temp_Delta_t_0"<<block_name.str()<<".out";
        filename= os.str();
        //outputGF_new(latt.Delta_,filename);



        os.str("");
        os.clear();
        os <<direct_out<<"temp_Sigma_t_0"<<block_name.str()<<".out";
        filename= os.str();
        //outputGF_new(latt.Sigma_,filename);



        os.str("");
        os.clear();
        os <<direct_out<<"temp_F_t_0"<<block_name.str()<<".out";
        filename= os.str();
        //outputGF_new(ltmp.Floc_,filename);




        os.str("");
        os.clear();
        os <<direct_out<<"temp_Sigma_Phonon_t_0"<<block_name.str()<<".out";
        filename= os.str();
        //outputGF_new(latt.Sigma_Phonon_,filename);


        os.str("");
        os.clear();
        os <<direct_out<<"temp_D0_Phonon_t_0"<<block_name.str()<<".out";
        filename= os.str();
        //outputGF_new(latt.Phonon_,filename);


	latt.Sigma_plus_Delta_.clear();
	latt.Sigma_plus_Delta_.incr(latt.Sigma_);
	latt.Sigma_plus_Delta_.incr(latt.Delta_);
        os.str("");
        os.clear();
        os <<direct_out<<"temp_Sigma_plus_Delta_t_0"<<block_name.str()<<".out";
        filename= os.str();
        //outputGF_new(latt.Sigma_plus_Delta_,filename);

        os.str("");
        os.clear();
        os <<direct_out<<"temp_Bath_fermion_bump_t_0"<<block_name.str()<<".out";
        filename= os.str();
        //outputGF_new(latt.Bath_fermion_bump_,filename);




	time[it]=it;
	occ_number[it]=R(0,0).real();	 
	kin_energy[it]=R_Ekin(0,0).real(); 
	pot_energy[it]=R_Epot(0,0).real();//+U*(R(0,0).real()-0.5); 
        fill_matrix(matrix_G,latt.Gloc_,it);
        //fill_matrix(matrix_Delta,latt.Delta_,it);
        //fill_matrix(matrix_Sigma,latt.Sigma_,it);
        latt.Sigma_plus_Delta_.clear();
        latt.Sigma_plus_Delta_=latt.Sigma_;
        latt.Sigma_plus_Delta_.incr(latt.Delta_);
        //fill_matrix(matrix_Sigma_plus_Delta,latt.Sigma_plus_Delta_,it);
        fill_matrix_der_F(matrix_der_F,latt.Gloc_,it,freq_to_plot);




	/*
        os.str("");
        os.clear();
        os<<direct_out<<"temp_occupation_plus_energy_equilibrium_"<<block_name.str()<<".out";
        filename= os.str();
        print_to_file(time,occ_number,kin_energy,pot_energy,filename,interval,it,first_time);//print in append mode
 	*/

     //setBath_diag(0.0,chemPot,beta_ph,latt.Bath_); //I do not need anymore the fermionic bath

     for(it=1;it<=nt;it++){

        if (it >= shift_Vt/h) gamma_final_dmft = gamma_final;
        else gamma_final_dmft = gamma_init;
        cout<<"gamma final dmft "<<gamma_final_dmft<<endl;


        cdmatrix R(1,1);
        double t0=(it-1)*h;
        // a Runge-Kutta step for
        // ii * d G</dt = Delta_R G< - G< Delta_A + Delta< GA -GR Delta<
        // where GR and Delta are time-dependent
        //set_les(ltmp,latt);
        
       //I calculate k1: 
         ltmp=latt;
       // U=sin_ramp(t0,tramp,U_i,U_f);
       // dFdt(lk1,ltmp,imp,U,lerr,NiterR,ipt,"R1_",mixing,ret);
         get_ISigma_F(lk1,ltmp);
	 smul_les_F(lk1,mult);


        //I calculate k2:
        ltmp=latt;
        incr_les_F(ltmp,lk1,h*0.5);
        V=sin_square(t0+0.5*h,Amp_Vt,T_Vt,shift_Vt);
	    //set_fermion_bath_bump(latt.Bath_fermion_bump_,V);
	    set_fermion_bath_bump(ltmp.Bath_fermion_bump_,V);
        U=sin_ramp(t0+0.5*h,tramp,U_i,U_f);
        dFdt(lk2,ltmp,imp,U,lerr,NiterR,ipt,"R2_",mixing,ret,it,g_ph_bath,gamma_final_dmft);
        smul_les_F(lk2,mult);

        //I calculate k3:       
        ltmp=latt;
        incr_les_F(ltmp,lk2,h*0.5);
        V=sin_square(t0+0.5*h,Amp_Vt,T_Vt,shift_Vt);
        //set_fermion_bath_bump(latt.Bath_fermion_bump_,V);
	set_fermion_bath_bump(ltmp.Bath_fermion_bump_,V);
        U=sin_ramp(t0+0.5*h,tramp,U_i,U_f);
        dFdt(lk3,ltmp,imp,U,lerr,NiterR,ipt,"R3_",mixing,ret,it,g_ph_bath,gamma_final_dmft);
        smul_les_F(lk3,mult);

	//I calculate k4:
        ltmp=latt;
        incr_les_F(ltmp,lk3,h);
        V=sin_square(t0+h,Amp_Vt,T_Vt,shift_Vt);
       	//set_fermion_bath_bump(latt.Bath_fermion_bump_,V);
       	set_fermion_bath_bump(ltmp.Bath_fermion_bump_,V);
        U=sin_ramp(t0+h,tramp,U_i,U_f);
        dFdt(lk4,ltmp,imp,U,lerr,NiterR,ipt,"R4_",mixing,ret,it,g_ph_bath,gamma_final_dmft);
        smul_les_F(lk4,mult);

        set_les_F(ltmp,lk1);
        incr_les_F(ltmp,lk2,2.0);
        incr_les_F(ltmp,lk3,2.0);
        incr_les_F(ltmp,lk4,1.0);
        incr_les_F(latt,ltmp,h/6.0);

/*    this part of the code is used just to calculate the error
        ltmp_F = latt;
        ltmp_F.get_floc();
        incr_les_F(latt,ltmp,h/6.0);
        latt.get_floc();
        err_F = GF2norm(latt.Floc_,ltmp_F.Floc_);
        cout<<"error on F = "<<err_F<<endl<<endl; 
       
*/        
        // dmft iteration at the new F:
        i=0;
        GF_err=10.0;

        V=sin_square(t0+h,Amp_Vt,T_Vt,shift_Vt);
       	set_fermion_bath_bump(latt.Bath_fermion_bump_,V);
        U=sin_ramp(t0+h,tramp,U_i,U_f);

               while(i<NiterR and GF_err>lerr){
	        GF_err=dmft_iteration(i,latt,imp,ipt,U,ret,mixing,it,"R5_",g_ph_bath,gamma_final_dmft);
          cout << "time -- iter " << i << " err " << GF_err << endl;
            i++;
        }
        ///

        DensityMatrix(latt.Gloc_,R);
        cout << "it = " << it << " t0 = "  << t0 << " U = " << U<<"  nt = "<<nt<<"  h = "<<h<<" V fermionic bump bath = "<<V<<endl;
        cout << "nA = " << R(0,0).real() << " " << R(0,0).imag()<<endl;
	//Calculate Ekin
	KinEnergy(latt.Gloc_,latt.Delta_,R_Ekin,it);
	cout << "EkinA = " << R_Ekin(0,0).real() << " " << R_Ekin(0,0).imag()<<endl;
	//Calculate Epot
	PotEnergy(latt.Gloc_,latt.Sigma_2U_,R_Epot,it);
	cout << "EpotA = " << R_Epot(0,0).real() << " " << R_Epot(0,0).imag();
	cout<<endl<<endl;

	time[it]=it;
	occ_number[it]=R(0,0).real();	 
	kin_energy[it]=R_Ekin(0,0).real(); 
	pot_energy[it]=R_Epot(0,0).real();//+U*(R(0,0).real()-0.5); 
        fill_matrix(matrix_G,latt.Gloc_,it);
        //fill_matrix(matrix_Delta,latt.Delta_,it);
        //fill_matrix(matrix_Sigma,latt.Sigma_,it);
        latt.Sigma_plus_Delta_.clear();
        latt.Sigma_plus_Delta_=latt.Sigma_;
        latt.Sigma_plus_Delta_.incr(latt.Delta_);
        //fill_matrix(matrix_Sigma_plus_Delta,latt.Sigma_plus_Delta_,it);
        fill_matrix_der_F(matrix_der_F,latt.Gloc_,it,freq_to_plot);



        /*
        if(V>0){   

        os.str("");
        os.clear();
        os <<direct_out<<"temp_Bath_fermion_bump_t_" << it<<block_name.str()<<".out";
        filename= os.str();
        outputGF_new(latt.Bath_fermion_bump_,filename);

	}
    */
		    
        //these are my output functions, that at the moment I do not want to print
        if(it%interval==0){   

	    //I want to print the scattering integral I=ltmp.Floc_ (real lesser) /6
        os.str("");
        os.clear();
        os <<direct_out<<"temp_F_t_" << it<<block_name.str()<<".out";
        filename= os.str();
        //outputGF_new(ltmp.Floc_,filename);


        os.str("");
        os.clear();
        os <<direct_out<<"temp_G_t_" << it<<block_name.str()<<".out";
        filename= os.str();
        outputGF_new(latt.Gloc_,filename);

        os.str("");
        os.clear();
        os <<direct_out<<"temp_Sigma_t_"<<it<<block_name.str()<<".out";
        filename= os.str();
        //outputGF_new(latt.Sigma_,filename);

        os.str("");
        os.clear();
        os <<direct_out<<"temp_Delta_t_"<<it<<block_name.str()<<".out";
        filename= os.str();
        //outputGF_new(latt.Delta_,filename);


        os.str("");
        os.clear();
        os <<direct_out<<"temp_Sigma_Phonon_t_"<<it<<block_name.str()<<".out";
        filename= os.str();
        //outputGF_new(latt.Sigma_Phonon_,filename);

 

        latt.Sigma_plus_Delta_.clear();
        latt.Sigma_plus_Delta_.incr(latt.Sigma_);
        latt.Sigma_plus_Delta_.incr(latt.Delta_);
        os.str("");
        os.clear();
        os <<direct_out<<"temp_Sigma_plus_Delta_t_"<<it<<block_name.str()<<".out";
        filename= os.str();
        //outputGF_new(latt.Sigma_plus_Delta_,filename);


 
        os.str("");
        os.clear();
        os<<direct_out<<"temp_occupation_plus_energy"<<block_name.str()<<".out";
        filename= os.str();
        print_to_file(time,occ_number,kin_energy,pot_energy,filename,interval,it,first_time,matrix_der_F);//print in append mode

        os.str("");
        os.clear();
        os<<direct_out<<"temp_matrix_G"<<block_name.str()<<".out";
        filename= os.str();
        //print_to_file_matrix_G(time,matrix_G,filename,interval,it,first_time_matrix_G);//print in append mode

        os.str("");
        os.clear();
        os<<direct_out<<"temp_matrix_der_F"<<block_name.str()<<".out";
        filename= os.str();
        //print_to_file_matrix_der_F(time,matrix_der_F,filename,interval,it,first_time_matrix_der_F);//print in append mode



        os.str("");
        os.clear();
        os<<direct_out<<"temp_matrix_Sigma"<<block_name.str()<<".out";
        filename= os.str();
        //print_to_file_matrix_G(time,matrix_Sigma,filename,interval,it,first_time_matrix_Sigma);//print in append mode


        os.str("");
        os.clear();
        os<<direct_out<<"temp_matrix_Delta"<<block_name.str()<<".out";
        filename= os.str();
        //print_to_file_matrix_G(time,matrix_Delta,filename,interval,it,first_time_matrix_Delta);


        os.str("");
        os.clear();
        os<<direct_out<<"temp_matrix_Sigma_plus_Delta"<<block_name.str()<<".out";
        filename= os.str();
        //print_to_file_matrix_G(time,matrix_Sigma_plus_Delta,filename,interval,it,first_time_matrix_Sigma_plus_Delta);//print in append mode



        }

		
        
    }	
	
        fftw_cleanup();
	return 0;
}

void print_to_file(vector <int>& time,vector <double>& occ_number,vector <double>& kin_energy,vector <double> & pot_energy,std::string filename,int interval,int it, int & first_time, vector <vector <double>> & matrix_G   ){
	int start,stop;
	std::ofstream outputFile;
	//outputFile.open(filename.c_str(),std::ios_base::app);
	if(first_time==1){
	outputFile.open(filename.c_str(),std::ios_base::out);
	outputFile << "# time" << "\t" << "Occupation" << "\t"<< "Kin_energy" <<"\t"<<"Pot Energy"<<"\t"<<"beta eff"<<std::endl;
	}
	else {
	outputFile.open(filename.c_str(),std::ios_base::app);
	}
	start=it-interval;
	stop=it;
	for (int i = start; i<stop; i++)
	{
		outputFile <<std::setprecision(10) << time[i] << "\t" << occ_number[i] << "\t"<< kin_energy[i] <<"\t"<<pot_energy[i]<<"\t";//<<std::endl;
                for (size_t j=0;j<matrix_G.size();j++){
                        outputFile<<std::setprecision(10)<<-4.*matrix_G[j][i]<<"\t";
                }
        outputFile<<std::endl;
	}	
	outputFile.close();
	first_time=0;
}

   void print_to_file_matrix_G(vector <int>& time,vector <vector <double>> & matrix_G,std::string filename,int interval,int it, int & first_time){
        int start,stop;
        std::ofstream outputFile;
        if(first_time==1){
        outputFile.open(filename.c_str(),std::ios_base::out);
        outputFile << "# time" << "\t" << "Im_G_ret_omega0" << "\t"<< "Im_G_les_omega0" <<"\t"<<"F_G_omega0"<<"\t"<<"der_omega0_Im_G_ret" << "\t"<< "der_omega0_Im_G_les" <<"\t"<<"der_omega0_F_G"<<std::endl;
        }
        else {
        outputFile.open(filename.c_str(),std::ios_base::app);
        }
        start=it-interval;
        stop=it;
        for (int i = start; i<stop; i++)
        {
                outputFile <<std::setprecision(10)<<time[i]<<"\t"<< matrix_G[0][i] << "\t"<< matrix_G[1][i] <<"\t"<<matrix_G[2][i]<<"\t"<< matrix_G[3][i]<<"\t"<< matrix_G[4][i] <<"\t"<<matrix_G[5][i]<<std::endl;
        }
        outputFile.close();
        first_time=0;
}


void print_to_file_matrix_der_F(vector <int>& time,vector <vector <double>> & matrix_G,std::string filename,int interval,int it, int & first_time){
        int start,stop;
        std::ofstream outputFile;
        if(first_time==1){
        outputFile.open(filename.c_str(),std::ios_base::out);
        }
        else {
        outputFile.open(filename.c_str(),std::ios_base::app);
        }
        start=it-interval;
        stop=it;
        for (int i = start; i<stop; i++)
        {
                outputFile <<std::setprecision(10)<<time[i]<<"\t";
                for (size_t j=0;j<matrix_G.size();j++){
                        outputFile<<std::setprecision(10)<<matrix_G[j][i]<<"\t";
                        //cout<<"matrix j,i = "<<matrix_G[j][i]<<endl;
                }
        outputFile<<std::endl;
        }
outputFile.close();
first_time=0;
}




/*
void print_to_file_equilibrium(vector <int>& time,vector <double>& occ_number,vector <double>& kin_energy,vector <double> & pot_energy,std::string filename,int interval,int it, int first_time,int numOfFreq){
	int start,stop;
	std::ofstream outputFile;
	outputFile.open(filename.c_str(),std::ios_base::app);

	 outputFile << "# Number of Frequency points" << "\t" << "Occupation" << "\t"<< "Kin_energy" <<"\t"<<"Pot Energy"<<std::endl;
	
	outputFile  << numOfFreq << "\t" << occ_number[0] << "\t"<< kin_energy[0] <<"\t"<<pot_energy[0]<<std::endl;
		
	outputFile.close();

}
*/

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
	std::cout << i << " lines read from data file" << std::endl;
	
	readFile.close();
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
