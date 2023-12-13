//Author: Antonio Picano

#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>
#include "cntr.hpp"
#include "moving_herm.hpp"
#include "read_inputfile.hpp"
#include <iomanip>
#include <mpi.h>

using namespace std;

#define CPLX complex<double>
#define GREEN cntr::herm_matrix<double>
#define CFUNC cntr::function<double>
#define GTRUNC cntr::moving_herm<double>
#define CFUNC_t cntr::moving_function<double>


void print_energy(int tstp,int kt, GREEN &G, GREEN &Sigma, double beta,double h){ 
       	ostringstream os;
     	string filename;
     	os.str("");
     	os.clear();
     	os << "temp_energy.out";
     	filename= os.str();
     	const char *filename_c = filename.c_str();
	std::ofstream f;
	f.open(filename_c);
	f.precision(10);

            if(tstp==-1 || tstp>=kt){
            CPLX rho1,rho2,rho3;
            for(int n=-1;n<=tstp;n++){
                cntr::convolution_density_matrix(n,&rho1,G,Sigma,integration::I<double>(kt),beta,h);
                cntr::convolution_density_matrix(n,&rho2,G,G,integration::I<double>(kt),beta,h);
                rho3=G.density_matrix(n);
                    f << " OBS: t= " << n;
                    f << " dens= " << rho3.real();
                    f << " eint= " << rho1.real();
                    f << " ekin= " << rho2.real();
                    f << endl;
            }
	f.close();
        }

}
void print_G_tstp_Sigma_tstp(GTRUNC& G, GTRUNC & Sigma,int tstp){
	//here I print G
       	ostringstream os;
     	string filename;
     	os.str("");
     	os.clear();
     	os << "temp_G_tstp_" <<tstp<<".out";
     	filename= os.str();
     	const char *filename_c = filename.c_str();
	std::ofstream f;
	f.open(filename_c);
	f.precision(10);
	int tc=G.tc();
	f<<"# les 00 Re Im | les 01 Re Im | les 10 Re Im | les 11 Re Im |ret 00 Re Im | ret 01 Re Im | ret 10 Re Im | ret 11 Re Im |"<<endl;
	  for(int t=0;t<=tc;t++){
	    f<<t;
	    f<<" "<<G.lesptr(0,t)[0].real();//les
	    f<<" "<<G.lesptr(0,t)[0].imag();
	    f<<" "<<G.lesptr(0,t)[1].real();
	    f<<" "<<G.lesptr(0,t)[1].imag();
	    f<<" "<<G.lesptr(0,t)[2].real();
	    f<<" "<<G.lesptr(0,t)[2].imag();
	    f<<" "<<G.lesptr(0,t)[3].real();
	    f<<" "<<G.lesptr(0,t)[3].imag();

	    f<<" "<<G.retptr(0,t)[0].real(); //ret
	    f<<" "<<G.retptr(0,t)[0].imag();
	    f<<" "<<G.retptr(0,t)[1].real();
	    f<<" "<<G.retptr(0,t)[1].imag();
	    f<<" "<<G.retptr(0,t)[2].real();
	    f<<" "<<G.retptr(0,t)[2].imag();
	    f<<" "<<G.retptr(0,t)[3].real();
	    f<<" "<<G.retptr(0,t)[3].imag();

	    f<<endl;
	  }
	 f.close();
	 
	//here I print Sigma
       	ostringstream os1;
     	string filename1;
     	os1.str("");
     	os1.clear();
     	os1 << "temp_Sigma_tstp_" <<tstp<<".out";
     	filename1= os1.str();
     	const char *filename_c1 = filename1.c_str();


	std::ofstream f1;
	f1.open(filename_c1);
	f1.precision(10);
	//int tc=G.tc();
	f1<<"# les 00 Re Im | les 01 Re Im | les 10 Re Im | les 11 Re Im |ret 00 Re Im | ret 01 Re Im | ret 10 Re Im | ret 11 Re Im |"<<endl;
	  for(int t=0;t<=tc;t++){
	    f1<<t;
	    f1<<" "<<Sigma.lesptr(0,t)[0].real();//les
	    f1<<" "<<Sigma.lesptr(0,t)[0].imag();
	    f1<<" "<<Sigma.lesptr(0,t)[1].real();
	    f1<<" "<<Sigma.lesptr(0,t)[1].imag();
	    f1<<" "<<Sigma.lesptr(0,t)[2].real();
	    f1<<" "<<Sigma.lesptr(0,t)[2].imag();
	    f1<<" "<<Sigma.lesptr(0,t)[3].real();
	    f1<<" "<<Sigma.lesptr(0,t)[3].imag();


	    f1<<" "<<Sigma.retptr(0,t)[0].real(); //ret
	    f1<<" "<<Sigma.retptr(0,t)[0].imag();
	    f1<<" "<<Sigma.retptr(0,t)[1].real();
	    f1<<" "<<Sigma.retptr(0,t)[1].imag();
	    f1<<" "<<Sigma.retptr(0,t)[2].real();
	    f1<<" "<<Sigma.retptr(0,t)[2].imag();
	    f1<<" "<<Sigma.retptr(0,t)[3].real();
	    f1<<" "<<Sigma.retptr(0,t)[3].imag();

	    f1<<endl;
	  }
	 f1.close();
}
void print_all2x2components_G_t_Sigma_t(GTRUNC& G_t, GTRUNC & Sigma_t,int tstp,double u0,double u1){

	std::stringstream stream_u0;
	stream_u0 << std::fixed << std::setprecision(1) << u0;
	std::string s_u0 = stream_u0.str();    
	std::stringstream stream_u1;
	stream_u1 << std::fixed << std::setprecision(1) << u1;
	std::string s_u1 = stream_u1.str();


       	ostringstream os;
     	string filename;
	string folder =  "/home/vault/mptf/mptf007h/libcntr/output/ui_"+s_u0+"/uf_"+s_u1+"/data/";
     	os.str("");
     	os.clear();
     	os <<folder<<"temp_G_t_" <<tstp<<".out";
	filename= os.str();
     	const char *filename_c = filename.c_str();
	G_t.print_to_file(filename_c); 
}
void print_magnetization(double u0, double u1,double beta,int tc,int nt,int ntheta, int tmax, int tstp_max, vector<double>& rho_up, vector<double>& rho_down, vector<double>& rho_up_t, vector<double>& rho_do_t,int neg_const,int time_interval,bool & prima_volta){
 
    

    std::ofstream output_n;
    //I neead a way to fix the precision of u0 and u1
    
    std::stringstream stream_u0;
    stream_u0 << std::fixed << std::setprecision(1) << u0;
    std::string s_u0 = stream_u0.str();
    
    std::stringstream stream_u1;
    stream_u1 << std::fixed << std::setprecision(1) << u1;
    std::string s_u1 = stream_u1.str();

    std::stringstream stream_tmax;
    stream_tmax << std::fixed << std::setprecision(1) << tmax;
    std::string s_tmax = stream_tmax.str();

    string folder =  "/home/vault/mptf/mptf007h/libcntr/output/ui_"+s_u0+"/uf_"+s_u1+"/data/";

    string os_n = folder+"temp_n_tc_"+ std::to_string(tc)+ "_ui_"+ s_u0+ "_uf_" +s_u1+ "_k_"+std::to_string(ntheta)+"_tmax_"+s_tmax + ".out";
    //local quantitities
    output_n.open(os_n,std::ios_base::app); 
    output_n.precision(10);
    if (prima_volta == false){ 
    for(int tstp=-1;tstp<=tstp_max;tstp++){
	output_n<<" t "<<tstp;
	if(tstp<=nt) {output_n<<" n_up "<<rho_up[tstp+1]; output_n<<" n_down "<<rho_down[tstp+1];}
	else {output_n<<" n_up "<<neg_const; output_n<<" n_down "<<neg_const;}
	output_n<<" n_up_t "<<rho_up_t[tstp+1]; output_n<<" n_down_t "<<rho_do_t[tstp+1];
	output_n<<" n_up_t + n_do_t"<<rho_up_t[tstp+1]+rho_do_t[tstp+1];
	output_n<<" n_up_t - n_do_t"<<rho_up_t[tstp+1]-rho_do_t[tstp+1];
	output_n<<endl;
	}
    output_n.close();
    prima_volta=true; 
   }

   else{
	 for(int tstp=tstp_max-time_interval+1;tstp<=tstp_max;tstp++){
		output_n<<" t "<<tstp;	
		output_n<<" n_up "<<neg_const; output_n<<" n_down "<<neg_const;
		output_n<<" n_up_t "<<rho_up_t[tstp+1]; output_n<<" n_down_t "<<rho_do_t[tstp+1];
		output_n<<" n_up_t + n_do_t"<<rho_up_t[tstp+1]+rho_do_t[tstp+1];
		output_n<<" n_up_t - n_do_t"<<rho_up_t[tstp+1]-rho_do_t[tstp+1];
		output_n<<endl;
		}	
	output_n.close();
	}
	
}
void print_Gk_hk(CFUNC_t & hk_t,GTRUNC & Gk_t,int & tstp,int &tid, int & k, int num_threads_per_process  ){

ostringstream os;
string filename;
os.str("");
os.clear();
os <<"temp_hk_eps_"<<k<<"_tstp_"<<tstp<<".out";
filename= os.str();
const char *filename_c = filename.c_str();
hk_t.print_to_file(filename_c); 

ostringstream os1;
string filename1;
os1.str("");
os1.clear();
os1 << "temp_Gk_eps_" <<k-(tid*num_threads_per_process)<<"_"<<tid<<"_tstp_"<<tstp<<".out";
filename1= os1.str();
const char *filename_c1 = filename1.c_str();
Gk_t.print_to_file(filename_c1); 

}//end of the function






void get_sigma(int n,GREEN &Sigma,CFUNC &ufunc,GREEN &G,GREEN &Gdo){                   //for the non-truncated version 1x1

    cntr::herm_matrix_timestep<double> chi(n,G.ntau(),1,1);
    cntr::Bubble1(n,chi,Gdo,Gdo);
    chi.left_multiply(ufunc,1.0);
    chi.right_multiply(ufunc,-1.0);
    cntr::Bubble2(n,Sigma,G,chi);

}

void get_sigma(int n,GREEN &Sigma,CFUNC &ufunc,GREEN &G){//for the non-truncated version 2x2
  
    cntr::herm_matrix_timestep<double> chi(n,G.ntau(),1,1);    
    cntr::herm_matrix_timestep<double> G00(n,G.ntau(),1,-1);    
    cntr::herm_matrix_timestep<double> G11(n,G.ntau(),1,-1);    
    cntr::herm_matrix_timestep<double> S(n,G.ntau(),1,-1);   //all these functions I am considering are timestep functions    


    G00.set_matrixelement(0,0,G,0,0);  //the first two indexes are the component of the herm_matrix_timestep in which I copy the data. The second two indexes are the ones of the herm_matrix from  which I copy the data
    G11.set_matrixelement(0,0,G,1,1);  //this set_matrixelement does not specify the tstp because G11 is already a tstp function

    
    cntr::Bubble1(n,chi,G00,G00);   // at timestep n I define chi as the product between G00 and G00
    chi.left_multiply(ufunc,1.0);
    chi.right_multiply(ufunc,-1.0);
    cntr::Bubble2(n,S,G11,chi);
    Sigma.set_matrixelement(n,1,1,S); //Sigma is a 2x2 matrix. I am setting the component 11 of Sigma to S (that is a 1x1 matrix).
                                      //In this case, since Sigma is a herm_matrix function (and not a herm_matrix_timestep function), I have to specify that S (that is defined at tstp) is copied into Sigma at tstp

    cntr::Bubble1(n,chi,G11,G11);
    chi.left_multiply(ufunc,1.0);
    chi.right_multiply(ufunc,-1.0);
    cntr::Bubble2(n,S,G00,chi);
    Sigma.set_matrixelement(n,0,0,S);//Sigma is a 2x2 matrix. I am setting the component 00 of Sigma to S (that is a 1x1 matrix)
    
}


void get_sigma(GTRUNC &Sigma,cntr::moving_function<double> &ufunc,GTRUNC &G,GTRUNC &chi,GTRUNC &G_down){   //for the truncated version 1x1
    cntr::MovBubble1(chi,G_down,G_down);
    chi.left_multiply(ufunc,1.0);
    chi.right_multiply(ufunc,-1.0);
    cntr::MovBubble2(Sigma,G,chi);
}

void get_sigma(GTRUNC &Sigma,cntr::moving_function<double> &ufunc,GTRUNC &G){    //truncated version 2x2--- only one timeslice

  //these matrixes are useful when I will calculate chi and S only at the current timestep (and not for all the tstp from 0 to tc)
    cntr::moving_herm_timestep<double> chi(G.tc(),1,1);
    cntr::moving_herm_timestep<double> G00(G.tc(),1,-1);
    cntr::moving_herm_timestep<double> G11(G.tc(),1,-1);
    cntr::moving_herm_timestep<double> S(G.tc(),1,-1);
    cntr::moving_herm_timestep<double> G_tstp(G.tc(),2,-1);


    G.get_timestep(0,G_tstp); //I put in G_tstp 2x2 the last timestep of the GTRUNC G 2x2
    G00.set_matrixelement(0,0,G_tstp,0,0); //   void  set_matrixelement(...) come method della classe moving_herm_timestep : OK
    G11.set_matrixelement(0,0,G_tstp,1,1);
    cntr::MovBubble1(chi,G00,G00);    //   void  MovBubble1(moving_herm_timestep)
    chi.left_multiply(ufunc,1.0);    // left_multiply come method della classe moving_herm_timestep
    chi.right_multiply(ufunc,-1.0);
    cntr::MovBubble2(S,G11,chi);       // void MovBubble2 (moving_herm_timestep)
    Sigma.set_matrixelement(1,1,S,0,0);


    cntr::MovBubble1(chi,G11,G11);
    chi.left_multiply(ufunc,1.0);
    chi.right_multiply(ufunc,-1.0);
    cntr::MovBubble2(S,G00,chi);
    Sigma.set_matrixelement(0,0,S,0,0);
    
    

}

void set_hk(int tstp, vector<CFUNC>& hk_, CFUNC mu_vector, CFUNC mu_vector_down,vector<double> eps_k,int ntheta_i,int tc,int nt){ //for the non-truncated code

     cdmatrix hk_matrix(2,2); //2x2 matrix

     hk_matrix(0,0)= mu_vector[tstp]; //sto convertendo una CFUNC in una cdmatrix
     hk_matrix(0,1)= eps_k[ntheta_i];
     hk_matrix(1,0)= eps_k[ntheta_i];
     hk_matrix(1,1)= mu_vector_down[tstp];

     hk_[ntheta_i].set_value(tstp,hk_matrix); //at timestep tstp I set all the 2x2 components of hk_ to the value given by hk_matrix
}

void set_hk_t(int t1,vector<CFUNC_t>& hk_t, CPLX mu_vector, CPLX mu_vector_down,vector<double> eps_k,int ntheta_i){ //for the truncated code

     cdmatrix hk_matrix(2,2); //2x2 matrix

     hk_matrix(0,0)= mu_vector; //sto convertendo una CFUNC_t in una cdmatrix
     hk_matrix(0,1)= eps_k[ntheta_i];
     hk_matrix(1,0)= eps_k[ntheta_i];
     hk_matrix(1,1)= mu_vector_down;
     // cout << "set hk at " << t1 << " " << hk_matrix << endl; 
     hk_t[ntheta_i].set_value(t1,hk_matrix); //at timestep tstp I set all the 2x2 components of hk_t to the value given by hk_matrix 
}


void print_hkt(vector<CFUNC_t>& hk_t, int ntheta, int tc, int size,int iter ){

  std::ofstream output_hk_t; 

  string os = "hk_t_iter_"+to_string(iter)+"mpi.out";   
output_hk_t.open(os);
output_hk_t.precision(10);

 int element_size=size*size;

 for (int tstp=0; tstp<=tc; tstp++){
 output_hk_t<<"tstp "<<tstp;
 for(int i=0;i<=ntheta;i++){
   for(int matr_elem=0; matr_elem< element_size; matr_elem++){
       output_hk_t<<" "<<hk_t[i].ptr(tc-tstp)[matr_elem].real();
     }
 output_hk_t<<"      ";
   }
   output_hk_t<<endl<<endl;
 }
 output_hk_t.close();
}


void print_hk( vector <CFUNC>& hk, int ntheta, int tstp, int size,int iter ){

   

  string os = "hk_iter_"+ to_string(iter)+".out";  
std::ofstream output_hk;
  output_hk.open(os);
output_hk.precision(10);  
		   

 
 int element_size=size*size; 
 for (int time=0; time<=tstp; time++){

   output_hk<<"tstp "<<time;

   
 for(int i=0;i<=ntheta;i++){

   for(int matr_elem=0; matr_elem<element_size; matr_elem++){
     
       output_hk<<" "<<hk[i].ptr(time)[matr_elem].real();


     } // end of the cycle over matrix elements

     output_hk<<"      "; 
   }//end of the cycle over k-vectors
   
   output_hk<<endl;
 }//end of the cycle over tstp
  output_hk.close();


 
}

void print_tstp(GTRUNC &G,char const *file,int tstp,int iter){    //  1x1 truncated
  char name[200];
  std::ofstream f;  
  sprintf(name,"temp_%s_%d_%d.out",file,tstp,iter);
  f.open(name);
  f.precision(10);
  int tc=G.tc();
  for(int t=0;t<=tc;t++){
    f<<t;
    f<<" "<<G.lesptr(0,t)[0].real();
    f<<" "<<G.lesptr(0,t)[0].imag();
    f<<" "<<G.retptr(0,t)[0].real();
    f<<" "<<G.retptr(0,t)[0].imag();
    f<<endl;
  }
  f.close();
}


void print_tstp(GREEN &G,char const *file,int tstp,int iter){      // 1x1 non-truncated 
  char name[200];
  sprintf(name,"temp_%s_%d_%d.out",file,tstp,iter);
  std::ofstream f;
  f.open(name);
  f.precision(10);
  
  for(int t=tstp;t>=0;t--){
    f<<t;
    CPLX re,le;
    G.get_ret(tstp,tstp-t,re);
    G.get_les(tstp,tstp-t,le);
    f<<" "<<le.real();
    f<<" "<<le.imag();
    f<<" "<<re.real();
    f<<" "<<re.imag();
    f<<" "<<endl;
  }
  f.close();
}


void print_tstp(GTRUNC &G,char const *file,int tstp,int iter,int up_down){    //  2x2 truncated
  char name[200];
  std::ofstream f;  
  sprintf(name,"temp_%s_%d_%d.out",file,tstp,iter);
  f.open(name);
  f.precision(10);
  int tc=G.tc();
  for(int t=0;t<=tc;t++){
    f<<t;
    f<<" "<<G.lesptr(0,t)[up_down].real();
    f<<" "<<G.lesptr(0,t)[up_down].imag();
    f<<" "<<G.retptr(0,t)[up_down].real();
    f<<" "<<G.retptr(0,t)[up_down].imag();
    f<<endl;
  }
  f.close();
}

void print_tstp_G_theta_les_trunc(GTRUNC &G,char const *file,int tstp,int theta){
char name[200];
std::ofstream f;  
sprintf(name,"temp_%s_%d_les.out",file,theta);
f.open(name, std::ios::app); 
f.precision(10);

for (int up_down = 0;up_down< G.element_size();up_down ++ ){
	f<<tstp;
	f<<" "<<G.lesptr(0,0)[up_down].real();
	f<<" "<<G.lesptr(0,0)[up_down].imag();
	}

f<<endl;
f.close(); 
}

void print_tstp(GREEN &G,char const *file,int tstp,int iter,int up_down){      // 2x2 NON-truncated 
  char name[200];
  sprintf(name,"temp_%s_%d_%d.out",file,tstp,iter);
  std::ofstream f;
  f.open(name);
  f.precision(10);
 
  for(int t=tstp;t>=0;t--){
    f<<t;
    f<<" "<<-G.lesptr(tstp-t,tstp)[up_down].real();  //N.B the indexes of lesser are switched with respect to the one of retarded!!!
                                                     //the minus sign is to pass from G_less taken from a row to G_less in the column at tstp .
                                                     //In this way I can compare the output with the one of the truncated code (in which also the lesser component is taken from the vertical)
    f<<" "<<G.lesptr(tstp-t,tstp)[up_down].imag();
    f<<" "<<G.retptr(tstp,tstp-t)[up_down].real();
    f<<" "<<G.retptr(tstp,tstp-t)[up_down].imag();
    f<<" "<<endl;
  }
  f.close();
  
}

void print_tstp_G_theta_les_nontrunc(GREEN &G,char const *file,int tstp,int theta){
char name[200];
std::ofstream f;  
sprintf(name,"temp_%s_%d_les.out",file,theta);
f.open(name, std::ios::app);
f.precision(10);

for (int up_down = 0;up_down< G.element_size();up_down ++ ){
	f<<tstp;
	f<<" "<<G.lesptr(tstp,tstp)[up_down].real();
	f<<" "<<G.lesptr(tstp,tstp)[up_down].imag();
	//f<<endl;
	}

f<<endl;
f.close(); 


}


//vector<double> calculate_weights(int ntheta,vector<double> wk_){
void calculate_weights(int ntheta,vector<double>& wk_,vector<double>& eps_k){
assert(ntheta>1 && (ntheta%2)==0);
double norm=0.0;
double Delta_theta=(M_PI/2.0)/(ntheta); //this is because eps in the AFM goes from -2 to 0  
for(int i=0;i<=ntheta;i++){
  eps_k[i]=-2.0*cos(i /ntheta);
    double theta=i*Delta_theta;
    double simpson=(i==0 || i==ntheta ? 1.0 : (1+(i%2))*2.0);   // integral(f) = h/3 * Somma_su_tutte_le_i di  (simpson[i]*f[i])
    // wk_[i]=sin(theta)*sin(theta)*dtheta*(simpson/3.0);
    wk_[i]= Delta_theta/3.0 *simpson* (2.0/M_PI *sin(theta) * sin(theta) );
    norm+=wk_[i];
   } 
// renormalize weights to sum_k wk =1
for(int i=0;i<=ntheta;i++) wk_[i]=wk_[i]/norm;
}

void calculate_weights_4th_order(int ntheta,vector<double>& wk_, vector<double>& eps_k){
assert(ntheta>1 && (ntheta%2)==0);
double norm=0.0;
double Delta_theta=(M_PI/2.0)/(ntheta); //this is because eps in the AFM goes from -2 to 0  
for(int i=0;i<=ntheta;i++){
    double theta=i*Delta_theta;
    double simpson;
    eps_k[i]=-2.0*cos(i*Delta_theta);
    if (i==0 or i== ntheta) simpson=0.0;
    else if (i==1 or i==ntheta-1) simpson=55.0/24.0;
    else if (i==2 or i==ntheta-2) simpson=-1.0/6.0;
    else if (i==3 or i==ntheta-3) simpson=11.0/8.0;
    else simpson=1.0;
    wk_[i]= Delta_theta *simpson* (2.0/M_PI *sin(theta) * sin(theta) );  // integral(f) = h/3 * Somma_su_tutte_le_i di  (simpson[i]*f[i])
    norm+=wk_[i];
   } 
// renormalize weights to sum_k wk =1
for(int i=0;i<=ntheta;i++) wk_[i]=wk_[i]/norm;
}



void calculate_weights_sin_square(int ntheta,vector<double>& wk_,vector<double>& eps_k ){
assert(ntheta>1 && (ntheta%2)==0);
double norm=0.0;
double Delta_theta=(M_PI/2.0)/(ntheta); //this is because eps in the AFM goes from -2 to 0  
for(int i=0;i<=ntheta;i++){
    double theta=i*Delta_theta;
    
    double simpson;    
    if (i==0 or i== ntheta) simpson=0.0;
    else if (i==1 or i==ntheta-1) simpson=55.0/24.0;
    else if (i==2 or i==ntheta-2) simpson=-1.0/6.0;
    else if (i==3 or i==ntheta-3) simpson=11.0/8.0;
    else simpson=1.0;

    eps_k[i]=8.0/M_PI*(i* (M_PI/2.0)/ntheta*0.5 - 0.125* sin(4.0*  i* (M_PI/2.0)/ntheta )) -2.0;    //1) case one: max of derivative is 8/M_PI
    wk_[i]= Delta_theta *simpson* (    4.0/(M_PI*M_PI) * sin(2*theta)*sin(2*theta) *   sqrt ( 4.0 -  eps_k[i]*eps_k[i])         );  // integral(f) = h * Somma_su_tutte_le_i di  (simpson[i]*f[i])

    norm+=wk_[i];
   } 
// renormalize weights to sum_k wk =1
 for(int i=0;i<=ntheta;i++) {wk_[i]=wk_[i]/norm; 
 //cout<<"wk "<<i<<" = "<<wk_[i]<<endl; 
	}
}


void calculate_weights_log(int ntheta,vector<double>& wk_,vector<double>& eps_k ){
assert(ntheta>1 && (ntheta%2)==0);
double norm=0.0;
 double a = 11.0;
 double Delta_s = (log(M_PI/2.0) + a)/(ntheta); //this is because eps in the AFM goes from -2 to 0  
for(int i=0;i<=ntheta;i++){
    double s=-a + Delta_s*i;
    
    double simpson;    
    if (i==0 or i== ntheta) simpson=0.0;
    else if (i==1 or i==ntheta-1) simpson=55.0/24.0;
    else if (i==2 or i==ntheta-2) simpson=-1.0/6.0;
    else if (i==3 or i==ntheta-3) simpson=11.0/8.0;
    else simpson=1.0;

    eps_k[i]= -2.0* cos(M_PI/2.0 - exp(s));
    wk_[i]= Delta_s *simpson* ((-2.0/M_PI) *  exp(s) * sin(M_PI/2.0 - exp(s))  *  sin(M_PI/2.0 - exp(s)) );  // integral(f) = h * Somma_su_tutte_le_i di  (simpson[i]*f[i])

    norm+=wk_[i];
   }
 
// renormalize weights to sum_k wk =1
  for(int i=0;i<=ntheta;i++) {
	wk_[i]=wk_[i]/norm; 
	//cout<<"wk "<<i<<" = "<<wk_[i]<<endl;
 }
}





int main(int argc,char** argv){
  //************NB*******************fundamental parameters******************************//
  //int k = 1; //this is the multiple of 40 Gk per node. The number of nodes should be 64/k 
  int num_threads_per_process =4;   //this is the number of Gk calculated per process
  int nodes =16;  // total number of nodes. Prima 8
  int ppn = 40;  //processes per node
  // I expect world_size being nodes*ppn 

  //int max_threads = 40; 
 //Initialize the MPI environment
  int ret = MPI_Init(NULL,NULL);
  if (ret != MPI_SUCCESS){
	cout<<"error: could not initialize MPI";
	MPI_Abort(MPI_COMM_WORLD,ret);
  }
  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  
  
  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // Get the name of the processor
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  int root=0;
  int tid=world_rank;
  if (tid == root) cout<<"the world_size, i.e. the total amount of processes is "<<world_size<<". I expect it is nodes x ppn: "<<nodes*ppn<<endl;

  //I declare the variables in all the ranks such that if I want to pass them from the root process, there is room for them in the other ranks.
  //However, I initialize Sigma and Gloc only in the root process, because only there I need them.
    
   
    GREEN Gloc,Sigma;
    vector<GREEN> Gk; 
    CFUNC ufunc,mu_vector,mu_down_vector;
    cntr::moving_function<double> ufunc_t;
    GTRUNC G_t,Sigma_t,Sigma_do_t,chi_t,G_do_t; 
    std::vector<GTRUNC> Gk_t;
    std::vector<double> rho_up,rho_down,rho_up_t,rho_do_t;//ek
       
    
    int tc,nt,size,tstp,ntau,kt,tmax,neps,itermax=200,itermax_trunc=10,ntheta,print=100,restart, time_restart_init, time_restart_end;
    double beta,h,u0,u1,errmax=1e-2,Delta_ionic,mu,mu_down,neg_const(-0.1),nup,ndown,tout;
    CPLX tmp,eps0_t,eps0_do_t;
    int multiple=1,time_interval=500;
    bool prima_volta = false, prima_volta_nk = false; 
    {
		if (argv[1] == NULL) 
		{
			std::cout << "no input file!" << std::endl;
			return 1;
		}
	
	find_param(argv[1], "__tc=",tc);
        find_param(argv[1], "__nt=",nt);
        find_param(argv[1], "__beta=",beta);
        find_param(argv[1], "__ntau=",ntau);
        find_param(argv[1], "__ntheta=",ntheta);
        find_param(argv[1], "__h=",h);
        find_param(argv[1], "__tmax=",tmax);
        find_param(argv[1], "__u0=",u0);
 	find_param(argv[1], "__u1=",u1);
	find_param(argv[1], "__neps=",neps);
        find_param(argv[1], "__Delta_ionic=",Delta_ionic);
	find_param(argv[1], "__errmax=",errmax);
	find_param(argv[1], "__restart=",restart);
	find_param(argv[1], "__time_restart_init=",time_restart_init);
	find_param(argv[1], "__time_restart_end=",time_restart_end);
	}
     
          if(tid==root){
            cout<<"ntau = "<<ntau<<endl;
            cout<<"beta = "<<beta<<endl;
            cout<<"nt = "<<nt<<endl;
            cout<<"tc = "<<tc<<endl;
            cout<<"ntheta = "<<ntheta<<endl;
            cout<<"h = "<<h<<endl;
            cout<<"tmax = "<<tmax<<endl;
            cout<<"u0 = "<<u0<<endl;
            cout<<"Delta_ionic = "<<Delta_ionic<<endl;
            cout<<"neps = "<<neps<<endl;
            cout<<"neps = "<<neps<<endl;
            cout<<"errmax = "<<errmax<<endl;
            cout<<"u1 = "<<u1<<endl;
            cout<<"restart = "<<restart<<endl;
          }


    
   if(tmax<nt){
        cerr << "tmax<nt, why would you want this??" << endl;
        return 0;
    }

   //I try to allocate these variables only on the root process
    //if(tid==root){
    rho_up.resize(nt+2);
    rho_down.resize(nt+2);
    rho_up_t.resize(tmax+2);
    rho_do_t.resize(tmax+2);
    //}
    size=2;
    kt=5;
    mu=0.0;

    ntheta =  world_size * num_threads_per_process-1;   //aggiungo -1 dato che poi i vettori vanno da 0 a ntheta- 1, ntheta=319
    //cout<<"ntheta that I have allocated is "<<ntheta<<endl;


    //I need only one wk and one eps_k per node 
    vector<double> wk(ntheta+1,0.0),eps_k(ntheta+1,0.0);  //ntheta elements, each of which is initialized to 0.0. Ho 320 elementi: da 0 a 319    
    calculate_weights_sin_square(ntheta,wk,eps_k);   // this part will be present in each process : I want each process having all the wk and epsk 
       double t_quench = 8.0/h;




    if (restart == 0) {
      Gloc=GREEN(nt,ntau,size,-1);
      Gloc.clear();
      Sigma=GREEN(nt,ntau,size,-1);   
      Sigma.clear(); 

    
    //Here I am saying that each rank will have num_threads_per_process = 40 Gk     
    Gk.resize(num_threads_per_process);   //one GREEN Gk for each k-point

    #pragma omp parallel for 
    for(int i=0;i<num_threads_per_process;i++) {
	Gk[i]=GREEN(nt,ntau,size,-1); 
	Gk[i].clear(); 
	}
    
    //Gloc_MPI will be present in each rank. Then it will take the value of Sum (wk*Gk). Then all the Gloc_MPI will be reduced into Gloc.
    //GREEN Gloc_MPI;        
    //Gloc_MPI=GREEN(nt,ntau,size,-1); Gloc_MPI.clear();
    
    //for me it is fine that these local variables are defined in all the processes, since then all the processes need them in order to solve the dyson equation
    ufunc=CFUNC(nt,1);
    mu_vector=CFUNC(nt,1);
    mu_down_vector=CFUNC(nt,1);
    
     //I initialize hk and hk_t for all the processes. So at the beginning I have a copy of them in all the processes
     vector<CFUNC> hk;
     hk.resize(ntheta+1);  //320 elementi: da 0 a 319
     #pragma omp parallel for
     for(int k = tid * num_threads_per_process; k < (tid+1)*num_threads_per_process; k++  ) hk[k]=CFUNC(nt,size);
    
    
   
     vector<CFUNC_t> hk_t;
     hk_t.resize(ntheta+1);
     #pragma omp parallel for 
     for(int k = tid * num_threads_per_process; k < (tid+1)*num_threads_per_process; k++)  hk_t[k]=CFUNC_t(tc,size);    
          
    
   // double t_quench = 8.0/h;
    for (tstp=0;tstp<=nt;tstp++){
      if(tstp>=0 and tstp<=t_quench)    ufunc[tstp]=u0+(u1-u0)/t_quench*tstp;
      else ufunc[tstp]=u1;
      //ufunc[tstp]=u1;
    }
    ufunc[-1]=u0;
    // if(tid==root){
    cntr::green_equilibrium_mat_bethe(Gloc,beta);     //these two functions are used to calculate Sigma with second order perturbation theory      
    get_sigma(-1,Sigma,ufunc,Gloc);   //tstp is equal to -1
    //}
    
    //************************************************* NON truncated****************************************************************************************************************************************//
    
	/// DMFT SOLUTION 2nd order self-consistent selfenergy up to timestep nt:
	for(tstp=-1;tstp<=nt;tstp++){//loop over all the tstp

		int nt1=(tstp>=0 && tstp<=kt ? 0 : tstp);  // 0, if tstp between 0 and 5; tstp, otherwise
		int nt2=(tstp>=0 && tstp<=kt ? kt : tstp);// 5, if tstp between 0 and 5; tstp, otherwise
		if(tstp>kt) {
		  cntr::extrapolate_timestep(tstp-1,Gloc,integration::I<double>(kt));
			get_sigma(tstp,Sigma,ufunc,Gloc);
		}

		cntr::herm_matrix_timestep<double> gloc_temp(tstp,ntau,size); // to store the last timestep tstp. 2x2 matrix
		
		for(int iter=0;iter<=itermax;iter++){                 // DMFT loop at fixed tstp	  

		Gloc.get_timestep(tstp,gloc_temp); //gloc_temp is used for the convergence in DMFT 
                                      
			           for(int n=nt1;n<=tstp;n++){  //The calculation is repeated from 0 to tstp (when tstp<=kt), because Gloc at tstp=0 for example changes during the starting procedure
				    
				    if (n==-1){
				    mu_vector[n]= ufunc[n]* (-1.0)* Gloc.matptr(ntau)[3].real() - ufunc[n]*0.5; //  - U / 2 + U < n >
				    mu_down_vector[n]= ufunc[n]* (-1.0)* Gloc.matptr(ntau)[0].real() - ufunc[n]*0.5;
					  }
                                   else{
                                    mu_vector[n]= ufunc[n]* Gloc.lesptr(n,n)[3].imag() - ufunc[n]*0.5; //  - U / 2 + U < n >
                                    mu_down_vector[n]= ufunc[n]*Gloc.lesptr(n,n)[0].imag() - ufunc[n]*0.5;
                                        }


				    if(n==-1 and iter<=2){
				       mu_vector[n]-=Delta_ionic;
				       mu_down_vector[n]+=Delta_ionic;
				   	
				      }
				    #pragma omp parallel for
				    for(int k = tid * num_threads_per_process; k < (tid+1) * num_threads_per_process; k++){                                     
				  	set_hk(n,hk,mu_vector,mu_down_vector,eps_k,k,tc,nt);
						}
				    }//end of the nt1 cycle

		
		        //here I start a cycle over all the Gk[i]			
				   	
                                   #pragma omp parallel for
				   for(int k = tid * num_threads_per_process; k < (tid+1)*num_threads_per_process; k++){ 

					  if(tstp==-1){
				    cntr::dyson_mat_fixpoint(Gk[k-(tid * num_threads_per_process )],Sigma,mu,hk[k],integration::I<double>(kt),beta,3,6);  //hk is a CFUNC(nt,2)  
				    }else if(tstp==0){
					cntr::set_t0_from_mat(Gk [ k- (tid * num_threads_per_process ) ]);  
				    }else if(tstp<=kt){
					 cntr::dyson_start(Gk[k-(tid * num_threads_per_process)],mu,hk[k],Sigma,integration::I<double>(tstp),beta,h); 
				    }else{
					 cntr::dyson_timestep(tstp,Gk[k-(tid * num_threads_per_process)],mu,hk[k],Sigma, integration::I<double>(kt),beta,h); 
				    	}
				  
				  }// end of the cycle over all the Gk[i]

			//now that I have all the Gk I can sum them up to obtain - in EACH PROCESS - Gloc MPI.
				   
 		  	for(int n=nt1;n<=tstp;n++){
			  
			  
			    Gloc.set_timestep_zero(n);	 
			 	
		    	    for(int k = tid * num_threads_per_process;k< (tid +1) *num_threads_per_process ;k++) {Gloc.incr_timestep(n,Gk[k-(tid * num_threads_per_process)],wk[k]); }
			   
			    // Now I simply have to copy Gloc_MPI from all the processes into Gloc_MPI_tstp ...
			    // ... and then reduce these functions into Gloc -at timestep tstp(or n) - of the ROOT process.
			    
			    cntr::herm_matrix_timestep<double> Gloc_MPI_tstp(n,ntau,size); 
			  
			    Gloc.get_timestep(n,Gloc_MPI_tstp);
			    Gloc_MPI_tstp.MPI_Reduce(root);
			    Gloc_MPI_tstp.Bcast_timestep(n,ntau,size,root);
			    Gloc.set_timestep(n,Gloc_MPI_tstp);	
			    get_sigma(n,Sigma,ufunc,Gloc);
		}			        	
		        

		       double err = cntr::distance_norm2(tstp,gloc_temp,Gloc);

			
		       if(tid==root and (tstp % multiple) ==0){
		       cout << "t = " << tstp;
		       cout << " iter = " << iter;
		       cout << " err = " << err;
			if(tstp == -1){	
		       		cout << " n_up "<< -1.0*Gloc.matptr(ntau)[0].real();
		       		cout << " n_do "<< -1.0*Gloc.matptr(ntau)[3].real();
		       		cout << " |n_up - n_do| " << abs(Gloc.matptr(ntau)[0].real() -  Gloc.matptr(ntau)[3].real()) ;
		       		cout << " n_up + n_do " << -1.0*Gloc.matptr(ntau)[0].real() +(-1.0)* Gloc.matptr(ntau)[3].real();
				}
			else{

				cout << " n_up "<< Gloc.lesptr(tstp,tstp)[0].imag();
                       		cout << " n_do "<<Gloc.lesptr(tstp,tstp)[3].imag();
                       		cout << " |n_up - n_do| " << abs(Gloc.lesptr(tstp,tstp)[0].imag() - Gloc.lesptr(tstp,tstp)[3].imag()) ;
                       		cout << " n_up + n_do " <<Gloc.lesptr(tstp,tstp)[0].imag() + Gloc.lesptr(tstp,tstp)[3].imag();

				}

		       		cout << endl;
			}

			



		       if(err<errmax) break;    //and iter > 3(I want to try to comment this part)
					       

	                }//end of the loop over all DMFT-iter
            
			if (tstp==-1){
                		rho_up[tstp+1] = -1.0*Gloc.matptr(ntau)[0].real();
                		rho_down[tstp+1]=-1.0*Gloc.matptr(ntau)[3].real();
               			 }

                else {
                		rho_up[tstp+1]=Gloc.lesptr(tstp,tstp)[0].imag();
                		rho_down[tstp+1]= Gloc.lesptr(tstp,tstp)[3].imag();
                }


        }//end of the loop over all tstp

        if (tid==root) print_energy(tstp,kt,Gloc, Sigma, beta, h);   //end of all the loops -----just print the outputs	

	for(tstp=-1;tstp<=tc;tstp++){
	rho_up_t[tstp+1]=rho_up[tstp+1];
	rho_do_t[tstp+1]=rho_down[tstp+1];
	}  //the remaining tmax-tc points are stored after they have been calculated
          
	if (tid==root){                            
	 Gloc.print_to_file("G.out");
	 Sigma.print_to_file("Sigma.out");
	}
										   
    G_t=GTRUNC(tc,tc,size,-1);
    G_t.set_from_G_backward(tc,Gloc);

    Sigma_t=GTRUNC(tc,tc,size,-1);        
    Sigma_t.set_from_G_backward(tc,Sigma);
    ufunc_t=CFUNC_t(tc,1);

    for(tstp=0;tstp<=tc;tstp++){
      
      if(tstp<=t_quench) {ufunc_t[tc-tstp] = u0 + (u1 - u0)/t_quench*tstp;}
      else {ufunc_t[tc-tstp]=u1;}

    } 
    
   
    Gk_t.resize(num_threads_per_process);
   #pragma omp parallel for
    for(int i=0; i< num_threads_per_process; i++){
      Gk_t[i] = GTRUNC(tc,tc,size,-1);
      Gk_t[i].set_from_G_backward(tc,Gk[i]);
    }
    
  //  cout<<"process "<<tid<<endl; 
    cntr::moving_herm_timestep<double> gtmp(tc,size,-1); //this is a time slice of length tc+1
  //  cntr::moving_herm_timestep<double> G_t_MPI_tstp(tc,size,-1); 
    //cntr::moving_herm_timestep<double> gloc_tc_00(tc,1,-1),gloc_tc_11(tc,1,-1);
    double err(10.0);

    #pragma omp parallel for
    for (int k = tid* num_threads_per_process; k< (tid+1)*num_threads_per_process;k++){
		for(int t1=0;t1<=tc;t1++){ 
			set_hk_t(t1,hk_t, mu_vector[tc-t1], mu_down_vector[tc-t1],eps_k,k);
				}
		}    
    //******************************************************************************* now starts the TRUNCATED time-evolution over all the tstp *****************************************************************//
    // cout<<"pp "<<tid<<endl; 

    for (tstp=tc+1;tstp<=tmax;tstp++){              // 1) iteration over all the tstp
      //cout<<"ppp "<<tid<<endl; 
   
	//every time I start an iteration at the next tstp I make room in the last column of the quantities of interest        
        Sigma_t.forward();
        G_t.forward();
       
	#pragma omp parallel for
	for (int i = 0; i< num_threads_per_process ;i++){

	  Gk_t[i].forward();

	}  
        #pragma omp parallel for
       	for (int k = tid * num_threads_per_process; k< (tid +1)*num_threads_per_process; k++){
         hk_t[k].forward();
	}
        ufunc_t.forward();
	if(tstp<=t_quench) ufunc_t[0] = u0 + (u1 - u0)/t_quench*tstp;
        else ufunc_t[0]=u1;

	cntr::extrapolate_timestep(G_t,integration::I<double>(kt));
	get_sigma(Sigma_t,ufunc_t,G_t);
        //*************************************************************************TRUNCATED DMFT **********************************************************************************//
	
        for(int iter=0;iter<=itermax_trunc;iter++){    // 2) DMFT iteration at tstp fixed
	        // if(tstp==tc+1 and iter==0) cout<<"******************* TRUNCATED time evolution *************************"<<endl;
	         G_t.get_timestep(0,gtmp);  // gtmp is used to calculate the error at each iteration
	 
		 eps0_t =ufunc_t[0] * (G_t.lesptr(0,0)[3].imag()-0.5);    
		 eps0_do_t= ufunc_t[0] * (G_t.lesptr(0,0)[0].imag()-0.5);





        //*************************************************************************loop over all the Gk **********************************************************************************//
		 #pragma omp parallel for 
		 for (int k = tid * num_threads_per_process; k< (tid +1) * num_threads_per_process ;k++){

		       set_hk_t(0,hk_t,eps0_t,eps0_do_t,eps_k,k);

		       cntr::dyson_timestep( Gk_t[k - ( tid * num_threads_per_process )],Sigma_t,hk_t[k],mu,integration::I<double>(kt),h);


		 }
		// if (tid==0 and tstp==tc+1 and iter==0) print_hkt(hk_t,ntheta,tc,size,iter);

		 G_t.clear_timestep(0);

		 for(int k = tid * num_threads_per_process; k< (tid +1)*num_threads_per_process ;k++){

		        G_t.incr_timestep(0,Gk_t[ k- (tid * num_threads_per_process) ],0,wk[k]);
		 		                                                                  
		 }                                                                //il primo argomento di incr_timestep mi fissa la riga della mia matrice (G_t in questo caso)
		                                                                  //il terzo argomento fissa la riga della matrice da cui prendo i valori(Gk_t in questo caso)
		                                                                  //Se le matrici sono 2x2 eseguo questa operazioe per tutte e 4 le submatrici
		                                                                  //Fissando la riga zero in entrambe le matrici, sono all'ultimo timestep
		 
		 G_t.MPI_Reduce_timestep(root);
		 G_t.Bcast_timestep(root);		 
		 get_sigma(Sigma_t,ufunc_t,G_t);

		 err=cntr::distance_norm2(0,G_t,0,gtmp);
		                  
                 if(tid==root and (tstp % multiple) ==0){		 
		 cout << "t= " << tstp;
                 cout << " iter= " << iter; cout << " err= " << err;
		 cout <<" n_up_t "<<G_t.lesptr(0,0)[0].imag()<<" n_do_t "<<G_t.lesptr(0,0)[3].imag();
	         cout<< " |n_up_t - n_down_t|= "<<abs(G_t.lesptr(0,0)[0].imag() - G_t.lesptr(0,0)[3].imag());
	         cout<<" n_up_t + n_down_t = "<< G_t.lesptr(0,0)[0].imag() + G_t.lesptr(0,0)[3].imag();
                 cout << endl;
		 }
		 

                 if(err<errmax ) break;   //and iter>3

	    
        }//end of DMFT iteration at tstp

	rho_up_t[tstp+1] = G_t.lesptr(0,0)[0].imag();
	rho_do_t[tstp+1] = G_t.lesptr(0,0)[3].imag();

	if (tid==root and tstp%time_interval ==0){
	print_magnetization(u0,u1,beta,tc,nt,ntheta,tmax,tstp,rho_up,rho_down,rho_up_t,rho_do_t,neg_const,time_interval,prima_volta);
	 print_all2x2components_G_t_Sigma_t(G_t,Sigma_t,tstp,u0,u1); 
	}


}//end of the loop over all tstp

}//end of the restart == 0 branch
else if (restart == 1){
    
    std::vector <double> rho_up_restart, rho_do_restart;
    rho_up_restart.resize(time_restart_end-time_restart_init);
    rho_do_restart.resize(time_restart_end-time_restart_init);											   

    G_t=GTRUNC(tc,tc,size,-1);
    ostringstream os2;
    string filename2;
    os2.str("");
    os2.clear();
    os2 << "temp_G_t_"<<time_restart_init<<".out";
    filename2= os2.str();
    const char *filename_c2 = filename2.c_str();
    G_t.read_from_file(filename_c2);

    Sigma_t=GTRUNC(tc,tc,size,-1);    
    ostringstream os3;
    string filename3;
    os3.str("");
    os3.clear();
    os3 << "temp_Sigma_t_"<<time_restart_init<<".out";
    filename3= os3.str();
    const char *filename_c3 = filename3.c_str();
    Sigma_t.read_from_file(filename_c3);

    ufunc_t=CFUNC_t(tc,1);
    ufunc_t.read_from_file("ufunc_t.out");
     vector<CFUNC_t> hk_t;
     hk_t.resize(ntheta+1);
     Gk_t.resize(num_threads_per_process);

     #pragma omp parallel for
     for(int k = tid * num_threads_per_process; k < (tid+1)*num_threads_per_process; k++){   
        
	hk_t[k]=CFUNC_t(tc,size);   
       	ostringstream os;
     	string filename;
     	os.str("");
     	os.clear();
     	os << "temp_hk_eps_"<<k<<"_tstp_"<<time_restart_init<<".out";
     	filename= os.str();
     	const char *filename_c = filename.c_str();
	hk_t[k].read_from_file(filename_c); 
        


  //now I have to read in the proper way all the Gk_t
	
	Gk_t[k - tid * num_threads_per_process] = GTRUNC(tc,tc,size,-1);	
       	ostringstream os1;
     	string filename1;
     	os1.str("");
     	os1.clear();
     	os1 << "temp_Gk_eps_" <<k-(tid*num_threads_per_process)<<"_"<<tid<<"_tstp_"<<time_restart_init<<".out";
     	filename1= os1.str();
     	const char *filename_c1 = filename1.c_str();
	Gk_t[k - tid * num_threads_per_process].read_from_file(filename_c1); 
	}

    cntr::moving_herm_timestep<double> gtmp(tc,size,-1); //this is a time slice of length tc+1
    double err(10.0);

    for (tstp=time_restart_init+1;tstp<=time_restart_end;tstp++){              // 1) iteration over all the tstp  
	//every time I start an iteration at the next tstp I make room in the last column of the quantities of interest
        
        Sigma_t.forward();
        G_t.forward();
    
	#pragma omp parallel for
	for (int i = 0; i< num_threads_per_process ;i++){

	  Gk_t[i].forward();

	}  
        #pragma omp parallel for
       	for (int k = tid * num_threads_per_process; k< (tid +1)*num_threads_per_process; k++){
         hk_t[k].forward();
	}
        ufunc_t.forward();
	if(tstp<=t_quench) ufunc_t[0] = u0 + (u1 - u0)/t_quench*tstp;
        else ufunc_t[0]=u1;

	cntr::extrapolate_timestep(G_t,integration::I<double>(kt));
	get_sigma(Sigma_t,ufunc_t,G_t);
        //*************************************************************************TRUNCATED DMFT **********************************************************************************//
	
        for(int iter=0;iter<=itermax_trunc;iter++){    // 2) DMFT iteration at tstp fixed
	        // if(tstp==tc+1 and iter==0) cout<<"******************* TRUNCATED time evolution *************************"<<endl;
	         G_t.get_timestep(0,gtmp);  // gtmp is used to calculate the error at each iteration
	 
		 eps0_t =ufunc_t[0] * (G_t.lesptr(0,0)[3].imag()-0.5);    
		 eps0_do_t= ufunc_t[0] * (G_t.lesptr(0,0)[0].imag()-0.5);

        //*************************************************************************loop over all the Gk **********************************************************************************//
		 #pragma omp parallel for
		 for (int k = tid * num_threads_per_process; k< (tid +1) * num_threads_per_process ;k++){

		       set_hk_t(0,hk_t,eps0_t,eps0_do_t,eps_k,k);

		       cntr::dyson_timestep( Gk_t[k - ( tid * num_threads_per_process )],Sigma_t,hk_t[k],mu,integration::I<double>(kt),h);


		 }

		 G_t.clear_timestep(0);

		 for(int k = tid * num_threads_per_process; k< (tid +1)*num_threads_per_process ;k++){

		        G_t.incr_timestep(0,Gk_t[ k- (tid * num_threads_per_process) ],0,wk[k]);
		 		                                                                  
		 }                                                                //il primo argomento di incr_timestep mi fissa la riga della mia matrice (G_t in questo caso)
		                                                                  //il terzo argomento fissa la riga della matrice da cui prendo i valori(Gk_t in questo caso)
		                                                                  //Se le matrici sono 2x2 eseguo questa operazioe per tutte e 4 le submatrici
		                                                                  //Fissando la riga zero in entrambe le matrici, sono all'ultimo timestep

		 
		 G_t.MPI_Reduce_timestep(root);
		 G_t.Bcast_timestep(root);		 
		 get_sigma(Sigma_t,ufunc_t,G_t);
		 err=cntr::distance_norm2(0,G_t,0,gtmp);
                 if(tid==root){		 
		 cout << "t= " << tstp;
                 cout << " iter= " << iter; cout << " err= " << err;
		 cout <<" n_up_t "<<G_t.lesptr(0,0)[0].imag()<<" n_do_t "<<G_t.lesptr(0,0)[3].imag();
	         cout<< " |n_up_t - n_down_t|= "<<abs(G_t.lesptr(0,0)[0].imag() - G_t.lesptr(0,0)[3].imag());
	         cout<<" n_up_t + n_down_t = "<< G_t.lesptr(0,0)[0].imag() + G_t.lesptr(0,0)[3].imag();
                 cout << endl;
		 }
                 if(err<errmax ) break;   //and iter>3

	    
        }//end of DMFT iteration at tstp

	rho_up_restart[tstp - (time_restart_init+1)] = G_t.lesptr(0,0)[0].imag();
	rho_do_restart[tstp - (time_restart_init+1)] = G_t.lesptr(0,0)[3].imag();

		
	if (tid==root and tstp%time_interval ==0){
		 print_all2x2components_G_t_Sigma_t(G_t,Sigma_t,tstp,u0,u1); 
    }

tout=tstp;
	}//end of the loop over all the tstp for restart == 1 branch 

     if (tid==root){
       
    ostringstream os2;
    string filename2;
    os2.str("");
    os2.clear();
    os2 << "temp_G_t_"<<time_restart_end<<".out";
    filename2= os2.str();
    const char * filename_c2 = filename2.c_str();
    G_t.print_to_file(filename_c2);

        
    ostringstream os3;
    string filename3;
    os3.str("");
    os3.clear();
    os3 << "temp_Sigma_t_"<<time_restart_end<<".out";
    filename3= os3.str();
    const char * filename_c3 = filename3.c_str();
    Sigma_t.print_to_file(filename_c3);
    }
	#pragma omp parallel for
	for(int k = tid * num_threads_per_process; k < (tid + 1 ) * num_threads_per_process; k++){
	print_Gk_hk(hk_t[k], Gk_t[k- tid*num_threads_per_process],time_restart_end,tid,k,num_threads_per_process);
	}

if (tid==root) cout<<"tout = "<<tout<<endl;

  }//end of the restart==1 branch	

    MPI_Finalize();
    return 0;
}




