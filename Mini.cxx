#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <time.h>

#include "TGraph.h"
#include "TMarker.h"
#include "TCanvas.h"
#include "TMultiGraph.h"

#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "TRandom2.h"
//========================================================================
std::vector<int> Temperature;
std::vector<double> Magneticfield;
std::vector<double> Frequency;
std::vector<double> Phi_h;
std::map<int, std::vector<int> > Group;

double templist[4];
double chi2list[4];
double k1list[4];

//========================================================================
double H(const double *xx)
{
  const double Ms     = xx[0];
  const double K1     = xx[1];
  const double K2     = xx[2];
  const double Phi_H  = xx[3];
  const double Phi_eq = xx[4];

  const double Eq1   = 1/(Ms*sin(Phi_H-Phi_eq));
  const double Eq2   = (4*M_PI*Ms*Ms-2*K1-4*K2)*sin(Phi_eq)*cos(Phi_eq);
  const double Eq3   = 4*K2*cos(Phi_eq)*sin(Phi_eq)*sin(Phi_eq)*sin(Phi_eq);

  return Eq1*(Eq2+Eq3);
}

double H1(const double *xx)
{
  const double Ms     = xx[0];
  const double K1     = xx[1];
  const double K2     = xx[2];
  const double Phi_H  = xx[3];
  const double Phi_eq = xx[4];

  const double H_a1 = 2*K1/Ms;
  const double H_a2 = 4*K2/Ms;
  const double H_a  = H_a1+H_a2;

  const double Eq1   = H(xx)*cos(Phi_H-Phi_eq);
  const double Eq2   = (4*M_PI*Ms-H_a)*cos(2*Phi_eq);
  const double Eq3   = H_a2*(3*sin(Phi_eq)*sin(Phi_eq)*cos(Phi_eq)*cos(Phi_eq)-sin(Phi_eq)*sin(Phi_eq)*sin(Phi_eq)*sin(Phi_eq));

  return Eq1+Eq2+Eq3;
}

double H2(const double *xx)
{
  const double Ms      = xx[0];
  const double K1      = xx[1];
  const double K2      = xx[2];
  const double Phi_H   = xx[3];
  const double Phi_eq  = xx[4];

  const double H_a1 = 2*K1/Ms;
  const double H_a2 = 4*K2/Ms;
  const double H_a  = H_a1+H_a2;

  const double Eq1   = H(xx)*cos(Phi_H-Phi_eq);
  const double Eq2   = (4*M_PI*Ms-H_a)*sin(Phi_eq)*sin(Phi_eq);
  const double Eq3   = H_a2*sin(Phi_eq)*sin(Phi_eq)*sin(Phi_eq)*sin(Phi_eq);

  return Eq1-Eq2-Eq3;
}

double F(const double *xx)
{
  const double gamma = 1.75878e7; // 2 * 8.7939E6 rad/(s.Oe)
  //std::cout << "F H1:" << H1(xx) << " H2: " << H2(xx) << std::endl;
  return gamma/(2*M_PI)*sqrt(H1(xx)*H2(xx));	
}

void HvsF(const double *xx, double *h, double *f){
  *h = H(xx);
  *f = F(xx);
}

//========================================================================
double AbsDiffToH(const double *pars)
{
  double xx[5] = {pars[0], pars[1], pars[2], pars[3], pars[4]};
  double hpoint = pars[5];

  return abs(H(xx)-hpoint);
}

double GetPhiFromH(const double *pars, const double hpoint){
  double Ms     = pars[0];  
  double K1     = pars[1]; 
  double K2     = pars[2]; 
  double Phi_H  = pars[3]; 
  double Phi_eq = pars[4]; 


  ROOT::Math::Functor funH(&AbsDiffToH,6); 

  // initialize minimizer
  const std::string minName = "Minuit2";
  const std::string algName = "";

  ROOT::Math::Minimizer* phimin = 
    ROOT::Math::Factory::CreateMinimizer(minName.c_str(), algName.c_str());

  phimin->SetMaxFunctionCalls(1000000); 
  phimin->SetMaxIterations(10000);  
  phimin->SetTolerance(0.1);
  phimin->SetPrintLevel(0);
  phimin->SetFunction(funH);

  // min parameters
  //double low = M_PI/4;
  double low = 0.0001;
  double up  = M_PI/2;
  double minstep[6]    = {0.0,  0.0,  0.0,  0.0,    0.1,     0.0};
  double startpoint[6] = {Ms,   K1,   K2,   Phi_H,  Phi_eq,  hpoint};
  int    randomSeed = time(NULL);

  TRandom2 r(randomSeed);

  while(true){
	  startpoint[4] = r.Uniform(low,up);

	  //std::cout << "startpoint of Phi_eq " << startpoint[4] << std::endl;
	  //std::cout << "want to fit to H:    " << hpoint << std::endl;

	  // Set the free variables to be minimized!
	  phimin->SetVariable(0,"Ms",              startpoint[0], minstep[0]);
	  phimin->SetVariable(1,"K1",              startpoint[1], minstep[1]);
	  phimin->SetVariable(2,"K2",              startpoint[2], minstep[2]);
	  phimin->SetVariable(3,"Phi_H",           startpoint[3], minstep[3]);
	  phimin->SetLimitedVariable(4,"Phi_eq",   startpoint[4], minstep[4], low, up);
	  phimin->SetVariable(5,"hpoint",          startpoint[5], minstep[5]);


	  // do the minimization
	  phimin->Minimize();

	  const double *xs = phimin->X();

	  if(abs(H(xs)-hpoint)<0.1)return xs[4];
  }

  // return xs[4];
}

//========================================================================
void readfile(std::string& datafilename){
	std::ifstream datafile(datafilename.c_str());

	std::string databuff;
	getline(datafile,databuff);
	std::cout << "read line: " << databuff << std::endl; 

	while(getline(datafile,databuff)){
		std::vector<std::string> substring;
		std::cout << "read line: " << databuff << std::endl; 
		boost::split( substring, databuff, boost::is_any_of( ";" ), boost::token_compress_on );
		Temperature.push_back(boost::lexical_cast<int>(substring[0]));
		Magneticfield.push_back(boost::lexical_cast<double>(substring[1]));
		Frequency.push_back(boost::lexical_cast<double>(substring[2])*1e9);
		Phi_h.push_back(boost::lexical_cast<double>(substring[3]));
	}
}

void printdata(){
	int size = Temperature.size();
	for(int i = 0; i < size; i++){
		std::cout << i << " " << Temperature[i] << ";" << Magneticfield[i] << ";" << Frequency[i] << ";" << Phi_h[i] << std::endl;
	}
}

void group(){

	int t = Temperature[0];
	std::vector<int> *pVec = new std::vector<int>;
	for(uint32_t i = 0; i < Temperature.size(); i++){
		if(Temperature[i] != t){
			Group.insert(std::pair<int,std::vector<int> >(t, *pVec));
			pVec = new std::vector<int>;
			t = Temperature[i];
		}
		pVec->push_back(i);
	}
	Group.insert(std::pair<int,std::vector<int> >(t, *pVec));

}

void printgroup(){
	std::map<int, std::vector<int> >::iterator it ;
	std::vector<int>::iterator vit;
	for(it = Group.begin(); it != Group.end(); it++){
		std::cout << it->first << std::endl;
		for(vit = it->second.begin(); vit != it->second.end(); vit++){
			std::cout << "    " << *vit << std::endl;
		}
	}
}

void plotdata(){
	TCanvas *cav  = new TCanvas("c1");

	std::map<int, std::vector<int> >::iterator mit = Group.begin();


	std::vector<TGraph*> gVect;
	std::vector<TGraph*>::iterator git;
	int color = 1;
	int style = 20;
	for(mit = Group.begin(); mit != Group.end(); mit++){
		int size = mit->second.size();
		double *hx  = new double[size];
		double *fy  = new double[size];

		std::cout << "Temperature: " << mit->first << " data size: " << size << std::endl;	

		for(int i = 0; i < size; i++){
			int idx = mit->second.at(i);
			//std::cout << idx << std::endl;
			hx[i] = Magneticfield[idx];
			fy[i] = Frequency[idx];
		}


		TGraph  *hvsf = new TGraph(size, hx, fy);
		hvsf->SetMarkerStyle(style);
		hvsf->SetMarkerColor(color);
		hvsf->SetMarkerSize(2);
		gVect.push_back(hvsf);

		color++;
		style++;
	}

	TMultiGraph *mgraph = new TMultiGraph();
	for(git = gVect.begin(); git != gVect.end(); git++){
		mgraph->Add(*git);
	}
	mgraph->Draw("AP");
	cav->Print("FunGraph.eps");

}

//========================================================================

double GetChi2ByTemp(const double *pars){
	const double P        = pars[0];
	const double Ms300    = pars[1];
	const double K1       = pars[2];
	const double K2       = pars[3];
	const int temperature = int(pars[4]);

	if(Group.find(temperature)==Group.end()){
		std::cout << "there is no this temperature " << temperature << " in data" << std::endl;
		return 0.;
	}
	std::vector<int> idx = Group[temperature];
	std::vector<int>::iterator it;


	//const double Ms300  = 1258.;
	const double Ms     = Ms300*(1-P*sqrt(temperature*temperature*temperature))/(1-P*5196.1524);

	const double Phi_eq = M_PI/4;

	std::vector<double> hv;
	std::vector<double> fv;
	std::vector<double> pv;
	for(it=idx.begin();it!=idx.end();it++){
		const double HPoint = Magneticfield[*it];
	        const double Phi_H  = Phi_h[*it];

		const double phipars[5] = {Ms,K1,K2,Phi_H,Phi_eq};
		double phi = GetPhiFromH(phipars,HPoint);

		double H,F;

		const double newpars[5] = {Ms,K1,K2,Phi_H,phi};
		HvsF(newpars, &H,&F);
		hv.push_back(H);
		fv.push_back(F);
		pv.push_back(phi);

	}

	double chi2 = 0;
	double diff = 0;
	it=idx.begin();
	for(uint32_t i = 0; i < fv.size(); i++){
		diff = Frequency[*it] - fv[i];
		chi2 += (diff*diff)/fv[i];
		//std::cout << "data point "<< *it << " H: " << Magneticfield[*it] << " F: " << Frequency[*it] << std::endl;
		//std::cout << "fitting point H: " << hv[i] << " F: " << fv[i] << " @phi:" << pv[i] << " diff: " << diff << " chi2: " << (diff*diff)/fv[i] << std::endl;
		it++;
	}
	//std::cout << "Temperature: " << temperature << " Chi2: " << chi2 << std::endl;
	return chi2;

}



double Getk1k2FromFixedTemperature(const double *pars, double temperature, double *kpars){
  double P           = pars[0];  
  double Ms0         = pars[1];  
  double K1          = pars[2]; 
  double K2          = pars[3]; 
  double Temperature = temperature; 

  ROOT::Math::Functor funH(&GetChi2ByTemp,5); 

  // initialize minimizer
  const std::string minName = "Minuit2";
  const std::string algName = "";

  ROOT::Math::Minimizer* k1k2min = 
    ROOT::Math::Factory::CreateMinimizer(minName.c_str(), algName.c_str());

  k1k2min->SetMaxFunctionCalls(1000000); 
  k1k2min->SetMaxIterations(10000);  
  k1k2min->SetTolerance(100);
  k1k2min->SetPrintLevel(0);
  k1k2min->SetFunction(funH);

  // min parameters
  double k1low = 1e7;
  double k1up  = 3e7;

  double minstep[5]    = {0.0,  0.0,  1e5,  0.0, 0.0};
  double startpoint[5] = {P,    Ms0,  K1,   K2, Temperature};
  int    randomSeed = time(NULL);

  TRandom2 r(randomSeed);

  double chi2min = -100.;

  int loop = 0;

  while(true){
	  loop++;
	  startpoint[2] = r.Uniform(k1low,k1up);
	  //startpoint[3] = r.Uniform(k2low,k2up);

	  std::cout << "startpoint of K1 " << startpoint[2] << std::endl;
	  std::cout << "startpoint of K2 " << startpoint[3] << std::endl;

	  // Set the free variables to be minimized!
	  k1k2min->SetVariable(0,"P",               startpoint[0], minstep[0]);
	  k1k2min->SetVariable(1,"Ms",              startpoint[1], minstep[1]);
	  k1k2min->SetLimitedVariable(2,"K1",       startpoint[2], minstep[2], k1low, k1up);
	  k1k2min->SetVariable(3,"K2",              startpoint[3], minstep[3]);
	  k1k2min->SetVariable(4,"Temperature",     startpoint[4], minstep[4]);


	  // do the minimization
	  k1k2min->Minimize();

	  const double *xs = k1k2min->X();
	  kpars[0] = xs[2];
	  kpars[1] = xs[3];

	  double newchi2 = GetChi2ByTemp(xs);
	  double oldchi2 = 0.0;
	  if(chi2min<0) chi2min=newchi2; 
	  else{

		  if(newchi2<chi2min){
			  oldchi2 = chi2min;
			  chi2min = newchi2;
		  }
		  if((abs(chi2min-oldchi2)<(oldchi2*0.01))and(loop>20)) break;
	  }
  }
  return chi2min;
}

void GetChi2ByAll(const double *pars){
	std::map<int, std::vector<int> >::iterator mit;
	double chi2 =0;
	double kpars[2];
	int idx=0;
	for(mit=Group.begin();mit != Group.end();mit++){
		const double temperature   = mit->first;
		chi2+=Getk1k2FromFixedTemperature(pars, temperature, kpars);
	        std::cout << "temperature " << temperature << " k1: " << kpars[0] << " chi2: " << chi2 << std::endl;
	        templist[idx]=temperature;
	        chi2list[idx]=chi2;
	        k1list[idx]=kpars[0];	
		idx++;
	}
}

int main(int argc, char *argv[]){
	if(argc != 2){
		std::cout << "Data file can not be loaded" << std::endl;
		return 0;
	}
	std::string datafilename(argv[1]);

	readfile(datafilename);
	//printdata();
	group();
	//printgroup();
	plotdata();

	//const int temperature      = 300;

	const double P      = 2.6679e-5;
	const double Ms300  = 1258.;
	const double Ms0    = Ms300;
	const double K1     = 13.2e6; 
	const double K2     = 150; 
	const double pars[4] = {P,Ms0,K1,K2};
	//double kpars[2];
	//double chi2 = Getk1k2FromFixedTemperature(pars, temperature, kpars);
	//std::cout << "temperature " << temperature << " k1: " << kpars[0] << " chi2: " << chi2 << std::endl;
        GetChi2ByAll(pars);
	for(int i = 0; i< 4; i++){
		std::cout <<  " temperature " << templist[i] << " chi2: " << chi2list[i] << " k1: " << k1list[i] << std::endl;
	}

}

