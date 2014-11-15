#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>

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
std::map<int, std::vector<int> > Group;

//========================================================================
double H(const double *xx)
{
  const double Ms     = xx[0];
  const double K1     = xx[1];
  const double K2     = xx[2];
  const double Phi_H  = xx[3];
  const double Phi_eq = xx[4];

  const double Eq1   = 1/(sin(Phi_H-Phi_eq));
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
  return Eq1+Eq2+Eq3;
}

double F(const double *xx)
{
  const double gamma = 1.75878e7; // 2 * 8.7939E6 rad/(s.Oe)
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
  phimin->SetTolerance(0.001);
  phimin->SetPrintLevel(0);
  phimin->SetFunction(funH);

  // min parameters
  double minstep[6]    = {0.0,  0.0,  0.0,  0.0,    0.0,     0.0};
  double startpoint[6] = {Ms,   K1,   K2,   Phi_H,  Phi_eq,  hpoint};
  int    randomSeed = 100000;

  TRandom2 r(randomSeed);
  startpoint[4] = r.Uniform(0.000001,M_PI/2);

  std::cout << "startpoint of Phi_eq " << startpoint[4] << std::endl;

  // Set the free variables to be minimized!
  phimin->SetVariable(0,"Ms",              startpoint[0], minstep[0]);
  phimin->SetVariable(1,"K1",              startpoint[1], minstep[1]);
  phimin->SetVariable(2,"K2",              startpoint[2], minstep[2]);
  phimin->SetVariable(3,"Phi_H",           startpoint[3], minstep[3]);
  phimin->SetLimitedVariable(4,"Phi_eq",   startpoint[4], minstep[4], 0.000001, M_PI/2);
  phimin->SetVariable(5,"hpoint",          startpoint[5], minstep[5]);


  // do the minimization
  phimin->Minimize();

  const double *xs = phimin->X();

  return xs[4];
}

//========================================================================
void readfile(std::string& datafilename){
	std::ifstream datafile(datafilename.c_str());

	std::string databuff;
	getline(datafile,databuff);

	while(getline(datafile,databuff)){
		std::vector<std::string> substring;
		boost::split( substring, databuff, boost::is_any_of( ";" ), boost::token_compress_on );
		Temperature.push_back(boost::lexical_cast<int>(substring[0]));
		Magneticfield.push_back(boost::lexical_cast<double>(substring[1]));
		Frequency.push_back(boost::lexical_cast<double>(substring[2]));
	}
}

void printdata(){
	int size = Temperature.size();
	for(int i = 0; i < size; i++){
		std::cout << i << " " << Temperature[i] << ";" << Magneticfield[i] << ";" << Frequency[i] << std::endl;
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

	const double Ms     = 1258.;
	const double K1     = 1.3e6;
	const double K2     = 20.1e3;
	const double Phi_H  = 0.1;
	const double Phi_eq = 0.1;
	const double HPoint = 4000.;
        const double pars[5] = {Ms,K1,K2,Phi_H,Phi_eq};
        double phi = GetPhiFromH(pars,HPoint);
        std::cout << "GetPhiFromH: " << phi << std::endl;
	std::cout << "H: " << H(pars) << std::endl;

	double H,F;

        const double newpars[5] = {Ms,K1,K2,Phi_H,phi};
	HvsF(newpars, &H,&F);
        std::cout << "HvsF H: " << H << " F: " << F << std::endl;


}

