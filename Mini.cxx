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
}

