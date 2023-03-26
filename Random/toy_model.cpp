//This script is a small simulation to test the per triger yields

//C++ libraries
#include <iostream>
#include <cmath>

using namespace std;

double factorial(double i){
	if (i < 0){
	 cout<<"No factorial of negatve"<<endl;
	 return 0;
	}
	if (i == 0){
		double fact = 1;
		return fact;
	}
	double fact = 1;
	for(int j = 1; j<=i;j++){
		fact = fact*j;
	}
	
	return fact;	 
}

void toy_model(){

	int Event = 1000000;
	double positive, negative;
	double pos_array[Event];
	double neg_array[Event];
	double SS_pairs, OS_pairs;
	//One can change these two variables
	double dom_prod = 20;//Number of particles 90% of the time
	double sub_prod = 21;//10% of the time
	
	for(int i = 0; i < Event;i++){
		pos_array[i] = dom_prod;
		neg_array[i] = dom_prod;
		if(i > 899999){
			pos_array[i] = sub_prod;
			neg_array[i] = sub_prod;
		}
		if(i%50000==0){
			cout<<"Positive particle produced is: "<<pos_array[i]<<endl;
		}
	}
	SS_pairs = 0;
	OS_pairs = 0;
	positive = 0;
	//Creating the pairs
	for(int i = 0; i < Event; i++){
	positive = positive + pos_array[i];
		if(pos_array[i] == dom_prod){
			OS_pairs = OS_pairs+pow(pos_array[i],2);
			if(pos_array[i] > 1){
				SS_pairs = SS_pairs+((factorial(pos_array[i])/factorial(pos_array[i]-2)));
			}
		}
		if(pos_array[i] == sub_prod){
			SS_pairs = SS_pairs+((factorial(pos_array[i])/factorial(pos_array[i]-2)));
			OS_pairs = OS_pairs+pow(pos_array[i],2); 
		}
		
		
	}

	cout<<"The number of pairs for the SS is: "<<SS_pairs<<" for the OS is:"<<OS_pairs<<endl;
	cout<<"The number of positive particles is: "<<positive<<endl;
	cout<<"The per triger yield for the SS is: "<<SS_pairs/positive<<" for the OS is: "<<OS_pairs/positive<<endl;
	
}
