#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

using namespace std;

//prototipi
struct retta
{
    double a;
    double b;
};

int main() {
    //dichiarazione delle variabili
	string filename;
    vector<double> x, y;
    vector<double> x2, y2;
	double xval,yval,sigma1,sigma2,n1,n2, x_cent1,x_cent2;
    int cov=0;
    retta r1;
    retta r2;
    double sxyres=0.0,syres=0.0, sxyres2=0.0,syres2=0.0;
    string filename2;

	cout<<"___________STIMA SIGNIFICANZA SPOSTAMENTO PICCO_____ "<<endl;
    //operazioni dal file e presa dati
	cout << "Inserire il nome del file di input numero 1: ";
    cin >> filename;
    cout << "Inserire il nome del file di input numero 2: ";
    cin >> filename2;
    ifstream fin(filename);
    ifstream fin2(filename);

    if(!fin){
    	cout<<"Il file non e' stato caricato correttamente"<<endl;
    	return -1;
	}
    if(!fin2){
    	cout<<"Il file non e' stato caricato correttamente"<<endl;
    	return -1;
	}
    cout<<"Inserisci i due valori x min e x max dell'intervallo in cui è possibile apprezzare un picco "<<endl;
    cout<<" x_min: ";
    cin>>n1;
    cout<<" x_max: "; 
    cin>>n2;

	while(fin>>xval>>yval){
		if(xval>=n1 && xval<=n2){
			if(yval!=0){
				x.push_back(xval);    		
          	    y.push_back(yval);
			}
        }
    }
    while(fin2>>xval>>yval){
		if(xval>=n1 && xval<=n2){
			if(yval!=0){
				x2.push_back(xval);    		
          	    y2.push_back(yval);
			}
        }
    }
    
    cout<<"-------------Ora inserisci le equazioni delle rette dei due spettri ( del tipo y=a+bx)------------ "<<endl;
    cout<<"a: ";
    cin>>r1.a;
    cout<<"b: ";
    cin>>r1.b;
    cout<<"Mentre la seconda: "<<endl;
    cout<<"a': ";
    cin>>r2.a;cout<<"b': ";
    cin>>r2.b;

    double diffx=0.0;
    double diffx2=0.0;

    //sommatorie intermedie
    for(int i=0;i<y.size();i++){
        sxyres+=x.at(i)*(y.at(i)-(r1.a+r1.b*x.at(i)));
        syres+=y.at(i)-(r1.a+r1.b*x.at(i));
        diffx+=pow((x.at(i)-x_cent1),2)*y.at(i);
    }
    x_cent1=sxyres/syres;

    for(int i=0;i<y2.size();i++){
        sxyres2+=x2.at(i)*(y2.at(i)-(r2.a+r2.b*x2.at(i)));
        syres2+=y2.at(i)-(r2.a+r2.b*x2.at(i));
        diffx2+=pow((x2.at(i)-x_cent2),2)*y2.at(i);
    }
    x_cent2=sxyres2/syres2;

    //calcolo incertezze contando che la distribuzione è di poisson
    sigma1=sqrt(diffx)/syres;
    sigma2=sqrt(diffx2)/syres2;

    //generiamo un file che contiene i vari residui per poterli rappresentare graficamente

    ofstream fout("residui2aprile.txt");
    ofstream fout2("residui3aprile.txt");

    for(int i=0;i<y.size();i++){
            fout << x.at(i) << " " << (y.at(i)-(r1.a+r1.b*x.at(i)))<< endl;
    }
    for(int i=0;i<y2.size();i++){
            fout2<< x2.at(i) << " " << (y2.at(i)-(r2.a+r2.b*x2.at(i)))<< endl;
    }

    cout<<"Il centroide del primo spettro : x_centroide = "<<x_cent1<<" con sigma: "<<sigma1<<endl;
    cout<<"Il centroide del secondo spettro : x_centroide = "<<x_cent2<<" con sigma: "<<sigma2<<endl;
    cout<<"Lo spostamento delta x dei picchi e' pari a: "<<abs(x_cent1-x_cent2)<<" con una incertezza sigma = "<<sqrt(pow(sigma1,2)+pow(sigma2,2))<<endl;

    return 0;    
}
