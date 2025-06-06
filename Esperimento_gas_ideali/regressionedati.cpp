#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
using namespace std;

void interpolazionesemplice(vector<double> x, vector<double> y, double &a, double&b, double &sigma_a, double &sigma_b);

void errorevolume(double b, vector<double> p, vector<double> y);


int main (){
	
	string filename;
	vector<double> x, y, p;
	double xval,yval,sigma;
    double zval; 
    double a,b;
    double sigma_a, sigma_b;

	cout<< "___________REGRESSIONE LINEARE V-1/P______ "<<endl;
    //operazioni dal file e presa dati
	cout << "Inserire il nome del file di input: ";
    cin >> filename;
    ifstream fin(filename);

    if(!fin){
    	cout<<"Il file non e' stato caricato correttamente"<<endl;
    	return -1;
    }
    while(fin>>xval>>yval>>zval){
	    x.push_back(1/xval);    		
    	y.push_back(yval);
    	p.push_back(xval);
    }  
    interpolazionesemplice(x,y,a,b,sigma_a, sigma_b);
    errorevolume(b,p,y);

double errcasb;
errcasb = sigma_b/b; 
cout << "L'errore casuale su b e'" << errcasb << endl; 



return 0; 
}


//INTERPOLAZIONE SEMPLICE CON INCERTEZZE SULLE Y TUTTE UGUALI
void interpolazionesemplice(vector<double> x, vector<double> y,double &a, double&b, double &sigma_a, double &sigma_b){  
    double sx = 0.0;   
    double sy = 0.0;  
    double sxx = 0.0;  
    double sxy = 0.0; 
    double n=x.size(); 
    double residui=0; 
    vector<double> res;
    string outputFileName;

	//calcolo somme intermedie
    for (int i = 0; i < n; i++) {
        sx  += x.at(i); 
        sy  += y.at(i);
        sxx += x.at(i) * x.at(i);
        sxy += x.at(i) * y.at(i);
    }

    // Calcolo del DELTA
    double delta = n*(sxx) - (sx * sx);

    // Stima dei parametri a e b
    a = ( (sxx * sy) - (sx * sxy) ) / delta;
    b = ( (n* sxy) - (sx * sy ) ) / delta;
        

 
    //calcolo sigma_post
    for(int i=0; i<n; i++){   
        residui+=pow(y.at(i)- (a+b*(x.at(i))),2);
        res.push_back(y.at(i)- (a+b*(x.at(i))));
    }
    //calcolosigma quadrato
    double sigma_postquadrato=residui/(n-2);

	double sigma_post= sqrt(sigma_postquadrato); 
	double xsomma=0.0;
	double varx=0.0;
    
    //calcolo varianza(x)
    for(int i=0; i< n; i++){
	    xsomma+= x.at(i);
	}
	double xmedio= xsomma/n;
	double scarti=0;
	for(int i=0; i<n; i++){
		scarti += pow(((x.at(i)-xmedio)),2);
    }
    varx=scarti/(n-1);
    
    //Calcolo delle incertezze relative ai parametri   
    sigma_a = sqrt ((sigma_postquadrato/n)*(1+((xmedio*xmedio)/(varx))));
    sigma_b = sqrt (sigma_postquadrato/(n*varx));

    //stampa risultati interpolazione
    cout << "=== RISULTATI DEL FIT LINEARE SEMPLICE per una retta y=a+bx ===" << endl;
    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "sigma a posteriori=" << sigma_post <<endl;
    cout << "sigma_a = " << sigma_a <<endl;
    cout << "sigma_b = " << sigma_b << endl;

    cout<<"Inserisci il nome del file di output con i residui: ";
    cin>>outputFileName;
    ofstream outputFile(outputFileName);
    for(int i=0;i<res.size();i++){
        outputFile<<res.at(i)<< '\n';
    }
    cout << "Dati dei residui salvati in '" << outputFileName << "'.\n";
}

void errorevolume(double b, vector<double> p, vector<double> y){
//calcolo deltav
double c;
double miny = y.at(0);
double maxy = y.at(0);
for(auto c:y){
    if(c>maxy){
        maxy = c;
    }
	if(c<miny){
        miny = c;
    }
}
	

double deltav = maxy - miny;

//calcolo deltap
double d;
double minp = p.at(0);
double maxp = p.at(0);
for (auto d:p)
	if(d>maxp)
	maxp = d;
	if(c<minp)
	minp = c;
double deltap = maxp - minp;

double bmax = b*(1+0.6/deltav);
double bmin = b*(1-0.6/deltav);
double deltab = (0.6/deltav)*b;
double errmaxrelativo = deltab/b;
double devstdv = 0.6/(sqrt(3)*deltav);
double devstdp = 0.06/(sqrt(3)*deltap);
cout << "errorerelativo= " << errmaxrelativo <<endl;
cout << "devstandardvolume = " << devstdv <<endl;
cout << "devstandardpressione = " << devstdp <<endl;
}



