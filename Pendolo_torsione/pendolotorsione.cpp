#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
using namespace std;

struct ParabolaFitResult {
    double a, b, c;
};
struct Punto {
    double x;
    double y;
};
Punto trovaVertice(ParabolaFitResult& parabola);
double sogliamin(vector<double> x, vector<double> y);
double sogliamax(vector<double> x, vector<double> y);
ParabolaFitResult fitParabola(vector<double> x, vector<double> y);
double media(vector<double> v);
double contadatipos(vector<double> p);

int main() {





string nomeFile;
vector<double> t, f, p, a, fa;
double tval, fval,pval,aval,faval; 
ParabolaFitResult par;
ofstream output("fitmassimi.txt");


    // cout<<"Inserisci il nome del file di una precisa frequenza impostata della forzante : ";
    // cin>>nomeFile;
    ifstream fin("f0.980.txt");
    if (!fin) {
        cout << "Errore nell'apertura del file in lettura.\n";
        return 1;
    }

    while(fin>>tval>>fval>>pval>>aval>>faval){
	    t.push_back(tval);    		
        f.push_back(fval);
    	p.push_back(pval);
        a.push_back(aval);
        fa.push_back(faval);
    }  
    double stimam;
    vector<double> stimamax;
    double stimax;
    vector<double> stimax1;
    //calcolo fit 
    vector<double> datifity,datifitx;
    vector<ParabolaFitResult> parabole;
    ofstream pic("picchicandidati.txt");
    // double soglia=sogliamax(t,p);
    // cout<<soglia<<endl;
    double soglia;

    for(int i=0; i<p.size();i++){
            soglia=p.at((contadatipos(p)/2)-1);
            cout<<"La soglia in questo particolare intervallo: "<<soglia<<endl;
        if(p.at(i)>=soglia){
           for (int j = i; j < p.size() && p.at(j) >= soglia; j++) {
                datifity.push_back(p.at(j));
                datifitx.push_back(t.at(j));
            }
            
            if (datifitx.size() >= 3) {
                // parabole.push_back(fitParabola(datifitx, datifity));
                stimax=media(datifitx)/datifitx.size();
                stimax1.push_back(stimax);
                stimam=media(datifity)/datifity.size();
                stimamax.push_back(stimam);

            }           
            i=i+datifity.size()-1;
        }    
        datifitx.clear();
        datifity.clear();
    }

    //stampa risultati
    cout<<"La soglia min: "<<sogliamin(t,p)<<" "<<endl;
    cout<<"La soglia max: "<<sogliamax(t,p)<<" "<<endl;
    
    //stampa parabole fittate, e inserisce in un file tutti i massimi trovati
    // for(auto a:parabole){
    //     Punto pt=trovaVertice(a);
    //     cout<<"Delle parabole "<<a.a<<" "<<a.b<<" "<<a.c<<"il vertice: "<<pt.x<<" l'ordinata: "<<pt.y<<endl;
    //     output<<pt.x<<" "<<pt.y<<endl; 
        
    // }

    for(int i=0; i<stimax1.size(); i++){
        pic<<stimax1.at(i)<<" "<<stimamax.at(i)<<endl;
        cout<<""<<stimax1.at(i)<<" "<<stimamax.at(i);
    }
    //sistemo considerando non covarianti
    vector<double> a1;
    for(int k=0; k<stimax1.size();k++){
        a1.push_back(k);
        a1.push_back(k+1);
        k=k+2;
    }

    //Stima media massimo
    double mediamassimo=media(a1);

    //stampa risultati


    return 0;
}

ParabolaFitResult fitParabola(vector<double> x, vector<double> y) {
    size_t n = x.size();
   
    double Sx = 0, Sx2 = 0, Sx3 = 0, Sx4 = 0;
    double Sy = 0, Sxy = 0, Sx2y = 0;

    for (size_t i = 0; i < n; ++i) {
        double xi = x[i];
        double yi = y[i];
        double xi2 = xi * xi;

        Sx   += xi;
        Sx2  += xi2;
        Sx3  += xi2 * xi;
        Sx4  += xi2 * xi2;
        Sy   += yi;
        Sxy  += xi * yi;
        Sx2y += xi2 * yi;
    }

    // Sistema lineare 3x3 da risolvere
    // | Sx4  Sx3  Sx2 |   | a |   = | Sx2y |
    // | Sx3  Sx2  Sx  | * | b |   = | Sxy  |
    // | Sx2  Sx   n   |   | c |   = | Sy   |

    double D = Sx4*(Sx2*n - Sx*Sx) - Sx3*(Sx3*n - Sx*Sx2) + Sx2*(Sx3*Sx - Sx2*Sx2);
    double Da = Sx2y*(Sx2*n - Sx*Sx) - Sxy*(Sx3*n - Sx*Sx2) + Sy*(Sx3*Sx - Sx2*Sx2);
    double Db = Sx4*(Sxy*n - Sx*Sy) - Sx3*(Sx2y*n - Sx*Sy) + Sx2*(Sx2y*Sx - Sxy*Sx2);
    double Dc = Sx4*(Sx2*Sy - Sx*Sx2y) - Sx3*(Sx3*Sy - Sx*Sxy) + Sx2*(Sx3*Sxy - Sx2*Sx2y);

    

    ParabolaFitResult result;
    result.a = Da / D;
    result.b = Db / D;
    result.c = Dc / D;
    return result;
}

Punto trovaVertice(ParabolaFitResult& parabola) {
    double x_vertice = -parabola.b / (2.0 * parabola.a);
    double y_vertice = parabola.a * x_vertice * x_vertice + parabola.b * x_vertice + parabola.c;

    return {x_vertice, y_vertice};
}
 
double sogliamax(vector<double> x, vector<double> y){
    double max=y.at(0);
    for(int i=0; i<y.size();i++){
        if(max<y.at(i)){
            max=y.at(i);
        }
    }

    return max*0.70;
}
double sogliamin(vector<double> x, vector<double> y){

     double min=y.at(0);
    for(int i=0; i<y.size(); i++){
        if(min>y.at(i)){
            min=y.at(i);
        }
    }
    return min*0.70;
}
double media(vector<double> v){
    double somma=0.0;
    for(int i=0;i<v.size();i++){
        somma+=v.at(i);
    }

    return somma/v.size();
}
//attenzione bisogna saltare il primo 
double contadatipos(vector<double> p){
    double i=0;
    if(p.at(0)>0){
         while(p.at(i+1)>=0){
            i++;
        }
    }else{

    }
    return i;
}