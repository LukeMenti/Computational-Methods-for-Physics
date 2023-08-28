#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
using namespace std;

int main()
{

    //tutte le grandezze sono in unità astronomiche

    //Terra
    double N,  x0, y0, vx0, vy0, Fx0, Fy0, Fx1, Fy1, r, deltat, ET,vt ;

    //inserisco numeri intervallo pianeti 
  cout << "inserire N numero di giri:" << endl;
    cin >> N;
    


    //Creo vettori N+1 dimensioni
    vector<double> x;
    vector<double> y;
    vector<double> vx;
    vector<double> vy;
    vector<double> et;


    //inserisco file dove salvare dati terra
    ofstream a;
    string nome;
    cout << "inserire nome dove salvare dati Terra :" << endl;
    cin >> nome;
    a.open(nome);


    // ho calcolato anche enrgia terra per dimostrare corettezza propagatore
    ofstream b;
    string nomeb;
    cout << "inserire nome dove salvare energia terra :" << endl;
    cin >> nomeb;
    b.open(nomeb);

    //inserisco dati iniziali al perielio (in unita astronomiche UA)
   
    x0 = 0.98;
    y0=0.0;
    vx0=0.0;
    vy0=2*M_PI;
  

    // imposto al primo elemento di ciascun vettore la posizione iniziale
    x.push_back(x0);
    y.push_back(y0);
   vx.push_back( vx0);
   vy.push_back(vy0) ;


    //costanti scritte in Unita Astronomiche 
    double GM =4* pow(M_PI,2);

    //massa in kg
    double m = 6E24;


    //definisco raggio
    r = sqrt(pow(x0, 2) + pow(y0, 2));

    //definisco modulo velocità terra iniziale
    vt = sqrt(pow(vx0, 2) + pow(vy0, 2));
    
    //energia condizioni iniziali
    ET = (0.5 * pow(vt, 2) * m) - (GM * m / r);

    //intervallo di tempo
    deltat = 0.001;

    //butto dentro primo elemento vettore energia
    et.push_back(vt);


    for (int i = 1; i < N; i++) {
       
        //Valore forza iniziale
        Fx0 = -GM * x0 / pow(r, 3);
        Fy0 = -GM * y0 / pow(r, 3);

        //VERLET PER POSIZIONI
        x0 += vx0 * deltat + 0.5 * pow(deltat, 2) * Fx0;
        y0 += vy0 * deltat + 0.5 * pow(deltat, 2) * Fy0;

        //Calcolo energia
        ET = (0.5 * pow(vt, 2) * m) - (GM * m / r);

        //definisco modulo velocità terra
        vt = sqrt(pow(vx0, 2) + pow(vy0, 2));
          
        //definisco raggio per nuove posizioni
        r = sqrt(pow(x0, 2) + pow(y0, 2));

        // Definisco F1 per usarlo nel propagatore velocità
        Fx1 = -GM * x0 / pow(r, 3);
        Fy1 = -GM * y0 / pow(r, 3);

        //Velocity Verlet velocità
        vx0 += 0.5 * deltat * (Fx0 + Fx1);
        vy0 += 0.5 * deltat * (Fy0 + Fy1);

        //butto ogni valore trovato nel corrispettivo vettore
        x.push_back(x0);
        y.push_back(y0);
        vx.push_back(vx0);
        vy.push_back(vy0);
        et.push_back(ET);
    }

    //Salvo nel file di testo 
    for (int p = 0; p < N; p++)
    {
      
        a << x[p] << " " << y[p] <<  endl;
        a << "" << "" << " " <<  endl;
    }




    a.close();


    //Salvo nel file di testo 
    for (int p = 0; p < N; p++)
    {

        b << p << " " << et[p] << endl;
       
    }

    b.close();


    //Marte
    double  x0M, y0M, vx0M, vy0M, Fx0M, Fy0M, Fx1M, Fy1M, rM;

    //Creo vettori N+1 dimensioni
    vector<double> xM;
    vector<double> yM;
    vector<double> vxM;
    vector<double> vyM;


    ofstream aM;
    string nomeM;
    cout << "inserire nome dove salvare dati Marte :" << endl;
    cin >> nomeM;
    aM.open(nomeM);


    //inserisco posizioni iniziali 
   
     x0M= 1.31;
     y0M=0.0;
    vx0M=0.0;
    vy0M= 1.8 *M_PI;


    // imposto al primo elemento di ciascun vettore la posizione iniziale
    xM.push_back(x0M);
    yM.push_back(y0M);
    vxM.push_back(vx0M);
    vyM.push_back(vy0M);


 
    double MassaMarte = 6E23;


    //definisco raggio
    rM = sqrt(pow(x0M, 2) + pow(y0M, 2));




    for (int i = 1; i < N; i++) {

        //Valore forza iniziale
        Fx0M = -GM * x0M / pow(rM, 3);
        Fy0M = -GM * y0M / pow(rM, 3);

        //VERLET PER POSIZIONI
        x0M += vx0M * deltat + 0.5 * pow(deltat, 2) * Fx0M;
        y0M += vy0M * deltat + 0.5 * pow(deltat, 2) * Fy0M;

        //definisco raggio per nuove posizioni
        rM = sqrt(pow(x0M, 2) + pow(y0M, 2));

        // Definisco F1 per usarlo nel propagatore velocità
        Fx1M = -GM * x0M / pow(rM, 3);
        Fy1M = -GM * y0M / pow(rM, 3);

        //Velocity Verlet velocità
        vx0M += 0.5 * deltat * (Fx0M + Fx1M);
        vy0M += 0.5 * deltat * (Fy0M + Fy1M);

        //butto ogni valore trovato nel corrispettivo vettore
        xM.push_back(x0M);
        yM.push_back(y0M);
        vxM.push_back(vx0M);
        vyM.push_back(vy0M);
    }

    //Salvo nel file di testo 
    for (int p = 0; p < N; p++)
    {
      
        aM << xM[p] << " " << yM[p] <<  endl;
        aM << "" << "" <<  endl;
    }
    aM.close();


    //Mercurio
    double  x0m, y0m, vx0m, vy0m, Fx0m, Fy0m, Fx1m, Fy1m, rm;
//Creo vettori N+1 dimensioni
    vector<double> xm;
    vector<double> ym;
    vector<double> vxm;
    vector<double> vym;


    ofstream am;
    string nomem;
    cout << "inserire nome dove salvare dati Mercurio :" << endl;
    cin >> nomem;
    am.open(nomem);


    //inserisco posizioni iniziali 
  
    x0m=0.313;
    y0m=0.0;
    vx0m=0.0;
    vy0m=3.9*M_PI;


    // imposto al primo elemento di ciascun vettore la posizione iniziale
    xm.push_back(x0m);
    ym.push_back(y0m);
    vxm.push_back(vx0m);
    vym.push_back(vy0m);


    double MassaMercurio = 3E23;
   

    //definisco raggio
    rm = sqrt(pow(x0m, 2) + pow(y0m, 2));




    for (int i = 1; i < N; i++) {

        //Valore forza iniziale
        Fx0m = -GM * x0m / pow(rm, 3);
        Fy0m = -GM * y0m / pow(rm, 3);

        //VERLET PER POSIZIONI
        x0m += vx0m * deltat + 0.5 * pow(deltat, 2) * Fx0m;
        y0m += vy0m * deltat + 0.5 * pow(deltat, 2) * Fy0m;

        //definisco raggio per nuove posizioni
        rm = sqrt(pow(x0m, 2) + pow(y0m, 2));

        // Definisco F1 per usarlo nel propagatore velocità
        Fx1m = -GM * x0m / pow(rm, 3);
        Fy1m = -GM * y0m / pow(rm, 3);

        //Velocity Verlet velocità
        vx0m += 0.5 * deltat * (Fx0m + Fx1m);
        vy0m += 0.5 * deltat * (Fy0m + Fy1m);

        //butto ogni valore trovato nel corrispettivo vettore
        xm.push_back(x0m);
        ym.push_back(y0m);
        vxm.push_back(vx0m);
        vym.push_back(vy0m);
    }

    //Salvo nel file di testo 
    for (int p = 0; p < N; p++)
    {
   
        am << xm[p] << " " << ym[p] << endl;
        am << " " << " " << endl;
    }
    am.close();


    //Venere 
    double  x0V, y0V, vx0V, vy0V, Fx0V, Fy0V, Fx1V, Fy1V, rV;

    //Creo vettori N+1 dimensioni
    vector<double> xV;
    vector<double> yV;
    vector<double> vxV;
    vector<double> vyV;


    ofstream aV;
    string nomeV;
    cout << "inserire nome dove salvare dati Venere :" << endl;
    cin >> nomeV;
    aV.open(nomeV);


    //inserisco posizioni iniziali 
  
     x0V=0.72;
    y0V=0.0;
     vx0V=0.0;
     vy0V=2.33*M_PI;


    // imposto al primo elemento di ciascun vettore la posizione iniziale
    xV.push_back(x0V);
    yV.push_back(y0V);
    vxV.push_back(vx0V);
    vyV.push_back(vy0V);



    double MassaVenere = 5E24;


    //definisco raggio
    rV = sqrt(pow(x0V, 2) + pow(y0V, 2));




    for (int i = 1; i < N; i++) {

        //Valore forza iniziale
        Fx0V = -GM * x0V / pow(rV, 3);
        Fy0V = -GM * y0V / pow(rV, 3);

        //VERLET PER POSIZIONI
        x0V += vx0V * deltat + 0.5 * pow(deltat, 2) * Fx0V;
        y0V += vy0V * deltat + 0.5 * pow(deltat, 2) * Fy0V;

        //definisco raggio per nuove posizioni
        rV = sqrt(pow(x0V, 2) + pow(y0V, 2));

        // Definisco F1 per usarlo nel propagatore velocità
        Fx1V = -GM * x0V / pow(rV, 3);
        Fy1V = -GM * y0V / pow(rV, 3);

        //Velocity Verlet velocità
        vx0V += 0.5 * deltat * (Fx0V + Fx1V);
        vy0V += 0.5 * deltat * (Fy0V + Fy1V);

        //butto ogni valore trovato nel corrispettivo vettore
        xV.push_back(x0V);
        yV.push_back(y0V);
        vxV.push_back(vx0V);
        vyV.push_back(vy0V);
    }

    //Salvo nel file di testo 
    for (int p = 0; p < N; p++)
    {

        aV << xV[p] << " " << yV[p] <<  endl;
        aV << "" << "" <<  endl;
    }
    aV.close();


    //Giove
    double  x0G, y0G, vx0G, vy0G, Fx0G, Fy0G, Fx1G, Fy1G, rG;

    //Creo vettori N+1 dimensioni
    vector<double> xG;
    vector<double> yG;
    vector<double> vxG;
    vector<double> vyG;


    ofstream aG;
    string nomeG;
    cout << "inserire nome dove salvare dati Giove :" << endl;
    cin >> nomeG;
    aG.open(nomeG);


    //inserisco posizioni iniziali 

     x0G=4.95;
     y0G=0.0;
    vx0G=0.0;
     vy0G=0.9*M_PI;


    // imposto al primo elemento di ciascun vettore la posizione iniziale
    xG.push_back(x0G);
    yG.push_back(y0G);
    vxG.push_back(vx0G);
    vyG.push_back(vy0G);


    double MassaGiove = 2E27;


    //definisco raggio
    rG = sqrt(pow(x0G, 2) + pow(y0G, 2));




    for (int i = 1; i < N; i++) {

        //Valore forza iniziale
        Fx0G = -GM * x0G / pow(rG, 3);
        Fy0G = -GM * y0G / pow(rG, 3);

        //VERLET PER POSIZIONI
        x0G += vx0G * deltat + 0.5 * pow(deltat, 2) * Fx0G;
        y0G += vy0G * deltat + 0.5 * pow(deltat, 2) * Fy0G;

        //definisco raggio per nuove posizioni
        rG = sqrt(pow(x0G, 2) + pow(y0G, 2));

        // Definisco F1 per usarlo nel propagatore velocità
        Fx1G = -GM * x0G / pow(rG, 3);
        Fy1G = -GM * y0G / pow(rG, 3);

        //Velocity Verlet velocità
        vx0G += 0.5 * deltat * (Fx0G + Fx1G);
        vy0G += 0.5 * deltat * (Fy0G + Fy1G);

        //butto ogni valore trovato nel corrispettivo vettore
        xG.push_back(x0G);
        yG.push_back(y0G);
        vxG.push_back(vx0G);
        vyG.push_back(vy0G);
    }

    //Salvo nel file di testo 
    for (int p = 0; p < N; p++)
    {
      
        aG<< xG[p] << " " << yG[p] << endl;
        aG << "" << "" <<  endl;
    }
    aG.close();

    //Saturno
    double  x0S, y0S, vx0S, vy0S, Fx0S, Fy0S, Fx1S, Fy1S, rS;

    //Creo vettori N+1 dimensioni
    vector<double> xS;
    vector<double> yS;
    vector<double> vxS;
    vector<double> vyS;


    ofstream aS;
    string nomeS;
    cout << "inserire nome dove salvare dati Saturno :" << endl;
    cin >> nomeS;
    aS.open(nomeS);


    //inserisco posizioni iniziali 
 
     x0S=9.02;
     y0S=0.0;
     vx0S=0.0;
     vy0S=0.7*M_PI;


    // imposto al primo elemento di ciascun vettore la posizione iniziale
    xS.push_back(x0S);
    yS.push_back(y0S);
    vxS.push_back(vx0S);
    vyS.push_back(vy0S);


 
    double MassaSaturno = 6E26;


    //definisco raggio
    rS = sqrt(pow(x0S, 2) + pow(y0S, 2));




    for (int i = 1; i < N; i++) {

        //Valore forza iniziale
        Fx0S = -GM * x0S / pow(rS, 3);
        Fy0S = -GM * y0S / pow(rS, 3);

        //VERLET PER POSIZIONI
        x0S += vx0S * deltat + 0.5 * pow(deltat, 2) * Fx0S;
        y0S += vy0S * deltat + 0.5 * pow(deltat, 2) * Fy0S;

        //definisco raggio per nuove posizioni
        rS = sqrt(pow(x0S, 2) + pow(y0S, 2));

        // Definisco F1 per usarlo nel propagatore velocità
        Fx1S = -GM * x0S / pow(rS, 3);
        Fy1S = -GM * y0S / pow(rS, 3);

        //Velocity Verlet velocità
        vx0S += 0.5 * deltat * (Fx0S + Fx1S);
        vy0S += 0.5 * deltat * (Fy0S + Fy1S);

        //butto ogni valore trovato nel corrispettivo vettore
        xS.push_back(x0S);
        yS.push_back(y0S);
        vxS.push_back(vx0S);
        vyS.push_back(vy0S);
    }

    //Salvo nel file di testo 
    for (int p = 0; p < N; p++)
    {
   
        aS << xS[p] << " " << yS[p] <<  endl;
        aS << "" << "" << endl;
    }
    aS.close();

    //Urano
    double  x0U, y0U, vx0U, vy0U, Fx0U, Fy0U, Fx1U, Fy1U, rU;

    //Creo vettori N+1 dimensioni
    vector<double> xU;
    vector<double> yU;
    vector<double> vxU;
    vector<double> vyU;


    ofstream aU;
    string nomeU;
    cout << "inserire nome dove salvare dati Urano :" << endl;
    cin >> nomeU;
    aU.open(nomeU);


    //inserisco posizioni iniziali 
    
     x0U=18.29;
     y0U=0.0;
     vx0U=0.0;
     vy0U=0.5*M_PI;


    // imposto al primo elemento di ciascun vettore la posizione iniziale
    xU.push_back(x0U);
    yU.push_back(y0U);
    vxU.push_back(vx0U);
    vyU.push_back(vy0U);


   
    double MassaUrano = 9E25;


    //definisco raggio
    rU = sqrt(pow(x0U, 2) + pow(y0U, 2));




    for (int i = 1; i < N; i++) {

        //Valore forza iniziale
        Fx0U = -GM * x0U / pow(rU, 3);
        Fy0U = -GM * y0U / pow(rU, 3);

        //VERLET PER POSIZIONI
        x0U += vx0U * deltat + 0.5 * pow(deltat, 2) * Fx0U;
        y0U += vy0U * deltat + 0.5 * pow(deltat, 2) * Fy0U;

        //definisco raggio per nuove posizioni
        rU = sqrt(pow(x0U, 2) + pow(y0U, 2));

        // Definisco F1 per usarlo nel propagatore velocità
        Fx1U = -GM * x0U / pow(rU, 3);
        Fy1U = -GM * y0U / pow(rU, 3);

        //Velocity Verlet velocità
        vx0U += 0.5 * deltat * (Fx0U + Fx1U);
        vy0U += 0.5 * deltat * (Fy0U + Fy1U);

        //butto ogni valore trovato nel corrispettivo vettore
        xU.push_back(x0U);
        yU.push_back(y0U);
        vxU.push_back(vx0U);
        vyU.push_back(vy0U);
    }

    //Salvo nel file di testo 
    for (int p = 0; p < N; p++)
    {
     
        aU << xU[p] << " " << yU[p] <<  endl;
        aU << "" << "" << endl;
    }
    aU.close();

    //Nettuno
    double  x0N, y0N, vx0N, vy0N, Fx0N, Fy0N, Fx1N, Fy1N, rN;

    //Creo vettori N+1 dimensioni
    vector<double> xN;
    vector<double> yN;
    vector<double> vxN;
    vector<double> vyN;


    ofstream aN;
    string nomeN;
    cout << "inserire nome dove salvare dati Nettuno :" << endl;
    cin >> nomeN;
    aN.open(nomeN);


    //inserisco posizioni iniziali 

     x0N=29.8;
     y0N=0.0;
     vx0N=0.0;
     vy0N=0.4*M_PI;


    // imposto al primo elemento di ciascun vettore la posizione iniziale
    xN.push_back(x0N);
    yN.push_back(y0N);
    vxN.push_back(vx0N);
    vyN.push_back(vy0N);


   
    double MassaNettuno = 1E26;


    //definisco raggio
    rN= sqrt(pow(x0N, 2) + pow(y0N, 2));




    for (int i = 1; i < N; i++) {

        //Valore forza iniziale
        Fx0N = -GM * x0N / pow(rN, 3);
        Fy0N = -GM * y0N / pow(rN, 3);

        //VERLET PER POSIZIONI
        x0N += vx0N * deltat + 0.5 * pow(deltat, 2) * Fx0N;
        y0N += vy0N * deltat + 0.5 * pow(deltat, 2) * Fy0N;

        //definisco raggio per nuove posizioni
        rN = sqrt(pow(x0N, 2) + pow(y0N, 2));

        // Definisco F1 per usarlo nel propagatore velocità
        Fx1N = -GM * x0N / pow(rN, 3);
        Fy1N = -GM * y0N / pow(rN, 3);

        //Velocity Verlet velocità
        vx0N += 0.5 * deltat * (Fx0N + Fx1N);
        vy0N += 0.5 * deltat * (Fy0N + Fy1N);

        //butto ogni valore trovato nel corrispettivo vettore
        xN.push_back(x0N);
        yN.push_back(y0N);
        vxN.push_back(vx0N);
        vyN.push_back(vy0N);
    }

    //Salvo nel file di testo 
    for (int p = 0; p < N; p++)
    {
  
        aN << xN[p] << " " << yN[p] <<  endl;
        aN << "" << "" << endl;
    }
    aN.close();


     //file con 0 che sarebbe Sole
    ofstream asole;
    string nomesole;
    cout << "inserire nome dove salvare Sole:" << endl;
    cin >> nomesole;
    asole.open(nomesole);

    //Salvo nel file di testo 
    for (int p = 0; p < N; p++)
    {

        asole << 0 << " " << 0 << endl;
        asole << "" << "" << endl;
    }
    asole.close();
    return 0;
}