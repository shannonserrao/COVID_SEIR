/* requirements for the containment of COVID-19 
 * compilation: g++ -std=c++11 SEAIRDM -O3 -o seairdm 
 * usage: ./seairdm L_i TT_i TTstart_i TestP_i DelayTest_i QuarantineP_i DelayQT_i d_f phi_f pop_f I0_f r_f a_f b_f outtype_i step_i Quarantine_i af_f sid_i 
 * with _i for integer numbers and _f for float numbers 
 */
#include <vector>
#include <iostream>
#include <string>
#include <cmath>
#include <ctime>
#include <random>
#include <iterator>
#include <algorithm>

using namespace std;

/* count adjacent neighbors and times it by the infection rate r */
double Aiv(vector<int> &ss, vector<int> &adi, double r) {
    int tot = 0;
    for (int j=0; j<4; j++) 
        tot += ((ss[adi[j]]==2)||(ss[adi[j]]==4));
    return r*tot;
}

/* parameter descriptions.
 * L: lateral size;
 * TT: simulation time;
 * TTstart: mitigation starting time;
 * TestP: testing period;
 * DelayTest: delay in time for getting test results (I->T);
 * QuarantineP: quarantine period;
 * DelayQT: delay time for quarantine;
 * d: diffusion rate;
 * phi: probability for forming a long-range link;
 * pop: population density;
 * I0: fraction of initial infectious population; pop*L*L*I0 gives the number of initial infected individuals;
 * r: infection rate;
 * a: recovery rate;
 * b: incubation rate;
 * outtype: control output data to state vectors (==0) or to the numbers of each species vs. time;
 * step: recording results every step steps;
 * Quarantine: nearest nbors without hopping=1; nearest nbors with hopping=2; long-range nbors included without hopping=3; long-range nbors included with hopping, no quarantine=0;
 * af: fraction of unidentifiable (asymptomatic) population
 * sid: seed id to avoid the same random number seed when multiple runs are invoked at the same time;
 * (A as unidentified individuals, I as identified individuals).
 */
void sir(int L, long TT, int TTstart, int TestP, int DelayTest, int QuarantineP, int DelayQT, double d, double phi, double pop, double I0, double r, double a, double b, int outtype, int step, int Quarantine, double af, int sid) 
{
    long N = L * L; // lattice size
    int left, right; // left and right of a link
    /* variables */
    int output[10]; // {T:0, NS:1, NI:2, NR:3, NA:4, NT:5, NQS:6, NQA:7, NE:8, NQE:9, NQI: 10}
    vector<int> S(N,0); // current states. S:1, I:2, R:3, A:4, T:5, QS:6, QA:7, E:8, QE:9, QI: 10.
    vector<int> RT(N,0); // recording releasing time for quarantined individuals.
    vector<int> DT(N,0); // recording delay time for quarantined individuals.
    vector<int> IDT(N,0); // recording delay time for tested individuals to get identified.
    // s: previous state; 
    // ini: initial infection indices
    vector<int> s, ini;
    vector<vector<int>> ad; // adjacent list for each node
    vector<int> nb; // neighbors, to be pushed back to ad

    int nbors=4; // number of nearest neighbors

    /* counts for each species */
    int NS, NI, NR, NA, NT, NQS, NQA, NE, NQE, NQI;
    NS=0; NI=0; NR=0; NA=0; NT=0; NQS=0; NQA=0; NE=0; NQE=0; NQI=0;

    // random number generators
    default_random_engine e;
    e.seed(time(NULL)+5*sid);
    uniform_int_distribution<unsigned> ui(0, N-1);
    uniform_real_distribution<double> ur(0, 1);

    // initial infection (set I0=1/(N*pop) if we only want one initial infection)
    int NI0 = ceil(N*pop*I0);
    double ps = pop * (1-I0);

    for(int i = 0; i < NI0; ++i) {
        int ri;
        do {
            ri = ui(e);
        }while( find(ini.begin(), ini.end(), ri)!=ini.end() );
        ini.push_back(ri);
    }

    // initially identifiables and unidentifiables split into 50-50
    for(auto i: ini) {
        if (ur(e)<0.5) {
            S[i] = 2; NI++;
        } else {
            S[i] = 4; NA++;
	  }
    }
    
    // initialization for susceptibles
    for(int i = 0; i < N; i++) {
	if (ur(e)<ps && S[i]!=2 && S[i]!=4) {
	    S[i] = 1; NS++;
	}
    }
    s = S;
    
    // constructing the 2d lattice
    for (int i=0; i < L; i++) {
        for (int j=0; j < L; j++) {
	   ad.push_back({((i-1+L)%L)*L+j, ((i+1)%L)*L+j, i*L+((j-1+L)%L), i*L+((j+1)%L)}); // up, down, left, right
	 }
    }
    // creating long-range links
    int NL = round(2*phi*N);
    for (int i=0; i < NL; i++) {
        left = ui(e);
	right = ui(e);
	ad[left].push_back(right);
	ad[right].push_back(left);
    }
    
    // main Monte Carlo loop
    for (int T = 0; T < TT; T++)
    {
      if (outtype == 1 && T%step==0) {// output state vectors; TT is usually short
            for (long i=0; i<N; i++)
                    cout << S[i] << " ";
            cout << endl;
       }
       else if (outtype == 2 && T%step==0) {// output NS, total NI, NR, NT and total NQ; TT can be long
            /* out put */
            cout << T << " " << NS << " " << NI+NA+NT+NQA+NQI << " " << NR << " " << NT<< " " << NQA+NQE+NQI+NQS << endl;
       }


        /* testing and quarantine */
        if (T>=TTstart && (T-TTstart)%TestP==0) {
            for (int i=0; i<N; i++) {
                if (s[i]==2 && IDT[i]==0) {
                    IDT[i]=T+DelayTest;
				
                if (1<=Quarantine && Quarantine<=4) {
                        if (Quarantine==1 || Quarantine==2) nbors=4;
                        else if (Quarantine==3 || Quarantine==4) nbors=ad[i].size();
			    		for (int j=0; j<nbors; j++) {
			    			if (s[ad[i][j]]==1) {DT[ad[i][j]] = T+DelayQT;} // S->Delayed S
			    			else if (s[ad[i][j]]==4) {DT[ad[i][j]] = T+DelayQT;} // A->Delayed A
			    			else if (s[ad[i][j]]==8) {DT[ad[i][j]] = T+DelayQT;} // E->Delayed E
			    		}
                }
			}

			if (s[i]==10 && IDT[i]==0) IDT[i]=T+DelayTest;
            }
        }

       

        for (int i = 0; i < N; i++)
        {
            /* Qurantine */
            if (DT[i]>0 && T>=DT[i]) {
            if (s[i]==1) {S[i]=6; NS--; NQS++; DT[i]=0; RT[i] = T+QuarantineP;} // Delayed S->QS
            else if (s[i]==2) {S[i]=10; NI--; NQI++; DT[i]=0; RT[i]=T+QuarantineP;}//Delayed E->Delayed I->QI
            else if (s[i]==4) {S[i]=7; NA--; NQA++; DT[i]=0; RT[i] = T+QuarantineP;} // Delayed A->QA
            else if (s[i]==8) {S[i]=9; NE--; NQE++; DT[i]=0; RT[i] = T+QuarantineP;} // Delayed E->QE
            s[i]=S[i];
            }
            /* I/QT to T */
            if (IDT[i]>0 && T>=IDT[i] && (s[i]==2||s[i]==10)) {
                S[i] = 5; NT++;
                if(s[i]==2) NI--;
                else if (s[i]==10) NQI--;
                IDT[i]=0;
                RT[i]=1; //Qurantine T guys
                s[i]=S[i];//I/QI->T
            }

            /* releasing after quarantine */
            if (1<=Quarantine && Quarantine<=4) {
                if (RT[i]>0 && T>=RT[i]) {
                    if (s[i]==6) {S[i]=1; RT[i]=0; NQS--; NS++;}
                    else if (s[i]==7) {S[i]=4; RT[i]=0; NQA--; NA++;}
                    else if (s[i]==9) {S[i]=8; RT[i]=0; NQE--; NE++;}

                s[i]=S[i];
                }
            }

            if ( (s[i]==1) && ( ur(e)<=Aiv(s, ad[i], r) ) ) { //exposed
            S[i] = 8; NS--; NE++; //S->E               
            }
            else if ( (s[i]==8) && (ur(e)<b) ) {
                if (ur(e)<af) {S[i] = 4; NE--; NA++;} // E->A 
                    else {S[i] = 2; NE--; NI++;}  //E->I (not I for the moment)
            }
            else if ( (s[i]==9) && (ur(e)<b) ) {
                if (ur(e)<af) {S[i] = 7; NQE--; NQA++;} // QE->QA 
                else {S[i] = 10; NQE--; NQI++;}  //QE->QI
            }
            else if ( (s[i]==2||s[i]==4||s[i]==5||s[i]==7||s[i]==10) && ( ur(e)<=a ) ) { //recovery
                S[i] = 3; NR++; // I/T/A/QA/QI->R with rate a
                if (s[i]==2) NI--;
                else if (s[i]==4) NA--;
                else if (s[i]==5) { NT--; RT[i]=0; }
                else if (s[i]==7) { NQA--; RT[i]=0; }
                else if (s[i]==10) { NQI--; RT[i]=0; }
            }
        
            if(S[i] && ur(e)<d) { //diffusion, no hopping if quarantined (RT[i]>0)
                int j = ad[i][rand() % ad[i].size()];
                if (Quarantine==0 || ((Quarantine==1 || Quarantine==3) && RT[i]==0) ){
                    if (S[j]==0)	
                    {S[j] = S[i]; S[i] = 0; DT[j]=DT[i]; DT[i]=0; IDT[j]=IDT[i]; IDT[i]=0;}
                }
                else if (Quarantine==2 || Quarantine==4){
                    if (S[j]==0)	
                    {S[j] = S[i]; S[i] = 0; RT[j]=RT[i]; RT[i]=0; DT[j]=DT[i]; DT[i]=0; IDT[j]=IDT[i]; IDT[i]=0;}
                }
            }
        }
        s=S;
    }
    return;
}

int main(int argc, char * argv[])
{
    int L=atoi(argv[1]); long TT=atoi(argv[2]); double TTstart=strtod(argv[3], NULL); int TestP=atoi(argv[4]); 
    int DelayTest=atoi(argv[5]); int QuarantineP=atoi(argv[6]); int DelayQT=atoi(argv[7]);
    double d=strtod(argv[8], NULL); double phi=strtod(argv[9], NULL); double pop=strtod(argv[10], NULL); 
    double I0=strtod(argv[11], NULL); double r=strtod(argv[12], NULL); double a=strtod(argv[13], NULL); 
    double b=strtod(argv[14], NULL); int outtype=atoi(argv[15]); int step=atoi(argv[16]); 
    int Quarantine=atoi(argv[17]); double af=strtod(argv[18], NULL); int sid=atoi(argv[19]);

    sir(L/*1*/, TT/*2*/, TTstart/*3*/, TestP/*4*/, DelayTest/*5*/, QuarantineP/*6*/, DelayQT/*7*/, d/*8*/, phi/*9*/, pop/*10*/, I0/*11*/, r/*12*/, a/*13*/, b/*14*/, outtype/*15*/, step/*16*/, Quarantine/*17*/, af/*18*/, sid/*19*/);
    return 0;
}
