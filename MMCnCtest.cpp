#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/timeb.h>
#include <iostream>
#include <pthread.h>
#include <math.h>
#include <immintrin.h>
#include <chrono>
#include <ctime>
#include <mutex>
#include <cnc/cnc.h>
#include <cnc/debug.h>

using namespace std;
#ifndef TYPE
#define TYPE int
#endif

#ifndef ALIGNMENT
#define ALIGNMENT 64
#endif

struct A_context;


struct A_Astep
{
    // declaration of execute method goes here
    int execute( const int & tag, A_context & c ) const;
};

struct A1_Astep
{
    // declaration of execute method goes here
    int execute( const int & tag, A_context & c ) const;
};

struct B1_Astep
{
    // declaration of execute method goes here
    int execute( const int & tag, A_context & c ) const;
};

struct C1_Astep
{
    // declaration of execute method goes here
    int execute( const int & tag, A_context & c ) const;
};

struct D1_Astep
{
    // declaration of execute method goes here
    int execute( const int & tag, A_context & c ) const;
};
struct A2_Astep
{
    // declaration of execute method goes here
    int execute( const int & tag, A_context & c ) const;
};

struct B2_Astep
{
    // declaration of execute method goes here
    int execute( const int & tag, A_context & c ) const;
};

struct C2_Astep
{
    // declaration of execute method goes here
    int execute( const int & tag, A_context & c ) const;
};

struct D2_Astep
{
    // declaration of execute method goes here
    int execute( const int & tag, A_context & c ) const;
};

struct A_context : public CnC::context< A_context >
{
    TYPE* Z; int zi; int zj; TYPE* X; int xi; int xj; TYPE* Y; int yi; int yj; int size; int N; int B;

    // step collections
    CnC::step_collection< A_Astep > m_Asteps;
    CnC::step_collection< A1_Astep > m_A1steps;
    CnC::step_collection< B1_Astep > m_B1steps;
    CnC::step_collection< C1_Astep > m_C1steps;
    CnC::step_collection< D1_Astep > m_D1steps;
    CnC::step_collection< A2_Astep > m_A2steps;
    CnC::step_collection< B2_Astep > m_B2steps;
    CnC::step_collection< C2_Astep > m_C2steps;
    CnC::step_collection< D2_Astep > m_D2steps;

    // item collection
    CnC::item_collection< int, int >  m_ABCD1syncitem;
    //CnC::item_collection< int, int >  m_BC2syncitem;

    // Tag collections
    CnC::tag_collection< int > m_Atags;
    CnC::tag_collection< int > m_A1tags;
    CnC::tag_collection< int > m_B1tags;
    CnC::tag_collection< int > m_C1tags;
    CnC::tag_collection< int > m_D1tags;
    CnC::tag_collection< int > m_A2tags;
    CnC::tag_collection< int > m_B2tags;
    CnC::tag_collection< int > m_C2tags;
    CnC::tag_collection< int > m_D2tags;

    // The context class constructor
    A_context()
            : CnC::context< A_context >(),
            // Initialize each step collection
              m_Asteps( *this ),
              m_A1steps( *this ),
              m_B1steps( *this ),
              m_C1steps( *this ),
              m_D1steps( *this ),
              m_A2steps( *this ),
              m_B2steps( *this ),
              m_C2steps( *this ),
              m_D2steps( *this ),

              m_ABCD1syncitem(*this),
              //m_BC2syncitem(*this),
            // Initialize each tag collection
              m_Atags( *this ),
              m_A1tags( *this ),
              m_B1tags( *this ),
              m_C1tags( *this ),
              m_D1tags( *this ),
              m_A2tags( *this ),
              m_B2tags( *this ),
              m_C2tags( *this ),
              m_D2tags( *this )
    {
        //this.X = x;
        //this.size = n;
        // Prescriptive relations
        //X = x;
        //size = n;
        m_Atags.prescribes( m_Asteps, *this );
        m_A1tags.prescribes( m_A1steps, *this );
        m_B1tags.prescribes( m_B1steps, *this );
        m_C1tags.prescribes( m_C1steps, *this );
        m_D1tags.prescribes( m_D1steps, *this );
        m_A2tags.prescribes( m_A2steps, *this );
        m_B2tags.prescribes( m_B2steps, *this );
        m_C2tags.prescribes( m_C2steps, *this );
        m_D2tags.prescribes( m_D2steps, *this );


        m_A1steps.produces( m_ABCD1syncitem );
        m_B1steps.produces( m_ABCD1syncitem );
        m_C1steps.produces( m_ABCD1syncitem );
        m_D1steps.produces( m_ABCD1syncitem );
        m_A2steps.consumes( m_ABCD1syncitem );
        m_B2steps.consumes( m_ABCD1syncitem );
        m_C2steps.consumes( m_ABCD1syncitem );
        m_D2steps.consumes( m_ABCD1syncitem );

    }
};


int A_Astep::execute( const int & tag, A_context & Actxt ) const
{
    TYPE * Z = Actxt.Z;
    TYPE * X = Actxt.X;
    TYPE * Y = Actxt.Y;
    int zi = Actxt.zi;
    int zj = Actxt.zj;
    int xi = Actxt.xi;
    int xj = Actxt.xj;
    int yi = Actxt.yi;
    int yj = Actxt.yj;
    int B = Actxt.B;
    int N = Actxt.N;
    int size = Actxt.size;

    if(size <= B) {
		for(int i =0; i<size; i++) {
                for(int k =0; k<size;k++) {
                        for(int j =0; j<size; j++) {
                                Z[(zi+i)*N+ zj+j] = Z[(zi+i)*N+ zj+j] + ( X[(xi+i)*N+ xj+k] * Y[(yi+k)*N+yj+j] );
                        }
                }
        }

	}
    else {
        Actxt.m_A1tags.put(0);
        Actxt.m_B1tags.put(0);
        Actxt.m_C1tags.put(0);
        Actxt.m_D1tags.put(0);
    }
    return 0;
}


int A1_Astep::execute( const int & tag, A_context & Actxt ) const {
    A_context actxt;
    actxt.Z = Actxt.Z;
    actxt.zi = Actxt.zi;
    actxt.zj = Actxt.zj;
    actxt.X = Actxt.X;
    actxt.xi = Actxt.xi;
    actxt.xj = Actxt.xj;
    actxt.Y = Actxt.Y;
    actxt.yi = Actxt.yi;
    actxt.yj = Actxt.yj;
    actxt.size = Actxt.size/2;
    actxt.N = Actxt.N;
    actxt.B = Actxt.B;

//ParallelRecursiveMM(Z,zi,zj, X, xi,xj, Y,yi,yj, n, N, B);
    actxt.m_Atags.put(0);
    actxt.wait();

    Actxt.m_ABCD1syncitem.put(0,1);
    Actxt.m_A2tags.put(0);
    Actxt.m_B2tags.put(0);
    Actxt.m_C2tags.put(0);
    Actxt.m_D2tags.put(0);

    return 0;
}


int B1_Astep::execute( const int & tag, A_context & Actxt ) const {
	int n = Actxt.size/2;
    A_context actxt;
    actxt.Z = Actxt.Z;
    actxt.zi = Actxt.zi;
    actxt.zj = Actxt.zj+n;
    actxt.X = Actxt.X;
    actxt.xi = Actxt.xi;
    actxt.xj = Actxt.xj;
    actxt.Y = Actxt.Y;
    actxt.yi = Actxt.yi;
    actxt.yj = Actxt.yj+n;
    actxt.size = n;
    actxt.N = Actxt.N;
    actxt.B = Actxt.B;
//ParallelRecursiveMM(Z,zi,zj+n, X, xi,xj, Y,yi,yj+n, n, N, B);
    actxt.m_Atags.put(0);
    actxt.wait();

    Actxt.m_ABCD1syncitem.put(1,1);
    Actxt.m_A2tags.put(0);
    Actxt.m_B2tags.put(0);
    Actxt.m_C2tags.put(0);
    Actxt.m_D2tags.put(0);

    return 0;
}

int C1_Astep::execute( const int & tag, A_context & Actxt ) const {
	int n = Actxt.size/2;
    A_context actxt;
    actxt.Z = Actxt.Z;
    actxt.zi = Actxt.zi+n;
    actxt.zj = Actxt.zj;
    actxt.X = Actxt.X;
    actxt.xi = Actxt.xi+n;
    actxt.xj = Actxt.xj;
    actxt.Y = Actxt.Y;
    actxt.yi = Actxt.yi;
    actxt.yj = Actxt.yj;
    actxt.size = n;
    actxt.N = Actxt.N;
    actxt.B = Actxt.B;

//ParallelRecursiveMM(Z,zi+n,zj, X, xi+n,xj, Y,yi,yj, n, N, B);
    actxt.m_Atags.put(0);
    actxt.wait();

    Actxt.m_ABCD1syncitem.put(2,1);
    Actxt.m_A2tags.put(0);
    Actxt.m_B2tags.put(0);
    Actxt.m_C2tags.put(0);
    Actxt.m_D2tags.put(0);

    return 0;
}

int D1_Astep::execute( const int & tag, A_context & Actxt ) const {
	int n = Actxt.size/2;
    A_context actxt;
    actxt.Z = Actxt.Z;
    actxt.zi = Actxt.zi+n;
    actxt.zj = Actxt.zj+n;
    actxt.X = Actxt.X;
    actxt.xi = Actxt.xi+n;
    actxt.xj = Actxt.xj;
    actxt.Y = Actxt.Y;
    actxt.yi = Actxt.yi;
    actxt.yj = Actxt.yj+n;
    actxt.size = n;
    actxt.N = Actxt.N;
    actxt.B = Actxt.B;

//ParallelRecursiveMM(Z,zi+n,zj+n, X, xi+n,xj, Y,yi,yj+n, n, N, B);
    actxt.m_Atags.put(0);
    actxt.wait();

    Actxt.m_ABCD1syncitem.put(3,1);
    Actxt.m_A2tags.put(0);
    Actxt.m_B2tags.put(0);
    Actxt.m_C2tags.put(0);
    Actxt.m_D2tags.put(0);

    return 0;
}


int A2_Astep::execute( const int & tag, A_context & Actxt ) const {

	int a,b,c,d;
    Actxt.m_ABCD1syncitem.get(0,a);
    Actxt.m_ABCD1syncitem.get(1,b);
    Actxt.m_ABCD1syncitem.get(2,c);
    Actxt.m_ABCD1syncitem.get(3,d);

    if(a==1 && b==1 && c==1 && d==1){
    	int n = Actxt.size/2;
    	A_context actxt;
    	actxt.Z = Actxt.Z;
    	actxt.zi = Actxt.zi;
    	actxt.zj = Actxt.zj;
    	actxt.X = Actxt.X;
    	actxt.xi = Actxt.xi;
    	actxt.xj = Actxt.xj+n;
    	actxt.Y = Actxt.Y;
    	actxt.yi = Actxt.yi+n;
    	actxt.yj = Actxt.yj;
    	actxt.size = n;
    	actxt.N = Actxt.N;
   		actxt.B = Actxt.B;

//ParallelRecursiveMM(Z,zi,zj, X, xi,xj+n, Y,yi+n,yj, n, N, B);
    	actxt.m_Atags.put(0);
    	actxt.wait();
    }
    return 0;
}


int B2_Astep::execute( const int & tag, A_context & Actxt ) const {

	int a,b,c,d;
    Actxt.m_ABCD1syncitem.get(0,a);
    Actxt.m_ABCD1syncitem.get(1,b);
    Actxt.m_ABCD1syncitem.get(2,c);
    Actxt.m_ABCD1syncitem.get(3,d);

    if(a==1 && b==1 && c==1 && d==1){
    	int n = Actxt.size/2;
    	A_context actxt;
    	actxt.Z = Actxt.Z;
    	actxt.zi = Actxt.zi;
    	actxt.zj = Actxt.zj+n;
    	actxt.X = Actxt.X;
    	actxt.xi = Actxt.xi;
    	actxt.xj = Actxt.xj+n;
    	actxt.Y = Actxt.Y;
    	actxt.yi = Actxt.yi+n;
    	actxt.yj = Actxt.yj+n;
    	actxt.size = n;
    	actxt.N = Actxt.N;
   		actxt.B = Actxt.B;

//ParallelRecursiveMM(Z,zi,zj+n, X, xi,xj+n, Y,yi+n,yj+n, n, N, B);
    	actxt.m_Atags.put(0);
    	actxt.wait();
    }
    return 0;
}

int C2_Astep::execute( const int & tag, A_context & Actxt ) const {

	int a,b,c,d;
    Actxt.m_ABCD1syncitem.get(0,a);
    Actxt.m_ABCD1syncitem.get(1,b);
    Actxt.m_ABCD1syncitem.get(2,c);
    Actxt.m_ABCD1syncitem.get(3,d);

    if(a==1 && b==1 && c==1 && d==1){
    	int n = Actxt.size/2;
    	A_context actxt;
    	actxt.Z = Actxt.Z;
    	actxt.zi = Actxt.zi+n;
    	actxt.zj = Actxt.zj;
    	actxt.X = Actxt.X;
    	actxt.xi = Actxt.xi+n;
    	actxt.xj = Actxt.xj+n;
    	actxt.Y = Actxt.Y;
    	actxt.yi = Actxt.yi+n;
    	actxt.yj = Actxt.yj;
    	actxt.size = n;
    	actxt.N = Actxt.N;
   		actxt.B = Actxt.B;

//ParallelRecursiveMM(Z,zi+n,zj, X, xi+n,xj+n, Y,yi+n,yj, n, N, B);
    	actxt.m_Atags.put(0);
    	actxt.wait();
    }
    return 0;
}

int D2_Astep::execute( const int & tag, A_context & Actxt ) const {

	int a,b,c,d;
    Actxt.m_ABCD1syncitem.get(0,a);
    Actxt.m_ABCD1syncitem.get(1,b);
    Actxt.m_ABCD1syncitem.get(2,c);
    Actxt.m_ABCD1syncitem.get(3,d);

    if(a==1 && b==1 && c==1 && d==1){
    	int n = Actxt.size/2;
    	A_context actxt;
    	actxt.Z = Actxt.Z;
    	actxt.zi = Actxt.zi+n;
    	actxt.zj = Actxt.zj+n;
    	actxt.X = Actxt.X;
    	actxt.xi = Actxt.xi+n;
    	actxt.xj = Actxt.xj+n;
    	actxt.Y = Actxt.Y;
    	actxt.yi = Actxt.yi+n;
    	actxt.yj = Actxt.yj+n;
    	actxt.size = n;
    	actxt.N = Actxt.N;
   		actxt.B = Actxt.B;

//ParallelRecursiveMM(Z,zi+n,zj+n, X, xi+n,xj+n, Y,yi+n,yj+n, n, N, B);
    	actxt.m_Atags.put(0);
    	actxt.wait();
    }
    return 0;
}

void IterMMBase(TYPE* Z, TYPE* X, TYPE* Y, int N) {

	for(int i =0; i<N; i++) {
                for(int k =0; k<N; k++) {
                        for(int j =0; j<N; j++) {
                                Z[(i)*N+ j] = Z[(i)*N+ j] + ( X[(i)*N+ k] * Y[(k)*N+ j] );
                        }
                }
        }
}
/*
void ParallelRecursiveMM(TYPE* Z, int zi, int zj, TYPE* X, int xi, int xj, TYPE* Y, int yi, int yj, int size, int N, int B) {
	
	if(size <= B) {

		for(int i =0; i<size; i++) {
                for(int k =0; k<size;k++) {
                        for(int j =0; j<size; j++) {
                                Z[(zi+i)*N+ zj+j] = Z[(zi+i)*N+ zj+j] + ( X[(xi+i)*N+ xj+k] * Y[(yi+k)*N+yj+j] );
                        }
                }
        }

	}

	else {

		int n = size/2;

		ParallelRecursiveMM(Z,zi,zj, X, xi,xj, Y,yi,yj, n, N, B);
		ParallelRecursiveMM(Z,zi,zj+n, X, xi,xj, Y,yi,yj+n, n, N, B);
		ParallelRecursiveMM(Z,zi+n,zj, X, xi+n,xj, Y,yi,yj, n, N, B);
		ParallelRecursiveMM(Z,zi+n,zj+n, X, xi+n,xj, Y,yi,yj+n, n, N, B);
		//cilk_sync;

		ParallelRecursiveMM(Z,zi,zj, X, xi,xj+n, Y,yi+n,yj, n, N, B);
		ParallelRecursiveMM(Z,zi,zj+n, X, xi,xj+n, Y,yi+n,yj+n, n, N, B);
		ParallelRecursiveMM(Z,zi+n,zj, X, xi+n,xj+n, Y,yi+n,yj, n, N, B);
		ParallelRecursiveMM(Z,zi+n,zj+n, X, xi+n,xj+n, Y,yi+n,yj+n, n, N, B);
		//cilk_sync;
		
	}
}
*/

int main(int argc, char *argv[]) {

for (int a=1; a<11; a++) {	

	int R = 10;
	int B = pow(2, a);

	if (argc > 1)
		R = atoi(argv[1]);

	int N = pow(2, R);


	TYPE * X = (TYPE *)_mm_malloc(N * N * sizeof(TYPE), ALIGNMENT);

	TYPE * Y= (TYPE *)_mm_malloc(N * N * sizeof(TYPE), ALIGNMENT);

	TYPE * Z = (TYPE *)_mm_malloc(N * N * sizeof(TYPE), ALIGNMENT);

	TYPE * MM = (TYPE *)_mm_malloc(N * N * sizeof(TYPE), ALIGNMENT);

	for(int i = 0; i < N; i++) {
			for(int j = 0; j < N; j++) {
				
				X[i * N + j] = abs((j + i) % 4) + 1;
				Y[i * N + j] = abs((j +i) % 7) + 1;
				Z[i * N + j] = 0;
				MM[i * N + j] = 0;	
			}
	}


	cout << "input size " << N << "\n";

	auto start = std::chrono::system_clock::now();

	A_context actxt;
    	actxt.Z = Z;
    	actxt.zi = 0;
    	actxt.zj = 0;
    	actxt.X = X;
    	actxt.xi = 0;
    	actxt.xj = 0;
    	actxt.Y = Y;
    	actxt.yi = 0;
    	actxt.yj = 0;
    	actxt.size = N;
    	actxt.N = N;
   		actxt.B = B;

   		actxt.m_Atags.put(0);
   		actxt.wait();
	

	//ParallelRecursiveMM(Z,0,0,X,0,0,Y,0,0, N, N, B);

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds= end-start;
	cout<<"Total time for B=:"<< B << " is "<<elapsed_seconds.count()<<"s\n";

	
	IterMMBase(MM, X, Y, N);
			
	cout << "If there is a mistake\n";
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++)  //
				{
					if (MM[i * N + j] != Z[i * N + j]) {
						cout << "Wrong at " << i << j << endl;
						break;
					}
				}
			}

	_mm_free(X);
	_mm_free(Y);
	_mm_free(Z);
	_mm_free(MM);

}
	return 0;
}