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

#define min(a, b) (a < b ? a : b)
#define max(a, b) (a > b ? a : b)

int N,B;
TYPE * dist;

struct FW_Control{
    TYPE* X;
    TYPE* U;
    TYPE* V;
    int cx,cu,cv;
    int size;
};

struct FW_context;

struct A_FWstep
{
    // declaration of execute method goes here
    int execute( const FW_Control & tag, FW_context & c ) const;
};

struct B_FWstep
{
    // declaration of execute method goes here
    int execute( const FW_Control & tag, FW_context & c ) const;
};

struct C_FWstep
{
    // declaration of execute method goes here
    int execute( const FW_Control & tag, FW_context & c ) const;
};

struct D_FWstep
{
    // declaration of execute method goes here
    int execute( const FW_Control & tag, FW_context & c ) const;
};


struct FW_context : public CnC::context< FW_context >
{

	CnC::step_collection< A_FWstep > m_Asteps;
    CnC::step_collection< B_FWstep > m_Bsteps;
    CnC::step_collection< C_FWstep > m_Csteps;
    CnC::step_collection< D_FWstep > m_Dsteps;

    CnC::item_collection< int, int, int >  m_FWsyncitem;

    CnC::tag_collection< FW_Control > m_Atags;
    CnC::tag_collection< FW_Control > m_Btags;
    CnC::tag_collection< FW_Control > m_Ctags;
    CnC::tag_collection< FW_Control > m_Dtags;

    FW_context()
            : CnC::context< FW_context >(),
            // Initialize each step collection
              m_Asteps( *this ),
              m_Bsteps( *this ),
              m_Csteps( *this ),
              m_Dsteps( *this ),

              m_FWsyncitem(*this),
              
            // Initialize each tag collection
              m_Atags( *this ),
              m_Btags( *this ),
              m_Ctags( *this ),
              m_Dtags( *this )
    {
        
        m_Atags.prescribes( m_Asteps, *this );
        m_Btags.prescribes( m_Bsteps, *this );
        m_Ctags.prescribes( m_Csteps, *this );
        m_Dtags.prescribes( m_Dsteps, *this );
        

        m_Asteps.produces( m_FWsyncitem );
        m_Bsteps.produces( m_FWsyncitem );
        m_Csteps.produces( m_FWsyncitem );
        m_Dsteps.produces( m_FWsyncitem );
        
        m_Asteps.consumes( m_FWsyncitem );
        m_Bsteps.consumes( m_FWsyncitem );
        m_Csteps.consumes( m_FWsyncitem );
        m_Dsteps.consumes( m_FWsyncitem );
        
    }
};


int A_FWstep::execute( const FW_Control & tag, FW_context & FW_ctxt ) const
{
    TYPE * x = tag.X;
    TYPE * u = tag.X;
    TYPE * v = tag.X;
    int cx = tag.cx;
    int cu = tag.cu;
    int cv = tag.cv;
    int size = tag.size;

    if (size <= B) {
        TYPE uv;
        TYPE *xx, *uu, *vv;

        xx = x;
        uu = u;
        vv = v;

        for (int k = 0; k < size; k++) {
            for (int i = 0; i < size; i++) {
                xx = x + i * size;
                for (int j = 0; j < size; j++) {
                    uv = (*uu) + (*vv++);
                    *xx = min(uv, *xx);
                    xx++;
                }
                vv = vv - size;
                uu = uu + size;
            }
            xx = x;
            vv = vv + size;
            uu = u + k + 1;
        }

        FW_ctxt.m_FWsyncitem.put(cx,0);
        return CnC::CNC_Success;
    }
    else {
        int n = size / 2;
        int nn = n * n;

        int m11 = 0;
        int mj = nn;
        int mi = nn * 2;
        int mij = nn * 2 + nn;

        FW_Control tagval;

        tagval.size = n;

        tagval.X = x+m11;
        tagval.U = x+m11;
        tagval.V = x+m11;
        tagval.cx = cx+m11;
        tagval.cu = cu+m11;
        tagval.cv = cv+m11; 
        FW_ctxt.m_Atags.put(tagval);
       
        tagval.X = x+mj;
        tagval.U = x+m11;
        tagval.V = x+mj;
        tagval.cx = cx+mj;
        tagval.cu = cu+m11;
        tagval.cv = cv+mj;
        FW_ctxt.m_Btags.put(tagval);

        tagval.X = x+mi;
        tagval.U = x+mi;
        tagval.V = x+m11;
        tagval.cx = cx+mi;
        tagval.cu = cu+mi;
        tagval.cv = cv+m11;
        FW_ctxt.m_Ctags.put(tagval);

        tagval.X = x+mij;
        tagval.U = x+mi;
        tagval.V = x+mj;
        tagval.cx = cx+mij;
        tagval.cu = cu+mi;
        tagval.cv = cv+mj;
        FW_ctxt.m_Dtags.put(tagval);

        m11 = nn * 2 + nn;
        mj = nn*2;
        mi = nn;
        mij = 0;

        tagval.X = x+m11;
        tagval.U = x+m11;
        tagval.V = x+m11;
        tagval.cx = cx+m11;
        tagval.cu = cu+m11;
        tagval.cv = cv+m11; 
        FW_ctxt.m_Atags.put(tagval);
       
        tagval.X = x+mj;
        tagval.U = x+m11;
        tagval.V = x+mj;
        tagval.cx = cx+mj;
        tagval.cu = cu+m11;
        tagval.cv = cv+mj;
        FW_ctxt.m_Btags.put(tagval);

        tagval.X = x+mi;
        tagval.U = x+mi;
        tagval.V = x+m11;
        tagval.cx = cx+mi;
        tagval.cu = cu+mi;
        tagval.cv = cv+m11;
        FW_ctxt.m_Ctags.put(tagval);

        tagval.X = x+mij;
        tagval.U = x+mi;
        tagval.V = x+mj;
        tagval.cx = cx+mij;
        tagval.cu = cu+mi;
        tagval.cv = cv+mj;
        FW_ctxt.m_Dtags.put(tagval);

    }
    return 0;
}


int B_FWstep::execute( const FW_Control & tag, FW_context & FW_ctxt ) const
{
	TYPE * x = tag.X;
    TYPE * u = tag.U;
    TYPE * v = tag.V;
    int cx = tag.cx;
    int cu = tag.cu;
    int cv = tag.cv;
    int size = tag.size;

    if (size <= B) {
    	int a;
    	FW_ctxt.m_FWsyncitem.get(cu,0,a);
    	if(a==1){
    		TYPE uv;
    		TYPE *xx, *uu, *vv;

    		xx = x;
    		uu = u;
    		vv = v;

    		for (int k = 0; k < size; k++) {

    		    for (int i = 0; i < size; i++) {

    		        for (int j = 0; j < size; j++) {
    		            uv = (*uu) + (*vv++);

    		            *xx = min(uv, *xx);
    		            xx++;
    		        }
    		        vv = vv - B;
    		        uu = uu + B;
    		    }
    		    xx = x;
    		    vv = vv + B;
    		    uu = u + k + 1;
    		}
    		FW_ctxt.m_FWsyncitem.put(cx,1);
    		return CnC::CNC_Success;
	    }
    }
    else{

    	int n = size / 2;
        int nn = n * n;

        int m11 = 0;
        int mj = nn;
        int mi = nn * 2;
        int mij = nn * 2 + nn;

        FW_Control tagval;

        tagval.size = n;

        tagval.X = x+m11;
        tagval.U = u+m11;
        tagval.V = v+m11;
        tagval.cx = cx+m11;
        tagval.cu = cu+m11;
        tagval.cv = cv+m11;
        FW_ctxt.m_Btags.put(tagval);

        tagval.X = x+mj;
        tagval.U = u+m11;
        tagval.V = v+mj;
        tagval.cx = cx+mj;
        tagval.cu = cu+m11;
        tagval.cv = cv+mj;
        FW_ctxt.m_Btags.put(tagval);

        tagval.X = x+mi;
        tagval.U = u+mi;
        tagval.V = v+m11;
        tagval.cx = cx+mi;
        tagval.cu = cu+mi;
        tagval.cv = cv+m11;
        FW_ctxt.m_Dtags.put(tagval);

        tagval.X = x+mij;
        tagval.U = u+mi;
        tagval.V = v+mj;
        tagval.cx = cx+mij;
        tagval.cu = cu+mi;
        tagval.cv = cv+mj;
        FW_ctxt.m_Dtags.put(tagval);

        tagval.X = x+mi;
        tagval.U = u+mij;
        tagval.V = v+mi;
        tagval.cx = cx+mi;
        tagval.cu = cu+mij;
        tagval.cv = cv+mi;
        FW_ctxt.m_Btags.put(tagval);

        tagval.X = x+mij;
        tagval.U = u+mij;
        tagval.V = v+mij;
        tagval.cx = cx+mij;
        tagval.cu = cu+mij;
        tagval.cv = cv+mij;
        FW_ctxt.m_Btags.put(tagval);

        tagval.X = x+m11;
        tagval.U = u+mj;
        tagval.V = v+mi;
        tagval.cx = cx+m11;
        tagval.cu = cu+mj;
        tagval.cv = cv+mi;
        FW_ctxt.m_Dtags.put(tagval);

        tagval.X = x+mj;
        tagval.U = u+mj;
        tagval.V = v+mij;
        tagval.cx = cx+mj;
        tagval.cu = cu+mj;
        tagval.cv = cv+mij;
        FW_ctxt.m_Dtags.put(tagval);

        
    }
}


int C_FWstep::execute( const FW_Control & tag, FW_context & FW_ctxt ) const
{
	TYPE * x = tag.X;
    TYPE * u = tag.U;
    TYPE * v = tag.V;
    int cx = tag.cx;
    int cu = tag.cu;
    int cv = tag.cv;
    int size = tag.size;

    if (size <= B) {
    	int a;
    	FW_ctxt.m_FWsyncitem.get(cv,0,a);
    	if(a==1){
    		TYPE uv;
    		TYPE *xx, *uu, *vv;

    		xx = x;
    		uu = u;
    		vv = v;

    		for (int k = 0; k < size; k++) {

    		    for (int i = 0; i < size; i++) {

    		        for (int j = 0; j < size; j++) {
    		            uv = (*uu) + (*vv++);

    		            *xx = min(uv, *xx);
    		            xx++;
    		        }
    		        vv = vv - B;
    		        uu = uu + B;
    		    }
    		    xx = x;
    		    vv = vv + B;
    		    uu = u + k + 1;
    		}
    		FW_ctxt.m_FWsyncitem.put(cx,2);
    		return CnC::CNC_Success;
    	}
    }
    else{
        int n = size / 2;
        int nn = n * n;

        int m11 = 0;
        int mj = nn;
        int mi = nn * 2;
        int mij = nn * 2 + nn;

        FW_Control tagval;

        tagval.size = n;

        tagval.X = x+m11;
        tagval.U = u+m11;
        tagval.V = v+m11;
        tagval.cx = cx+m11;
        tagval.cu = cu+m11;
        tagval.cv = cv+m11;
        FW_ctxt.m_Ctags.put(tagval);

        tagval.X = x+mi;
        tagval.U = u+mi;
        tagval.V = v+m11;
        tagval.cx = cx+mi;
        tagval.cu = cu+mi;
        tagval.cv = cv+m11;
        FW_ctxt.m_Ctags.put(tagval);

        tagval.X = x+mj;
        tagval.U = u+m11;
        tagval.V = v+mj;
        tagval.cx = cx+mj;
        tagval.cu = cu+m11;
        tagval.cv = cv+mj;
        FW_ctxt.m_Dtags.put(tagval);

        tagval.X = x+mij;
        tagval.U = u+mi;
        tagval.V = v+mj;
        tagval.cx = cx+mij;
        tagval.cu = cu+mi;
        tagval.cv = cv+mj;
        FW_ctxt.m_Dtags.put(tagval);

        tagval.X = x+mj;
        tagval.U = u+mj;
        tagval.V = v+mij;
        tagval.cx = cx+mj;
        tagval.cu = cu+mj;
        tagval.cv = cv+mij;
        FW_ctxt.m_Ctags.put(tagval);

        tagval.X = x+mij;
        tagval.U = u+mij;
        tagval.V = v+mij;
        tagval.cx = cx+mij;
        tagval.cu = cu+mij;
        tagval.cv = cv+mij;
        FW_ctxt.m_Ctags.put(tagval);

        tagval.X = x+m11;
        tagval.U = u+mj;
        tagval.V = v+mi;
        tagval.cx = cx+m11;
        tagval.cu = cu+mj;
        tagval.cv = cv+mi;
        FW_ctxt.m_Dtags.put(tagval);

        tagval.X = x+mi;
        tagval.U = u+mij;
        tagval.V = v+mi;
        tagval.cx = cx+mi;
        tagval.cu = cu+mij;
        tagval.cv = cv+mi;
        FW_ctxt.m_Dtags.put(tagval);
    }

}

int D_FWstep::execute( const FW_Control & tag, FW_context & FW_ctxt ) const
{
	TYPE * x = tag.X;
    TYPE * u = tag.U;
    TYPE * v = tag.V;
    int cx = tag.cx;
    int cu = tag.cu;
    int cv = tag.cv;
    int size = tag.size;

    if (size <= B) {
    	int a,b;
    	FW_ctxt.m_FWsyncitem.get(cu,1);
    	FW_ctxt.m_FWsyncitem.get(cv,2);
    	if(a==1 && b==1){
    		int n = size;
    		TYPE V[n * n];

    		int uin = 0;


    		for (int i = 0; i < n; i++) {
    		    int in = i * n;

    		    for (int j = 0; j < n; j++) {
    		        V[j * n + i] = v[in + j];
    		    }
    		}

    		TYPE *uu, *vv;
    		int in = 0;
    		uin = 0;

    		for (int i = 0; i < n; i++) {
    		    for (int j = 0; j < n; j++) {
    		        TYPE x_ij = x[in + j];

    		        uu = u + uin;
    		        vv = V + j * n;

    		        for (int k = 0; k < n; k++) {
    		            x_ij = min(x_ij, (*uu + *vv));
    		            uu++;
    		            vv++;
    		        }
    		        x[in + j] = x_ij;
    		    }
    		    in = in + n;
    		    uin = uin + n;
    		}
    		FW_ctxt.m_FWsyncitem.put(cx,3);
    	}
    }
    else{
    	int n = size / 2;
        int nn = n * n;

        int m11 = 0;
        int mj = nn;
        int mi = nn * 2;
        int mij = nn * 2 + nn;

        FW_Control tagval;

        tagval.size = n;

    	tagval.X = x+m11;
        tagval.U = u+m11;
        tagval.V = v+m11;
        tagval.cx = cx+m11;
        tagval.cu = cu+m11;
        tagval.cv = cv+m11;
        FW_ctxt.m_Dtags.put(tagval);

        tagval.X = x+mj;
        tagval.U = u+m11;
        tagval.V = v+mj;
        tagval.cx = cx+mj;
        tagval.cu = cu+m11;
        tagval.cv = cv+mj;
        FW_ctxt.m_Dtags.put(tagval);

        tagval.X = x+mi;
        tagval.U = u+mi;
        tagval.V = v+m11;
        tagval.cx = cx+mi;
        tagval.cu = cu+mi;
        tagval.cv = cv+m11;
        FW_ctxt.m_Dtags.put(tagval);

        tagval.X = x+mij;
        tagval.U = u+mi;
        tagval.V = v+mj;
        tagval.cx = cx+mij;
        tagval.cu = cu+mi;
        tagval.cv = cv+mj;
        FW_ctxt.m_Dtags.put(tagval);

        tagval.X = x+m11;
        tagval.U = u+mj;
        tagval.V = v+mi;
        tagval.cx = cx+m11;
        tagval.cu = cu+mj;
        tagval.cv = cv+mi;
        FW_ctxt.m_Dtags.put(tagval);

        tagval.X = x+mj;
        tagval.U = u+mj;
        tagval.V = v+mij;
        tagval.cx = cx+mj;
        tagval.cu = cu+mj;
        tagval.cv = cv+mij;
        FW_ctxt.m_Dtags.put(tagval);

        tagval.X = x+mj;
        tagval.U = u+mij;
        tagval.V = v+mi;
        tagval.cx = cx+mi;
        tagval.cu = cu+mij;
        tagval.cv = cv+mi;
        FW_ctxt.m_Dtags.put(tagval);

        tagval.X = x+mij;
        tagval.U = u+mij;
        tagval.V = v+mij;
        tagval.cx = cx+mij;
        tagval.cu = cu+mij;
        tagval.cv = cv+mij;
        FW_ctxt.m_Dtags.put(tagval);

    }


}

void conv_RM_2_ZM_RM(TYPE *x, int ix, int jx, int len, int R) {
    if (len <= 0)
        return;
    if (len <= B) {
        for (int i = ix; i < ix + len; i++) {
            for (int j = jx; j < jx + len; j++) {
                (*x++) = dist[(i) * (N) + j];
            }
        }
    }
    else {
        int n = len/R;
        int nn = n*n;

        for(int i=0; i<R; i++) {
            for(int j=0; j<R; j++) {
                int xx =  nn*R*i + nn*j;
                int xi =  n*i;
                int xj =  n*j;
                conv_RM_2_ZM_RM(x+xx, ix+xi, jx+xj, n, R);
            }
        }
    }
}

void conv_ZM_2_RM_RM(TYPE *x, TYPE *V, int ix, int jx, int len, int R) {
    if (len <= 0)
        return;
    if (len <= B) {
        for (int i = ix; i < ix + len; i++) {
            for (int j = jx; j < jx + len; j++) {
                V[(i) * (N) + j] = (*x++);
            }
        }
    }
    else {
        int n = len/R;
        int nn = n*n;

        for(int i=0; i<R; i++) {
            for(int j=0; j<R; j++) {
                int xx =  nn*R*i + nn*j;
                int xi =  n*i;
                int xj =  n*j;
                conv_ZM_2_RM_RM(x+xx, V, ix+xi, jx+xj, n, R);
            }
        }
    }
}


void FloydWarshall2(TYPE* D, int n) {
    for (int k = 0; k < n; k++) {
        int ink = k * N;
        for(int i = 0; i < n; i++) {
            int in = i * N;
            int *d = D + in;
            int v = D[in + k];
            int *dk = D + ink;

            for (int j = 0; j < n; j++) {
                {
                    *d = min(*d, v + *dk);
                    dk++;
                    d++;
                }
            }
        }
    }
}


int main(int argc, char *argv[]) {
    N = 32;
    int R = 2;
    B = 16;
    if (argc > 1)
        N = atoi(argv[1]);
    if (argc > 2)
        B = atoi(argv[2]);
    if (argc > 3)
        R = atoi(argv[3]);

    if (B > N)
        B = N;

    dist = (TYPE *)_mm_malloc(N * N * sizeof(TYPE), ALIGNMENT);

    TYPE* X = (TYPE *)_mm_malloc(N * N * sizeof(TYPE), ALIGNMENT);

    TYPE* D = (TYPE *)_mm_malloc(N * N * sizeof(TYPE), ALIGNMENT);

    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            if (i == j) {
                dist[i * N + j] = 0;
                D[i * N + j] = 0;
            }
            else {
                dist[i * N + j] = abs((j - i) % 4) + 1;
                D[i * N + j] = abs((j - i) % 4) + 1;
            }
        }
    }

    conv_RM_2_ZM_RM(X, 0, 0, N, R);

    FW_context FW_ctxt;
    FW_Control FW_tag;

    FW_tag.X = X;
    FW_tag.U = X;
    FW_tag.V = X;
    FW_tag.cx = 0;
    FW_tag.cu = 0;
    FW_tag.cv = 0;
    FW_tag.size = N;

    FW_ctxt.m_Atags.put(FW_tag);
    FW_ctxt.wait();

    //func_A(X, X, X, N, R);

    conv_ZM_2_RM_RM(X, dist, 0, 0, N, R);

    _mm_free(X);

    FloydWarshall2(D, N);

    /*
    cout << "loop version\n";
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    cout << D[i * N + j] << " ";

                }
                cout << "\n";
            }

    cout << "recursive version\n";
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    cout << dist[i * N + j] << " ";
                }
                cout << "\n";
            }
    */

    cout << "If there is a mistake\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++)  //
        {
            if (D[i * N + j] != dist[i * N + j]) {
                cout << "Wrong at" << i << j << endl;
                break;
            }
        }
    }

    cout << "There is no mistake\n";
    _mm_free(dist);
    _mm_free(D);
    return 0;
}