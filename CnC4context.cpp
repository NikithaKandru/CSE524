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

struct A_context;
struct B_context;
struct C_context;
struct D_context;

struct FW_Control;

struct FW_Control{
    int X;
    int U;
    int V;
    int xi,xj,ui,uj,vi,vj;
    int N,size,B;
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
    TYPE * X;
    //TYPE * U;
    //TYPE * V;
    int B;
    int size;

    // step collections
    CnC::step_collection< A1_Astep > m_A1steps;
    CnC::step_collection< B1_Astep > m_B1steps;
    CnC::step_collection< C1_Astep > m_C1steps;
    CnC::step_collection< D1_Astep > m_D1steps;
    CnC::step_collection< A2_Astep > m_A2steps;
    CnC::step_collection< B2_Astep > m_B2steps;
    CnC::step_collection< C2_Astep > m_C2steps;
    CnC::step_collection< D2_Astep > m_D2steps;

    // item collection
    CnC::item_collection< int, int >  m_BC1syncitem;
    CnC::item_collection< int, int >  m_BC2syncitem;

    // Tag collections
    CnC::tag_collection< FW_Control > m_A1tags;
    CnC::tag_collection< FW_Control > m_B1tags;
    CnC::tag_collection< FW_Control > m_C1tags;
    CnC::tag_collection< FW_Control > m_D1tags;
    CnC::tag_collection< FW_Control > m_A2tags;
    CnC::tag_collection< FW_Control > m_B2tags;
    CnC::tag_collection< FW_Control > m_C2tags;
    CnC::tag_collection< FW_Control > m_D2tags;

    // The context class constructor
    A_context()
            : CnC::context< A_context >(),
            // Initialize each step collection
              m_A1steps( *this ),
              m_B1steps( *this ),
              m_C1steps( *this ),
              m_D1steps( *this ),
              m_A2steps( *this ),
              m_B2steps( *this ),
              m_C2steps( *this ),
              m_D2steps( *this ),

              m_BC1syncitem(*this),
              m_BC2syncitem(*this),
            // Initialize each tag collection
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
        m_A1tags.prescribes( m_A1steps, *this );
        m_B1tags.prescribes( m_B1steps, *this );
        m_C1tags.prescribes( m_C1steps, *this );
        m_D1tags.prescribes( m_D1steps, *this );
        m_A2tags.prescribes( m_A2steps, *this );
        m_B2tags.prescribes( m_B2steps, *this );
        m_C2tags.prescribes( m_C2steps, *this );
        m_D2tags.prescribes( m_D2steps, *this );


        m_B1steps.produces( m_BC1syncitem );
        m_C1steps.produces( m_BC1syncitem );
        m_D1steps.consumes( m_BC1syncitem );
        m_B2steps.produces( m_BC2syncitem );
        m_C2steps.produces( m_BC2syncitem );
        m_D2steps.consumes( m_BC2syncitem );

    }
};

struct B_Bstep
{
    // declaration of execute method goes here
    int execute( const int & tag, B_context & c ) const;
};

struct B1_Bstep
{
    // declaration of execute method goes here
    int execute( const int & tag, B_context & c ) const;
};

struct B2_Bstep
{
    // declaration of execute method goes here
    int execute( const int & tag, B_context & c ) const;
};

struct B3_Bstep
{
    // declaration of execute method goes here
    int execute( const int & tag, B_context & c ) const;
};

struct B4_Bstep
{
    // declaration of execute method goes here
    int execute( const int & tag, B_context & c ) const;
};

struct D1_Bstep
{
    // declaration of execute method goes here
    int execute( const int & tag, B_context & c ) const;
};

struct D2_Bstep
{
    // declaration of execute method goes here
    int execute( const int & tag, B_context & c ) const;
};

struct D3_Bstep
{
    // declaration of execute method goes here
    int execute( const int & tag, B_context & c ) const;
};

struct D4_Bstep
{
    // declaration of execute method goes here
    int execute( const int & tag, B_context & c ) const;
};


struct B_context : public CnC::context< B_context >
{
    TYPE * X;
    TYPE * U;
    int B;
    int size;

    // step collections
    CnC::step_collection< B_Bstep > m_Bsteps;
    CnC::step_collection< B1_Bstep > m_B1steps;
    CnC::step_collection< B2_Bstep > m_B2steps;
    CnC::step_collection< B3_Bstep > m_B3steps;
    CnC::step_collection< B4_Bstep > m_B4steps;
    CnC::step_collection< D1_Bstep > m_D1steps;
    CnC::step_collection< D2_Bstep > m_D2steps;
    CnC::step_collection< D3_Bstep > m_D3steps;
    CnC::step_collection< D4_Bstep > m_D4steps;

    // item collection
    CnC::item_collection< int, int >  m_BD1syncitem;
    CnC::item_collection< int, int >  m_DB1syncitem;
    CnC::item_collection< int, int >  m_BD2syncitem;
    //CnC::item_collection< int, int >  m_DB2syncitem;
    //CnC::item_collection< int, int >  m_Ditem;

    // Tag collections
    CnC::tag_collection< int > m_Btags;
    CnC::tag_collection< int > m_B1tags;
    CnC::tag_collection< int > m_B2tags;
    CnC::tag_collection< int > m_B3tags;
    CnC::tag_collection< int > m_B4tags;
    CnC::tag_collection< int > m_D1tags;
    CnC::tag_collection< int > m_D2tags;
    CnC::tag_collection< int > m_D3tags;
    CnC::tag_collection< int > m_D4tags;

    // The context class constructor
    B_context()
            : CnC::context< B_context >(),
            // Initialize each step collection
              m_Bsteps( *this ),
              m_B1steps( *this ),
              m_B2steps( *this ),
              m_B3steps( *this ),
              m_B4steps( *this ),
              m_D1steps( *this ),
              m_D2steps( *this ),
              m_D3steps( *this ),
              m_D4steps( *this ),

              m_BD1syncitem(*this),
              m_DB1syncitem(*this),
              m_BD2syncitem(*this),
            //m_DB2syncitem(*this),
            // Initialize each tag collection
              m_Btags( *this ),
              m_B1tags( *this ),
              m_B2tags( *this ),
              m_B3tags( *this ),
              m_B4tags( *this ),
              m_D1tags( *this ),
              m_D2tags( *this ),
              m_D3tags( *this ),
              m_D4tags( *this )
    {
        // Prescriptive relations
        m_Btags.prescribes( m_Bsteps, *this );
        m_B1tags.prescribes( m_B1steps, *this );
        m_B2tags.prescribes( m_B2steps, *this );
        m_B3tags.prescribes( m_B3steps, *this );
        m_B4tags.prescribes( m_B4steps, *this );
        m_D1tags.prescribes( m_D1steps, *this );
        m_D2tags.prescribes( m_D2steps, *this );
        m_D3tags.prescribes( m_D3steps, *this );
        m_D4tags.prescribes( m_D4steps, *this );

        m_B1steps.produces( m_BD1syncitem );
        m_B2steps.produces( m_BD1syncitem );
        m_D1steps.consumes( m_BD1syncitem );
        m_D2steps.consumes( m_BD1syncitem );
        m_D1steps.produces( m_DB1syncitem );
        m_D2steps.produces( m_DB1syncitem );
        m_B3steps.consumes( m_DB1syncitem );
        m_B4steps.consumes( m_DB1syncitem );
        m_B3steps.produces( m_BD2syncitem );
        m_B4steps.produces( m_BD2syncitem );
        m_D3steps.consumes( m_BD2syncitem );
        m_D4steps.consumes( m_BD2syncitem );

    }
};

struct C_Cstep
{
    // declaration of execute method goes here
    int execute( const int & tag, C_context & c ) const;
};

struct C1_Cstep
{
    // declaration of execute method goes here
    int execute( const int & tag, C_context & c ) const;
};

struct C2_Cstep
{
    // declaration of execute method goes here
    int execute( const int & tag, C_context & c ) const;
};

struct C3_Cstep
{
    // declaration of execute method goes here
    int execute( const int & tag, C_context & c ) const;
};

struct C4_Cstep
{
    // declaration of execute method goes here
    int execute( const int & tag, C_context & c ) const;
};

struct D1_Cstep
{
    // declaration of execute method goes here
    int execute( const int & tag, C_context & c ) const;
};

struct D2_Cstep
{
    // declaration of execute method goes here
    int execute( const int & tag, C_context & c ) const;
};

struct D3_Cstep
{
    // declaration of execute method goes here
    int execute( const int & tag, C_context & c ) const;
};

struct D4_Cstep
{
    // declaration of execute method goes here
    int execute( const int & tag, C_context & c ) const;
};



struct C_context : public CnC::context< C_context >
{
    TYPE * X;
    TYPE * V;
    int B;
    int size;

    // step collections
    CnC::step_collection< C_Cstep > m_Csteps;
    CnC::step_collection< C1_Cstep > m_C1steps;
    CnC::step_collection< C2_Cstep > m_C2steps;
    CnC::step_collection< C3_Cstep > m_C3steps;
    CnC::step_collection< C4_Cstep > m_C4steps;
    CnC::step_collection< D1_Cstep > m_D1steps;
    CnC::step_collection< D2_Cstep > m_D2steps;
    CnC::step_collection< D3_Cstep > m_D3steps;
    CnC::step_collection< D4_Cstep > m_D4steps;

    // item collection
    CnC::item_collection< int, int >  m_CD1syncitem;
    CnC::item_collection< int, int >  m_DC1syncitem;
    CnC::item_collection< int, int >  m_CD2syncitem;
    //CnC::item_collection< int, int >  m_DB2syncitem;
    //CnC::item_collection< int, int >  m_Ditem;

    // Tag collections
    CnC::tag_collection< int > m_Ctags;
    CnC::tag_collection< int > m_C1tags;
    CnC::tag_collection< int > m_C2tags;
    CnC::tag_collection< int > m_C3tags;
    CnC::tag_collection< int > m_C4tags;
    CnC::tag_collection< int > m_D1tags;
    CnC::tag_collection< int > m_D2tags;
    CnC::tag_collection< int > m_D3tags;
    CnC::tag_collection< int > m_D4tags;

    // The context class constructor
    C_context()
            : CnC::context< C_context >(),
            // Initialize each step collection
              m_Csteps( *this ),
              m_C1steps( *this ),
              m_C2steps( *this ),
              m_C3steps( *this ),
              m_C4steps( *this ),
              m_D1steps( *this ),
              m_D2steps( *this ),
              m_D3steps( *this ),
              m_D4steps( *this ),

              m_CD1syncitem(*this),
              m_DC1syncitem(*this),
              m_CD2syncitem(*this),
            //m_DB2syncitem(*this),
            // Initialize each tag collection
              m_Ctags( *this ),
              m_C1tags( *this ),
              m_C2tags( *this ),
              m_C3tags( *this ),
              m_C4tags( *this ),
              m_D1tags( *this ),
              m_D2tags( *this ),
              m_D3tags( *this ),
              m_D4tags( *this )
    {
        // Prescriptive relations
        m_Ctags.prescribes( m_Csteps, *this );
        m_C1tags.prescribes( m_C1steps, *this );
        m_C2tags.prescribes( m_C2steps, *this );
        m_C3tags.prescribes( m_C3steps, *this );
        m_C4tags.prescribes( m_C4steps, *this );
        m_D1tags.prescribes( m_D1steps, *this );
        m_D2tags.prescribes( m_D2steps, *this );
        m_D3tags.prescribes( m_D3steps, *this );
        m_D4tags.prescribes( m_D4steps, *this );

        m_C1steps.produces( m_CD1syncitem );
        m_C2steps.produces( m_CD1syncitem );
        m_D1steps.consumes( m_CD1syncitem );
        m_D2steps.consumes( m_CD1syncitem );
        m_D1steps.produces( m_DC1syncitem );
        m_D2steps.produces( m_DC1syncitem );
        m_C3steps.consumes( m_DC1syncitem );
        m_C4steps.consumes( m_DC1syncitem );
        m_C3steps.produces( m_CD2syncitem );
        m_C4steps.produces( m_CD2syncitem );
        m_D3steps.consumes( m_CD2syncitem );
        m_D4steps.consumes( m_CD2syncitem );

    }
};


struct D_Dstep
{
    // declaration of execute method goes here
    int execute( const int & tag, D_context & c ) const;
};

struct D1_Dstep
{
    // declaration of execute method goes here
    int execute( const int & tag, D_context & c ) const;
};
struct D2_Dstep
{
    // declaration of execute method goes here
    int execute( const int & tag, D_context & c ) const;
};
struct D3_Dstep
{
    // declaration of execute method goes here
    int execute( const int & tag, D_context & c ) const;
};
struct D4_Dstep
{
    // declaration of execute method goes here
    int execute( const int & tag, D_context & c ) const;
};
struct D5_Dstep
{
    // declaration of execute method goes here
    int execute( const int & tag, D_context & c ) const;
};
struct D6_Dstep
{
    // declaration of execute method goes here
    int execute( const int & tag, D_context & c ) const;
};
struct D7_Dstep
{
    // declaration of execute method goes here
    int execute( const int & tag, D_context & c ) const;
};
struct D8_Dstep
{
    // declaration of execute method goes here
    int execute( const int & tag, D_context & c ) const;
};


struct D_context : public CnC::context< D_context >
{
    TYPE * X;
    TYPE * U;
    TYPE * V;
    int B;
    int size;

    // step collections
    CnC::step_collection< D_Dstep > m_Dsteps;
    CnC::step_collection< D1_Dstep > m_D1steps;
    CnC::step_collection< D2_Dstep > m_D2steps;
    CnC::step_collection< D3_Dstep > m_D3steps;
    CnC::step_collection< D4_Dstep > m_D4steps;
    CnC::step_collection< D5_Dstep > m_D5steps;
    CnC::step_collection< D6_Dstep > m_D6steps;
    CnC::step_collection< D7_Dstep > m_D7steps;
    CnC::step_collection< D8_Dstep > m_D8steps;

    //CnC::step_collection< D_Cstep > m_Dsteps;

    // item collection
    CnC::item_collection< int, int >  m_DDsyncitem;
    //CnC::item_collection< int, int >  m_Ditem;

    // Tag collections

    CnC::tag_collection< int > m_Dtags;
    CnC::tag_collection< int > m_D1tags;
    CnC::tag_collection< int > m_D2tags;
    CnC::tag_collection< int > m_D3tags;
    CnC::tag_collection< int > m_D4tags;
    CnC::tag_collection< int > m_D5tags;
    CnC::tag_collection< int > m_D6tags;
    CnC::tag_collection< int > m_D7tags;
    CnC::tag_collection< int > m_D8tags;

    //CnC::tag_collection< int > m_Dtags;

    // The context class constructor
    D_context()
            : CnC::context< D_context >(),
            // Initialize each step collection
              m_Dsteps( *this ),
              m_D1steps( *this ),
              m_D2steps( *this ),
              m_D3steps( *this ),
              m_D4steps( *this ),
              m_D5steps( *this ),
              m_D6steps( *this ),
              m_D7steps( *this ),
              m_D8steps( *this ),

              m_DDsyncitem(*this),
            // Initialize each tag collection
            //m_Ctags( *this ),
              m_Dtags( *this ),
              m_D1tags( *this ),
              m_D2tags( *this ),
              m_D3tags( *this ),
              m_D4tags( *this ),
              m_D5tags( *this ),
              m_D6tags( *this ),
              m_D7tags( *this ),
              m_D8tags( *this )
    {

        // Prescriptive relations
        m_Dtags.prescribes( m_Dsteps, *this );
        m_D1tags.prescribes( m_D1steps, *this );
        m_D2tags.prescribes( m_D2steps, *this );
        m_D3tags.prescribes( m_D3steps, *this );
        m_D4tags.prescribes( m_D4steps, *this );
        m_D5tags.prescribes( m_D5steps, *this );
        m_D6tags.prescribes( m_D6steps, *this );
        m_D7tags.prescribes( m_D7steps, *this );
        m_D8tags.prescribes( m_D8steps, *this );

        m_D1steps.produces( m_DDsyncitem );
        m_D2steps.produces( m_DDsyncitem );
        m_D3steps.produces( m_DDsyncitem );
        m_D4steps.produces( m_DDsyncitem );
        m_D5steps.consumes( m_DDsyncitem );
        m_D6steps.consumes( m_DDsyncitem );
        m_D7steps.consumes( m_DDsyncitem );
        m_D8steps.consumes( m_DDsyncitem );
    }
};



int A1_Astep::execute( const int & tag, A_context & Actxt ) const
{
    TYPE * x = Actxt.X;
    TYPE * u = Actxt.X;
    TYPE * v = Actxt.X;
    int B = Actxt.B;
    int size = Actxt.size;

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
        return CnC::CNC_Success;
    }
    else {
        int n = size / 2;
        int nn = n * n;

        int m11 = 0;
        int mj = nn;
        int mi = nn * 2;
        int mij = nn * 2 + nn;

        A_context actxt;
        //func_A(x + m11, u + m11, v + m11, n, R);
        actxt.X = x+m11;
        actxt.B = B;
        actxt.size = n;

        actxt.m_A1tags.put(0);
        actxt.wait();

        Actxt.m_B1tags.put(0);
        Actxt.m_C1tags.put(0);

    }
    return 0;
}



int B1_Astep::execute( const int & tag, A_context & Actxt ) const {

    TYPE * x = Actxt.X;
    TYPE * u = Actxt.X;
    TYPE * v = Actxt.X;
    int size = Actxt.size;

    int n = size / 2;
    int nn = n * n;

    int m11 = 0;
    int mj = nn;
    int mi = nn * 2;
    int mij = nn * 2 + nn;

    //func_B(x+mj, u+m11 , v+mj, n, R);

    B_context bcxtx;
    bcxtx.X = x+mj;
    bcxtx.U = u+m11;
    //bcxtx.V = v+mj;
    bcxtx.B = Actxt.B;
    bcxtx.size = size;

    bcxtx.m_Btags.put(0);
    bcxtx.wait();

    Actxt.m_BC1syncitem.put(0,1);
    Actxt.m_D1tags.put(0);

    return CnC::CNC_Success;
}

int C1_Astep::execute( const int & tag, A_context & Actxt ) const {

    TYPE * x = Actxt.X;
    TYPE * u = Actxt.X;
    TYPE * v = Actxt.X;
    int B = Actxt.B;
    int size = Actxt.size;

    int n = size / 2;
    int nn = n * n;

    int m11 = 0;
    int mj = nn;
    int mi = nn * 2;
    int mij = nn * 2 + nn;

    //func_C(x+mi, u+mi , v+m11, n, R);

    C_context cctxt;
    cctxt.X = x+mi;
    //cctxt.U = u+mi;
    cctxt.V = v+m11;
    cctxt.B = B;
    cctxt.size = size;

    cctxt.m_Ctags.put(0);
    cctxt.wait();

    Actxt.m_BC1syncitem.put(1,1);
    Actxt.m_D1tags.put(1);

    return CnC::CNC_Success;
}


int D1_Astep::execute( const int & tag, A_context & Actxt ) const {

    int a,b;
    Actxt.m_BC1syncitem.get(0,a);
    Actxt.m_BC1syncitem.get(1,b);

    if(a==1 && b==1) {

        TYPE * x = Actxt.X;
        TYPE * u = Actxt.X;
        TYPE * v = Actxt.X;
        int B = Actxt.B;
        int size = Actxt.size;

        int n = size / 2;
        int nn = n * n;

        int m11 = 0;
        int mj = nn;
        int mi = nn * 2;
        int mij = nn * 2 + nn;

        //func_D(x+mij, u+mi, v+mj, n, R);

        D_context dctxt;
        dctxt.X = x+mij;
        dctxt.U = u+mi;
        dctxt.V = v+mj;
        dctxt.B = B;
        dctxt.size = size;

        dctxt.m_Dtags.put(0);
        dctxt.wait();

        Actxt.m_A2tags.put(0);

    }

    return CnC::CNC_Success;
}

int A2_Astep::execute( const int & tag, A_context & Actxt ) const {

    TYPE *x = Actxt.X;
    TYPE *u = Actxt.X;
    TYPE *v = Actxt.X;
    int B = Actxt.B;
    int size = Actxt.size;

    int n = size / 2;
    int nn = n * n;
    int m11 = nn*2 + nn;
    int mj = nn*2;
    int mi = nn;
    int mij = 0;

    //func_A(x + m11, u + m11, v + m11, n, R);
    A_context actxt;
    actxt.X = x+m11;
    actxt.B = B;
    actxt.size = n;
    actxt.m_A1tags.put(0);
    actxt.wait();

    Actxt.m_B2tags.put(0);
    Actxt.m_C2tags.put(0);

}

int B2_Astep::execute( const int & tag, A_context & Actxt ) const {

    TYPE * x = Actxt.X;
    TYPE * u = Actxt.X;
    TYPE * v = Actxt.X;
    int B = Actxt.B;
    int size = Actxt.size;

    int n = size / 2;
    int nn = n * n;
    int m=0;
    int m11 = m + nn*2 + nn;
    int mj = nn*2;
    int mi = nn;
    int mij = m;

    //func_B(x+mj, u+m11 , v+mj, n, R);

    B_context bctxt;
    bctxt.X = x+mj;
    bctxt.U = u+m11;
    bctxt.B = B;
    bctxt.size = size;

    bctxt.m_Btags.put(0);
    bctxt.wait();

    Actxt.m_BC2syncitem.put(0,1);
    Actxt.m_D2tags.put(0);

    return CnC::CNC_Success;
}

int C2_Astep::execute( const int & tag, A_context & Actxt ) const {

    TYPE * x = Actxt.X;
    TYPE * u = Actxt.X;
    TYPE * v = Actxt.X;
    int B = Actxt.B;
    int size = Actxt.size;

    int n = size / 2;
    int nn = n * n;
    int m=0;
    int m11 = m + nn*2 + nn;
    int mj = nn*2;
    int mi = nn;
    int mij = m;
    //func_C(x+mi, u+mi , v+m11, n, R);

    C_context cctxt;
    cctxt.X = x+mi;
    //cctxt.U = u+mi;
    cctxt.V = v+m11;
    cctxt.B = B;
    cctxt.size = size;

    cctxt.m_Ctags.put(0);
    cctxt.wait();

    Actxt.m_BC2syncitem.put(1,1);
    Actxt.m_D2tags.put(0);

    return CnC::CNC_Success;
}

int D2_Astep::execute( const int & tag, A_context & Actxt ) const {

    int a,b;
    Actxt.m_BC2syncitem.get(0,a);
    Actxt.m_BC2syncitem.get(1,b);

    if(a==1 && b==1) {

        TYPE * x = Actxt.X;
        TYPE * u = Actxt.X;
        TYPE * v = Actxt.X;
        int B = Actxt.B;
        int size = Actxt.size;

        int n = size / 2;
        int nn = n * n;

        int m11 = 0;
        int mj = nn;
        int mi = nn * 2;
        int mij = nn * 2 + nn;

        //func_D(x+m11, u+mj, v+mi, n, R);

        D_context dctxt;
        dctxt.X = x+m11;
        dctxt.U = u+mj;
        dctxt.V = v+mi;
        dctxt.B = B;
        dctxt.size = size;

        dctxt.m_Dtags.put(0);
        dctxt.wait();
    }

    return CnC::CNC_Success;
}



int B_Bstep::execute( const int & tag, B_context & Bctxt ) const {

    TYPE * x = Bctxt.X;
    TYPE * u = Bctxt.U;
    TYPE * v = Bctxt.X;
    int B = Bctxt.B;
    int size = Bctxt.size;

    if (size <= B) {
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
        return CnC::CNC_Success;
    }
    else{
        Bctxt.m_B1tags.put(0);
        Bctxt.m_B2tags.put(0);
    }
}

int B1_Bstep::execute( const int & tag, B_context & Bctxt ) const {
    TYPE * x = Bctxt.X;
    TYPE * u = Bctxt.U;
    TYPE * v = Bctxt.X;
    int size = Bctxt.size;

    int n = size/2;
    int nn = n*n;
    int m = 0;
    int m11 = m;
    int mj = nn;
    int mi = nn*2;
    int mij = nn*2 + nn;

    //func_B(x+m11, u+m11 , v+m11, n, R);
    //func_B(x+mj, u+m11 , v+mj, n, R);

    B_context bctxt;
    bctxt.X = x+m11;
    bctxt.U = u+m11;
    bctxt.B = B;
    //bctxt.V = v+m11;
    bctxt.size = size;

    bctxt.m_Btags.put(0);
    bctxt.wait();

    Bctxt.m_BD1syncitem.put(0,1);
    Bctxt.m_D1tags.put(0);
    Bctxt.m_D2tags.put(0);

}

int B2_Bstep::execute( const int & tag, B_context & Bctxt ) const {
    TYPE * x = Bctxt.X;
    TYPE * u = Bctxt.U;
    TYPE * v = Bctxt.X;
    int B = Bctxt.B;
    int size = Bctxt.size;

    int n = size/2;
    int nn = n*n;
    int m = 0;
    int m11 = m;
    int mj = nn;
    int mi = nn*2;
    int mij = nn*2 + nn;

    //func_B(x+m11, u+m11 , v+m11, n, R);
    //func_B(x+mj, u+m11 , v+mj, n, R);

    B_context bctxt;
    bctxt.X = x+mj;
    bctxt.U = u+m11;
    bctxt.B = B;
    //bctxt.V = v+m11;
    bctxt.size = size;

    bctxt.m_Btags.put(0);
    bctxt.wait();

    Bctxt.m_BD1syncitem.put(1,1);
    Bctxt.m_D1tags.put(0);
    Bctxt.m_D2tags.put(0);

}

int D1_Bstep::execute( const int & tag, B_context & Bctxt ) const {

    int a,b;
    Bctxt.m_BD1syncitem.get(0,a);
    Bctxt.m_BD1syncitem.get(1,b);

    if(a==1 && b==1) {
        TYPE * x = Bctxt.X;
        TYPE * u = Bctxt.U;
        TYPE * v = Bctxt.X;
        int B = Bctxt.B;
        int size = Bctxt.size;

        int n = size/2;
        int nn = n*n;
        int m = 0;
        int m11 = m;
        int mj = nn;
        int mi = nn*2;
        int mij = nn*2 + nn;

        //func_D(x+mi, u+mi, v+m11, n, R);
        //func_D(x+mij, u+mi, v+mj, n, R);

        D_context dctxt;
        dctxt.X = x+mi;
        dctxt.U = u+mi;
        dctxt.V = v+m11;
        dctxt.B = B;
        dctxt.size = size;

        dctxt.m_Dtags.put(0);
        dctxt.wait();

        Bctxt.m_DB1syncitem.put(0,1);
        Bctxt.m_B3tags.put(0);
        Bctxt.m_B4tags.put(0);
    }

}

int D2_Bstep::execute( const int & tag, B_context & Bctxt ) const {

    int a,b;
    Bctxt.m_BD1syncitem.get(0,a);
    Bctxt.m_BD1syncitem.get(1,b);

    if(a==1 && b==1) {
        TYPE * x = Bctxt.X;
        TYPE * u = Bctxt.U;
        TYPE * v = Bctxt.X;
        int B = Bctxt.B;
        int size = Bctxt.size;

        int n = size/2;
        int nn = n*n;
        int m = 0;
        int m11 = m;
        int mj = nn;
        int mi = nn*2;
        int mij = nn*2 + nn;

        D_context dctxt;
        dctxt.X = x+mij;
        dctxt.U = u+mi;
        dctxt.V = v+mj;
        dctxt.B = B;
        dctxt.size = size;

        dctxt.m_Dtags.put(0);
        dctxt.wait();

        Bctxt.m_DB1syncitem.put(1,1);
        Bctxt.m_B3tags.put(0);
        Bctxt.m_B4tags.put(0);
    }
}

int B3_Bstep::execute( const int & tag, B_context & Bctxt ) const {

    int a,b;
    Bctxt.m_DB1syncitem.get(0,a);
    Bctxt.m_DB1syncitem.get(1,b);

    if(a==1 && b==1) {
        TYPE * x = Bctxt.X;
        TYPE * u = Bctxt.U;
        TYPE * v = Bctxt.X;
        int B = Bctxt.B;
        int size = Bctxt.size;

        int n = size/2;
        int nn = n*n;
        int m = 0;
        int m11 = m;
        int mj = nn;
        int mi = nn*2;
        int mij = nn*2 + nn;

        //func_D(x+mi, u+mi, v+m11, n, R);
        //func_D(x+mij, u+mi, v+mj, n, R);

        B_context bctxt;
        bctxt.X = x+mi;
        bctxt.U = u+mij;
        bctxt.B = B;
        //bctxt.V = v+m11;
        bctxt.size = size;

        bctxt.m_Btags.put(0);
        bctxt.wait();

        Bctxt.m_BD2syncitem.put(0,1);
        Bctxt.m_D3tags.put(0);
        Bctxt.m_D4tags.put(0);
    }

}

int B4_Bstep::execute( const int & tag, B_context & Bctxt ) const {

    int a,b;
    Bctxt.m_BD1syncitem.get(0,a);
    Bctxt.m_BD1syncitem.get(1,b);

    if(a==1 && b==1) {
        TYPE * x = Bctxt.X;
        TYPE * u = Bctxt.U;
        TYPE * v = Bctxt.X;
        int B = Bctxt.B;
        int size = Bctxt.size;

        int n = size/2;
        int nn = n*n;
        int m = 0;
        int m11 = m;
        int mj = nn;
        int mi = nn*2;
        int mij = nn*2 + nn;

        //func_B(x+m11, u+m11 , v+m11, n, R);
        //func_B(x+mj, u+m11 , v+mj, n, R);

        B_context bctxt;
        bctxt.X = x+mij;
        bctxt.U = u+mij;
        bctxt.B = B;
        //bctxt.V = v+m11;
        bctxt.size = size;

        bctxt.m_Btags.put(0);
        bctxt.wait();

        Bctxt.m_BD2syncitem.put(1,1);
        Bctxt.m_D3tags.put(0);
        Bctxt.m_D4tags.put(0);
    }
}

int D3_Bstep::execute( const int & tag, B_context & Bctxt ) const {

    int a,b;
    Bctxt.m_BD2syncitem.get(0,a);
    Bctxt.m_BD2syncitem.get(1,b);

    if(a==1 && b==1) {
        TYPE * x = Bctxt.X;
        TYPE * u = Bctxt.U;
        TYPE * v = Bctxt.X;
        int B = Bctxt.B;
        int size = Bctxt.size;

        int n = size/2;
        int nn = n*n;
        int m = 0;
        int m11 = m;
        int mj = nn;
        int mi = nn*2;
        int mij = nn*2 + nn;

        //func_D(x+mi, u+mi, v+m11, n, R);
        //func_D(x+mij, u+mi, v+mj, n, R);

        D_context dctxt;
        dctxt.X = x+m11;
        dctxt.U = u+mj;
        dctxt.V = v+mi;
        dctxt.B = B;
        dctxt.size = size;

        dctxt.m_Dtags.put(0);
        dctxt.wait();
    }

}

int D4_Bstep::execute( const int & tag, B_context & Bctxt ) const {

    int a,b;
    Bctxt.m_BD1syncitem.get(0,a);
    Bctxt.m_BD1syncitem.get(1,b);

    if(a==1 && b==1) {
        TYPE * x = Bctxt.X;
        TYPE * u = Bctxt.U;
        TYPE * v = Bctxt.X;
        int B = Bctxt.B;
        int size = Bctxt.size;

        int n = size/2;
        int nn = n*n;
        int m = 0;
        int m11 = m;
        int mj = nn;
        int mi = nn*2;
        int mij = nn*2 + nn;

        //func_B(x+m11, u+m11 , v+m11, n, R);
        //func_B(x+mj, u+m11 , v+mj, n, R);

        D_context dctxt;
        dctxt.X = x+mj;
        dctxt.U = u+mj;
        dctxt.V = v+mij;
        dctxt.B = B;
        dctxt.size = size;

        dctxt.m_Dtags.put(0);
        dctxt.wait();
    }
}


int C_Cstep::execute( const int & tag, C_context & Cctxt ) const {

    TYPE * x = Cctxt.X;
    TYPE * u = Cctxt.X;
    TYPE * v = Cctxt.V;
    int B = Cctxt.B;
    int size = Cctxt.size;

    if (size <= B) {
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
        return CnC::CNC_Success;
    }
    else{
        Cctxt.m_C1tags.put(0);
        Cctxt.m_C2tags.put(0);
    }
    return CnC::CNC_Success;
}

int C1_Cstep::execute( const int & tag, C_context & Cctxt ) const {
    TYPE * x = Cctxt.X;
    TYPE * u = Cctxt.X;
    TYPE * v = Cctxt.V;
    int B = Cctxt.B;
    int size = Cctxt.size;

    int n = size/2;
    int nn = n*n;
    int m = 0;
    int m11 = m;
    int mj = nn;
    int mi = nn*2;
    int mij = nn*2 + nn;

    //func_B(x+m11, u+m11 , v+m11, n, R);
    //func_B(x+mj, u+m11 , v+mj, n, R);

    C_context cctxt;
    cctxt.X = x+m11;
    cctxt.V = v+m11;
    cctxt.B = B;
    //bctxt.V = v+m11;
    cctxt.size = size;

    cctxt.m_Ctags.put(0);
    cctxt.wait();

    Cctxt.m_CD1syncitem.put(0,1);
    Cctxt.m_D1tags.put(0);
    Cctxt.m_D2tags.put(0);


}

int C2_Cstep::execute( const int & tag, C_context & Cctxt ) const {
    TYPE * x = Cctxt.X;
    TYPE * u = Cctxt.X;
    TYPE * v = Cctxt.V;
    int B = Cctxt.B;
    int size = Cctxt.size;

    int n = size/2;
    int nn = n*n;
    int m = 0;
    int m11 = m;
    int mj = nn;
    int mi = nn*2;
    int mij = nn*2 + nn;

    //func_B(x+m11, u+m11 , v+m11, n, R);
    //func_B(x+mj, u+m11 , v+mj, n, R);

    C_context cctxt;
    cctxt.X = x+mi;
    cctxt.V = v+m11;
    cctxt.B = B;
    //bctxt.V = v+m11;
    cctxt.size = size;

    cctxt.m_Ctags.put(0);
    cctxt.wait();

    Cctxt.m_CD1syncitem.put(1,1);
    Cctxt.m_D1tags.put(0);
    Cctxt.m_D2tags.put(0);

}

int D1_Cstep::execute( const int & tag, C_context & Cctxt ) const {

    int a,b;
    Cctxt.m_CD1syncitem.get(0,a);
    Cctxt.m_CD1syncitem.get(1,b);

    if(a==1 && b==1) {
        TYPE * x = Cctxt.X;
        TYPE * u = Cctxt.X;
        TYPE * v = Cctxt.V;
        int B = Cctxt.B;
        int size = Cctxt.size;

        int n = size/2;
        int nn = n*n;
        int m = 0;
        int m11 = m;
        int mj = nn;
        int mi = nn*2;
        int mij = nn*2 + nn;

        //func_D(x+mi, u+mi, v+m11, n, R);
        //func_D(x+mij, u+mi, v+mj, n, R);

        D_context dctxt;
        dctxt.X = x+mj;
        dctxt.U = u+m11;
        dctxt.V = v+mj;
        dctxt.B = B;
        dctxt.size = size;

        dctxt.m_Dtags.put(0);
        dctxt.wait();

        Cctxt.m_DC1syncitem.put(0,1);
        Cctxt.m_C3tags.put(0);
        Cctxt.m_C4tags.put(0);
    }

}

int D2_Cstep::execute( const int & tag, C_context & Cctxt ) const {

    int a,b;
    Cctxt.m_CD1syncitem.get(0,a);
    Cctxt.m_CD1syncitem.get(1,b);

    if(a==1 && b==1) {
        TYPE * x = Cctxt.X;
        TYPE * u = Cctxt.X;
        TYPE * v = Cctxt.V;
        int B = Cctxt.B;
        int size = Cctxt.size;

        int n = size/2;
        int nn = n*n;
        int m = 0;
        int m11 = m;
        int mj = nn;
        int mi = nn*2;
        int mij = nn*2 + nn;

        //func_D(x+mi, u+mi, v+m11, n, R);
        //func_D(x+mij, u+mi, v+mj, n, R);

        D_context dctxt;
        dctxt.X = x+mij;
        dctxt.U = u+mi;
        dctxt.V = v+mj;
        dctxt.B = B;
        dctxt.size = size;

        dctxt.m_Dtags.put(0);
        dctxt.wait();

        Cctxt.m_DC1syncitem.put(1,1);
        Cctxt.m_C3tags.put(0);
        Cctxt.m_C4tags.put(0);
    }
}

int C3_Cstep::execute( const int & tag, C_context & Cctxt ) const {

    int a,b;
    Cctxt.m_DC1syncitem.get(0,a);
    Cctxt.m_DC1syncitem.get(1,b);

    if(a==1 && b==1) {
        TYPE * x = Cctxt.X;
        TYPE * u = Cctxt.X;
        TYPE * v = Cctxt.X;
        int B = Cctxt.B;
        int size = Cctxt.size;

        int n = size/2;
        int nn = n*n;
        int m = 0;
        int m11 = m;
        int mj = nn;
        int mi = nn*2;
        int mij = nn*2 + nn;

        //func_D(x+mi, u+mi, v+m11, n, R);
        //func_D(x+mij, u+mi, v+mj, n, R);

        C_context cctxt;
        cctxt.X = x+mj;
        cctxt.V = v+mij;
        cctxt.B = B;
        //bctxt.V = v+m11;
        cctxt.size = size;

        cctxt.m_Ctags.put(0);
        cctxt.wait();

        Cctxt.m_CD2syncitem.put(0,1);
        Cctxt.m_D3tags.put(0);
        Cctxt.m_D4tags.put(0);
    }

}

int C4_Cstep::execute( const int & tag, C_context & Cctxt ) const {

    int a,b;
    Cctxt.m_DC1syncitem.get(0,a);
    Cctxt.m_DC1syncitem.get(1,b);

    if(a==1 && b==1) {
        TYPE * x = Cctxt.X;
        TYPE * u = Cctxt.X;
        TYPE * v = Cctxt.X;
        int B = Cctxt.B;
        int size = Cctxt.size;

        int n = size/2;
        int nn = n*n;
        int m = 0;
        int m11 = m;
        int mj = nn;
        int mi = nn*2;
        int mij = nn*2 + nn;

        //func_D(x+mi, u+mi, v+m11, n, R);
        //func_D(x+mij, u+mi, v+mj, n, R);

        C_context cctxt;
        cctxt.X = x+mij;
        cctxt.V = v+mij;
        cctxt.B = B;
        //bctxt.V = v+m11;
        cctxt.size = size;

        cctxt.m_Ctags.put(0);
        cctxt.wait();

        Cctxt.m_CD2syncitem.put(1,1);
        Cctxt.m_D3tags.put(0);
        Cctxt.m_D4tags.put(0);
    }
}

int D3_Cstep::execute( const int & tag, C_context & Cctxt ) const {

    int a,b;
    Cctxt.m_CD2syncitem.get(0,a);
    Cctxt.m_CD2syncitem.get(1,b);

    if(a==1 && b==1) {
        TYPE * x = Cctxt.X;
        TYPE * u = Cctxt.X;
        TYPE * v = Cctxt.V;
        int B = Cctxt.B;
        int size = Cctxt.size;

        int n = size/2;
        int nn = n*n;
        int m = 0;
        int m11 = m;
        int mj = nn;
        int mi = nn*2;
        int mij = nn*2 + nn;

        //func_D(x+mi, u+mi, v+m11, n, R);
        //func_D(x+mij, u+mi, v+mj, n, R);

        D_context dctxt;
        dctxt.X = x+m11;
        dctxt.U = u+mj;
        dctxt.V = v+mi;
        dctxt.B = B;
        dctxt.size = size;

        dctxt.m_Dtags.put(0);
        dctxt.wait();
    }

}

int D4_Cstep::execute( const int & tag, C_context & Cctxt ) const {

    int a,b;
    Cctxt.m_CD1syncitem.get(0,a);
    Cctxt.m_CD1syncitem.get(1,b);

    if(a==1 && b==1) {
        TYPE * x = Cctxt.X;
        TYPE * u = Cctxt.X;
        TYPE * v = Cctxt.V;
        int B = Cctxt.B;
        int size = Cctxt.size;

        int n = size/2;
        int nn = n*n;
        int m = 0;
        int m11 = m;
        int mj = nn;
        int mi = nn*2;
        int mij = nn*2 + nn;

        //func_D(x+mi, u+mi, v+m11, n, R);
        //func_D(x+mij, u+mi, v+mj, n, R);

        D_context dctxt;
        dctxt.X = x+mi;
        dctxt.U = u+mij;
        dctxt.V = v+mi;
        dctxt.B = B;
        dctxt.size = size;

        dctxt.m_Dtags.put(0);
        dctxt.wait();
    }
}


int D_Dstep::execute( const int & tag, D_context & Dctxt ) const {

    TYPE * x = Dctxt.X;
    TYPE * u = Dctxt.U;
    TYPE * v = Dctxt.V;
    int B = Dctxt.B;
    int size = Dctxt.size;

    if (size <= B) {
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
        return CnC::CNC_Success;
    }
    else{
        Dctxt.m_D1tags.put(0);
        Dctxt.m_D2tags.put(0);
        Dctxt.m_D3tags.put(0);
        Dctxt.m_D4tags.put(0);

    }
}

int D1_Dstep::execute( const int & tag, D_context & Dctxt ) const {

    TYPE * x = Dctxt.X;
    TYPE * u = Dctxt.U;
    TYPE * v = Dctxt.V;
    int B = Dctxt.B;
    int size = Dctxt.size;

    int n = size/2;
    int nn = n*n;
    int m = 0;
    int m11 = m;
    int mj = nn;
    int mi = nn*2;
    int mij = nn*2 + nn;

    //func_D(x+m11, u+m11, v+m11, n, R);
    D_context dctxt;
    dctxt.X = x+m11;
    dctxt.U = u+m11;
    dctxt.V = v+m11;
    dctxt.B = B;
    dctxt.size = size;

    dctxt.m_Dtags.put(0);
    dctxt.wait();

    Dctxt.m_DDsyncitem.put(0,1);
    Dctxt.m_D5tags.put(0);
    Dctxt.m_D6tags.put(0);
    Dctxt.m_D7tags.put(0);
    Dctxt.m_D8tags.put(0);
}

int D2_Dstep::execute( const int & tag, D_context & Dctxt ) const {

    TYPE * x = Dctxt.X;
    TYPE * u = Dctxt.U;
    TYPE * v = Dctxt.V;
    int B = Dctxt.B;
    int size = Dctxt.size;

    int n = size/2;
    int nn = n*n;
    int m = 0;
    int m11 = m;
    int mj = nn;
    int mi = nn*2;
    int mij = nn*2 + nn;

    //func_D(x+m11, u+m11, v+m11, n, R);
    D_context dctxt;
    dctxt.X = x+mj;
    dctxt.U = u+m11;
    dctxt.V = v+mj;
    dctxt.B = B;
    dctxt.size = size;

    dctxt.m_Dtags.put(0);
    dctxt.wait();

    Dctxt.m_DDsyncitem.put(1,1);
    Dctxt.m_D5tags.put(0);
    Dctxt.m_D6tags.put(0);
    Dctxt.m_D7tags.put(0);
    Dctxt.m_D8tags.put(0);
}

int D3_Dstep::execute( const int & tag, D_context & Dctxt ) const {

    TYPE * x = Dctxt.X;
    TYPE * u = Dctxt.U;
    TYPE * v = Dctxt.V;
    int B = Dctxt.B;
    int size = Dctxt.size;

    int n = size/2;
    int nn = n*n;
    int m = 0;
    int m11 = m;
    int mj = nn;
    int mi = nn*2;
    int mij = nn*2 + nn;

    //func_D(x+m11, u+m11, v+m11, n, R);
    D_context dctxt;
    dctxt.X = x+mi;
    dctxt.U = u+mi;
    dctxt.V = v+m11;
    dctxt.B = B;
    dctxt.size = size;

    dctxt.m_Dtags.put(0);
    dctxt.wait();

    Dctxt.m_DDsyncitem.put(2,1);
    Dctxt.m_D5tags.put(0);
    Dctxt.m_D6tags.put(0);
    Dctxt.m_D7tags.put(0);
    Dctxt.m_D8tags.put(0);
}

int D4_Dstep::execute( const int & tag, D_context & Dctxt ) const {

    TYPE * x = Dctxt.X;
    TYPE * u = Dctxt.U;
    TYPE * v = Dctxt.V;
    int B = Dctxt.B;
    int size = Dctxt.size;

    int n = size/2;
    int nn = n*n;
    int m = 0;
    int m11 = m;
    int mj = nn;
    int mi = nn*2;
    int mij = nn*2 + nn;

    //func_D(x+m11, u+m11, v+m11, n, R);
    D_context dctxt;
    dctxt.X = x+mij;
    dctxt.U = u+mi;
    dctxt.V = v+mj;
    dctxt.B = B;
    dctxt.size = size;

    dctxt.m_Dtags.put(0);
    dctxt.wait();

    Dctxt.m_DDsyncitem.put(3,1);
    Dctxt.m_D5tags.put(0);
    Dctxt.m_D6tags.put(0);
    Dctxt.m_D7tags.put(0);
    Dctxt.m_D8tags.put(0);
}


int D5_Dstep::execute( const int & tag, D_context & Dctxt ) const {
    int a,b,c,d;
    Dctxt.m_DDsyncitem.get(0,a);
    Dctxt.m_DDsyncitem.get(1,b);
    Dctxt.m_DDsyncitem.get(2,b);
    Dctxt.m_DDsyncitem.get(3,b);

    if(a==1 && b==1 && c==1 && d==1) {
        TYPE * x = Dctxt.X;
        TYPE * u = Dctxt.U;
        TYPE * v = Dctxt.V;
        int B = Dctxt.B;
        int size = Dctxt.size;

        int n = size/2;
        int nn = n*n;
        int m = 0;
        int m11 = m;
        int mj = nn;
        int mi = nn*2;
        int mij = nn*2 + nn;

        //func_D(x+m11, u+m11, v+m11, n, R);
        D_context dctxt;
        dctxt.X = x+m11;
        dctxt.U = u+mj;
        dctxt.V = v+mi;
        dctxt.B = B;
        dctxt.size = size;

        dctxt.m_Dtags.put(0);
        dctxt.wait();
    }
}

int D6_Dstep::execute( const int & tag, D_context & Dctxt ) const {
    int a,b,c,d;
    Dctxt.m_DDsyncitem.get(0,a);
    Dctxt.m_DDsyncitem.get(1,b);
    Dctxt.m_DDsyncitem.get(2,b);
    Dctxt.m_DDsyncitem.get(3,b);

    if(a==1 && b==1 && c==1 && d==1) {
        TYPE * x = Dctxt.X;
        TYPE * u = Dctxt.U;
        TYPE * v = Dctxt.V;
        int B = Dctxt.B;
        int size = Dctxt.size;

        int n = size/2;
        int nn = n*n;
        int m = 0;
        int m11 = m;
        int mj = nn;
        int mi = nn*2;
        int mij = nn*2 + nn;

        //func_D(x+m11, u+m11, v+m11, n, R);
        D_context dctxt;
        dctxt.X = x+mj;
        dctxt.U = u+mj;
        dctxt.V = v+mij;
        dctxt.B = B;
        dctxt.size = size;

        dctxt.m_Dtags.put(0);
        dctxt.wait();
    }
}

int D7_Dstep::execute( const int & tag, D_context & Dctxt ) const {
    int a,b,c,d;
    Dctxt.m_DDsyncitem.get(0,a);
    Dctxt.m_DDsyncitem.get(1,b);
    Dctxt.m_DDsyncitem.get(2,b);
    Dctxt.m_DDsyncitem.get(3,b);

    if(a==1 && b==1 && c==1 && d==1) {
        TYPE * x = Dctxt.X;
        TYPE * u = Dctxt.U;
        TYPE * v = Dctxt.V;
        int B = Dctxt.B;
        int size = Dctxt.size;

        int n = size/2;
        int nn = n*n;
        int m = 0;
        int m11 = m;
        int mj = nn;
        int mi = nn*2;
        int mij = nn*2 + nn;

        //func_D(x+m11, u+m11, v+m11, n, R);
        D_context dctxt;
        dctxt.X = x+mi;
        dctxt.U = u+mij;
        dctxt.V = v+mi;
        dctxt.B = B;
        dctxt.size = size;

        dctxt.m_Dtags.put(0);
        dctxt.wait();
    }
}

int D8_Dstep::execute( const int & tag, D_context & Dctxt ) const {
    int a,b,c,d;
    Dctxt.m_DDsyncitem.get(0,a);
    Dctxt.m_DDsyncitem.get(1,b);
    Dctxt.m_DDsyncitem.get(2,b);
    Dctxt.m_DDsyncitem.get(3,b);

    if(a==1 && b==1 && c==1 && d==1) {
        TYPE * x = Dctxt.X;
        TYPE * u = Dctxt.U;
        TYPE * v = Dctxt.V;
        int B = Dctxt.B;
        int size = Dctxt.size;

        int n = size/2;
        int nn = n*n;
        int m = 0;
        int m11 = m;
        int mj = nn;
        int mi = nn*2;
        int mij = nn*2 + nn;

        //func_D(x+m11, u+m11, v+m11, n, R);
        D_context dctxt;
        dctxt.X = x+mij;
        dctxt.U = u+mij;
        dctxt.V = v+mij;
        dctxt.B = B;
        dctxt.size = size;

        dctxt.m_Dtags.put(0);
        dctxt.wait();
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

    A_context actxt;
    B_context bctxt;
    C_context cctxt;
    D_context dctxt;

    actxt.X = X;
    actxt.B = B;
    actxt.size = N;

    actxt.m_A1tags.put(0);
    //actxt.wait();

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