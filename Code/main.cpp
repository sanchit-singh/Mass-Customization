#include <ilcplex/ilocplex.h>
ILOSTLBEGIN
#include <stdio.h>
#include <algorithm>    // std::sort
#include <vector>       // std::vector
#include <set>
#include <typeinfo>

#include <boost/unordered_map.hpp>
using namespace boost;
//#include <unordered_map> //uncomment this only when working with stdc++ support (doesnt work with stdlibc++) and c++11 onwards compiler;

char* foldername = new char[200];
char* foldername_LP_files = new char[200];
char* foldername_Dataset_files = new char[200];
char* filename_rd = new char[200];
char* filename_wr = new char[200];

int return_current_multiplicity_of_kj(int k, int j);
void read_data(IloEnv & env);

class CLS_DATA_PRODS;
class CLS_DATA_SUBS;
class CLS_KJ;
class CLS_KJM;
class CLS_VARS_J;
class CLS_VARS_KJ;
class CLS_VARS_KJM;
class CLS_RNGS_K;
class CLS_RNGS_J;
class CLS_RNGS_KJ;
class CLS_RNGS_KJM;

set<int> set_prods;
set<int> set_subs;
set<class CLS_KJ> set_kj;
set<class CLS_KJM> set_kjm;

map<int,CLS_DATA_PRODS> map_k_to_data;
map<int,CLS_DATA_SUBS> map_j_to_data;
map<int, int> map_j_to_M_j;
map<CLS_KJ, int> map_kj_to_multiplicity;
map<CLS_KJ, int> map_kj_to_M_kj;
map<CLS_KJM, int> map_kjm_to_M_kjm;
map<CLS_KJM,vector<CLS_KJM>> map_kjm_to_vec_klm;

class CLS_DATA_PRODS
{
    int root_sub;
    int demand_qty;
    set<int> all_subs;
public:
    CLS_DATA_PRODS(int a, int b)
    {
        root_sub=a;
        demand_qty=b;
    }
    ~CLS_DATA_PRODS()
    {}
    int return_root_sub() const
    {
        return root_sub;
    }
    int return_demand_qty() const
    {
        return demand_qty;
    }
    set<int> return_all_subs() const
    {
        return all_subs;
    }
    void fill_all_subs(int a)
    {
        all_subs.insert(a);
    }
    void print_this() const
    {
        cout<<"demand:"<<demand_qty<<" ,root sub:"<<root_sub<<endl;
    }
};

class CLS_DATA_SUBS
{
    vector<int> vec_child;
    set<int> all_prods;
    float ass_setup_cost;
    float ass_unit_cost;
    float total_setup_cost;
    float total_unit_cost;
public:
    CLS_DATA_SUBS()
    {
    }
    ~CLS_DATA_SUBS()
    {
        //vec_child.clear();
    }
    vector<int> return_vec_child() const
    {
        return vec_child;
    }
    set<int> return_all_prods() const
    {
        return all_prods;
    }
    float return_ass_setup_cost() const
    {
        return ass_setup_cost;
    }
    float return_ass_unit_cost() const
    {
        return ass_unit_cost;
    }
    float return_total_setup_cost() const
    {
        return total_setup_cost;
    }
    float return_total_unit_cost() const
    {
        return total_unit_cost;
    }
    void fill_vec_child(int a)
    {
        vec_child.push_back(a);
    }
    void fill_all_prods(int a)
    {
        all_prods.insert(a);
    }
    void initialize_ass_costs(float a, float b)
    {
        ass_setup_cost=a;
        ass_unit_cost=b;
    }
    void initialize_total_costs(float a,float b)
    {
        total_setup_cost=a;
        total_unit_cost=b;
    }
    void print_this(bool if_print_vec_child)
    {
        if(if_print_vec_child)
        {
            cout<<"child subs:";
            for(vector<int>::iterator itr=vec_child.begin();itr!=vec_child.end();itr++)
            {
                cout<<*itr<<" ,";
            }
        }
        cout<<"ass_setup_cost:"<<ass_setup_cost<<", ass_unit_cost:"<<ass_unit_cost<<endl;
        cout<<"total_setup_cost:"<<total_setup_cost<<" ,total_unit_cost:"<<total_unit_cost<<endl;
    }
};

class CLS_KJ
{
    int k,j;
public:
    CLS_KJ(int k1,int j1)
    {
        k=k1; j=j1;
    }
    ~CLS_KJ()
    {}
    bool operator<( const CLS_KJ & other ) const
    {
        if(k<other.k) return true;
        else if (k==other.k)
        {
            if(j<other.j) return true;
            else return false;
        }
        else return false;
    }
    int return_k() const
    {
        return k;
    }
    int return_j() const
    {
        return j;
    }
    void print_this() const
    {
        cout<<"k:"<<k<<" ,j:"<<j<<endl;
        cout<<"multiplicity:"<<map_kj_to_multiplicity[CLS_KJ(k,j)]<<" ,M_kj:"<<map_kj_to_M_kj[CLS_KJ(k,j)]<<endl;
    }
};

class CLS_KJM
{
    int k,j,m;
public:
    CLS_KJM(int k1,int j1,int m1)
    {
        k=k1; j=j1; m=m1;
    }
    ~CLS_KJM()
    {}
    bool operator<( const CLS_KJM & other ) const
    {
        if(k<other.k) return true;
        else if (k==other.k)
        {
            if(j<other.j) return true;
            else if (j==other.j)
            {
                if(m<other.m) return true;
                else return false;
            }
            else return false;
        }
        else return false;
    }
    int return_k() const
    {
        return k;
    }
    int return_j() const
    {
        return j;
    }
    int return_m() const
    {
        return m;
    }
    void print_this() const
    {
        cout<<"k:"<<k<<" ,j:"<<j<<" ,m:"<<m<<endl;
    }
};

class CLS_VARS_J
{
public: //[sanchit] make variables as private
    IloNumVar w2;
    IloBoolVar delta;
    IloNumVar f; //corresponds to cost of making w2
    CLS_VARS_J()
    {
    }
    ~CLS_VARS_J()
    {
    }
    void initialize_w2(IloEnv & env,float lb, float ub,char * buffer)
    {
        w2=IloNumVar(env,lb,ub,buffer);
    }
    void initialize_f(IloEnv & env,float lb, float ub, char * buffer)
    {
        f=IloNumVar(env,lb,ub,buffer);
    }
    void initialize_delta(IloEnv & env,char * buffer)
    {
        delta=IloBoolVar(env, buffer);
    }
};

class CLS_VARS_KJ
{
public: //[sanchit] make variables as private
    IloIntVar w1;
    IloBoolVar delta;
    IloNumVar f; //corresponds to cost of making w1
    CLS_VARS_KJ()
    {
    }
    ~CLS_VARS_KJ()
    {
    }
    void initialize_w1(IloEnv & env,float lb, float ub, char * buffer)
    {
        w1=IloIntVar(env,lb,ub,buffer);
    }
    void initialize_f(IloEnv & env,float lb, float ub, char * buffer)
    {
        f=IloNumVar(env,lb,ub,buffer);
    }
    void initialize_delta(IloEnv & env,char * buffer)
    {
        delta=IloBoolVar(env, buffer);
    }
};

class CLS_VARS_KJM
{
public: //[sanchit] make variables as private
    IloIntVar w1;
    IloIntVar w2;
    IloIntVar w2_1;
    IloIntVar w2_2;
    CLS_VARS_KJM()
    {
    }
    ~CLS_VARS_KJM()
    {
    }
    void initialize_w1(IloEnv & env,float lb, float ub, char * buffer)
    {
        w1=IloIntVar(env,lb,ub,buffer);
    }
    void initialize_w2(IloEnv & env,float lb, float ub, char * buffer)
    {
        w2=IloIntVar(env,lb,ub,buffer);
    }
    void initialize_w2_1(IloEnv & env,float lb, float ub, char * buffer)
    {
        w2_1=IloIntVar(env,lb,ub,buffer);
    }
    void initialize_w2_2(IloEnv & env,float lb, float ub, char * buffer)
    {
        w2_2=IloIntVar(env,lb,ub,buffer);
    }
};

class CLS_RNGS_K
{
public: //[sanchit] make variables as private
    IloRange c1;
    CLS_RNGS_K()
    {
    }
    ~CLS_RNGS_K()
    {
    }
    void initialize_c1(IloRange & rng)
    {
        c1 = rng;
    }
};

class CLS_RNGS_J
{
public: //[sanchit] make variables as private
    IloRange c6,c7,c8;
    CLS_RNGS_J()
    {
    }
    ~CLS_RNGS_J()
    {
    }
    void initialize_c6(IloRange & rng)
    {
        c6=rng;
    }
    void initialize_c7(IloRange & rng)
    {
        c7=rng;
    }
    void initialize_c8(IloRange & rng)
    {
        c8=rng;
    }
};

class CLS_RNGS_KJ
{
public: //[sanchit] make variables as private
    IloRange c3,c4,c5;
    CLS_RNGS_KJ()
    {
    }
    ~CLS_RNGS_KJ()
    {
    }
    void initialize_c3(IloRange & rng)
    {
        c3=rng;
    }
    void initialize_c4(IloRange & rng)
    {
        c4=rng;
    }
    void initialize_c5(IloRange & rng)
    {
        c5=rng;
    }
};

class CLS_RNGS_KJM
{
public: //[sanchit] make variables as private
    IloRangeArray c2;
    IloRange c9;
    IloRangeArray c10;
    CLS_RNGS_KJM(IloEnv & env,IloInt n_child)
    {
        c2=IloRangeArray(env, n_child); //indexed from 0ton_child
        c10=IloRangeArray(env, n_child);
    }
    ~CLS_RNGS_KJM()
    {
        //[sanchit] figure out problem and uncomment this later
        //c2.clear();
        //c10.clear();
    }
    void add_to_c2(int index,IloRange & rng)
    {
        c2[index]=rng;
    }
    void initialize_c9(IloRange & rng)
    {
        c9=rng;
    }
    void add_to_c10(int index,IloRange & rng)
    {
        c10[index]=rng;
    }
};

int return_current_multiplicity_of_kj(int k, int j)
{
    if(set_kjm.find(CLS_KJM(k,j,1))!=set_kjm.end())
    {
        int m=map_kj_to_multiplicity.find(CLS_KJ(k,j))->second;
        return m;
    }
    else
    {
        return 0;
    }
}

void read_data(IloEnv & env)
{
    char buffer[100];
    
    ifstream file_rd(filename_rd);
    if (!file_rd)
    {
        sprintf("ERROR: could not open file--%s", filename_rd);
        getchar();
    }
    int temp_c;
    file_rd>>buffer>>buffer>>temp_c;
    for(int ind=0;ind<temp_c;ind++)
    {
        file_rd>>buffer;
        set_prods.insert(atoi(buffer));
    }
    
    file_rd>>buffer>>buffer>>temp_c;
    for(int ind=0;ind<temp_c;ind++)
    {
        file_rd>>buffer;
        set_subs.insert(atoi(buffer));
    }
    
    int prod_k,root_sub,demand_qty;
    file_rd>>buffer>>buffer;
    for(int ind=0;ind<set_prods.size();ind++)
    {
        file_rd>>prod_k>>root_sub>>demand_qty;
        map_k_to_data.insert(make_pair(prod_k, CLS_DATA_PRODS(root_sub,demand_qty)));
        file_rd>>buffer;
    }
    
    for(auto itr=set_subs.begin();itr!=set_subs.end();itr++)
    {
        //map_j_to_data[*itr]=CLS_DATA_SUBS(); //an empty initialization
        map_j_to_data.insert(make_pair(*itr, CLS_DATA_SUBS())); //an empty initialization
    }
    int sub_j,child_total,sub_l;
    file_rd>>buffer>>buffer;
    for(int ind=0;ind<set_subs.size();ind++)
    {
        file_rd>>sub_j>>child_total;
        auto itr=map_j_to_data.find(sub_j);
        for(int ind2=0;ind2<child_total;ind2++)
        {
            file_rd>>sub_l;
            itr->second.fill_vec_child(sub_l);
        }
        file_rd>>buffer;
    }
    int sub_j2;
    float cost_temp_a,cost_temp_b;
    file_rd>>buffer>>buffer;
    for(int ind=0;ind<set_subs.size();ind++)
    {
        file_rd>>sub_j2;
        map<int,CLS_DATA_SUBS>::iterator itr=map_j_to_data.find(sub_j2);
        file_rd>>cost_temp_a>>cost_temp_b;
        itr->second.initialize_ass_costs(cost_temp_a, cost_temp_b);
        file_rd>>cost_temp_a>>cost_temp_b;
        itr->second.initialize_total_costs(cost_temp_a, cost_temp_b);
        file_rd>>buffer;
    }
}


int main()
{
    sprintf(foldername, "./subassemblies");
    
    sprintf(foldername_LP_files, "%s/LP_format_files", foldername);
    sprintf(foldername_Dataset_files, "%s/Dataset_files", foldername);
    sprintf(filename_rd, "%s/read_data2.dat", foldername_Dataset_files);
    sprintf(filename_wr, "%s/file.lp", foldername_LP_files);

    char* buffer = new char[200];
    IloEnv env;
    
    IloModel model(env);
    IloCplex cplex(model);
    
    read_data(env);
    
    vector<int> vec_jm_temp;
    for(auto itr=set_prods.begin();itr!=set_prods.end();itr++)
    {
        int k,j,l,m_j,m_l;
        deque<vector<int>> deq_jm_temp;
        k=*itr;
        auto itrk=map_k_to_data.find(k);
        j=itrk->second.return_root_sub();
        m_j=0;
        itrk->second.fill_all_subs(j);
        set_kj.insert(CLS_KJ(k,j));
        auto itrj0=map_j_to_data.find(j);
        itrj0->second.fill_all_prods(k);
        set_kjm.insert(CLS_KJM(k, j, m_j+1));
        map_kj_to_multiplicity[CLS_KJ(k, j)]=m_j+1; //either creates new entry or overwrites existing one
        //
        vec_jm_temp.push_back(j);
        vec_jm_temp.push_back(m_j+1);
        deq_jm_temp.push_back(vec_jm_temp);
        vec_jm_temp.clear();
        while(!deq_jm_temp.empty())
        {
            j=deq_jm_temp.front()[0];
            m_j=deq_jm_temp.front()[1];
            auto itrj=map_j_to_data.find(j);
            deq_jm_temp.pop_front();
            vector<int> vec_child_j_temp=itrj->second.return_vec_child();
            for(auto itr4=vec_child_j_temp.begin();itr4!=vec_child_j_temp.end();itr4++)
            {
                l=*itr4;
                m_l=return_current_multiplicity_of_kj(k, l);
                itrk->second.fill_all_subs(l);
                set_kj.insert(CLS_KJ(k,l));
                auto itrl=map_j_to_data.find(l);
                itrl->second.fill_all_prods(k);
                set_kjm.insert(CLS_KJM(k, l, m_l+1));
                map_kj_to_multiplicity[CLS_KJ(k, l)]=m_l+1;
                map_kjm_to_vec_klm[CLS_KJM(k,j,m_j)].push_back(CLS_KJM(k, l, m_l+1));
                //
                vec_jm_temp.push_back(l);
                vec_jm_temp.push_back(m_l+1);
                deq_jm_temp.push_back(vec_jm_temp);
                vec_jm_temp.clear();
            }
        }
        deq_jm_temp.clear();
    }
    
    for(map<int,CLS_DATA_SUBS>::iterator itr=map_j_to_data.begin();itr!=map_j_to_data.end();itr++)
    {
        int j=itr->first;
        set<int> set_prods=itr->second.return_all_prods();
        int M_j=0;
        for(set<int>::iterator itr2=set_prods.begin();itr2!=set_prods.end();itr2++)
        {
            int k=*itr2;
            int m_j=map_kj_to_multiplicity.find(CLS_KJ(k,j))->second;
            int D_k=map_k_to_data.find(k)->second.return_demand_qty();
            int M_kj=D_k*m_j;
            M_j+=M_kj;
            map_kj_to_M_kj[CLS_KJ(k,j)]=M_kj;
        }
        map_j_to_M_j[j]=M_j;
    }
    for(set<CLS_KJM>::iterator itr=set_kjm.begin();itr!=set_kjm.end();itr++)
    {
        int k=itr->return_k();
        int j=itr->return_j();
        int m=itr->return_m();
        int D_k=map_k_to_data.find(k)->second.return_demand_qty();
        map_kjm_to_M_kjm[CLS_KJM(k,j,m)]=D_k;
    }
    
    for(map<int,CLS_DATA_SUBS>::iterator itr=map_j_to_data.begin();itr!=map_j_to_data.end();itr++)
    {
        cout<<"j:"<<itr->first<<endl;
        itr->second.print_this(true);
        cout<<"M_j:"<<map_j_to_M_j[itr->first]<<endl;
    }

    cout<<endl;
    for(auto itr=set_kjm.begin();itr!=set_kjm.end();itr++)
    {
        itr->print_this();
    }
    
    cout<<endl;
    for(auto itr=map_kjm_to_vec_klm.begin();itr!=map_kjm_to_vec_klm.end();itr++)
    {
        itr->first.print_this();
        cout<<"child nodes:"<<endl;
        for (auto itr2=itr->second.begin();itr2!=itr->second.end();itr2++)
        {
            cout<<"  ";
            itr2->print_this();
        }
    }

    //VARS...
    map<CLS_KJM,CLS_VARS_KJM> map_kjm_to_vars;
    for(auto itr=set_kjm.begin();itr!=set_kjm.end();itr++)
    {
        int k=itr->return_k();
        int j=itr->return_j();
        int m=itr->return_m();
        int M_kjm=map_kjm_to_M_kjm[CLS_KJM(k,j,m)];
        auto pair=map_kjm_to_vars.insert(make_pair(CLS_KJM(k, j, m), CLS_VARS_KJM()));
        if(pair.second==true)
        {
            sprintf(buffer, "w1_kjm(%d,%d,%d)", k,j,m);
            pair.first->second.initialize_w1(env, 0, M_kjm, buffer);
            sprintf(buffer, "w2_kjm(%d,%d,%d)", k,j,m);
            pair.first->second.initialize_w2(env, 0, M_kjm, buffer);
            sprintf(buffer, "w2_1_kjm(%d,%d,%d)", k,j,m);
            pair.first->second.initialize_w2_1(env, 0, M_kjm, buffer);
            sprintf(buffer, "w2_2_kjm(%d,%d,%d)", k,j,m);
            pair.first->second.initialize_w2_2(env, 0, M_kjm, buffer);
        }
    }

    map<int,CLS_VARS_J> map_j_to_vars;
    for(auto itr=set_subs.begin();itr!=set_subs.end();itr++)
    {
        int j=*itr;
        int M_j=map_j_to_M_j[j];
        auto pair=map_j_to_vars.insert(make_pair(j, CLS_VARS_J()));
        if(pair.second==true)
        {
            sprintf(buffer, "w2_j(%d)", j);
            pair.first->second.initialize_w2(env, 0, M_j, buffer);
            sprintf(buffer, "delta_j(%d)", j);
            pair.first->second.initialize_delta(env, buffer);
            sprintf(buffer, "f_j(%d)", j);
            pair.first->second.initialize_f(env, 0, IloInfinity, buffer);
        }
    }
    
    map<CLS_KJ,CLS_VARS_KJ> map_kj_to_vars;
    for(auto itr=set_prods.begin();itr!=set_prods.end();itr++)
    {
        int k=*itr;
        auto itr2=map_k_to_data.find(k);
        set<int> set_subs=itr2->second.return_all_subs();
        for(auto itr2=set_subs.begin();itr2!=set_subs.end();itr2++)
        {
            int j=*itr2;
            int M_kj=map_kj_to_M_kj[CLS_KJ(k,j)];
            auto pair=map_kj_to_vars.insert(make_pair(CLS_KJ(k,j), CLS_VARS_KJ()));
            if(pair.second==true)
            {
                sprintf(buffer, "w1_kj(%d,%d)", k,j);
                pair.first->second.initialize_w1(env, 0, M_kj, buffer);
                sprintf(buffer, "delta_kj(%d,%d)", k,j);
                pair.first->second.initialize_delta(env, buffer);
                sprintf(buffer, "f_kj(%d,%d)", k,j);
                pair.first->second.initialize_f(env, 0, IloInfinity, buffer);
            }
        }
    }

    //OBJ...
    IloExpr cost_expr(env);
    for(auto itr=set_kj.begin();itr!=set_kj.end();itr++)
    {
        int k=itr->return_k();
        int j=itr->return_j();
        cost_expr+=map_kj_to_vars[CLS_KJ(k,j)].f;
    }
    for(auto itr=set_subs.begin();itr!=set_subs.end();itr++)
    {
        int j=*itr;
        cost_expr+=map_j_to_vars[j].f;
    }
    IloObjective obj_cost(env, cost_expr, IloObjective::Minimize);
    obj_cost.setName("obj_cost");
    model.add(obj_cost);

    //RANGES....
    map<int,CLS_RNGS_K> map_k_to_rngs;
    for(auto itr=set_prods.begin();itr!=set_prods.end();itr++)
    {
        int k=*itr;
        auto itr_find=map_k_to_data.find(k); //[sanchit] why is normal access method not working
        int j_root=itr_find->second.return_root_sub();
         auto pair=map_k_to_rngs.insert(make_pair(k, CLS_RNGS_K()));
        if(pair.second==true)
        {
            //c1...
            sprintf(buffer, "c1_k(%d)", k);
            float lb,ub;
            lb=itr_find->second.return_demand_qty();
            ub=lb;
            IloExpr expr(env);
            expr+=map_kjm_to_vars[CLS_KJM(k, j_root, 1)].w1+map_kjm_to_vars[CLS_KJM(k, j_root, 1)].w2_1;
            IloRange rng(env,lb, expr, ub, buffer);
            pair.first->second.initialize_c1(rng);
            model.add(pair.first->second.c1); //always reference this way when using range in model
        }
    }

    map<CLS_KJM,CLS_RNGS_KJM> map_kjm_to_rngs;
    for(auto itr=set_kjm.begin();itr!=set_kjm.end();itr++)
    {
        int k=itr->return_k();
        int j=itr->return_j();
        int m_j=itr->return_m();
        int n_child=map_kjm_to_vec_klm[CLS_KJM(k,j,m_j)].size();
        auto pair=map_kjm_to_rngs.insert(make_pair(CLS_KJM(k,j,m_j), CLS_RNGS_KJM(env,n_child)));
        if(pair.second==true)
        {
            //c9...
            sprintf(buffer, "c9_kjm(%d,%d,%d)",k,j,m_j);
            float lb0,ub0;
            lb0=0;
            ub0=0;
            IloExpr expr0(env);
            expr0+=map_kjm_to_vars[CLS_KJM(k, j, m_j)].w2_1+map_kjm_to_vars[CLS_KJM(k, j, m_j)].w2_2-map_kjm_to_vars[CLS_KJM(k, j, m_j)].w2;
            IloRange rng0(env,lb0, expr0, ub0, buffer);
            pair.first->second.initialize_c9(rng0);
            model.add(pair.first->second.c9); //always reference this way when using range in model

            int index=0;
            for(auto itr2=map_kjm_to_vec_klm[CLS_KJM(k,j,m_j)].begin();itr2!=map_kjm_to_vec_klm[CLS_KJM(k,j,m_j)].end();itr2++)
            {
                int l=itr2->return_j();
                int m_l=itr2->return_m();
                //c2...
                sprintf(buffer, "c2_k(%d)_jm(%d,%d)_lm(%d,%d)",k,j,m_j,l,m_l);
                float lb,ub;
                lb=0;
                ub=0;
                IloExpr expr(env);
                expr+=map_kjm_to_vars[CLS_KJM(k, l, m_l)].w1+map_kjm_to_vars[CLS_KJM(k, l, m_l)].w2_1-map_kjm_to_vars[CLS_KJM(k, j, m_j)].w1;
                IloRange rng(env,lb, expr, ub, buffer);
                pair.first->second.add_to_c2(index, rng);
                model.add(pair.first->second.c2[index]); //always reference this way when using range in model
                //c10...
                sprintf(buffer, "c10_k(%d)_jm(%d,%d)_lm(%d,%d)",k,j,m_j,l,m_l);
                lb=0;
                ub=0;
                IloExpr expr2(env);
                expr2+=map_kjm_to_vars[CLS_KJM(k, l, m_l)].w2_2-map_kjm_to_vars[CLS_KJM(k, j, m_j)].w2;
                IloRange rng2(env,lb, expr2, ub, buffer);
                pair.first->second.add_to_c10(index, rng2);
                model.add(pair.first->second.c10[index]); //always reference this way when using range in model
                //...
                index++;
            }
        }
    }
    
    map<CLS_KJ,CLS_RNGS_KJ> map_kj_to_rngs;
    for(auto itr=set_kj.begin();itr!=set_kj.end();itr++)
    {
        int k=itr->return_k();
        int j=itr->return_j();
        int m_j=map_kj_to_multiplicity[CLS_KJ(k,j)];
        auto pair=map_kj_to_rngs.insert(make_pair(CLS_KJ(k,j), CLS_RNGS_KJ()));
        if(pair.second==true)
        {
            //c3...
            sprintf(buffer, "c3_kj(%d,%d)",k,j);
            float lb,ub;
            lb=0;
            ub=0;
            IloExpr expr(env);
            for(int m=1;m<m_j+1;m++)
            {
                expr+=map_kjm_to_vars[CLS_KJM(k,j,m)].w1;
            }
            IloRange rng(env,lb, expr-map_kj_to_vars[CLS_KJ(k,j)].w1, ub, buffer);
            pair.first->second.initialize_c3(rng);
            model.add(pair.first->second.c3); //always reference this way when using range in model
            //c4...
            sprintf(buffer, "c4_kj(%d,%d)",k,j);
            float lb2,ub2;
            lb2=0;
            ub2=0;
            IloExpr expr2(env);
            expr2+=map_j_to_data[j].return_ass_unit_cost()*map_kj_to_vars[CLS_KJ(k,j)].w1+map_j_to_data[j].return_ass_setup_cost()*map_kj_to_vars[CLS_KJ(k,j)].delta-map_kj_to_vars[CLS_KJ(k,j)].f;
            IloRange rng2(env,lb2,expr2,ub2,buffer);
            pair.first->second.initialize_c4(rng2);
            model.add(pair.first->second.c4); //always reference this way when using range in model
            //c5...
            sprintf(buffer, "c5_kj(%d,%d)",k,j);
            float lb3,ub3;
            lb3=-IloInfinity;
            ub3=0;
            IloExpr expr3(env);
            expr3+=map_kj_to_vars[CLS_KJ(k,j)].w1-map_kj_to_M_kj[CLS_KJ(k,j)]*map_kj_to_vars[CLS_KJ(k,j)].delta;
            IloRange rng3(env,lb3,expr3,ub3,buffer);
            pair.first->second.initialize_c5(rng3);
            model.add(pair.first->second.c5); //always reference this way when using range in model
        }
    }
    
    map<int,CLS_RNGS_J> map_j_to_rngs;
    for(auto itr=set_subs.begin();itr!=set_subs.end();itr++)
    {
        int j=*itr;
        auto pair=map_j_to_rngs.insert(make_pair(j, CLS_RNGS_J()));
        if(pair.second==true)
        {
            //c6...
            sprintf(buffer, "c6_j(%d)",j);
            float lb,ub;
            lb=0;
            ub=0;
            IloExpr expr(env);
            set<int> set_prods=map_j_to_data[j].return_all_prods();
            for(set<int>::iterator itr2=set_prods.begin();itr2!=set_prods.end();itr2++)
            {
                int k=*itr2;
                int m_j=map_kj_to_multiplicity.find(CLS_KJ(k,j))->second;
                for(int m=1;m<m_j+1;m++)
                {
                    expr+=map_kjm_to_vars[CLS_KJM(k,j,m)].w2;
                }
            }
            IloRange rng(env,lb, expr-map_j_to_vars[j].w2, ub, buffer);
            pair.first->second.initialize_c6(rng);
            model.add(pair.first->second.c6); //always reference this way when using range in model
            //c7...
            sprintf(buffer, "c7_j(%d)",j);
            float lb2,ub2;
            lb2=0;
            ub2=0;
            IloExpr expr2(env);
            expr2+=map_j_to_data[j].return_ass_unit_cost()*map_j_to_vars[j].w2+map_j_to_data[j].return_ass_setup_cost()*map_j_to_vars[j].delta-map_j_to_vars[j].f;
            IloRange rng2(env,lb2,expr2,ub2,buffer);
            pair.first->second.initialize_c7(rng2);
            model.add(pair.first->second.c7); //always reference this way when using range in model
            //c8...
            sprintf(buffer, "c8_j(%d)",j);
            float lb3,ub3;
            lb3=-IloInfinity;
            ub3=0;
            IloExpr expr3(env);
            expr3+=map_j_to_vars[j].w2-map_j_to_M_j[j]*map_j_to_vars[j].delta;
            IloRange rng3(env,lb3,expr3,ub3,buffer);
            pair.first->second.initialize_c8(rng3);
            model.add(pair.first->second.c8); //always reference this way when using range in model
        }
    }

    cplex.exportModel(filename_wr);
    cplex.solve();
    
    for(auto itr=set_kjm.begin();itr!=set_kjm.end();itr++)
    {
        int k=itr->return_k();
        int j=itr->return_j();
        int m=itr->return_m();
        printf("w1_kjm(%d,%d,%d):%f \t w2:%f \n",k,j,m,cplex.getValue(map_kjm_to_vars.find(CLS_KJM(k, j, m))->second.w1),cplex.getValue(map_kjm_to_vars.find(CLS_KJM(k, j, m))->second.w2));
    }
    cout<<endl;
    
    for(auto itr=set_kj.begin();itr!=set_kj.end();itr++)
    {
        int k=itr->return_k();
        int j=itr->return_j();
        printf("w1_kj(%d,%d):%f \t f:%f \n",k,j,cplex.getValue(map_kj_to_vars.find(CLS_KJ(k, j))->second.w1),cplex.getValue(map_kj_to_vars.find(CLS_KJ(k, j))->second.f));
    }
    cout<<endl;

    for(auto itr=set_subs.begin();itr!=set_subs.end();itr++)
    {
        int j=*itr;
        printf("w2_j(%d):%f \t f:%f \n",j,cplex.getValue(map_j_to_vars.find(j)->second.w2),cplex.getValue(map_j_to_vars.find(j)->second.f));
    }
    cout<<endl;

    getchar();
}
