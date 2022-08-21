// Girvan-Newman algorithm (skg4078@snu.ac.kr)

// Input description.
// First line, integer n,m : each indicates # of nodes, # of edges.
// Following m lines, each line consists of i,j,w : edge with weight w connects node i, node j.
// i and j must follow 0-based ordering.

// Command example : ./Girvan_Newman_algorithm < test_18.txt

#include <bits/stdc++.h>
#define endl '\n'
#define cediv(a,b) ((a)%(b)==0?((a)/(b)):((a)/(b))+1)
#define fi first
#define se second
#define pb push_back

using namespace std;

typedef long long ll;
typedef unsigned int ui;
typedef unsigned long long ull;

template<typename T>
inline T umax(T &u, T v){return u = max(u, v);}
template<typename T>
inline T umin(T &u, T v){return u = min(u, v);}

#define N 100100
#define M 1001000
int n,m; // # of nodes, # of edges, respectively.


// For division
vector<int> edges[N];
map<pair<int,int>,double> adj_con; // adj와 adj_con은 같은 값이지만, adj_con은 알고리즘이 진행되는 중에 변경되지 않음. 
map<pair<int,int>,double> adj; // Newman algorithm은 weighted graph에 대해서 특별하게 처리하지는 않는다.
map<pair<int,int>,double> betweenness;
map<pair<int,int>,double> betweenness_temp; // 하나의 점에서 시작했을 때의 betweenness
bool visited[N];
ll depth[N];
ll path_cnt[N]; // 각 점까지의 최단경로의 개수
pair<int,int> disconnect_order[M];
double tot_weight;


void prt_map(map<pair<int,int>,double>& m){
    // For debug
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            cout<<m[{i,j}]<<' ';
        }
        cout<<endl;
    }
}

void BFS(int s){
    // update path_cnt
    queue<pair<int,int>> q;
    fill(depth, depth+n+1, 0);
    fill(path_cnt, path_cnt+n+1,0);
    path_cnt[s]=1;
    depth[s]=0;
    q.push({s,0});
    while(!q.empty()){
        pair<int,int> cur=q.front();
        q.pop();
        int idx=cur.fi, d=cur.se;
        if(visited[idx]) continue;
        visited[idx]=true;
        depth[idx]=d;
        for(int i=0; i<edges[idx].size(); i++){
            int next=edges[idx][i];
            if(adj[{idx,next}]==0) continue;
            if(depth[next]==depth[idx]-1) path_cnt[idx]+=path_cnt[next];
            if(!visited[next]) q.push({next,d+1});
        }
    }
}

void DFS(int idx){
    if(visited[idx]) return;
    visited[idx]=true;
    double p_weight=1;
    for(int i=0; i<edges[idx].size(); i++){
        int next=edges[idx][i];
        if(adj[{idx,next}]==0) continue;
        if(!visited[next] && depth[next]==depth[idx]+1) DFS(next);
        if(depth[next]==depth[idx]+1) p_weight+=betweenness_temp[{idx,next}];
    }
    for(int i=0; i<edges[idx].size(); i++){
        int next=edges[idx][i];
        if(adj[{idx,next}]==0) continue;
        if(depth[next]==depth[idx]-1){
            betweenness_temp[{idx,next}]=(p_weight*path_cnt[next])/path_cnt[idx];
            betweenness_temp[{next,idx}]=(p_weight*path_cnt[next])/path_cnt[idx];
        }
    }
}

void Recalculate(){
    // Recalculation. 제거한 edges는 weight matrix의 
    for(int i=0; i<n; i++){
        for(int j=0; j<edges[i].size(); j++){
            if(adj[{i,edges[i][j]}]==0) continue;
            betweenness[{i,edges[i][j]}]=0;
        }
    }
    for(int u=0; u<n; u++){
        for(int i=0; i<n; i++){
            for(int j=0; j<edges[i].size(); j++){
                if(adj[{i,edges[i][j]}]==0) continue;
                betweenness_temp[{i,edges[i][j]}]=0;
            }
        }
        fill(visited,visited+n+1,false);
        BFS(u);
        fill(visited,visited+n+1,false);
        path_cnt[u]=1;
        DFS(u);
        for(int i=0; i<n; i++){
            for(int j=0; j<edges[i].size(); j++){
                if(adj[{i,edges[i][j]}]==0) continue;
                betweenness[{i,edges[i][j]}]+=betweenness_temp[{i,edges[i][j]}];
            }
        }
    }
}

void prt_cut(){
    for(int i=0; i<4*n+20; i++) cout<<'=';
    cout<<endl;
    cout<<"Sequence of edge cut"<<endl;
    for(int i=0; i<m; i++){
        cout<<"Cut #";
        cout.width(2);
        cout<<left<<i+1<<" : "<<disconnect_order[i].fi<<' '<<disconnect_order[i].se<<endl;
    }
    cout<<endl;
}

// For dendrogram
class node{
public:
    node(){}
    node(int l_idx, int r_idx){
        this->l_idx=l_idx, this->r_idx=r_idx;
    }
    int l_idx=-1, r_idx=-1;
    int l,r,cut;
    double Q=0; // Modularity
};
int parent[2*N];
int tree_sz, visit_cnt;
int new_order[N]; // Dendrogram 만들기 위해 새롭게 정의된 node의 순서.
int rev_order[N]; // new_order의 역변환
int cut_tree_idx[N]; // 위에서부터의 cut순서대로, dendrogram에서의 node index
node dendro_tree[2*N];
int root(int idx){return parent[idx]=(parent[idx]==idx?idx:root(parent[idx]));}
void union_c(int a, int b){parent[root(a)]=parent[root(b)]=tree_sz;}
double mx_Q;
bool cut[N], mx_cut[N];
int cluster_cnt=0;
int cluster[N];
vector<int> cluster_element[N];

void dendro_DFS(int idx){
    if(idx<n){
        dendro_tree[idx].l=visit_cnt;
        dendro_tree[idx].r=visit_cnt;
        dendro_tree[idx].cut=-1;
        new_order[visit_cnt++]=idx;
    }
    else{
        dendro_DFS(dendro_tree[idx].l_idx);
        dendro_DFS(dendro_tree[idx].r_idx);
        dendro_tree[idx].l=dendro_tree[dendro_tree[idx].l_idx].l;
        dendro_tree[idx].r=dendro_tree[dendro_tree[idx].r_idx].r;
        dendro_tree[idx].cut=dendro_tree[dendro_tree[idx].l_idx].r;
    }
}

double cluster_Q(int l, int r){
    // new_order기준으로 modularity를 계산할 cluster의 왼쪽 끝과 오른쪽 끝의 idx를 입력.
    double e=0, a=0;
    for(int i_=l; i_<=r; i_++){
        int i=new_order[i_];
        for(int j_=0; j_<edges[i].size(); j_++){
            int j=edges[i][j_];
            a+=adj_con[{i,j}];
            if(l<=rev_order[j] && rev_order[j]<=r){
                e+=adj_con[{i,j}];   
            }
        }
    }
    e=e/(2*tot_weight);
    a=a/(2*tot_weight);
    return e-a*a;
}

void prt_dendrogram(int digit_sz=2, int star_len=20){
    for(int i=0; i<4*n+20; i++) cout<<'=';
    cout<<endl;
    cout<<"Dendrogram";
    for(int i=0; i<4*n-7; i++) cout<<' ';
    cout<<"Modularity"<<endl;
    
    // Print modularity.
    fill(cut,cut+n,false);
    for(int j=0; j<n; j++){
        cout.width(digit_sz);
        cout<<left<<new_order[j];
        cout<<(cut[j]?'|':' ');
        cout<<' ';
    }
    double cur_Q=0; // Modularity is 0 at initial state.
    cout<<":  ";
    if(cur_Q<0) cout<<"(Negative)";
    else for(int t=0; t<cur_Q*star_len; t++) cout<<'*';
    cout<<endl;
    for(int i=0; i<n-1; i++){
        node cur=dendro_tree[cut_tree_idx[i]];
        cut[cur.cut]=true;
        cur_Q=dendro_tree[cut_tree_idx[i]].Q;
        for(int j=0; j<n; j++){
            cout.width(digit_sz);
            cout<<new_order[j];
            cout<<(cut[j]?'|':' ');
            cout<<' ';
        }
        cout<<":  ";
        if(cur_Q<0) cout<<"(Negative)";
        else for(int t=0; t<cur_Q*star_len; t++) cout<<'*';
        cout<<endl;
    }
    cout<<endl;
}

void prt_cluster(){
    // Print clusters that maximize modularity.
    for(int i=0; i<4*n+20; i++) cout<<'=';
    cout<<endl;
    cout<<"Total of "<<cluster_cnt<<" clusters with a modularity "<<mx_Q<<endl;
    for(int i=0; i<cluster_cnt; i++){
        cout<<"Cluster #";
        cout.width(2);
        cout<<left<<i<<" : ";
        for(int j=0; j<cluster_element[i].size(); j++){
            cout.width(2);
            cout<<cluster_element[i][j]<<' ';
        }
        cout<<endl;
    }
}


int main(){
    ios_base::sync_with_stdio(false); cin.tie(NULL);
    cin>>n>>m;
    int a,b;
    double w;
    tot_weight=0;
    

    // Initialize : O(m)
    // First line, integer n,m : each indicates # of nodes, # of edges.
    // Following m lines, each line consists of i,j,w : edge with weight w connects node i, node j.
    // i and j must follow 0-based ordering.
    for(int i=0; i<n+1; i++) edges[i].clear();
    for(int i=0; i<m; i++){
        cin>>a>>b>>w;
        edges[a].push_back(b);
        edges[b].push_back(a);
        adj[{a,b}]=w;
        adj[{b,a}]=w;
        adj_con[{a,b}]=w;
        adj_con[{b,a}]=w;
        tot_weight+=w;
    }


    // Main iteration (Calculate betweenness, remove maximun-between edge.) : O(m^2 n log m)
    Recalculate();
    for(int cut_cnt=0; cut_cnt<m; cut_cnt++){
        Recalculate();
        //cout<<cut_cnt<<"th betweenness"<<endl;
        //prt_map(betweenness);
        
        pair<int,int> mx_edge;
        double mx_betweenness=0;
        for(int i=0; i<n; i++){
            for(int j=0; j<edges[i].size(); j++){
                if(adj[{i,edges[i][j]}]==0) continue;
                if(betweenness[{i,edges[i][j]}]>mx_betweenness){
                    mx_betweenness=betweenness[{i,edges[i][j]}];
                    mx_edge={i,edges[i][j]};
                }
            }
        }
        adj[{mx_edge.fi,mx_edge.se}]=0;
        adj[{mx_edge.se,mx_edge.fi}]=0;
        disconnect_order[cut_cnt]={min(mx_edge.fi,mx_edge.se),max(mx_edge.fi,mx_edge.se)};
    }


    // Build cutting tree using union-find : O(m log n)
    tree_sz=n; // Dendrogram을 구성하는 새로운 node의 index
    for(int i=0; i<2*n; i++) parent[i]=i;
    for(int i=m-1; i>=0; i--){
        int a=disconnect_order[i].fi, b=disconnect_order[i].se;
        if(root(a)!=root(b)){
            dendro_tree[tree_sz]=node(root(a),root(b));
            cut_tree_idx[2*n-2-tree_sz]=tree_sz;
            union_c(a,b);
            tree_sz++;
        }
    }
    visit_cnt=0;
    dendro_DFS(tree_sz-1);
    for(int i=0; i<n; i++) rev_order[new_order[i]]=i;


    // Build dendrogram : O(n m)
    double cur_Q=0;
    mx_Q=0;
    fill(cut,cut+n,false);
    for(int i=0; i<n-1; i++){
        node &cur=dendro_tree[cut_tree_idx[i]];
        cut[cur.cut]=true;
        int last_cut=-1;
        cur_Q=0;
        for(int j=0; j<n; j++){
            if(cut[j]){
                cur_Q+=cluster_Q(last_cut+1,j);
                last_cut=j;
            }
        }
        if(last_cut!=n-1) cur_Q+=cluster_Q(last_cut+1,n-1);
        if(cur_Q>mx_Q){
            mx_Q=cur_Q;
            for(int j=0; j<n; j++) mx_cut[j]=cut[j];
        }
        cur.Q=cur_Q;
    }


    // Assign clusters to maximize modularity.
    cluster_cnt=0;
    for(int i=0; i<n; i++) cluster_element[i].clear();
    for(int i=0; i<n; i++){
        cluster[new_order[i]]=cluster_cnt;
        cluster_element[cluster_cnt].push_back(new_order[i]);
        if(mx_cut[i]) cluster_cnt++;
    }
    cluster_cnt++;
    for(int i=0; i<cluster_cnt; i++) sort(cluster_element[i].begin(), cluster_element[i].end());



    // Print sequence of removed edges.
    prt_cut();
    // Print dendrogram.
    prt_dendrogram();
    // Print cluster.
    prt_cluster();
    return 0;
}