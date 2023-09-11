/*********************************************
 * OPL 20.1.0.0 Model
 * Author: Son Jeong Min
 * Creation Date: 2023. 3. 29. at 오후 6:12:14
 *********************************************/
int bsize = ...;
range B = 1..bsize; 
 
int tsize = ...;
range T = 1..tsize;

int isize = ...;
range I = 1..isize;

int nsize = ...;
range N = 1..nsize;


float p[T] = ...;
float s[T] = ...;
float Emax = ...;
float Emin = ...;
float E_0 = ...;
float C = ...;
float U = ...;
float cp[N] = ...; //충전 power(KW)
float dp[N] = ...; //방전 power(KW)
float ce = ...; //충전 효율(%)
float de = ...; //방전 효율(%)
float g[I] = ...; //average energy consumption per time slot(kWh)
float tstart[I] = ...;
float tend[I] = ...;

dvar boolean x[B][N][T]; //충전여부
dvar boolean y[B][N][T]; //방전여부
dvar boolean o[B][I][T]; //b가 t시점에 trip i를 운행여부
dvar float+ e[B][T]; //버스 b 의 t시점에서 전력량(kWh)
dvar float+ wbuy[T]; //t 시점의 구매량(kWh)
dvar float+ wsell[T]; //t 시점의 판매량(kWh)

maximize sum(t in T) wsell[t] * s[t] - sum(t in T) wbuy[t] * p[t];

subject to {
  forall(b in B, t in T) 
    sum(i in I) o[b][i][t] + sum(n in N) x[b][n][t] + sum(n in N) y[b][n][t] <= 1;
  forall(i in I,t in T) 
    if (tstart[i] <= t <= tend[i]){
      sum(b in B) o[b][i][t] == 1;
        } else{
        sum(b in B) o[b][i][t] == 0;
          }       
  forall(i in I, b in B,t in T) 
    if (tstart[i] <= t <= tend[i] - 1){
      o[b][i][t+1] >= o[b][i][t];
    } 
  forall(n in N, t in T) 
    sum(b in B) x[b][n][t] + sum(b in B) y[b][n][t] <=1;
  forall(b in B,t in 2..145)
     e[b][t] == e[b][t-1] + sum(n in N) ce * cp[n] * x[b][n][t] - sum(n in N) dp[n] * y[b][n][t] / de - sum(i in I) g[i] * o[b][i][t] ;
  forall(t in T) 
     sum(b in B,n in N) cp[n] * x[b][n][t] == wbuy[t];
  forall(t in T) 
     sum(b in B,n in N) dp[n] * y[b][n][t] == wsell[t];
  forall(b in B,t in T)	
	 C * Emin <= e[b][t] <= C * Emax ;
  forall(b in B) 
     e[b][1] == C * E_0 ;
  forall(t in T) 
    sum(b in B,n in N) cp[n] * x[b][n][t] + sum(b in B,n in N) dp[n] * y[b][n][t] <= U; 
     }
    
   