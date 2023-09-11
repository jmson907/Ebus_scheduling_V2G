/*********************************************
 * OPL 20.1.0.0 Model
 * Author: Son Jeong Min
 * Creation Date: 2023. 3. 29. at ���� 6:12:14
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
float cp[N] = ...; //���� power(KW)
float dp[N] = ...; //���� power(KW)
float ce = ...; //���� ȿ��(%)
float de = ...; //���� ȿ��(%)
float g[I] = ...; //average energy consumption per time slot(kWh)
float tstart[I] = ...;
float tend[I] = ...;

dvar boolean x[B][N][T]; //��������
dvar boolean y[B][N][T]; //��������
dvar boolean o[B][I][T]; //b�� t������ trip i�� ���࿩��
dvar float+ e[B][T]; //���� b �� t�������� ���·�(kWh)
dvar float+ wbuy[T]; //t ������ ���ŷ�(kWh)
dvar float+ wsell[T]; //t ������ �Ǹŷ�(kWh)

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
    
   