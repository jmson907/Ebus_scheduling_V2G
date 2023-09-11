/*********************************************
 * OPL 20.1.0.0 Model
 * Author: Son Jeong Min
 * Creation Date: 2023. 4. 12. at 오후 2:01:14
 *********************************************/
int tsize = ...;
range T = 3..tsize; //3~26

int bsize = ...;
range B = 1..bsize; 

float p[T] = ...; // time t의 구매가격
float s[T] = ...; // time t의 판매가격
float capa = ...; // 버스 한대배터리의 용량
float U[B] = ...; // 변압기 capa
float o[B][T] = ...; // time t의 out 대수
float i[B][T] = ...; // time t의 in 대수
float hb[B][T] = ...; // time t의  before 시점의 버스 대수
float alpha[B][T] = ...; //time t에서의 out 비율
float total[B] = ...; // 총주행 energy
float HA0[B] = ...; // 차고지 b별 초기에너지량
float ave[B] = ...;//차고지b별 노선당 1회운행 평균 연료소비량

dvar float+ D[B][T]; //time t의 방전량
dvar float+ C[B][T]; //time t의 충전량
dvar float+ I[B][T]; //time t에 들어오는 에너지량
dvar float+ O[B][T]; //time t에 나가는에너지량
dvar float+ HA[B][T]; //time t의 after시점 에너지량
dvar float+ HB[B][T]; //time t의 before시점 에너지량

maximize sum(b in B, t in T) s[t] * D[b][t] - sum(b in B, t in T) p[t] * C[b][t];


subject to {
  forall(b in B, t in 3..5)
     I[b][t] ==0;
  
  forall(b in B)   
     HB[b][3] == I[b][3] - O[b][3] + HA0[b];
  
  forall(b in B,t in 4..26) //1
    HB[b][t] == I[b][t] - O[b][t] + HA[b][t-1];
  forall(b in B,t in T) //2
    HA[b][t] == C[b][t] - D[b][t] + HB[b][t];
  forall(b in B,t in T) //3
    C[b][t] + D[b][t] <= U[b];      
  forall(b in B)
    O[b][3] == HA0[b] * alpha[b][3];
  forall(b in B, t in 4..26) //4
    O[b][t] == HA[b][t-1] * alpha[b][t];
  forall(b in B,t in T) //5
    0 <= C[b][t] && C[b][t] <= (0.8 * capa * hb[b][t]) - HB[b][t];
  forall(b in B,t in T) //6
    0 <= D[b][t] && D[b][t] <= HB[b][t] - (0.2 * capa * hb[b][t]) ;
  forall(b in B)
     HA0[b] + sum(t in T) C[b][t] - sum(t in T) D[b][t] >= total[b]; //7
  
  forall(b in B, t in 6..26) //9
   if(o[b][t-2] + o[b][t-3] != 0){
     I[b][t] == ((O[b][t-2]  + O[b][t-3]) - ((o[b][t-2] + o[b][t-3]) * ave[b])) * i[b][t] / (o[b][t-2] + o[b][t-3]);}
   else{I[b][t] == 0;}
  
    
 
  }
    