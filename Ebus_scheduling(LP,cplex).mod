/*********************************************
 * OPL 20.1.0.0 Model
 * Author: Son Jeong Min
 * Creation Date: 2023. 4. 12. at ���� 2:01:14
 *********************************************/
int tsize = ...;
range T = 3..tsize; //3~26

int bsize = ...;
range B = 1..bsize; 

float p[T] = ...; // time t�� ���Ű���
float s[T] = ...; // time t�� �ǸŰ���
float capa = ...; // ���� �Ѵ���͸��� �뷮
float U[B] = ...; // ���б� capa
float o[B][T] = ...; // time t�� out ���
float i[B][T] = ...; // time t�� in ���
float hb[B][T] = ...; // time t��  before ������ ���� ���
float alpha[B][T] = ...; //time t������ out ����
float total[B] = ...; // ������ energy
float HA0[B] = ...; // ������ b�� �ʱ⿡������
float ave[B] = ...;//������b�� �뼱�� 1ȸ���� ��� ����Һ�

dvar float+ D[B][T]; //time t�� ������
dvar float+ C[B][T]; //time t�� ������
dvar float+ I[B][T]; //time t�� ������ ��������
dvar float+ O[B][T]; //time t�� �����¿�������
dvar float+ HA[B][T]; //time t�� after���� ��������
dvar float+ HB[B][T]; //time t�� before���� ��������

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
    