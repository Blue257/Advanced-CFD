#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define len 129

//Code takes around 15-20 mins to converge

float Re=400,dx=0.0078125,dy=0.0078125,dt=0.001;

float erru (double a1[len+1][len],double a2[len+1][len]){
	float s1=0,s2=0,z;
	int i,j;
	for(i=1;i<len;i++){
		for(j=1;j<len-1;j++){
			s1=s1+fabs(a1[i][j]-a2[i][j]);
			//s2=s2+fabs(a2[i][j]);
		}
	}
	z=fabs((s1-s2));
	return z;
}

float errp (double a1[len+1][len+1],double a2[len+1][len+1]){
	float s1=0,s2=0,z=0;
	int i,j;
	for(i=1;i<len;i++){
		for(j=1;j<len;j++){
			s1=s1+fabs(a1[i][j]-a2[i][j]);
			//s2=s2+fabs(a2[i][j]);
		}
	}
	z=fabs(s1-s2);
	return z;
}

void print(double u[len+1][len+1]){
	int i,j;
	for(i=0;i<len+1;i++){
		for(j=0;j<len+1;j++){
			printf("%f ",u[i][j]);
		}
		printf("\n");
	}	
}

///////////////////////////////////////
//MAIN
///////////////////////////////////////

int main (){
	double u0[len+1][len],u1[len+1][len],v0[len][len+1],v1[len][len+1],p0[len+1][len+1],p1[len+1][len+1];
	double error,error_p;
	double RHSU[len][len]={0},RHSV[len][len]={0};
	int i,j,count=0;
	double t1,t2,t3,t4,pint;
	double uvel[len][len],vvel[len][len];
	
	for(j=0;j<len+1;j++){
		for(i=0;i<len;i++){				
			u0[j][i]=0;
			u1[j][i]=0;			
		}						
	}	
	
	for(i=0;i<len;i++){
		u0[len][i]=1;
		u0[len-1][i]=1;		
	}
	
	for(j=0;j<len;j++){
		for(i=0;i<len+1;i++){
			v0[j][i]=0;													
		}
	}
	
	for(j=0;j<len+1;j++){
		for(i=0;i<len+1;i++){
			p0[j][i]=0;
			p1[j][i]=0;													
		}
	}
	
	for(j=0;j<len;j++){
		for(i=0;i<len;i++){
			RHSU[j][i]=0;
			RHSV[j][i]=0;
		}
	}
	
	printf("Hello\n");
///////////////////////////////////////
//LOOPING
///////////////////////////////////////	

	do{
		for(j=1;j<=len-1;j++){
			for(i=1;i<=len-2;i++){
				t1=(dt/dx)*0.25*(pow(((u0[j][i+1]+u0[j][i])),2)-pow(((u0[j][i]+u0[j][i-1])),2));
				t2=(dt/dy)*(((u0[j][i]+u0[j+1][i])*(v0[j][i]+v0[j][i+1])*0.25)-((u0[j][i]+u0[j-1][i])*(v0[j-1][i]+v0[j-1][i+1])*0.25));
				t3=((dt/(dx*dx*Re))*(u0[j][i-1]-(2*u0[j][i])+u0[j][i+1]));
				t4=((dt/(dy*dy*Re))*(u0[j-1][i]-(2*u0[j][i])+u0[j+1][i]));
				RHSU[j][i]=u0[j][i]-t1-t2+t3+t4;
				//printf("%f\t\t ",RHSU[j][i]);
			}
		}	
	
		for(j=1;j<=len-2;j++){
			for(i=1;i<=len-1;i++){
				t1=(dt/dx)*(((u0[j][i]+u0[j+1][i])*(v0[j][i]+v0[j][i+1])*0.25)-((u0[j][i-1]+u0[j+1][i-1])*(v0[j][i]+v0[j][i-1])*0.25));
				t2=(dt/dy)*(pow(((v0[j+1][i]+v0[j][i])),2)-pow(((v0[j][i]+v0[j-1][i])),2))*0.25;
				t3=((dt/(dx*dx*Re))*(v0[j][i-1]-(2*v0[j][i])+v0[j][i+1]));
				t4=((dt/(dy*dy*Re))*(v0[j-1][i]-(2*v0[j][i])+v0[j+1][i]));
				RHSV[j][i]=v0[j][i]-t1-t2+t3+t4;
			}
		}				
		
		error=0;
		do{
			error_p=0;
			for(j=1;j<=len-1;j++){
				for(i=1;i<=len-1;i++){
						p0[j][i]=0.25*((p0[j][i+1]+p0[j][i-1]+p0[j+1][i]+p0[j-1][i])-((dx*dx/dt)*(((RHSU[j][i]-RHSU[j][i-1])/dx)+((RHSV[j][i]-RHSV[j-1][i])/dy))));						
					}
				}
			
		//print(p1);printf("\n\n");print(p0);
			error_p=errp(p1,p0);
			//printf("Error P: %f\n",error_p);
			for(j=0;j<len+1;j++){
				for(i=0;i<len+1;i++){
					p1[j][i]=p0[j][i];																	
				}
			}
			//count++;						
		}while((error_p>0.001)&&(count<=50));
		count=0;
		
		for(j=1;j<len;j++){
			for(i=1;i<len-1;i++){				
					u0[j][i]=((-dt/dx)*(p0[j][i+1]-p0[j][i]))+RHSU[j][i];
				}
			}
		
		
		for(j=1;j<len-1;j++){
			for(i=1;i<len;i++){
				v0[j][i]=((-dt/dx)*(p0[j+1][i]-p0[j][i]))+RHSV[j][i];
				}
			}
		
		for(i=0;i<=len-1;i++){
			u0[len][i]=2-u0[len-1][i];	//TW
			u0[0][i]=-u0[1][i];			//BW
		}	
		
		for(j=2;j<=len-2;j++){
			u0[j][0]=0;					//LW
			u0[j][len-1]=0;				//RW
		}
		
		for(i=0;i<=len;i++){
			v0[len-1][i]=0;				//TW
			v0[0][i]=0;					//BW
		}	
		
		for(j=2;j<=len-3;j++){
			v0[j][0]=-v0[j][1];			//LW
			v0[j][len]=-v0[j][len-1];	//RW
		}
		
		for(i=0;i<=len;i++){
			p0[len][i]=p0[len-1][i];
			p0[0][i]=p0[1][i];
		}
		
		for(j=0;j<=len-2;j++){
			p0[j][0]=p0[j][1];
			p0[j][len]=p0[j][len-1];
		}
		
		error=erru(u1,u0);
		printf("Error: %f\n",error);
		
		for(j=0;j<len+1;j++){
			for(i=0;i<len;i++){			
				u1[j][i]=u0[j][i];}						
			}
		
		
		for(j=0;j<len;j++){
			for(i=0;i<len+1;i++){
				v1[j][i]=v0[j][i];												
			}
		}			
		
	}while(error>0.0001);	
	
	printf("\n\n\nThe u mid velocity is:\n");
	
	double um[len];
	
	for(i=0;i<len;i++){
	    um[i]=(u0[i][64]+u0[i+1][64])*0.5;
	}
	
	for(i=0;i<129;i++){
	    printf("%lf\n",um[i]);
	}
	
	for(i=0;i<len;i++){
		for(j=0;j<len;j++){
			uvel[i][j]=(u0[i][j]+u0[i+1][j])*0.5;
			vvel[i][j]=(v0[i][j]+v0[i][j+1])*0.5;
		}
	}	
	
	double si[len][len]={0};
	
	si[0][0]=0;
	for(j=0;j<len-1;j++){
		si[j+1][0]=si[j][0]-(uvel[j+1][0]*dy);
	}
	
	for(j=0;j<len;j++){
		for(i=0;i<len-1;i++){
			si[j][i+1]=si[j][i]+((vvel[j][i+1])*dx);
		}		
	}
	
	////////////////////////////////////
	//PRINTING
	///////////////////////////////////
	
	FILE *fp;
	
	/*fp=fopen("u-v-vel.txt","w");
	for(i=0;i<len;i++){
		for(j=0;j<len;j++){
			fprintf(fp,"%lf\t%lf\t%lf\t%lf\n",(i*dx),(j*dy),uvel[i][j],vvel[i][j]);
		}	
	}
	fclose(fp);*/
	
	/*fp=fopen("si.txt","w");
	for(i=0;i<len;i++){
		for(j=0;j<len;j++){
			fprintf(fp,"%lf\t%lf\t%lf\n",(i*dx),(j*dy),si[i][j]);
		}		
	}
	fclose(fp);*/
	
	/*fp=fopen("pressure.txt","w");
	for(i=0;i<len;i++){
		for(j=0;j<len;j++){
			fprintf(fp,"%lf\t%lf\t%lf\n",(i*dx),(j*dy),p0[i][j]);		}
		
	}
	fclose(fp);*/
	
	/*double mid=0.5;
	fp=fopen("umid.txt","w");
	for(j=0;j<len;j++){
		fprintf(fp,"%lf\t%lf\t%lf\n",mid,(j*dy),uvel[j][64]);
	}
	fclose(fp);*/
	
	/*fp=fopen("vmid.txt","w");
	for(j=0;j<len;j++){
		fprintf(fp,"%lf\t%lf\t%lf\n",(j*dy),mid,vvel[64][j]);
	}
	fclose(fp);*/
	
	return 0;
}
