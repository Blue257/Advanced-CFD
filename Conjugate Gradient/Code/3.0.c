#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define len 101
#define pi 3.14159265

void print1(double u[len*len]){
	int i,j;
	for(j=0;j<len*len;j++){
			printf("%f\n",u[j]);
		}				
}

void print(double u[len][len]){
	int i,j;
	for(i=0;i<len;i++){
		for(j=0;j<len;j++){
			printf("%f ",u[i][j]);
		}
		printf("\n");
	}	
}

double err(double a1[len*len], double a0[len*len]){
	int i;
	double error=0;
	for(i=0;i<len*len;i++){
		error=error+fabs(a1[i]-a0[i]);
	}
	return error;
}

void printA(double a[len*len],double b[len*len],double c[len*len],double d[len*len],double e[len*len]){
	int i;
	printf("Ap\tAn\tAs\tAe\tAw\n\n");
	for(i=0;i<len*len;i++){
		printf("%f\t%f\t%f\t%f\t%f\n",a[i],b[i],c[i],d[i],e[i]);
	}
}

double dot(double a1[len*len],double a2[len*len]){
	int i;
	double z,sum=0;
	for(i=0;i<len*len;i++){
		sum=sum+(a1[i]*a2[i]);
	}
	z=sum;
	return z;
}

void mat_mult(double Ap[len*len],double An[len*len],double As[len*len],double Ae[len*len],double Aw[len*len],double x[len*len],double a[len*len]){
	double t1[len*len]={0},t2[len*len]={0},t3[len*len]={0},t4[len*len]={0},t5[len*len]={0};
	int i;
	
	for(i=0;i<len*len;i++){
		t1[i]=Ap[i]*x[i];
	}
	
	for(i=0;i<(len*len)-1;i++){
		t2[i]=An[i]*x[i+1];
	}
	
	for(i=1;i<(len*len);i++){
		t3[i]=As[i]*x[i-1];
	}
	
	for(i=0;i<(len*len)-len;i++){
		t4[i]=Ae[i]*x[i+len];
	}
	
	for(i=len;i<(len*len);i++){
		t5[i]=Aw[i]*x[i-len];
	}
	
	for(i=0;i<(len*len);i++){
		a[i]=(t1[i]+t2[i]+t3[i]+t4[i]+t5[i]);
	}
	
}

//////////////////////
//MAIN
/////////////////////

int main(){
	double length=len,dx=(1/(length-1)),dy=(1/(length-1));
	double Tan[len][len]={0},T[len][len]={0},x[len*len]={0},b[len*len]={0};
	double Ap[len*len]={0},An[len*len]={0},As[len*len]={0},Ae[len*len]={0},Aw[len*len]={0};
	int i,j,k;
	double sum2=0,error=100,alpha,beta;
	double r1[len*len]={0},r0[len*len]={0},p[len*len]={0};
	double temp[len*len]={0},t1=0,t2=0;	
			
	//////////////////////
	//Analytical
	/////////////////////
	
	for(j=1;j<len-1;j++){
		for(i=1;i<len-1;i++){
			Tan[j][i]=sin(2*pi*i*dx)*sin(2*pi*j*dy);		    
		}
	}	
	
	//////////////////////
	//A,b creation
	/////////////////////
	
	for(i=0;i<len*len;i++){
		Ap[i]=1;		
	}	
	for(i=len;i<((len*len)-len);i++){
		if(((i%len)!=0)&&((i%len)!=(len-1))){
			Ap[i]=-4;
			An[i]=1;
			As[i]=1;
			Ae[i]=1;
			Aw[i]=1;
		}
	}	
	
	k=0;
	for(i=0;i<((len));i++){
		for(j=0;j<len;j++){
			b[k]=-8*pi*pi*(sin(2*pi*i*dx)*sin(2*pi*j*dy))*(dx*dx);
			k++;
		}
	}	
	
	//////////////////////
	//Start
	/////////////////////
	
	mat_mult(Ap,An,As,Ae,Aw,x,temp);
	for(i=0;i<(len*len);i++){
		r1[i]=b[i]-temp[i];
	}
	for(i=0;i<(len*len);i++){
		p[i]=r1[i];
	}
	
	//print1(r1);
	//printf("\n\n");
	
	for(i=0;i<(len*len);i++){
		temp[i]=0;
	}	
	
	k=0;
	//////////////////////
	//Looping
	/////////////////////
	
	while(fabs(error)>1.0e-10){
	//while(k<=len-1){
	
		t1=dot(r1,r1);		
		mat_mult(Ap,An,As,Ae,Aw,p,temp);
	
		//print1(temp);
		//printf("\n\n");
		t2=dot(p,temp);	
		alpha=(t1/t2);
		for(i=0;i<(len*len);i++){
			temp[i]=0;
		}
		
		for(i=0;i<(len*len);i++){
			x[i]=x[i]+(alpha*p[i]);
		}
		
		for(i=0;i<(len*len);i++){
			r0[i]=r1[i];
		}
				
		mat_mult(Ap,An,As,Ae,Aw,p,temp);
		
		for(i=0;i<(len*len);i++){
			r1[i]=r1[i]-(alpha*temp[i]);
		}
		for(i=0;i<(len*len);i++){
			temp[i]=0;
		}
		
		//print1(r1);
		//printf("\n\n");
		
		t1=dot(r1,r1);
		t2=dot(r0,r0);
		beta=(t1/t2);
		
		for(i=0;i<(len*len);i++){
			p[i]=r1[i]+(beta*p[i]);
		}
		
		error=dot(r1,r1);
		//printf("Error: %lf\n",error);
		k++;
		
	}
	
	printf("\n\n");
	print(Tan);
	
	//////////////////////
	//Convert 1D to 2D
	/////////////////////
		
	k=0;
	
	for(i=0;i<len;i++){
	    for(j=0;j<len;j++){
	        T[j][i]=x[k];
	        k++;
	    }
	}
	
	printf("\n\n");
	print(T);
	
	////////////////////////////////////
	//PRINTING
	///////////////////////////////////
	
	FILE *fp;
	
	fp=fopen("Temp3.txt","w");
	for(i=0;i<len;i++){
		for(j=0;j<len;j++){
			fprintf(fp,"%lf\t%lf\t%lf\t%lf\n",(j*dx),(i*dy),T[i][j],Tan[i][j]);
		}	
	}
	fclose(fp);
	
	int mid=(len-1)*0.5;
	fp=fopen("T mid x y 3.txt","w");
	for(j=0;j<len;j++){
		fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",(j*dx),T[mid][j],T[j][mid],Tan[mid][j],Tan[j][mid]);
	}	
	fclose(fp);
	
	return 0;
}
