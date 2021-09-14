#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define len 101
#define pi 3.14159265

void print1(double u[len*len][len*len]){
	int i,j;
	for(i=0;i<len*len;i++){
		for(j=0;j<len*len;j++){
			printf("%f ",u[i][j]);
		}
		printf("\n");
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

void printval(double a[len*len],double b[len*len],double c[len*len],double d[len*len],double e[len*len],double f[len*len],double g[len*len],double h[len*len]){
	int i;
	printf("Lp\t\tLs\t\tLw\t\tUn\n\n");
	for(i=0;i<len*len;i++){
		printf("%f\t%f\t%f\t%f\n",a[i],b[i],c[i],d[i]);
	}
	printf("\n\nUe\t\tMnw\t\tMse\t\tb\n\n");
	for(i=0;i<len*len;i++){
		printf("%f\t%f\t%f\t%f\n",e[i],f[i],g[i],h[i]);
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

void crtL(double Lp[len*len], double Ls[len*len], double Lw[len*len]){
	double L[len*len][len*len]={0};
	int i;
	for(i=0;i<len*len;i++){
		L[i][i]=Lp[i];
		if(i>0){
			L[i][i-1]=Ls[i];
		}
		if(i>=len){
			L[i][i-len]=Lw[i];
		}
	}
	print1(L);
}

void printA(double a[len*len],double b[len*len],double c[len*len],double d[len*len],double e[len*len]){
	int i;
	printf("Ap\tAn\tAs\tAe\tAw\n\n");
	for(i=0;i<len*len;i++){
		printf("%f\t%f\t%f\t%f\t%f\n",a[i],b[i],c[i],d[i],e[i]);
	}
}

//////////////////////
//MAIN
/////////////////////

int main(){
	double length=len,dx=(1/(length-1)),dy=(1/(length-1));
	double Tan[len][len]={0},x[len*len]={0},b[len*len]={0};
	double Ap[len*len]={0},An[len*len]={0},As[len*len]={0},Ae[len*len]={0},Aw[len*len]={0};
	double Lw[len*len]={0},Ls[len*len]={0},Lp[len*len]={0},Un[len*len]={0},Ue[len*len]={0},Mnw[len*len]={0},Mse[len*len]={0};
	int i,j,k;
	double sum=0,stp1,stp2,stp3;
	double phi1[len*len]={0},phi0[len*len]={0};
	double b1[len*len]={0}, b2[len*len]={0}, bfinal[len*len]={0};
	double sum1=0,temp[len*len]={0},sum2=0,error=1;	
			
	//////////////////////
	//Analytical
	/////////////////////
	
	for(i=0;i<len;i++){
		Tan[len-1][i]=1;
	}
	
	for(j=1;j<len-1;j++){
		for(i=1;i<len-1;i++){
			for (k=1;k<111;k++)
		    {
		      stp1=(1+pow(-1,k+1))/(k);
		      stp2=(sinh(k*pi*(dy*j)))/(sinh(k*pi));
		      stp3=sin(k*pi*(i*dx));
		      sum=sum+(stp1*stp2*stp3);		      
		    }
		    Tan[j][i]=(2/pi)*sum;
		    sum=0;
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
	
	for(i=0;i<((len*len));i++){
		if(i%len==(len-1)){
			b[i]=1;
		}
	}
	
	//////////////////////
	//L,U,N creation
	/////////////////////
	
	//for(i=len;i<((len*len)-len);i++){
	//for(i=len;i<((len*len)-len);i++){
	for(i=0;i<((len*len));i++){
		Lw[i]=Aw[i];
		Ls[i]=As[i];
		if(i==0){
		    Lp[i]=Ap[i];
		}
	    else if((i>=1)&&(i<len)){
		    Lp[i]=Ap[i]-(Ls[i]*Un[i-1]);
	    }
		else{ 
		    Lp[i]=Ap[i]-(Lw[i]*Ue[i-len])-(Ls[i]*Un[i-1]);
		}
		Un[i]=An[i]/Lp[i];
		Ue[i]=Ae[i]/Lp[i];
	}
	
	for(i=0;i<((len*len));i++){
		if(i>=len){
		    Mnw[i]=Lw[i]*Un[i-len];
		}
		if(i>=1){
			Mse[i]=Ls[i]*Ue[i-1];
		}
	}	
	
	//printval(Lp,Ls,Lw,Un,Ue,Mnw,Mse,b);
	//crtL(Lp,Ls,Lw);
	
	while(error>0.0000001){	
	
		//////////////////////
		//b,Nphi0 creation
		/////////////////////
		
		for(i=0;i<(len*len)-len+1;i++){
			b1[i]=Mse[i]*phi0[i+len-1];
		}
		for(i=len-1;i<len*len;i++){
			b2[i]=Mnw[i]*phi0[i-len+1];
		}
		for(i=0;i<len*len;i++){
			bfinal[i]=b[i]+b1[i]+b2[i];
		}
		
		//////////////////////
		//Lower diagonal creation
		/////////////////////
		
		for(i=0;i<len*len;i++){
			if((i>=1)&&(i<len)){
				sum1=(Ls[i]*temp[i-1]);
			}
			else if(i>=len){
				sum1=(Ls[i]*temp[i-1])+(Lw[i]*temp[i-len]);
			}
			//sum1=(Ls[i]*temp[i-1])+(Lw[i]*temp[i-len]);
			temp[i]=(bfinal[i]-sum1)/Lp[i];
			sum1=0;
		}
		
		//////////////////////
		//Upper diagonal creation
		/////////////////////
		
		for(i=(len*len)-1;i>=0;i--){
			if((i<(len*len)-1)&&(i>(len*len)-len)){
				sum2=(Un[i]*phi1[i+1]);
			}
			else if(i<=((len*len)-len)){
				sum2=(Un[i]*phi1[i+1])+(Ue[i]*phi1[i+len]);
			}
			//sum2=(Un[i]*phi1[i+1])+(Ue[i]*phi1[i+len]);
			phi1[i]=(temp[i]-sum2)/1;
			sum2=0;
		}
		
		//////////////////////
		//Phi interchange
		/////////////////////
		
		error=err(phi1,phi0);
		//printf("Error: %f\n",error);
		
		for(i=0;i<len*len;i++){
			phi0[i]=phi1[i];
		}
		
	}
	
	//////////////////////
	//Convert 1D to 2D
	/////////////////////
	
	double T[len][len]={0};
	k=0;
	
	for(i=0;i<len;i++){
	    for(j=0;j<len;j++){
	        T[j][i]=phi1[k];
	        k++;
	    }
	}
	
	print(T);
	//printA(Ap,An,As,Ae,Aw);
	
	////////////////////////////////////
	//PRINTING
	///////////////////////////////////
	
	FILE *fp;
	
	/*fp=fopen("Temp1.txt","w");
	for(i=0;i<len;i++){
		for(j=0;j<len;j++){
			fprintf(fp,"%lf\t%lf\t%lf\t%lf\n",(j*dx),(i*dy),T[i][j],Tan[i][j]);
		}	
	}
	fclose(fp);
	
	int mid=(len-1)*0.5;
	fp=fopen("T mid x y 1.txt","w");
	for(j=0;j<len;j++){
		fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",(j*dx),T[mid][j],T[j][mid],Tan[mid][j],Tan[j][mid]);
	}	
	fclose(fp);*/
	
	return 0;	
}
