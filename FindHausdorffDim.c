// Compiling commands:
// (base) mariacameron@wireless-10-104-29-24 Rayyan % gcc FindHausdorffDim.c -lm -O3
// (base) mariacameron@wireless-10-104-29-24 Rayyan % ./a.out                       


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define KMIN 0
#define KMAX 10
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) <= (b) ? (a) : (b))

struct quart_interval {
	int level;
	int x_ind;
	int y_ind;
	struct quart_interval *child_1;
	struct quart_interval *child_2;
	struct quart_interval *child_3;
	struct quart_interval *child_4;
};

// Hausdorff measure function from Stewart's paper
double Hausdorff_measure(double d,struct quart_interval *node,int kmin,int kmax);

// Checks if the binary interval contains at least one point of the Cantor set
char set_test(struct quart_interval *node, double x0, double y0,double lx,double ly);

// Gets the left endpoint of the binary interval
double get_x_left(struct quart_interval *node);
double get_y_left(struct quart_interval *node);

// allocates memory for binary interval and outputs pointer to it
struct quart_interval *create_node(int level,int x_ind, int y_ind);

void print_node(struct quart_interval *node);

// the driver for the code, must be the last function in the code
int main(void);

// This code computed the Hausdorff dimension of the Cantor set
// The Hausdorff dimension of a set A is defined as 
//
// dim_H(A) = inf{d : H_d(A) = 0} = sup{d : H_d(A) = +\infinity}
//
// where H_d(A) is the Hausdorff measure of the set A defined as
//
// H_d(A) = inf_{delta-->0}inf_{{U_j}} \sum_j (diam U_j)^d,
//
// where the infimum is taken over all covers of A by sets U_j with diameter <= delta
//
// In this code, the test set is the Cantor set in 1D:
//
// C = [0,1] \ Union_{n=0}^{infinity} Union_{k=0}^{3^n-1} ((3k+1)/3^{n+1},(3k+2)/3^{n+1})

double get_x_left(struct quart_interval *node) {
	double left_x = node->x_ind/pow(2,node->level);
	return left_x;
}

double get_y_left(struct quart_interval *node) {
	double down_y = node->y_ind/pow(2,node->level);
	return down_y;
}
char set_test(struct quart_interval *node,double x0,double y0,double lx,double ly) {
// tests if the binary interval belongs or not to the Cantor set
	int n, NMAX;
	//print_node(node);
	double left_x = get_x_left(node);
	double right_x = left_x + 1.0/pow(2,node->level);
	double down_y = get_y_left(node);
	double up_y = down_y + 1.0/pow(2,node->level);
	double test_l = 1/pow(2, node->level);
	double test = 1;
	double t_x,t_y;
	double flag = 0;
	double ind = pow(2,16-node->level)*node->x_ind;
	double rind = ind + pow(2,16-node->level);
	int y;
	
	// to do for BM: create a function which outputs the y value of the BM at inputted x values
	// then we input both left_x and right_y into the BM and check if the given y values fall
	// in our range for y
	// imports x from matlab to an array BM
	FILE *fid;
	fid = fopen("brownMot.txt","r");
	int N = pow(2,16);
	double x;
	double BM[N];
	int j;
	for( j = 0; j < N; j++) {
	  fscanf(fid,"%le\n",&x);
	  BM[j] = x;
	}
	fclose(fid);
	for(y=ind;y<rind;y++){
		BM[y] = pow(BM[y]/300,2); //rescale to [0,1], squared to make positive *dependant on graph!*
		if(down_y > BM[y] || up_y < BM[y]){
			test = 0;}
	    else{
		 test = 1; 
		 break;
		 } 
	}
	return test;
}

double Hausdorff_measure(double d,struct quart_interval *node,int kmin,int kmax) {
	int level = node->level,x_ind = node->x_ind,y_ind = node->y_ind;
	double measure = 0.0, min_measure = 1.0/pow(2,d*kmax);
	char test_node;
	// return test_node = 0 is the interval contains no points of Cantor's set
	// return test_node = 1 otherwise

    // test if node intersects with the set
    test_node = set_test(node,0,0,1,1);
    if( test_node == 1 ) {	
		if( level == kmax ) {
			measure = min_measure;
			free(node);
			return measure;
		}
    	node->child_1 = create_node(level + 1,2*x_ind,2*y_ind);
		node->child_2 = create_node(level + 1,2*x_ind + 1,2*y_ind);
		node->child_3 = create_node(level + 1,2*x_ind,2*y_ind + 1);
		node->child_4 = create_node(level + 1,2*x_ind + 1,2*y_ind + 1);
    	
    	measure = Hausdorff_measure(d,node->child_1,kmin,kmax) + 
    			  Hausdorff_measure(d,node->child_2,kmin,kmax) + 
				  Hausdorff_measure(d,node->child_3,kmin,kmax) + 
    			  Hausdorff_measure(d,node->child_4,kmin,kmax);
    	// Masha is not sure about the need for this if	(needed for countable sets?)	  
    	 if(level >= kmin ) {
		//	if(measure < 1.0/pow(2,d*level)){
			measure = measure;		   
			}
		//	else{ 
		//		measure = 1.0/pow(2,d*level);
		//	}
    	//}
    }
	if( test_node == 0) {
		measure = 0;
	}
    free(node);
    return measure;
}

struct quart_interval *create_node(int level,int x_ind, int y_ind){
	struct quart_interval *node;
	
	// allocate memory for a single binary interval
	// required size of memory piece: sizeof(struct binary_interval)
	// request to allocate it: malloc(...)
	// request to structure it for recording the datatype struct binary_interval *
	node = (struct quart_interval *)malloc(sizeof(struct quart_interval));
	node->level = level;
	node->x_ind = x_ind;
	node->y_ind = y_ind;
	node->child_1 = NULL;
	node->child_2 = NULL;
	node->child_3 = NULL;
	node->child_4 = NULL;
	
	
	return node;
}

// function for debugging
 void print_node(struct quart_interval *node){
	printf("node: level = %i\n",node->level);
	printf("node: x_ind = %i\n",node->x_ind);
	printf("node: y_ind = %i\n",node->y_ind);
	//printf("node: child1 = %i\n",&(node->child_1));
	//printf("node: child2 = %i\n",node->child_2);
	//printf("node: child3 = %i\n",node->child_3);
	//printf("node: child4 = %i\n",node->child_4);
}

int main(void) {
	struct quart_interval *node;
	double d,dmin = 1.45,dmax = 1.55,dstep = 5.0e-2;	
	double Haus;  // Hausdorff measure for the current dimension
	int kmin = KMIN, kmax;
    clock_t CPUbegin; // for measuring CPU time
    double cpu; // for recording CPU time
    FILE *fid; // 
    
    fid = fopen("Cantor.txt","w");
	for( d = dmin; d <= dmax; d+=dstep ) {
		printf("d = %.4f\n",d);
		fprintf(fid,"d = %.4f\n",d);
		for( kmax = 5; kmax <= KMAX; kmax++ ) {
		    CPUbegin=clock(); // start time measurement
			node = create_node(0,0,0); // create node of level 0 x index 0 and y index 0
			 // Hausdorff measure of the Cantor set of dimension d
			Haus = Hausdorff_measure(d,node,kmin,kmax);
		    cpu = (clock()-CPUbegin)/((double)CLOCKS_PER_SEC);	// end time measurement		
			printf("kmax = %i\t Hmeasure = %.6e \t CPU time = %g\n",kmax,Haus,cpu);
			fprintf(fid,"%i\t %.6e\t %g\n",kmax,Haus,cpu);
		}			
	}
	fclose(fid);
	return 0;
}








