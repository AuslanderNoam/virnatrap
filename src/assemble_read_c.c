#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>

float SCORE_THR = 0.5;  /*lower threshold for a contig to be considered viral*/
int RUNS = 100; /*maximum number of rounds of extending a contig in one directorion*/
int SUBLEN = 24; 
int MAX_SIZE = 2560000000;
//int NUMEL = 0;
//float TOTSC = 0;


/*Takes as input string s, which is trusted to be of length len
 *allocates a new string of size 2* len +1 and copies s to the beginning of  
 *the new string*/
char *doubleString(char *s, int len)
{
    char *s2 = malloc(2 * len + 1);
    strcpy(s2, s); //Assumes accuracy of len
    return s2;
}

/*structure to return a string and some ancillary information
  c is the string; len is the length without the null terminator; numel is number of scores summed; totsc is the total score*/
struct ret {
    char *c;
    int len;
    int numel;
    float totsc;
};

/*find_sub_right  tries to extend a contig on the right
 * reads is a two-dimensional array of all reads
 * sb0 is the end piece to extend
 * read_used could be a bit array to indicate which reads have been used, but is currently ignored
 * num_reads is the number of reads, equivalently the size of the major dimension of **reads */ 
int * find_sub_right(char **reads, char *sb0, int *read_used, int num_reads) {


	int pos;
	char *ptr;

	static int cmin[2] = {100,100};
	
	cmin[1]=100;
	cmin[0]=100;	
	int i=0, cidx=0;
	while (i < num_reads)
	{	
		

		if (!(read_used[i]))
		{	
			ptr = strstr(reads[i], sb0);
			if (ptr != NULL) { 
			  read_used[i] = 1;
			  pos = ptr-reads[i];
			  if(pos<cmin[0] && pos>=0)
			    {	
			      cmin[0] = pos;
			      cmin[1] = i;
			    }
			}
						
		}
		i++;
	}



	return cmin;
}

/*find_sub_left  tries to extend a contig on the left
 *reads is a two-dimensional array of all reads
 *sb0 is the end piece to extend
 *read_used could be a bit array to indicate which reads have been used, but is currently ignored
 *num_reads is the number of reads, equivalently the size of the major dimension of **reads */ 
int * find_sub_left(char **reads, char *sb0, int *read_used, int num_reads) {
	int pos;
	char *ptr;

	static int cmax[2] = {-100,-100};

	cmax[1]=-100;
	cmax[0]=-100;	
	
	int i=0;
	while (i < num_reads)
	{		

		if (!(read_used[i]))
		{	
			
			ptr = strstr(reads[i], sb0);
			if (ptr != NULL) { 
				read_used[i] = 1;

			  	pos = ptr-reads[i];
			  	if(pos>cmax[0] && pos>=0 && pos>=0)
			    	{		
			      		cmax[0] = pos;
			      		cmax[1] = i;
			    	}
			}
		}
		i++;

	}
	
	return cmax;
}

/*assemble_right builds a contig starting from read and extending to the right; read_list is the two-dimensional array of reads
 *score is the total score; score_list is the score of each read; read)used is a bit array of which reads have been used
 *num_reads is the number of reads or the size of  upper dimension of read_list
 *overall logic is similar to assemble_right but the management of string memory to extend the contig is different*/ 
struct ret assemble_right(char read[], char *read_list[],float score, float *score_list, int *read_used, int SEGMENT_LENGTH, int num_reads){
	
	float TOTSC = score;
	int NUMEL = 1;
	TOTSC = score;
	NUMEL = 1;

	char* contig;
        int clen = strlen(read);
        int ccap = 2*clen; /*current upper bound on the length of the contig array*/ 
	contig = doubleString(&(read[0]), clen);
	
	char sb0[SUBLEN+1];
	strncpy(sb0,&read[SEGMENT_LENGTH-SUBLEN], SUBLEN);
        sb0[SUBLEN] = '\0';   //need to add a termination character because strncpy does not
	bool flag = true;
	int cnt = 0;
	int *poss;

	while (flag==true && cnt<RUNS) //keep trying to extend the contig as long as there is a possible extension and average score per element stays above 0.5
	{	

		poss = find_sub_right(read_list, sb0, read_used, num_reads);

		if (poss[1]<100 && poss[0]< (SEGMENT_LENGTH-SUBLEN) && poss[1]>=0 && TOTSC/NUMEL>0.5){	
			read_used[poss[1]] = 1;
				
                        //If the length of the contig is getting close to ccap, then double the available length and free the memory for the shorter array
			if (clen + SEGMENT_LENGTH +1 > ccap) {
				ccap = clen * 2;
				char *new_contig = doubleString(contig, clen);
                                free(contig);
                                contig = new_contig;
			}
			clen += SEGMENT_LENGTH;                        
			//Append first SEG_LEN chars of read_list[...] to contig
           		//strncat(contig, & read_list[poss[1]][0], SEGMENT_LENGTH); //Here was the error
			strncat(contig, & read_list[poss[1]][poss[0]+strlen(sb0)], SEGMENT_LENGTH);


			strncpy(sb0, &contig[clen-SUBLEN], SUBLEN);
			sb0[SUBLEN] = '\0'; //need to add a termination character because strncpy does not
                        NUMEL++;	
			TOTSC+=poss[1]; //update  total score
			cnt++;	//update number of loop iterations

		}
		else{
			flag = false;
		}	
	}
        struct ret rres = { contig, clen, NUMEL, TOTSC };
	return rres;
}

/*assemble_left builds a contig starting from read and extending to the right; read_list is the two-dimensional array of reads
 *score is the total score; score_list is the score of each read; read)used is a bit array of which reads have been used
 *num_reads is the number of reads or the size of  upper dimension of read_list
 *overall logic is similar to assemble_right but the management of string memory to extend the contig is different*/ 
struct ret assemble_left(char read[], char *read_list[],float score, float *score_list,int *read_used, int SEGMENT_LENGTH, int num_reads, int NUMEL, float TOTSC){

	TOTSC += score;
	NUMEL++;
	int clen = 0;
	char* contig = calloc(1, sizeof(char));
		
	char sb0[SUBLEN+1];
	strncpy(sb0,read, SUBLEN);
        sb0[SUBLEN] = '\0'; //need to add a termination character because strncpy does not	
	bool flag = true;
	int cnt = 0;
	int *poss;

	while (flag==true && cnt<RUNS) //keep trying to extend the contig as long as there is a possible extension and average score per element stays above 0.5
	{	

		poss = find_sub_left(read_list, sb0, read_used, num_reads);
		if (poss[1]>-100 && poss[0]> 0 && poss[0]< 100 && poss[1]>=0 && TOTSC/NUMEL>0.5){
			read_used[poss[1]] = 1;		
			
                        //Allocate new space for the contig because there is no easy way to prepend the added read on the left
			char *new_contig = malloc((poss[0] + clen + 1) * sizeof(char));			
                        //Copy the new piece on the left
                        strncpy(new_contig, & read_list[poss[1]][0], poss[0]);
                        new_contig[poss[0]] = '\0'; // add termination character because strncpy does not
                        //Assuming SEGMENT_LENGTH >= poss[0], this did nothing
                        //new_contig[SEGMENT_LENGTH] = '\0';
			          
			clen += poss[0]; 

                        //Copy the old contig on the right
			strcat(new_contig+poss[0], contig);
			
			free(contig);
                        contig = new_contig;
	
			strncpy(sb0,&read_list[poss[1]][0], SUBLEN);
			sb0[SUBLEN] = '\0'; // add termination character because strncpy does not
                        NUMEL++;	
			TOTSC+=poss[1]; //update  total score
			cnt++;	//update number of loop iterations

		}
		else{
			flag = false;
		}	
	}
	struct ret lres = { contig, clen, NUMEL, TOTSC };
        return lres;
}


char* assemble_read(float *f_arr, char *ch_arr[], char *read, float score, int SEGMENT_LENGTH, int num_reads, int *read_used){
	
	
	int i;
        char *fullcont;
	
	
	struct ret rres = assemble_right(read, ch_arr, score, f_arr, read_used, SEGMENT_LENGTH, num_reads);
	//The current TOTSC & NUMEL values are passed as fields in rres to match previous behavior
        struct ret lres = assemble_left(read, ch_arr, score, f_arr, read_used, SEGMENT_LENGTH, num_reads, rres.numel, rres.totsc);
	//free(read_used);
        //Some diagnostics now commented out
	// fprintf(stdout,"Got left contig %s with expected length %d and actual length %d\n",lres.c, lres.len, strlen(lres.c));
        // fprintf(stdout,"Got right contig %s with expected length %d and actual length %d\n",rres.c, rres.len, strlen(rres.c));
        fullcont = malloc((lres.len + rres.len + 1) * sizeof(char));	
        strcpy(fullcont, lres.c);
        strcpy(fullcont + lres.len, rres.c);
        free(lres.c);
	free(rres.c);
	
        float av_score = (lres.totsc)/(lres.numel);

        //Diagnostics now commented out
	//fprintf(stdout,"About to return %s with expected length %d and actual length %d\n",fullcont, lres.len + rres.len, strlen(fullcont));
	return fullcont;
}



void remove_used_reads(char *ch_arr[], char *ch_arr2[], char* contig, int num_reads, int nvr, int *read_used, int *read_used2){
	int i;
	char *ptr;

   	for(i = 0; i < num_reads; i++)
	{
		if (!(read_used[i])){
			ptr = strstr(contig, ch_arr[i]);
			if (ptr != NULL) { 
				read_used[i] = 1;
			}
		}
	}
	
	for(i = 0; i < nvr; i++)
	{
		if (!(read_used2[i])){
			ptr = strstr(contig, ch_arr2[i]);
			if (ptr != NULL) { 
				read_used2[i] = 1;
			}
		}
	}
    		
	
}

int assemble_read_loop(float *f_arr,float *f_arr2, char *ch_arr[], char *ch_arr2[], int SEGMENT_LENGTH, int num_reads,int nvr, char *fname){
	
	int * read_used;
	int * read_used2;
	read_used = (int *) malloc(num_reads * sizeof(int));
	int i;
	for(i = 0; i < num_reads; i++)
    		read_used[i] = false;

	read_used2 = (int *) malloc(nvr * sizeof(int));

	for(i = 0; i < nvr; i++)
    		read_used2[i] = false;

	FILE *fp;
	int wi = 0;
	fp = fopen(fname, "w+");
	char *contig;
	char vi[SEGMENT_LENGTH+1];
	

	for(i = 0; i < nvr ; i++){
		if (!(read_used2[i])){
			strncpy(vi, ch_arr2[i], SEGMENT_LENGTH);
			vi[SEGMENT_LENGTH] = '\0';
			contig = assemble_read(f_arr, ch_arr, vi, f_arr2[i], SEGMENT_LENGTH, num_reads, read_used);
			remove_used_reads(ch_arr, ch_arr2, contig, num_reads, nvr, read_used, read_used2);
			fprintf(fp,">contig_%d[]\n", i);
			fprintf(fp, contig);
			fprintf(fp,"\n");
		}


	}
	fclose(fp);
	return 1;
}

