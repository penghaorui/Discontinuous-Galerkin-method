
/*this file contains the function for reading and writing data to files*/

/*prototyping functions*/
void ReadDouble(char path[], double *matrix, long int Row, long int Col);

void ReadInt(char path[], long int *matrix, long int Row, long int Col);

void Save_Int(char path[], long *matrix, long int Row, long int Col);

void Save_Record(char path[], double *matrix, long int Row);

void epartinput(char path[], long int *matrix, long int Num);

void PrintDouble(double *matrix, long int Row, long int Col, long int start, long int end);

void PrintInt(long int *matrix, long int Row, long int Col, long int start, long int end);


/************************************************************************************************/
/*read double precision matrix*/
void ReadDouble(char path[], double *matrix, long int Row, long int Col)
{
  long int i;
  long int j;
  long int k;
  long int num;
  long int sym1;
  long int sym2;
  long int power;
  long int m;
  long int n;
  long int fd;
  long int size;
    
  static char s[100000000];
  extern int errno;
  
  fd = open(path,O_RDONLY);

  if(fd!=-1)
  {
    //printf("opened file %s.\n",path);
  }
  else
  {
    printf("can't open file %s.\n",path);
    printf("errno: %d\n",errno);
    //printf("ERR : %s\n",strerror(errno));   
  }

  size = read(fd,s,sizeof(s));
  close(fd);
  
  num = 0;
  m = 0;
  n = 0;
  
  i = 0;
  while(i<100000000)
  { 
       
    if(s[i]==101)
    { 
      if(s[i-10]==32)
      {
        sym1 = 1;
      }
      if(s[i-10]==45)
      {
        sym1 = -1;
      }
      
      if(s[i+1]==43)
      {
        sym2 = 1;
      }
      if(s[i+1]==45)
      {
        sym2 = -1;
      }
            
      *(matrix+m*Col+n) = (s[i-9]-48);        
          
      k = -7;
      for(j=i-1;j>=i-7;j--)
      {
        *(matrix+m*Col+n) = *(matrix+m*Col+n) + (s[j]-48)*pow(10,k);
        k++;        
      }
      
      power = sym2*((s[i+2]-48)*10 + (s[i+3]-48)*1);
      *(matrix+m*Col+n) = *(matrix+m*Col+n) * pow(10,power) * sym1;
      
      n = n + 1;
      num = num + 1;
    }
            
    if(num==Row*Col)
    {
      break;    
    }
    
    if(s[i]==10)
    {
      m = m + 1;
      n = 0;
    }
            
    i++;
  }
  
} 

/*read Int precision matrix*/
void ReadInt(char path[], long int *matrix, long int Row, long int Col)
{
  long int i;
  long int j;
  long int k;
  long int num;
  long int sym1;
  long int sym2;
  long int power;
  long int m;
  long int n;
  long int fd;
  long int size;
    
  static char s[100000000];
  extern int errno;
  
  fd = open(path,O_RDONLY);

  if(fd!=-1)
  {
    //printf("opened file %s.\n",path);
  }
  else
  {
    printf("can't open file %s.\n",path);
    printf("errno: %d\n",errno);
    //printf("ERR : %s\n",strerror(errno));   
  }

  size = read(fd,s,sizeof(s));
  close(fd);
  
  num = 0;
  m = 0;
  n = 0;
  
  i = 0;
  while(i<100000000)
  { 
       
    if(s[i]==101)
    { 
      if(s[i-10]==32)
      {
        sym1 = 1;
      }
      if(s[i-10]==45)
      {
        sym1 = -1;
      }
      
      if(s[i+1]==43)
      {
        sym2 = 1;
      }
      if(s[i+1]==45)
      {
        sym2 = -1;
      }
      
      power = sym2*((s[i+2]-48)*10 + (s[i+3]-48)*1);      
      *(matrix+m*Col+n) = (s[i-9]-48) * pow(10,power) * sym1;        
                 
      k = power - 1;
      for(j=i-7;j<=i-7+power-1;j++)
      {              
        *(matrix+m*Col+n) = *(matrix+m*Col+n) + (s[j]-48)*pow(10,k);                
        k--;               
      }                
      n = n + 1;
      num = num + 1;
    }
                
    if(num==Row*Col)
    {
      break;    
    }
    
    if(s[i]==10)
    {
      m = m + 1;
      n = 0;
    }
            
    i++;
  }
  
}

/*save Int precision matrix*/

void Save_Int(char path[], long *matrix, long int Row, long int Col)
{
   FILE *fp;
   long int i,j;

   if((fp = fopen(path,"wb"))==NULL)
   {
      printf("cannot open this file!");
      exit(0);
   }
   
   for(i = 0; i < Row; i++)
   {  
     for(j = 0; j < Col; j++)
     {   
        fprintf(fp,"%ld ", *(matrix+i*Col+j));
     }
     fprintf(fp,"\n");

   } 

  fclose(fp);
}

/*save double precision matrix*/
void Save_Data(char path[], double *matrix, long int Row, long int Col)
{
   FILE *fp;
   long int i,j;

   if((fp = fopen(path,"wb"))==NULL)
   {
      printf("cannot open this file!");
      exit(0);
   }
   
   for(i = 0; i < Row; i++)
   {  
     for(j = 0; j < Col; j++)
     {   
        fprintf(fp,"%0.14f ", *(matrix+i*Col+j));
     }
     fprintf(fp,"\n");

   } 

  fclose(fp);
}



/*save double precision matrix for records*/
void Save_Record(char path[], double *matrix, long int Row)
{
   FILE *fp;
   long int i,j;

   if((fp = fopen(path,"wb"))==NULL)
   {
      printf("cannot open this file!");
      exit(0);
   }
   
   for(i = 0; i < Row; i++)
   {  

        fprintf(fp,"%0.14f ", *(matrix+i));

     fprintf(fp,"\n");

   } 

  fclose(fp);
}

/*read partition file output from metis*/
void epartinput(char path[], long int *matrix, long int Num)
{
	long i,j;
	int temp;
	FILE *fp;
	fp=fopen(path,"r");

	for (i=0;i<Num;i++)
	{
		fscanf(fp,"%d\n",&temp);
		matrix[i]=(long int)temp;
	}
}

/*print double precision matrix*/
void PrintDouble(double *matrix, long int Row, long int Col, long int start, long int end)
{
  long int m;
  long int n;
  for(m=start-1;m<=end-1;m++)
  {
    for(n=0;n<Col;n++)
    {
      printf("%0.8e ",*(matrix+m*Col+n));
    }
    printf("  line  %ld",m+1);
    printf("\n");
  }
}

/*print Int precision matrix*/
void PrintInt(long int *matrix, long int Row, long int Col, long int start, long int end)
{
  long int m;
  long int n;
  for(m=start-1;m<=end-1;m++)
  {
    for(n=0;n<Col;n++)
    {
      printf("%ld ",*(matrix+m*Col+n));
    }
    printf("  line  %ld",m+1);
    printf("\n");
  }
}
